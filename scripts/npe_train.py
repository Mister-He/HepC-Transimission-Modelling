#!/usr/bin/env python3
"""Train a Neural Posterior Estimator (NPE_C / SNPE, MAF) for the HCV model.

Multi-round sequential NPE (--round r, final round = --n-rounds):
    --round 1  train on prior draws; pickle the inference object; write
               proposal_theta_round2{suffix}.csv (next proposal draws).
    --round r>1  load the pickled round-(r-1) inference; append round-r
               simulations with proposal = round-(r-1) posterior; train;
               if r == n_rounds, sample the final posterior + SBC; else
               pickle and write the next proposal.

Transforms (documented in bayes_methodology.md):
    z = theta / 2                 (regularized Normal prior representation)
    x_t = [logit(p_1..6), log(N_1..6)] standardised by round-1 stats,
          clipped to [-8, 8]

Outputs in --out:
    posterior_samples_npe{suffix}.csv   final posterior draws (raw theta)
    sbc_summary{suffix}.csv             per-parameter SBC coverage
    loss_curve{suffix}.csv              training/validation loss
    npe_config{suffix}.json             seeds, transforms, settings
    npe_inference_round1{suffix}.pkl    pickled round-1 inference
    proposal_theta_round2{suffix}.csv   (round 1 only)
"""
import argparse
import csv
import json
import math
import os
import pickle

import numpy as np
import torch

from sbi.inference import NPE_C
from sbi.neural_nets import posterior_nn
from sbi.utils import MultipleIndependent
from torch.distributions import Normal, StudentT


def read_numeric_csv(path):
    """Read the numeric CSVs emitted by R without requiring pandas."""
    return np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)


def write_numeric_csv(path, values, columns):
    values = np.asarray(values)
    if values.ndim == 1:
        values = values[:, None]
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(columns)
        writer.writerows(values.tolist())


def add_observation_noise(x_raw, n_prev):
    """Approximate the calibration observation model for NPE training.

    The R simulator returns noiseless summaries. Training directly on those
    summaries makes NPE learn an unrealistically sharp inverse map. Inject
    the same uncertainty represented by the calibration likelihood: Binomial
    noise for age-specific prevalence and log-normal noise for population counts.
    """
    x = np.asarray(x_raw, dtype=np.float64).copy()
    n_prev = np.asarray(n_prev, dtype=np.int64)
    p = np.clip(x[:, :6], 1e-8, 1 - 1e-8)
    x[:, :6] = np.random.binomial(n_prev[None, :], p) / n_prev[None, :]
    x[:, 6:12] = np.exp(np.log(np.maximum(x[:, 6:12], 1e-8)) +
                         np.random.normal(0.0, 0.10, size=x[:, 6:12].shape))
    return x


def parse_vector(value, name):
    values = [float(v) for v in value.split(",")]
    if len(values) != 6:
        raise ValueError(f"{name} must contain six comma-separated values")
    return values


def build_prior(contact_anchor, beta_anchor, contact_sd, beta_sd, beta_df):
    dists = []
    for a in contact_anchor:
        dists.append(Normal(
            torch.tensor([math.log(a) / 2.0]),
            torch.tensor([contact_sd / 2.0]),
        ))
    for a in beta_anchor:
        if math.isfinite(beta_df):
            dists.append(StudentT(
                df=torch.tensor([beta_df]),
                loc=torch.tensor([math.log(a) / 2.0]),
                scale=torch.tensor([beta_sd / 2.0]),
            ))
        else:
            dists.append(Normal(
                torch.tensor([math.log(a) / 2.0]),
                torch.tensor([beta_sd / 2.0]),
            ))
    return MultipleIndependent(dists)


def load_config(args, suffix):
    cfg_path = os.path.join(args.out, f"npe_config{suffix}.json")
    with open(cfg_path) as fh:
        return json.load(fh)


def make_x_transform(x_train, cfg=None):
    if cfg is None:
        p = np.clip(x_train[:, :6], 1e-6, 1 - 1e-6)
        xt = np.hstack([np.log(p / (1 - p)), np.log(x_train[:, 6:])])
        mu = xt.mean(axis=0)
        sd = xt.std(axis=0) + 1e-8
    else:
        mu = np.array(cfg["mu_x"])
        sd = np.array(cfg["sd_x"])
    def trans(mat):
        pm = np.clip(mat[:, :6], 1e-6, 1 - 1e-6)
        y = np.hstack([np.log(pm / (1 - pm)), np.log(mat[:, 6:])])
        return np.clip((y - mu) / sd, -8.0, 8.0)
    return trans, mu, sd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--theta", required=True)
    ap.add_argument("--x", required=True)
    ap.add_argument("--x-obs", required=True)
    ap.add_argument("--theta-sbc", required=True)
    ap.add_argument("--x-sbc", required=True)
    ap.add_argument("--theta2", default=None)
    ap.add_argument("--x2", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--n-draws", type=int, default=60000)
    ap.add_argument("--n-sbc", type=int, default=300)
    ap.add_argument("--sbc-draws", type=int, default=1000)
    ap.add_argument("--n-proposal", type=int, default=20000)
    ap.add_argument("--round", type=int, default=1)
    ap.add_argument("--n-rounds", type=int, default=2)
    ap.add_argument("--suffix", default="")
    ap.add_argument("--max-epochs", type=int, default=150)
    ap.add_argument("--contact-anchor", required=True)
    ap.add_argument("--beta-anchor", required=True)
    ap.add_argument("--contact-sd", type=float, required=True)
    ap.add_argument("--beta-sd", type=float, required=True)
    ap.add_argument("--beta-df", type=float, default=float("inf"))
    ap.add_argument("--hidden-features", type=int, default=64)
    ap.add_argument("--num-transforms", type=int, default=8)
    args = ap.parse_args()

    suffix = ("_" + args.suffix) if args.suffix else ""
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    contact_anchor = parse_vector(args.contact_anchor, "--contact-anchor")
    beta_anchor = parse_vector(args.beta_anchor, "--beta-anchor")
    prior = build_prior(contact_anchor, beta_anchor, args.contact_sd,
                        args.beta_sd, args.beta_df)
    device = "mps" if torch.backends.mps.is_available() else "cpu"
    print(f"[npe] round={args.round}, seed={args.seed}, device={device}")

    x_obs_raw = read_numeric_csv(args.x_obs).astype(np.float64)[0]
    n_prev = np.rint(x_obs_raw[6:12]).astype(np.int64)

    if args.round == 1:
        theta_raw = read_numeric_csv(args.theta).astype(np.float64)
        x_raw = read_numeric_csv(args.x).astype(np.float64)
        x_raw = add_observation_noise(x_raw, n_prev)
        z_train = theta_raw / 2.0
        trans_x, mu_x, sd_x = make_x_transform(x_raw)
        x_z = trans_x(x_raw)

        density_estimator = posterior_nn(
            "maf", hidden_features=args.hidden_features,
            num_transforms=args.num_transforms,
        )
        inference = NPE_C(prior=prior, density_estimator=density_estimator,
                          device=device)
        estimator = inference.append_simulations(
            torch.as_tensor(z_train, dtype=torch.float32),
            torch.as_tensor(x_z, dtype=torch.float32),
        ).train(
            training_batch_size=256,
            max_num_epochs=args.max_epochs,
            validation_fraction=0.1,
            show_train_summary=False,
        )
        pkl_path = os.path.join(args.out, f"npe_inference_round1{suffix}.pkl")
        with open(pkl_path, "wb") as fh:
            pickle.dump(inference, fh)

        posterior = inference.build_posterior(estimator)
        posterior.set_default_x(torch.as_tensor(trans_x(x_obs_raw[None, :])[0],
                                                dtype=torch.float32))
        if args.n_rounds == 1:
            # A one-round run is useful as a robustness check: unlike the
            # sequential proposal, it retains prior support in conflict
            # regions instead of concentrating round-2 simulations early.
            z_draws = posterior.sample((args.n_draws,)).numpy()
            write_numeric_csv(
                os.path.join(args.out,
                             f"posterior_samples_npe{suffix}.csv"),
                z_draws * 2.0, [f"theta{i}" for i in range(1, 13)])
            theta_sbc_raw = read_numeric_csv(args.theta_sbc).astype(np.float64)
            x_sbc_raw = add_observation_noise(
                read_numeric_csv(args.x_sbc).astype(np.float64), n_prev)
            z_sbc = theta_sbc_raw / 2.0
            x_z_sbc = trans_x(x_sbc_raw)
            n_sbc = min(args.n_sbc, len(z_sbc))
            ranks = np.zeros((n_sbc, 12))
            for j in range(n_sbc):
                draws = posterior.sample(
                    (args.sbc_draws,),
                    x=torch.as_tensor(x_z_sbc[j], dtype=torch.float32)).numpy()
                ranks[j] = (draws < z_sbc[j]).mean(axis=0)
            cov95 = ((ranks >= 0.025) & (ranks <= 0.975)).mean(axis=0)
            cov90 = ((ranks >= 0.05) & (ranks <= 0.95)).mean(axis=0)
            with open(os.path.join(
                    args.out, f"sbc_summary{suffix}.csv"), "w", newline="") as fh:
                writer = csv.writer(fh)
                writer.writerow(["parameter", "coverage_95", "coverage_90",
                                 "mean_rank", "sd_rank"])
                for i in range(12):
                    writer.writerow([f"theta{i + 1}", cov95[i], cov90[i],
                                     ranks[:, i].mean(), ranks[:, i].std()])
            np.save(os.path.join(args.out, f"sbc_ranks{suffix}.npy"), ranks)
            config = {
                "seed": args.seed, "round": 1, "n_train": len(z_train),
                "n_draws": args.n_draws, "device": device,
                "density_estimator": "maf",
                "hidden_features": args.hidden_features,
                "num_transforms": args.num_transforms,
                "transform": "theta/2; logit(p),log(N)",
                "prior_contact_anchor": contact_anchor,
                "prior_beta_anchor": beta_anchor,
                "prior_contact_sd": args.contact_sd,
                "prior_beta_sd": args.beta_sd,
                "prior_beta_df": args.beta_df,
                "observation_noise": "Binomial age p; logNormal N sigma=0.10",
                "mu_x": mu_x.tolist(), "sd_x": sd_x.tolist(),
                "max_epochs": args.max_epochs,
            }
            with open(os.path.join(args.out,
                                   f"npe_config{suffix}.json"), "w") as fh:
                json.dump(config, fh, indent=2)
            print(f"[npe] one-round final posterior -> "
                  f"posterior_samples_npe{suffix}.csv")
            return
        prop = posterior.sample((args.n_proposal,)).numpy()
        write_numeric_csv(
            os.path.join(args.out, f"proposal_theta_round2{suffix}.csv"),
            prop * 2.0, [f"theta{i}" for i in range(1, 13)])

        config = {
            "seed": args.seed, "round": 1, "n_train": len(z_train),
            "n_proposal": args.n_proposal, "device": device,
            "density_estimator": "maf",
            "hidden_features": args.hidden_features,
            "num_transforms": args.num_transforms,
            "transform": "theta/2; logit(p),log(N)",
            "prior_contact_anchor": contact_anchor,
            "prior_beta_anchor": beta_anchor,
            "prior_contact_sd": args.contact_sd,
            "prior_beta_sd": args.beta_sd,
            "prior_beta_df": args.beta_df,
            "observation_noise": "Binomial age p; logNormal N sigma=0.10",
            "mu_x": mu_x.tolist(), "sd_x": sd_x.tolist(),
            "max_epochs": args.max_epochs,
        }
        with open(os.path.join(args.out, f"npe_config{suffix}.json"), "w") as fh:
            json.dump(config, fh, indent=2)
        print(f"[npe] round 1 done; proposal -> proposal_theta_round2{suffix}.csv")
        return

    # ---------------- round r > 1 ----------------
    cfg = load_config(args, suffix)
    trans_x, _, _ = make_x_transform(x_obs_raw[None, :], cfg)
    prev_round = args.round - 1
    with open(os.path.join(args.out,
                          f"npe_inference_round{prev_round}{suffix}.pkl"),
              "rb") as fh:
        inference = pickle.load(fh)
    proposal = inference.build_posterior()
    proposal.set_default_x(torch.as_tensor(trans_x(x_obs_raw[None, :])[0],
                                           dtype=torch.float32))

    theta_r = read_numeric_csv(args.theta2).astype(np.float64) / 2.0
    x_r_raw = read_numeric_csv(args.x2).astype(np.float64)
    x_r_raw = add_observation_noise(x_r_raw, n_prev)
    x_r_z = trans_x(x_r_raw)

    inference.append_simulations(
        torch.as_tensor(theta_r, dtype=torch.float32),
        torch.as_tensor(x_r_z, dtype=torch.float32),
        proposal=proposal,
    )
    estimator = inference.train(
        training_batch_size=256,
        max_num_epochs=args.max_epochs,
        validation_fraction=0.1,
        force_first_round_loss=True,
        show_train_summary=False,
    )
    posterior = inference.build_posterior(estimator)
    posterior.set_default_x(torch.as_tensor(trans_x(x_obs_raw[None, :])[0],
                                            dtype=torch.float32))

    if args.round < args.n_rounds:
        pkl_path = os.path.join(args.out,
                                f"npe_inference_round{args.round}{suffix}.pkl")
        with open(pkl_path, "wb") as fh:
            pickle.dump(inference, fh)
        prop = posterior.sample((args.n_proposal,)).numpy()
        write_numeric_csv(
            os.path.join(args.out,
                         f"proposal_theta_round{args.round + 1}{suffix}.csv"),
            prop * 2.0, [f"theta{i}" for i in range(1, 13)])
        cfg["round"] = args.round
        with open(os.path.join(args.out, f"npe_config{suffix}.json"), "w") as fh:
            json.dump(cfg, fh, indent=2)
        print(f"[npe] round {args.round} done; proposal -> "
              f"proposal_theta_round{args.round + 1}{suffix}.csv")
        return

    z_draws = posterior.sample((args.n_draws,)).numpy()
    write_numeric_csv(
        os.path.join(args.out, f"posterior_samples_npe{suffix}.csv"),
        z_draws * 2.0, [f"theta{i}" for i in range(1, 13)])

    theta_sbc_raw = read_numeric_csv(args.theta_sbc).astype(np.float64)
    x_sbc_raw = read_numeric_csv(args.x_sbc).astype(np.float64)
    x_sbc_raw = add_observation_noise(x_sbc_raw, n_prev)
    z_sbc = theta_sbc_raw / 2.0
    x_z_sbc = trans_x(x_sbc_raw)
    n_sbc = min(args.n_sbc, len(z_sbc))
    ranks = np.zeros((n_sbc, 12))
    for j in range(n_sbc):
        draws = posterior.sample(
            (args.sbc_draws,),
            x=torch.as_tensor(x_z_sbc[j], dtype=torch.float32)).numpy()
        ranks[j] = (draws < z_sbc[j]).mean(axis=0)
    cov95 = ((ranks >= 0.025) & (ranks <= 0.975)).mean(axis=0)
    cov90 = ((ranks >= 0.05) & (ranks <= 0.95)).mean(axis=0)
    with open(os.path.join(args.out, f"sbc_summary{suffix}.csv"), "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["parameter", "coverage_95", "coverage_90",
                         "mean_rank", "sd_rank"])
        for i in range(12):
            writer.writerow([f"theta{i + 1}", cov95[i], cov90[i],
                             ranks[:, i].mean(), ranks[:, i].std()])
    np.save(os.path.join(args.out, f"sbc_ranks{suffix}.npy"), ranks)

    losses = getattr(inference, "_summary", None)
    if losses is not None and hasattr(losses, "train_log_probs"):
        with open(os.path.join(args.out, f"loss_curve{suffix}.csv"), "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["epoch", "train_loss", "val_loss"])
            train = [-v for v in losses.train_log_probs]
            val = ([-v for v in losses.val_log_probs]
                   if losses.val_log_probs is not None else [np.nan] * len(train))
            writer.writerows(zip(range(len(train)), train, val))

    cfg["round"] = args.round
    cfg[f"n_round{args.round}"] = len(theta_r)
    cfg["n_draws"] = args.n_draws
    cfg["observation_noise"] = "Binomial age p; logNormal N sigma=0.10"
    with open(os.path.join(args.out, f"npe_config{suffix}.json"), "w") as fh:
        json.dump(cfg, fh, indent=2)
    print(f"[npe] round {args.round} (final) done; posterior -> "
          f"posterior_samples_npe{suffix}.csv")


if __name__ == "__main__":
    main()
