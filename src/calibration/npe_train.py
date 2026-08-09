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
    z = theta / 2                 (prior N(0,1)^6 x StudentT(3)^6)
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
import json
import os
import pickle

import numpy as np
import pandas as pd
import torch

from sbi.inference import NPE_C
from sbi.utils import MultipleIndependent
from torch.distributions import Normal, StudentT


def build_prior():
    return MultipleIndependent(
        [Normal(torch.zeros(1), torch.ones(1))] * 6 +
        [StudentT(df=torch.tensor([3.0]), loc=torch.zeros(1),
                  scale=torch.ones(1))] * 6
    )


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
    args = ap.parse_args()

    suffix = ("_" + args.suffix) if args.suffix else ""
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    prior = build_prior()
    device = "mps" if torch.backends.mps.is_available() else "cpu"
    print(f"[npe] round={args.round}, seed={args.seed}, device={device}")

    x_obs_raw = pd.read_csv(args.x_obs).values.astype(np.float64)[0]

    if args.round == 1:
        theta_raw = pd.read_csv(args.theta).values.astype(np.float64)
        x_raw = pd.read_csv(args.x).values.astype(np.float64)
        z_train = theta_raw / 2.0
        trans_x, mu_x, sd_x = make_x_transform(x_raw)
        x_z = trans_x(x_raw)

        inference = NPE_C(prior=prior, density_estimator="maf", device=device)
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
        prop = posterior.sample((args.n_proposal,)).numpy()
        pd.DataFrame(prop * 2.0,
                     columns=[f"theta{i}" for i in range(1, 13)]).to_csv(
            os.path.join(args.out, f"proposal_theta_round2{suffix}.csv"),
            index=False)

        config = {
            "seed": args.seed, "round": 1, "n_train": len(z_train),
            "n_proposal": args.n_proposal, "device": device,
            "density_estimator": "maf", "transform": "theta/2; logit(p),log(N)",
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

    theta_r = pd.read_csv(args.theta2).values.astype(np.float64) / 2.0
    x_r_raw = pd.read_csv(args.x2).values.astype(np.float64)
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
        pd.DataFrame(prop * 2.0,
                     columns=[f"theta{i}" for i in range(1, 13)]).to_csv(
            os.path.join(args.out,
                         f"proposal_theta_round{args.round + 1}{suffix}.csv"),
            index=False)
        cfg["round"] = args.round
        with open(os.path.join(args.out, f"npe_config{suffix}.json"), "w") as fh:
            json.dump(cfg, fh, indent=2)
        print(f"[npe] round {args.round} done; proposal -> "
              f"proposal_theta_round{args.round + 1}{suffix}.csv")
        return

    z_draws = posterior.sample((args.n_draws,)).numpy()
    pd.DataFrame(z_draws * 2.0,
                 columns=[f"theta{i}" for i in range(1, 13)]).to_csv(
        os.path.join(args.out, f"posterior_samples_npe{suffix}.csv"), index=False)

    theta_sbc_raw = pd.read_csv(args.theta_sbc).values.astype(np.float64)
    x_sbc_raw = pd.read_csv(args.x_sbc).values.astype(np.float64)
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
    pd.DataFrame({
        "parameter": [f"theta{i}" for i in range(1, 13)],
        "coverage_95": cov95, "coverage_90": cov90,
        "mean_rank": ranks.mean(axis=0), "sd_rank": ranks.std(axis=0),
    }).to_csv(os.path.join(args.out, f"sbc_summary{suffix}.csv"), index=False)
    np.save(os.path.join(args.out, f"sbc_ranks{suffix}.npy"), ranks)

    losses = getattr(inference, "_summary", None)
    if losses is not None and hasattr(losses, "train_log_probs"):
        pd.DataFrame({
            "epoch": np.arange(len(losses.train_log_probs)),
            "train_loss": [-v for v in losses.train_log_probs],
            "val_loss": ([-v for v in losses.val_log_probs]
                         if losses.val_log_probs is not None else np.nan),
        }).to_csv(os.path.join(args.out, f"loss_curve{suffix}.csv"), index=False)

    cfg["round"] = args.round
    cfg[f"n_round{args.round}"] = len(theta_r)
    cfg["n_draws"] = args.n_draws
    with open(os.path.join(args.out, f"npe_config{suffix}.json"), "w") as fh:
        json.dump(cfg, fh, indent=2)
    print(f"[npe] round {args.round} (final) done; posterior -> "
          f"posterior_samples_npe{suffix}.csv")


if __name__ == "__main__":
    main()
