#!/usr/bin/env python3
# =============================================================================
# npe_train.py — trains single-round NPE (sbi NPE_C, standard NLL loss) on R-generated
# prior simulations and writes posterior draws for the HCV calibration parameters. The
# posterior is an equal-weight ensemble of --n-ensemble (default 5) networks that differ
# only in their torch seed.
#
# Network and training: MAF with 50 hidden features and 5 transforms, batch 256,
# Adam lr 5e-4, validation fraction 0.10, early stopping after 20 epochs, at most
# 250 epochs, gradient clipping 5.0, per-feature z-scoring by sbi, CPU only.
# Parameters are the 7 log-scale theta (six contact scales, first-arrest multiplier c) with an
# independent Normal(--prior-mean, --prior-sd) prior truncated to [--theta-lo, --theta-hi]
# (posterior draws outside are rejected). Data are x = logit(p_1..6) with p clipped to
# [1e-4, 1-1e-4]; Beta(p*phi, (1-p)*phi) noise (--beta-precision phi) is added here so
# training pairs follow the calibration observation model. Outputs in --out:
# posterior_samples_npe.csv, npe_config.json, npe_timing.json, loss_curve.csv,
# npe_inference.pkl and, when --n-sbc > 0, sbc_summary.csv / sbc_ranks.npy /
# sbc_rank_histograms.png.
# =============================================================================
import argparse
import csv
import json
import math
import os
import pickle
import time

T_IMPORT0 = time.time()
import numpy as np
import torch
import sbi
from sbi.inference import NPE_C
from sbi.neural_nets import posterior_nn
from torch.distributions import Independent, Normal
T_IMPORT = time.time() - T_IMPORT0
N_THETA = 7  # six log contact scales + log c


def read_numeric_csv(path):
    return np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)


def write_numeric_csv(path, values, columns):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(columns)
        writer.writerows(np.asarray(values).tolist())


def parse_vector(value):
    return [float(v) for v in value.split(",")]


def build_prior(mean, sd):
    if len(mean) != len(sd):
        raise ValueError("--prior-mean and --prior-sd must have the same length")
    return Independent(Normal(torch.tensor(mean, dtype=torch.float32),
                              torch.tensor(sd, dtype=torch.float32)), 1)


CLIP = 1e-4  # observed or simulated prevalence is clipped to [CLIP, 1 - CLIP] before logit


def add_observation_noise(x_raw, phi, rng):
    """Beta(p * phi, (1 - p) * phi) noise on each age group's prevalence (no counts)."""
    p = np.clip(np.asarray(x_raw, dtype=np.float64), CLIP, 1.0 - CLIP)
    return rng.beta(p * phi, (1.0 - p) * phi)


def transform_x(x):
    p = np.clip(x, CLIP, 1.0 - CLIP)
    return np.log(p / (1.0 - p))


class Ensemble:
    """Equal-weight mixture of independently trained posteriors (n/K draws from each), truncated
    to the parameter bounds by rejection (the prior is a Normal truncated to the bounds)."""
    def __init__(self, members, lo, hi):
        self.members, self.lo, self.hi = members, torch.tensor(lo), torch.tensor(hi)

    def _one(self, m, n, x, show_progress_bars):
        out = torch.empty((0, len(self.lo)))
        for _ in range(100):
            d = m.sample((2 * n,), x=x, show_progress_bars=show_progress_bars)
            out = torch.cat([out, d[((d >= self.lo) & (d <= self.hi)).all(dim=1)]])
            if len(out) >= n:
                return out[:n]
        raise RuntimeError(f"only {len(out)} of {200 * n} posterior draws fell inside the bounds")

    def sample(self, shape, x=None, show_progress_bars=False):
        k = len(self.members)
        if shape[0] % k:
            raise ValueError(f"draw count {shape[0]} must be a multiple of the ensemble size {k}")
        return torch.cat([self._one(m, shape[0] // k, x, show_progress_bars) for m in self.members])


def rank_shape(r, n, chisq_p):
    """Predeclared skew rule: outer-10% rank share or mean rank more than 3 binomial SE off,
    else any deviation with a Bonferroni-corrected (per parameter) chi-square p < 0.001.
    The two 3-SE checks are per parameter, so a calibrated posterior still fails the
    all-uniform gate about 4% of the time (7 parameters); a failure should be confirmed on fresh data."""
    outer = np.mean((r < 0.05) | (r > 0.95))
    se_outer, se_mean = math.sqrt(0.1 * 0.9 / n), math.sqrt(1 / 12 / n)
    if outer > 0.1 + 3 * se_outer:
        return "U-shaped"
    if outer < 0.1 - 3 * se_outer:
        return "modal"
    if abs(r.mean() - 0.5) > 3 * se_mean:
        return "biased_high" if r.mean() > 0.5 else "biased_low"
    return "nonuniform" if chisq_p < 0.001 / N_THETA else "uniform"


def run_sbc(posterior, args, rng):
    from scipy import stats
    theta = read_numeric_csv(args.theta_sbc)
    x = transform_x(add_observation_noise(read_numeric_csv(args.x_sbc), args.beta_precision, rng))
    n = min(args.n_sbc, len(theta))
    ranks = np.zeros((n, theta.shape[1]))
    for j in range(n):
        d = posterior.sample((args.sbc_draws,), x=torch.as_tensor(x[j], dtype=torch.float32),
                             show_progress_bars=False).numpy()
        ranks[j] = (d < theta[j]).mean(axis=0)
    cov95 = ((ranks >= 0.025) & (ranks <= 0.975)).mean(axis=0)
    cov90 = ((ranks >= 0.05) & (ranks <= 0.95)).mean(axis=0)
    n_bins = 20
    with open(os.path.join(args.out, "sbc_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["parameter", "n_sbc", "coverage_95", "coverage_95_se", "coverage_90", "coverage_90_se",
                    "mean_rank", "sd_rank", "outer10_share", "chisq20_p", "ks_p", "rank_shape"])
        for i in range(ranks.shape[1]):
            r = ranks[:, i]
            chisq_p = stats.chisquare(np.histogram(r, bins=n_bins, range=(0, 1))[0]).pvalue
            w.writerow([f"theta{i + 1}", n, cov95[i], math.sqrt(cov95[i] * (1 - cov95[i]) / n),
                        cov90[i], math.sqrt(cov90[i] * (1 - cov90[i]) / n), r.mean(), r.std(),
                        np.mean((r < 0.05) | (r > 0.95)), chisq_p,
                        stats.kstest(r, "uniform").pvalue, rank_shape(r, n, chisq_p)])
    np.save(os.path.join(args.out, "sbc_ranks.npy"), ranks)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 4, figsize=(10, 5), sharex=True)
    lo, hi = stats.binom.ppf([0.005, 0.995], n, 1 / n_bins)
    for ax in axes.flat[ranks.shape[1]:]:
        ax.set_visible(False)
    for i, ax in zip(range(ranks.shape[1]), axes.flat):
        ax.hist(ranks[:, i], bins=n_bins, range=(0, 1), color="0.6", edgecolor="white")
        ax.axhspan(lo, hi, color="tab:blue", alpha=0.15, lw=0)
        ax.set_title(f"theta{i + 1}", fontsize=9)
    fig.supxlabel("Normalised rank of the true parameter (band: 99% range under uniformity)")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "sbc_rank_histograms.png"), dpi=150)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--theta", required=True)
    ap.add_argument("--x", required=True)
    ap.add_argument("--x-obs", required=True)
    ap.add_argument("--theta-sbc", default=None)
    ap.add_argument("--x-sbc", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20261012)
    ap.add_argument("--n-draws", type=int, default=20000)
    ap.add_argument("--n-sbc", type=int, default=0)
    ap.add_argument("--n-ensemble", type=int, default=5)
    ap.add_argument("--sbc-draws", type=int, default=1000)
    ap.add_argument("--prior-mean", required=True)  # comma-separated, log scale
    ap.add_argument("--prior-sd", required=True)
    ap.add_argument("--beta-precision", type=float, required=True)  # phi of the Beta likelihood
    ap.add_argument("--theta-lo", required=True)  # parameter bounds: draws outside are rejected
    ap.add_argument("--theta-hi", required=True)
    ap.add_argument("--hidden-features", type=int, default=50)
    ap.add_argument("--num-transforms", type=int, default=5)
    ap.add_argument("--batch-size", type=int, default=256)
    ap.add_argument("--learning-rate", type=float, default=5e-4)
    ap.add_argument("--validation-fraction", type=float, default=0.10)
    ap.add_argument("--stop-after-epochs", type=int, default=20)
    ap.add_argument("--max-epochs", type=int, default=250)
    ap.add_argument("--clip-max-norm", type=float, default=5.0)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--sampling-label", default="unspecified")
    args = ap.parse_args()
    torch.set_num_threads(args.threads)
    timing = {"python_import_sec": T_IMPORT}
    t0 = time.time()

    # Separate reproducible streams: observation noise, network/training, posterior sampling.
    noise_rng = np.random.default_rng([args.seed, 11])
    torch_train_seed = args.seed * 100 + 21
    torch_sample_seed = args.seed * 100 + 31

    prior = build_prior(parse_vector(args.prior_mean), parse_vector(args.prior_sd))
    x_obs_raw = read_numeric_csv(args.x_obs)[0]
    x_obs = torch.as_tensor(transform_x(x_obs_raw[None, :])[0], dtype=torch.float32)

    theta = read_numeric_csv(args.theta)
    if theta.shape[1] != N_THETA:
        raise ValueError(f"expected {N_THETA} parameters, got {theta.shape[1]}")
    x = transform_x(add_observation_noise(read_numeric_csv(args.x), args.beta_precision, noise_rng))
    timing["load_sec"] = time.time() - t0
    t1 = time.time()
    # Single round on prior simulations (standard NLL loss). Each ensemble member differs only in
    # its network initialisation and batch order (torch seed); a single trained flow has small
    # seed-specific location errors that the equal-weight mixture averages out.
    members, summaries = [], []
    for k in range(args.n_ensemble):
        torch.manual_seed(torch_train_seed + k)
        inference = NPE_C(prior=prior, device=args.device,
                          density_estimator=posterior_nn("maf", hidden_features=args.hidden_features,
                                                         num_transforms=args.num_transforms,
                                                         z_score_theta="independent", z_score_x="independent"))
        inference.append_simulations(torch.as_tensor(theta, dtype=torch.float32),
                                     torch.as_tensor(x, dtype=torch.float32))
        estimator = inference.train(training_batch_size=args.batch_size, learning_rate=args.learning_rate,
                                    validation_fraction=args.validation_fraction,
                                    stop_after_epochs=args.stop_after_epochs, max_num_epochs=args.max_epochs,
                                    clip_max_norm=args.clip_max_norm, show_train_summary=False)
        members.append(inference.build_posterior(estimator).set_default_x(x_obs))
        summaries.append(inference._summary)
    timing["train_sec"] = time.time() - t1
    posterior = Ensemble(members, parse_vector(args.theta_lo), parse_vector(args.theta_hi))
    with open(os.path.join(args.out, "npe_inference.pkl"), "wb") as fh:
        pickle.dump(posterior, fh)

    t2 = time.time()
    torch.manual_seed(torch_sample_seed)
    draws = posterior.sample((args.n_draws,), x=x_obs, show_progress_bars=False).numpy()
    write_numeric_csv(os.path.join(args.out, "posterior_samples_npe.csv"), draws,
                      [f"theta{i}" for i in range(1, N_THETA + 1)])
    timing["sample_sec"] = time.time() - t2
    if args.n_sbc > 0 and args.theta_sbc:
        t3 = time.time()
        run_sbc(posterior, args, noise_rng)
        timing["sbc_sec"] = time.time() - t3
    timing["total_sec"] = time.time() - t0 + T_IMPORT

    with open(os.path.join(args.out, "loss_curve.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["member", "epoch", "train_loss", "val_loss", "epoch_seconds"])
        for k, summary in enumerate(summaries, 1):
            tl, vl, es = (summary.get("training_loss", []), summary.get("validation_loss", []),
                          summary.get("epoch_durations_sec", []))
            for i in range(len(tl)):
                w.writerow([k, i + 1, tl[i], vl[i] if i < len(vl) else "", es[i] if i < len(es) else ""])

    cfg = {
        "algorithm": "single-round NPE (sbi NPE_C, standard NLL loss)",
        "prior_sampling": args.sampling_label,
        "sbi_version": sbi.__version__, "torch_version": torch.__version__,
        "device": args.device, "torch_threads": torch.get_num_threads(),
        "density_estimator": "maf", "n_ensemble": args.n_ensemble, "hidden_features": args.hidden_features,
        "num_transforms": args.num_transforms, "z_score_theta": "independent", "z_score_x": "independent",
        "batch_size": args.batch_size, "learning_rate": args.learning_rate, "optimizer": "Adam",
        "validation_fraction": args.validation_fraction, "stop_after_epochs": args.stop_after_epochs,
        "max_epochs": args.max_epochs, "clip_max_norm": args.clip_max_norm,
        "x_transform": f"logit(p) with p clipped to [{CLIP}, 1-{CLIP}]",
        "observation_noise": f"Beta(p*phi, (1-p)*phi), phi={args.beta_precision}",
        "truncation": {"lo": parse_vector(args.theta_lo), "hi": parse_vector(args.theta_hi)},
        "prior": {"mean": parse_vector(args.prior_mean), "sd": parse_vector(args.prior_sd),
                  "form": "independent Normal on log scales"},
        "n_train": len(theta), "n_draws": args.n_draws, "seed": args.seed,
        "seeds": {"noise": [args.seed, 11], "torch_train": [torch_train_seed + k for k in range(args.n_ensemble)],
                  "torch_sample": torch_sample_seed},
        "epochs_trained": [sm.get("epochs_trained") for sm in summaries],
        "best_validation_loss": [sm.get("best_validation_loss") for sm in summaries],
    }
    with open(os.path.join(args.out, "npe_config.json"), "w") as fh:
        json.dump(cfg, fh, indent=2)
    with open(os.path.join(args.out, "npe_timing.json"), "w") as fh:
        json.dump(timing, fh, indent=2)
    print(f"[npe] trained on {len(theta)} simulations in {timing['train_sec']:.1f}s, "
          f"{args.n_ensemble} members, epochs {[sm.get('epochs_trained') for sm in summaries]}")


if __name__ == "__main__":
    main()
