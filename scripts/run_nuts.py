#!/usr/bin/env python3
# =============================================================================
# run_nuts.py — NumPyro NUTS 推断流水线
#
# 流程：
#   1. 读取 config/（模型参数、模拟设置、目标、先验）
#   2. NUTS 多链采样（以 legacy Nelder-Mead 点估计为初值）
#   3. 收敛诊断（R-hat / ESS / 发散数 / 接受率）
#   4. 后验预测检验（目标级 + 数据级）
#   5. 与 legacy 后验（legacy/posterior_samples_mcmc.csv）逐参数比较
#   6. 输出结果与图到 results/numpyro_nuts/
#
# 用法示例：
#   python scripts/run_nuts.py --num-warmup 300 --num-samples 500 \
#       --num-chains 4 --seed 2026 --out-dir results/numpyro_nuts
# =============================================================================
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcache")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import jax
import jax.numpy as jnp
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import arviz as az
import numpyro
from numpyro.infer import Predictive

from src.inference import (
    THETA_COLS,
    compare_with_legacy,
    ppc_data_summary,
    ppc_target_summary,
    summarize_posterior,
    theta_df_from_samples,
)
from src.model import make_mcmc
from src.simulator import (
    build_params_from_theta,
    expand_y0,
    load_model_parameters,
    load_priors,
    load_simulation,
    load_targets,
    make_sim_fn,
)


def parse_args():
    ap = argparse.ArgumentParser(description="NumPyro NUTS 推断")
    ap.add_argument("--out-dir", type=str, default=str(ROOT / "results" / "numpyro_nuts"))
    ap.add_argument("--num-warmup", type=int, default=300)
    ap.add_argument("--num-samples", type=int, default=500)
    ap.add_argument("--num-chains", type=int, default=4)
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--dt-days", type=float, default=1.0,
                    help="模拟步长（天）。1 = 与 legacy 完全一致；>1 为加速实验。")
    ap.add_argument("--chain-method", type=str, default="sequential",
                    choices=["sequential", "parallel", "vectorized"])
    ap.add_argument("--ppc-draws", type=int, default=400)
    ap.add_argument("--target-accept", type=float, default=0.8)
    ap.add_argument("--penalty-mode", type=str, default="none",
                    choices=["legacy", "smooth", "none"],
                    help="平衡态惩罚模式：legacy=1e6 硬惩罚（同 legacy）；smooth=平滑二次惩罚；none=无（后处理过滤，默认）")
    ap.add_argument("--parameterization", type=str, default="whiten",
                    choices=["theta", "whiten"],
                    help="采样参数化：theta=原始 12 参数空间；whiten=legacy 后验协方差白化空间（默认）")
    ap.add_argument("--steps", type=str, default="all",
                    help="all | fit | ppc | compare")
    return ap.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sim_cfg = load_simulation()
    if args.dt_days != 1.0:
        sim_cfg = dict(sim_cfg)
        sim_cfg["dt"] = args.dt_days / 365.0
    base_params = load_model_parameters()
    targets = load_targets()
    priors = load_priors()
    y0 = expand_y0(sim_cfg)
    sim_fn = make_sim_fn(sim_cfg)

    legacy_pt = json.loads((ROOT / "legacy" / "nelder_mead_point_estimate.json").read_text())
    theta_nm = jnp.asarray(legacy_pt["theta_log"], dtype=jnp.float64)
    legacy_df = pd.read_csv(ROOT / "legacy" / "posterior_samples_mcmc.csv")
    legacy_theta = legacy_df[THETA_COLS].values

    whiten_mean = None
    whiten_chol = None
    if args.parameterization == "whiten":
        cov = np.cov(legacy_theta, rowvar=False)
        whiten_mean = jnp.asarray(np.mean(legacy_theta, axis=0), dtype=jnp.float64)
        whiten_chol = jnp.asarray(
            np.linalg.cholesky(cov + 1e-8 * np.eye(12)), dtype=jnp.float64
        )

    run_meta = {
        "dt_days": args.dt_days,
        "num_warmup": args.num_warmup,
        "num_samples": args.num_samples,
        "num_chains": args.num_chains,
        "seed": args.seed,
        "chain_method": args.chain_method,
        "penalty_mode": args.penalty_mode,
        "parameterization": args.parameterization,
        "whiten_source": "legacy posterior mean + chol(cov + 1e-8 I)",
        "numpyro_version": numpyro.__version__,
        "jax_version": jax.__version__,
        "init": "legacy Nelder-Mead point estimate",
        "x64": str(jax.config.jax_enable_x64),
    }

    if args.steps in ("all", "fit"):
        print("== NUTS fitting ==", run_meta)
        mcmc = make_mcmc(
            sim_fn,
            y0,
            base_params,
            targets,
            priors,
            theta_nm,
            num_warmup=args.num_warmup,
            num_samples=args.num_samples,
            num_chains=args.num_chains,
            seed=args.seed,
            chain_method=args.chain_method,
            target_accept_prob=args.target_accept,
            penalty_mode=args.penalty_mode,
            parameterization=args.parameterization,
            whiten_mean=whiten_mean,
            whiten_chol=whiten_chol,
        )
        samples = mcmc.get_samples(group_by_chain=True)
        extra = mcmc.get_extra_fields(group_by_chain=True)

        # ---- theta 空间长表 + 摘要 ----
        theta_df = theta_df_from_samples(samples, priors)
        theta_df["logpost"] = -np.asarray(extra["potential_energy"]).reshape(-1)
        theta_df.to_csv(out_dir / "posterior_samples_theta.csv", index=False)
        summary = summarize_posterior(theta_df)
        summary.to_csv(out_dir / "posterior_summary.csv")
        with open(out_dir / "run_config.json", "w") as fh:
            json.dump(run_meta, fh, indent=2, ensure_ascii=False)

        # ---- 收敛诊断（arviz）----
        idata = az.from_numpyro(mcmc)
        summary = az.summary(idata, var_names=["theta"], kind="diagnostics")
        # 列名兼容不同 arviz 版本
        rhat_col = next(c for c in summary.columns if "r_hat" in c or "rhat" in c)
        ess_bulk_col = next(c for c in summary.columns if "ess_bulk" in c)
        ess_tail_col = next(c for c in summary.columns if "ess_tail" in c)
        rhat_flat = summary[rhat_col].values
        ess_flat = summary[ess_bulk_col].values
        ess_tail_flat = summary[ess_tail_col].values
        div = int(np.asarray(extra["diverging"]).sum())
        accept = float(np.asarray(extra["accept_prob"]).mean())
        num_steps = np.asarray(extra["num_steps"])
        tree_depth = float(np.mean(np.log2(np.maximum(num_steps, 1)) + 1.0))

        diag_rows = []
        for i, p in enumerate(THETA_COLS):
            diag_rows.append(
                {
                    "parameter": p,
                    "rhat": float(rhat_flat[i]),
                    "ess_bulk": float(ess_flat[i]),
                    "ess_tail": float(ess_tail_flat[i]),
                    "n_divergences": int(div),
                    "mean_accept_prob": float(accept),
                    "mean_tree_depth": float(tree_depth),
                }
            )
        diag = pd.DataFrame(diag_rows)
        diag.to_csv(out_dir / "diagnostics_numpyro.csv", index=False)
        print(diag.to_string(index=False))
        print(f"\nR-hat range: [{diag.rhat.min():.4f}, {diag.rhat.max():.4f}]")
        print(f"ESS bulk range: [{diag.ess_bulk.min():.1f}, {diag.ess_bulk.max():.1f}]")
        print(f"ESS tail range: [{diag.ess_tail.min():.1f}, {diag.ess_tail.max():.1f}]")
        print(f"total divergences: {div}")

        # ---- 图：trace + 边际密度（NUTS vs legacy）----
        legacy_df = pd.read_csv(ROOT / "legacy" / "posterior_samples_mcmc.csv")
        fig, axes = plt.subplots(4, 3, figsize=(16, 18))
        for i, p in enumerate(THETA_COLS):
            ax = axes.flat[i]
            for c in range(args.num_chains):
                ax.plot(theta_df.loc[theta_df.chain == c, p].values,
                        lw=0.4, alpha=0.7, label=f"chain {c}" if i == 0 else None)
            ax.set_title(p)
        axes.flat[0].legend(fontsize=8)
        fig.suptitle("NUTS trace (theta space)")
        fig.tight_layout()
        fig.savefig(out_dir / "trace_numpyro.png", dpi=110)
        plt.close(fig)

        fig, axes = plt.subplots(4, 3, figsize=(16, 18))
        for i, p in enumerate(THETA_COLS):
            ax = axes.flat[i]
            ax.hist(legacy_df[p], bins=60, density=True, alpha=0.5,
                    label="legacy AM-MCMC", color="tab:orange")
            ax.hist(theta_df[p], bins=60, density=True, alpha=0.5,
                    label="NumPyro NUTS", color="tab:blue")
            ax.set_title(p)
        axes.flat[0].legend(fontsize=8)
        fig.suptitle("Marginal posterior: legacy vs NumPyro")
        fig.tight_layout()
        fig.savefig(out_dir / "density_compare.png", dpi=110)
        plt.close(fig)

        # ---- 与 legacy 比较 ----
        cmp = compare_with_legacy(theta_df, legacy_df)
        cmp.to_csv(out_dir / "compare_legacy.csv", index=False)
        print("\n== Comparison with legacy posterior ==")
        print(cmp.round(4).to_string(index=False))

    if args.steps in ("all", "ppc"):
        post = pd.read_csv(out_dir / "posterior_samples_theta.csv")
        n_draws = post.groupby("chain").size().iloc[0]
        step = max(1, n_draws // args.ppc_draws)
        thin = post.iloc[::step]
        th_flat = np.concatenate(
            [thin.loc[thin.chain == c, THETA_COLS].values
             for c in sorted(thin.chain.unique())],
            axis=0,
        )
        if args.parameterization == "whiten":
            u_flat = (th_flat - np.asarray(whiten_mean)) @ np.linalg.inv(
                np.asarray(whiten_chol)
            )
            flat = {"u": u_flat}
        else:
            priors_cfg = load_priors()
            beta_anchor = np.log(np.asarray(priors_cfg["beta"]["anchor_scale"], dtype=float))
            beta_sd = priors_cfg["beta"]["sd"]
            flat = {
                "contact_log": th_flat[:, 0:6],
                "z_beta": (th_flat[:, 6:12] - beta_anchor[None, :]) / beta_sd,
            }
        predictive = Predictive(
            _ppc_model,
            posterior_samples=flat,
            return_sites=["p_hat", "N_hat", "eq_pass", "prev_obs", "pop_obs"],
        )
        ppc = predictive(
            jax.random.PRNGKey(args.seed + 1),
            sim_fn=sim_fn,
            y0=y0,
            base_params=base_params,
            targets=targets,
            priors=priors,
            parameterization=args.parameterization,
            whiten_mean=whiten_mean,
            whiten_chol=whiten_chol,
        )
        eq_frac = float(np.mean(np.asarray(ppc["eq_pass"]) == 1.0))
        print(f"\nPPC equilibrium feasibility fraction: {eq_frac:.4f}")
        with open(out_dir / "run_config.json") as fh:
            cfg = json.load(fh)
        cfg["ppc_eq_pass_fraction"] = eq_frac
        with open(out_dir / "run_config.json", "w") as fh:
            json.dump(cfg, fh, indent=2, ensure_ascii=False)
        ppc_target_summary(ppc, targets).to_csv(out_dir / "ppc_target_summary.csv", index=False)
        ppc_data_summary(ppc, targets).to_csv(out_dir / "ppc_data_summary.csv", index=False)
        print("\n== PPC target summary ==")
        print(pd.read_csv(out_dir / "ppc_target_summary.csv").round(4).to_string(index=False))
        print("\n== PPC data summary ==")
        print(pd.read_csv(out_dir / "ppc_data_summary.csv").round(4).to_string(index=False))

        _plot_ppc(ppc, targets, out_dir)

    print("\nDone. Results in", out_dir)


def _ppc_model(
    sim_fn,
    y0,
    base_params,
    targets,
    priors,
    parameterization="whiten",
    whiten_mean=None,
    whiten_chol=None,
):
    """PPC 专用模型：固定后验样本，重新采样观测与摘要（不带 obs）。"""
    from src.model import hcv_model

    return hcv_model(
        sim_fn, y0, base_params, targets, priors,
        parameterization=parameterization,
        whiten_mean=whiten_mean,
        whiten_chol=whiten_chol,
        obs_prev=None, obs_pop=None,
    )


def _plot_ppc(ppc, targets, out_dir):
    p_hat = np.asarray(ppc["p_hat"])
    n_hat = np.asarray(ppc["N_hat"])
    p_obs = np.asarray(targets["x_prev"]) / np.asarray(targets["prison_total"])
    n_obs = np.asarray(targets["prison_total"])
    ages = targets["age_groups"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    x = np.arange(6)
    ax = axes[0]
    lo, hi = np.quantile(p_hat, [0.025, 0.975], axis=0)
    ax.errorbar(x, np.median(p_hat, axis=0), yerr=[np.median(p_hat, axis=0) - lo,
                                                   hi - np.median(p_hat, axis=0)],
                fmt="o", capsize=4, label="PPC 95% CrI")
    ax.plot(x, p_obs, "rs", label="observed")
    ax.set_xticks(x, ages)
    ax.set_title("Prevalence (target-level PPC)")
    ax.legend()

    ax = axes[1]
    lo, hi = np.quantile(n_hat, [0.025, 0.975], axis=0)
    ax.errorbar(x, np.median(n_hat, axis=0), yerr=[np.median(n_hat, axis=0) - lo,
                                                   hi - np.median(n_hat, axis=0)],
                fmt="o", capsize=4, label="PPC 95% CrI")
    ax.plot(x, n_obs, "rs", label="observed")
    ax.set_xticks(x, ages)
    ax.set_title("Prison population (target-level PPC)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "ppc_targets.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
