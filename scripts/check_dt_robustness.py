#!/usr/bin/env python3
# =============================================================================
# check_dt_robustness.py — 量化模拟步长（dt）对模型摘要的影响
#
# 主推断采用 dt=7 天以控制计算成本；本脚本用精确模型（dt=1/365）与
# dt=7 天模型在同一些后验样本上比较目标摘要（p_hat、N_hat）与平衡态指标，
# 量化离散化近似对 PPC/后验比较的影响。
#
# 用法：
#   python scripts/check_dt_robustness.py \
#       --posterior results/numpyro_nuts/posterior_samples_theta.csv \
#       --n-draws 200 --out results/numpyro_nuts/dt_robustness.csv
# =============================================================================
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import jax
import jax.numpy as jnp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.simulator import (  # noqa: E402
    build_params_from_theta,
    equilibrium_metrics,
    expand_y0,
    j_summary,
    load_model_parameters,
    load_simulation,
    make_sim_fn,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--posterior", type=str, required=True)
    ap.add_argument("--n-draws", type=int, default=200)
    ap.add_argument("--out", type=str, default="results/numpyro_nuts/dt_robustness.csv")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    base = load_model_parameters()
    sim_cfg = load_simulation()
    y0 = expand_y0(sim_cfg)

    post = pd.read_csv(args.posterior)
    cols = [f"theta{i}" for i in range(1, 13)]
    rng = np.random.default_rng(args.seed)
    th = post[cols].values
    idx = rng.choice(len(th), size=min(args.n_draws, len(th)), replace=False)
    th = th[idx]

    cfg1 = dict(sim_cfg)
    cfg7 = dict(sim_cfg)
    cfg7["dt"] = 7.0 / 365.0
    sim1 = make_sim_fn(cfg1)
    sim7 = make_sim_fn(cfg7)

    def make_batch(sim_fn):
        @jax.jit
        def summaries_batch(thetas):
            def one(t):
                pm = build_params_from_theta(t, base)
                yT, yT5, ytgt = sim_fn(y0, pm)
                n, p, _ = j_summary(ytgt)
                eq = equilibrium_metrics(yT, yT5)
                return n, p, eq["pass"].astype(jnp.float64), eq["max_prev_change"]

            return jax.vmap(one)(thetas)

        return summaries_batch

    batch1 = make_batch(sim1)
    batch7 = make_batch(sim7)
    n1, p1, pass1, prev1 = batch1(jnp.asarray(th))
    n7, p7, pass7, prev7 = batch7(jnp.asarray(th))

    dN = np.abs(np.asarray(n1) - np.asarray(n7))
    dp = np.abs(np.asarray(p1) - np.asarray(p7))
    rows = []
    for i, age in enumerate(sim_cfg.get("age_groups", [f"age{i+1}" for i in range(6)])):
        rows.append(
            {
                "age_group": age,
                "N_hat_1d_median": np.median(np.asarray(n1)[:, i]),
                "N_hat_7d_median": np.median(np.asarray(n7)[:, i]),
                "N_abs_diff_mean": dN[:, i].mean(),
                "N_abs_diff_max": dN[:, i].max(),
                "N_rel_diff_pct": 100 * dN[:, i].mean() / np.median(np.asarray(n1)[:, i]),
                "p_hat_1d_median": np.median(np.asarray(p1)[:, i]),
                "p_hat_7d_median": np.median(np.asarray(p7)[:, i]),
                "p_abs_diff_mean": dp[:, i].mean(),
                "p_abs_diff_max": dp[:, i].max(),
                "eq_pass_1d": float(np.mean(np.asarray(pass1) == 1.0)),
                "eq_pass_7d": float(np.mean(np.asarray(pass7) == 1.0)),
            }
        )
    out = pd.DataFrame(rows)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(out.round(5).to_string(index=False))
    print("\nSaved to", out_path)


if __name__ == "__main__":
    main()
