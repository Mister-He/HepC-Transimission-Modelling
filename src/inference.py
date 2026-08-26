# =============================================================================
# inference.py — 推断结果整理：收敛诊断、后验预测检验、与 legacy 后验比较
# =============================================================================
from __future__ import annotations

from typing import Any, Dict, List

import jax.numpy as jnp
import numpy as np
import pandas as pd
import scipy.stats as st


THETA_COLS = [f"theta{i+1}" for i in range(12)]


def theta_df_from_samples(samples: Dict[str, Any], priors: Dict[str, Any]) -> pd.DataFrame:
    """把 numpyro 后验样本整理为 legacy 同构的 theta1..12 长表（chain, iteration）。"""
    if "theta" in samples:
        th = np.asarray(samples["theta"])          # (n_chains, n_draws, 12)
        n_chains, n_draws, _ = th.shape
        frames = []
        for c in range(n_chains):
            d = {
                "chain": c,
                "iteration": np.arange(1, n_draws + 1),
                **{f"theta{i+1}": th[c, :, i] for i in range(12)},
            }
            frames.append(pd.DataFrame(d))
        return pd.concat(frames, ignore_index=True)
    n_chains = samples["contact_log"].shape[0]
    n_draws = samples["contact_log"].shape[1]
    beta_anchor = np.log(np.asarray(priors["beta"]["anchor_scale"], dtype=float))
    beta_sd = priors["beta"]["sd"]
    frames = []
    for c in range(n_chains):
        contact = np.asarray(samples["contact_log"][c])          # (n,6)
        z_beta = np.asarray(samples["z_beta"][c])                # (n,6)
        beta_log = beta_anchor[None, :] + beta_sd * z_beta
        d = {
            "chain": c,
            "iteration": np.arange(1, n_draws + 1),
            **{f"theta{i+1}": contact[:, i] for i in range(6)},
            **{f"theta{i+7}": beta_log[:, i] for i in range(6)},
        }
        frames.append(pd.DataFrame(d))
    return pd.concat(frames, ignore_index=True)


def summarize_posterior(theta_df: pd.DataFrame, cols: List[str] = THETA_COLS) -> pd.DataFrame:
    s = pd.DataFrame(
        {
            "mean": theta_df[cols].mean(),
            "median": theta_df[cols].median(),
            "sd": theta_df[cols].std(),
            "q025": theta_df[cols].quantile(0.025),
            "q975": theta_df[cols].quantile(0.975),
        }
    )
    s.index.name = "parameter"
    return s


def compare_with_legacy(
    new_df: pd.DataFrame,
    legacy_df: pd.DataFrame,
    cols: List[str] = THETA_COLS,
) -> pd.DataFrame:
    """逐参数比较新后验（NUTS）与 legacy 后验（AM-MCMC）：
    中位数/95% 区间、区间重叠、2 样本 KS 检验（作用于全部样本）。"""
    rows = []
    for p in cols:
        new = new_df[p]
        leg = legacy_df[p]
        new_med, new_q025, new_q975 = new.median(), new.quantile(0.025), new.quantile(0.975)
        leg_med, leg_q025, leg_q975 = leg.median(), leg.quantile(0.025), leg.quantile(0.975)
        overlap = max(0.0, min(new_q975, leg_q975) - max(new_q025, leg_q025))
        ks = st.ks_2samp(new.values, leg.values)
        rows.append(
            {
                "parameter": p,
                "numpyro_median": new_med,
                "numpyro_q025": new_q025,
                "numpyro_q975": new_q975,
                "legacy_median": leg_med,
                "legacy_q025": leg_q025,
                "legacy_q975": leg_q975,
                "interval_overlap": overlap,
                "median_shift": new_med - leg_med,
                "ks_statistic": ks.statistic,
                "ks_pvalue": ks.pvalue,
            }
        )
    return pd.DataFrame(rows)


def ppc_target_summary(ppc: Dict[str, Any], targets: Dict[str, Any]) -> pd.DataFrame:
    """目标级后验预测摘要（p_hat、N_hat）+ 观测覆盖 + 平衡态通过率。"""
    p_hat = np.asarray(ppc["p_hat"])          # (n_ppc, 6)
    n_hat = np.asarray(ppc["N_hat"])
    eq_pass = np.asarray(ppc["eq_pass"]).ravel() == 1.0
    p_hat = p_hat[eq_pass]
    n_hat = n_hat[eq_pass]
    p_obs = np.asarray(targets["x_prev"]) / np.asarray(targets["prison_total"])
    n_obs = np.asarray(targets["prison_total"])
    rows = []
    for i, age in enumerate(targets["age_groups"]):
        p_lo, p_hi = np.quantile(p_hat[:, i], [0.025, 0.975])
        n_lo, n_hi = np.quantile(n_hat[:, i], [0.025, 0.975])
        rows.append(
            {
                "age_group": age,
                "p_obs": p_obs[i],
                "p_ppc_median": np.median(p_hat[:, i]),
                "p_ppc_lo": p_lo,
                "p_ppc_hi": p_hi,
                "p_obs_in_cri": bool(p_lo <= p_obs[i] <= p_hi),
                "N_obs": n_obs[i],
                "N_ppc_median": np.median(n_hat[:, i]),
                "N_ppc_lo": n_lo,
                "N_ppc_hi": n_hi,
                "N_obs_in_cri": bool(n_lo <= n_obs[i] <= n_hi),
            }
        )
    out = pd.DataFrame(rows)
    out["eq_pass_fraction"] = float(np.mean(eq_pass))
    out["n_draws_used"] = len(p_hat)
    return out


def ppc_data_summary(ppc: Dict[str, Any], targets: Dict[str, Any]) -> pd.DataFrame:
    """数据级后验预测：二项患病率与监狱人口观测覆盖 + Bayesian p 值。"""
    prev_obs = np.asarray(ppc["prev_obs"])    # (n_ppc, 6)
    pop_obs = np.asarray(ppc["pop_obs"])
    eq_pass = np.asarray(ppc["eq_pass"]).ravel() == 1.0
    prev_obs = prev_obs[eq_pass]
    pop_obs = pop_obs[eq_pass]
    n_obs = np.asarray(targets["prison_total"])
    p_obs = np.asarray(targets["x_prev"]) / n_obs
    rows = []
    for i, age in enumerate(targets["age_groups"]):
        p_sim = prev_obs[:, i] / n_obs[i]
        lo, hi = np.quantile(p_sim, [0.025, 0.975])
        pop_lo, pop_hi = np.quantile(pop_obs[:, i], [0.025, 0.975])
        rows.append(
            {
                "age_group": age,
                "prev_obs": p_obs[i],
                "prev_ppc_median": np.median(p_sim),
                "prev_ppc_lo": lo,
                "prev_ppc_hi": hi,
                "prev_obs_in_cri": bool(lo <= p_obs[i] <= hi),
                "bayes_p_prev": float(np.mean(p_sim <= p_obs[i])),
                "pop_obs": n_obs[i],
                "pop_ppc_median": np.median(pop_obs[:, i]),
                "pop_ppc_lo": pop_lo,
                "pop_ppc_hi": pop_hi,
                "pop_obs_in_cri": bool(pop_lo <= n_obs[i] <= pop_hi),
                "bayes_p_pop": float(np.mean(pop_obs[:, i] <= n_obs[i])),
            }
        )
    out = pd.DataFrame(rows)
    out["n_draws_used"] = len(prev_obs)
    return out
