"""推断工具函数测试。"""

import numpy as np
import pandas as pd

from src.inference import (
    compare_with_legacy,
    ppc_data_summary,
    ppc_target_summary,
    summarize_posterior,
    theta_df_from_samples,
)


def test_theta_df_from_samples_whiten(priors):
    samples = {
        "theta": np.random.default_rng(0).normal(size=(2, 10, 12)),
    }
    df = theta_df_from_samples(samples, priors)
    assert list(df.columns[:2]) == ["chain", "iteration"]
    assert len(df) == 20
    assert np.allclose(df.loc[df.chain == 0, "theta1"].values, samples["theta"][0, :, 0])


def test_summarize_and_compare(legacy_posterior):
    df = legacy_posterior.copy()
    cols = [f"theta{i}" for i in range(1, 13)]
    summary = summarize_posterior(df, cols)
    assert summary.loc["theta1", "median"] == df["theta1"].median()
    cmp = compare_with_legacy(df, df, cols)
    # 与自身比较：KS p 值应接近 1，重叠为整个区间
    assert (cmp["ks_pvalue"] > 0.05).all()
    assert (cmp["interval_overlap"] > 0).all()


def test_ppc_summaries(targets):
    rng = np.random.default_rng(1)
    n = 100
    ppc = {
        "p_hat": rng.uniform(0.1, 0.5, size=(n, 6)),
        "N_hat": rng.uniform(80, 1900, size=(n, 6)),
        "eq_pass": np.ones((n, 1)),
        "prev_obs": rng.binomial(
            np.asarray(targets["prison_total"]).astype(int), 0.3, size=(n, 6)
        ),
        "pop_obs": rng.normal(
            np.asarray(targets["prison_total"]), 100, size=(n, 6)
        ),
    }
    ts = ppc_target_summary(ppc, targets)
    ds = ppc_data_summary(ppc, targets)
    assert len(ts) == 6 and len(ds) == 6
    assert (ts["eq_pass_fraction"] == 1.0).all()
