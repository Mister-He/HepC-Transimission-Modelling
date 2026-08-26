"""似然/先验数值与 legacy 标定输出的对照测试。"""

import numpy as np
import jax.numpy as jnp

from src.model import legacy_log_prior, _prior_bounds_factor
from src.simulator import build_params_from_theta, j_summary


def test_statistical_nll_matches_legacy(sim_fn, y0, base_params, targets, theta_nm):
    """在 NM 点估计处，二项 + log-normal NLL 与 legacy solutions.csv 一致。"""
    th = jnp.asarray(theta_nm)
    pm = build_params_from_theta(th, base_params)
    _, _, y_tgt = sim_fn(y0, pm)
    n_hat, p_sero, _ = j_summary(y_tgt)
    n_obs = np.asarray(targets["prison_total"])
    x_obs = np.asarray(targets["x_prev"])
    p_safe = np.clip(np.asarray(p_sero), targets["eps_prev"], 1 - targets["eps_prev"])
    from scipy.stats import binom

    nll_prev = -np.sum(binom.logpmf(x_obs, n_obs, p_safe))
    nll_pop = 0.5 * np.sum(((np.log(np.asarray(n_hat)) - np.log(n_obs)) / targets["sigma_pop"]) ** 2)
    # legacy solutions.csv（best start）：nll_prev=20.5579308229303, nll_pop=0.860023710776245
    assert np.isclose(nll_prev, 20.5579308229303, rtol=1e-6, atol=1e-5)
    assert np.isclose(nll_pop, 0.860023710776245, rtol=1e-6, atol=1e-5)


def test_prior_log_density(priors):
    """legacy 先验对数密度与手工计算一致。"""
    import numpy as np

    theta = jnp.asarray([0.0] * 6 + [0.0] * 6)
    lp = float(legacy_log_prior(theta, priors))
    contact_sd = priors["contact"]["sd"]
    beta_anchor = np.asarray(priors["beta"]["anchor_scale"])
    beta_sd = priors["beta"]["sd"]
    df = priors["beta"]["df"]
    from scipy.stats import norm, t

    # theta[0:6] = 0：x=0, loc=log(anchor)
    expected = np.sum(
        norm.logpdf(0.0, loc=np.log(np.asarray(priors["contact"]["anchor_scale"])),
                    scale=contact_sd)
    )
    z = (0.0 - np.log(beta_anchor)) / beta_sd
    expected += np.sum(t.logpdf(z, df=df))
    assert np.isclose(lp, expected, rtol=1e-6)


def test_bounds_factor():
    ok = float(_prior_bounds_factor(jnp.zeros(6), jnp.zeros(6)))
    assert ok == 0.0
    bad = float(_prior_bounds_factor(jnp.zeros(6), jnp.full(6, 20.0)))
    assert bad == -np.inf
