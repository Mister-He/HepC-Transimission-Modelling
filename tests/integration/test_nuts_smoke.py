"""集成测试：短时域 NumPyro NUTS + 后验预测冒烟。

为保证 CI 可行，这里使用缩短的模拟时域（t_end=5 年、dt=7 天）与极少量样本，
只验证整条推断流水线可运行、输出有限且诊断可计算。完整推断见
scripts/run_nuts.py 与 results/numpyro_nuts/。
"""

import numpy as np
import pytest

import jax
import jax.numpy as jnp

from src.model import make_mcmc
from src.simulator import expand_y0, make_sim_fn

pytestmark = pytest.mark.integration


def _short_cfg(sim_cfg):
    cfg = dict(sim_cfg)
    cfg["t_end"] = 5.0
    cfg["target_time"] = 5.0
    cfg["equilibrium_t_lag"] = 1.0
    cfg["dt"] = 7.0 / 365.0
    return cfg


def test_nuts_runs_end_to_end(base_params, sim_cfg, targets, priors, theta_nm):
    cfg = _short_cfg(sim_cfg)
    y0 = expand_y0(cfg)
    sim_fn = make_sim_fn(cfg)
    th = jnp.asarray(theta_nm)

    mcmc = make_mcmc(
        sim_fn,
        y0,
        base_params,
        targets,
        priors,
        th,
        num_warmup=5,
        num_samples=5,
        num_chains=1,
        seed=7,
        chain_method="sequential",
        penalty_mode="none",
        parameterization="theta",
    )
    samples = mcmc.get_samples()
    assert samples["contact_log"].shape == (5, 6)
    assert np.all(np.isfinite(np.asarray(samples["contact_log"])))
    assert np.all(np.isfinite(np.asarray(samples["z_beta"])))

    import arviz as az

    idata = az.from_numpyro(mcmc)
    assert idata.posterior["contact_log"].shape[1] == 5


def test_ppc_smoke(base_params, sim_cfg, targets, priors, theta_nm):
    cfg = _short_cfg(sim_cfg)
    y0 = expand_y0(cfg)
    sim_fn = make_sim_fn(cfg)
    th = jnp.asarray(theta_nm)

    # 用 NM 点估计附近的少量伪后验样本跑 Predictive
    from numpyro.infer import Predictive
    from src.model import hcv_model

    rng = np.random.default_rng(0)
    contact_log = np.asarray(th)[0:6][None, :] + 0.01 * rng.normal(size=(4, 6))
    beta_log = np.asarray(th)[6:12][None, :] + 0.01 * rng.normal(size=(4, 6))
    priors_cfg = priors
    beta_anchor = np.log(np.asarray(priors_cfg["beta"]["anchor_scale"]))
    z_beta = (beta_log - beta_anchor[None, :]) / priors_cfg["beta"]["sd"]

    predictive = Predictive(
        hcv_model,
        posterior_samples={"contact_log": contact_log, "z_beta": z_beta},
        return_sites=["p_hat", "N_hat", "eq_pass", "prev_obs", "pop_obs"],
    )
    ppc = predictive(
        jax.random.PRNGKey(0),
        sim_fn=sim_fn,
        y0=y0,
        base_params=base_params,
        targets=targets,
        priors=priors,
        parameterization="theta",
        obs_prev=None,
        obs_pop=None,
    )
    assert np.asarray(ppc["p_hat"]).shape == (4, 6)
    assert np.all(np.isfinite(np.asarray(ppc["p_hat"])))
    assert np.all(np.isfinite(np.asarray(ppc["prev_obs"])))
