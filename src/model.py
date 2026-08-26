# =============================================================================
# model.py — NumPyro 概率模型
#
# 与 legacy（R 版）完全同构：
#   先验      contact_log ~ Normal(log(anchor), 0.5^2)（6 个）
#             z_beta ~ StudentT(3)；beta_log = log(anchor) + 0.8*z_beta
#   似然      p_hat[i] ~ Binomial(x_i | n_i)（血清学患病率）
#             log N_hat[i] ~ Normal(log N_obs[i], sigma_pop^2)
#   惩罚      equilibrium gate（与 legacy mcmc.R::log_posterior 一致，
#             仅 log_pop/0.01 + prev/0.005 + total/0.01 三项，系数 1e6）
#   边界      contact_scale ∈ [0.01, 1000]；beta_scale ∈ [0.01, 1e6]
# =============================================================================
from __future__ import annotations

from typing import Any, Dict, Optional

import jax
import jax.numpy as jnp
import numpyro
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS, init_to_value

from .simulator import (
    build_params_from_theta,
    equilibrium_metrics,
    j_summary,
)


def _prior_bounds_factor(contact_log: jax.Array, beta_log: jax.Array) -> jax.Array:
    contact_scale = jnp.exp(contact_log)
    beta_scale = jnp.exp(beta_log)
    in_bounds = jnp.all((contact_scale >= 0.01) & (contact_scale <= 1000.0)) & jnp.all(
        (beta_scale >= 0.01) & (beta_scale <= 1e6)
    )
    return jnp.where(in_bounds, 0.0, -jnp.inf)


def legacy_log_prior(theta: jax.Array, priors: Dict[str, Any]) -> jax.Array:
    """legacy 先验在 theta 空间的对数密度（mcmc.R::log_prior）。
    contact_log ~ Normal(log(anchor), sd^2)；z_beta ~ StudentT(df)，
    其中 z_beta = (beta_log - log(anchor)) / sd。"""
    contact_anchor = jnp.asarray(priors["contact"]["anchor_scale"], dtype=jnp.float64)
    beta_anchor = jnp.asarray(priors["beta"]["anchor_scale"], dtype=jnp.float64)
    z_beta = (theta[6:12] - jnp.log(beta_anchor)) / priors["beta"]["sd"]
    lp = dist.Normal(
        jnp.log(contact_anchor), priors["contact"]["sd"]
    ).log_prob(theta[0:6]).sum()
    lp = lp + dist.StudentT(priors["beta"]["df"]).log_prob(z_beta).sum()
    return lp


def _model_core(
    theta: jax.Array,
    sim_fn,
    y0: jax.Array,
    base_params: Dict[str, Any],
    targets: Dict[str, Any],
    eq_penalty_factor: float,
    penalty_mode: str,
    obs_prev: Optional[jax.Array],
    obs_pop: Optional[jax.Array],
):
    """共享的确定性部分：模拟 -> 摘要 -> 似然 -> 平衡态惩罚。"""
    numpyro.factor("bounds", _prior_bounds_factor(theta[0:6], theta[6:12]))

    pm = build_params_from_theta(theta, base_params)
    y_T, y_T5, y_tgt = sim_fn(y0, pm)
    n_hat, p_sero, _ = j_summary(y_tgt)
    eq = equilibrium_metrics(y_T, y_T5)

    numpyro.deterministic("N_hat", n_hat)
    numpyro.deterministic("p_hat", p_sero)
    numpyro.deterministic("eq_pass", eq["pass"].astype(jnp.float64))
    numpyro.deterministic("max_log_pop_ratio", eq["max_log_pop_ratio"])
    numpyro.deterministic("max_prev_change", eq["max_prev_change"])
    numpyro.deterministic("total_log_ratio", eq["total_log_ratio"])

    eps = targets["eps_prev"]
    p_safe = jnp.clip(p_sero, eps, 1.0 - eps)
    n_obs = jnp.asarray(targets["prison_total"], dtype=jnp.float64)
    n_obs_int = jnp.asarray(targets["prison_total"], dtype=jnp.int64)
    sigma_pop = targets["sigma_pop"]
    with numpyro.plate("age", 6):
        numpyro.sample(
            "prev_obs",
            dist.Binomial(total_count=n_obs_int, probs=p_safe),
            obs=None if obs_prev is None else jnp.asarray(obs_prev, dtype=jnp.int64),
        )
        numpyro.sample(
            "pop_obs",
            dist.LogNormal(loc=jnp.log(n_obs), scale=sigma_pop),
            obs=None if obs_pop is None else jnp.asarray(obs_pop, dtype=jnp.float64),
        )

    if penalty_mode in ("legacy", "smooth"):
        relu = lambda x: jnp.maximum(x, 0.0)  # noqa: E731
        if penalty_mode == "legacy":
            pen = eq_penalty_factor * (
                relu(eq["max_log_pop_ratio"] / 0.01)
                + relu(eq["max_prev_change"] / 0.005)
                + relu(eq["total_log_ratio"] / 0.01)
            )
        else:
            excess = jnp.stack(
                [
                    relu(eq["max_log_pop_ratio"] / 0.01 - 1.0),
                    relu(eq["max_prev_change"] / 0.005 - 1.0),
                    relu(eq["total_log_ratio"] / 0.01 - 1.0),
                ]
            )
            pen = eq_penalty_factor * jnp.sum(excess**2)
        numpyro.factor("equilibrium_penalty", -pen)


def hcv_model(
    sim_fn,
    y0: jax.Array,
    base_params: Dict[str, Any],
    targets: Dict[str, Any],
    priors: Dict[str, Any],
    eq_penalty_factor: float = 1e6,
    penalty_mode: str = "none",
    parameterization: str = "theta",
    whiten_mean: Optional[jax.Array] = None,
    whiten_chol: Optional[jax.Array] = None,
    obs_prev: Optional[jax.Array] = None,
    obs_pop: Optional[jax.Array] = None,
):
    """NumPyro 模型（12 个 log 参数）。

    parameterization:
      "theta"  —— 直接在 theta 空间采样（先验原生采样）；
      "whiten" —— 在单位尺度 u 空间采样，theta = whiten_mean + whiten_chol @ u，
                  先验以 factor 加入（白化可显著改善 NUTS 几何，u 空间近似 N(0,I)）。
    """
    if parameterization == "theta":
        contact_anchor = jnp.asarray(priors["contact"]["anchor_scale"], dtype=jnp.float64)
        beta_anchor = jnp.asarray(priors["beta"]["anchor_scale"], dtype=jnp.float64)
        with numpyro.plate("contact", 6):
            contact_log = numpyro.sample(
                "contact_log",
                dist.Normal(jnp.log(contact_anchor), priors["contact"]["sd"]),
            )
        with numpyro.plate("beta", 6):
            z_beta = numpyro.sample("z_beta", dist.StudentT(priors["beta"]["df"]))
        beta_log = jnp.log(beta_anchor) + priors["beta"]["sd"] * z_beta
        theta = jnp.concatenate([contact_log, beta_log], axis=-1)
        _model_core(
            theta, sim_fn, y0, base_params, targets,
            eq_penalty_factor, penalty_mode, obs_prev, obs_pop,
        )
    elif parameterization == "whiten":
        if whiten_mean is None or whiten_chol is None:
            raise ValueError("whiten parameterization requires whiten_mean/whiten_chol")
        with numpyro.plate("whiten_plate", 12):
            u = numpyro.sample("u", dist.Normal(0.0, 1.0))
        theta = jnp.asarray(whiten_mean, dtype=jnp.float64) + jnp.asarray(
            whiten_chol, dtype=jnp.float64
        ) @ u
        numpyro.deterministic("theta", theta)
        numpyro.factor("prior", legacy_log_prior(theta, priors))
        _model_core(
            theta, sim_fn, y0, base_params, targets,
            eq_penalty_factor, penalty_mode, obs_prev, obs_pop,
        )
    else:
        raise ValueError(f"unknown parameterization: {parameterization}")


def nm_init_values(
    theta_nm: jax.Array, priors: Dict[str, Any], parameterization: str = "theta"
) -> Dict[str, jax.Array]:
    """以 legacy Nelder-Mead 点估计构造 NUTS 初值（z_beta 为反变换）。"""
    if parameterization == "whiten":
        return {"u": jnp.zeros(12, dtype=jnp.float64)}
    contact_log = jnp.asarray(theta_nm[0:6], dtype=jnp.float64)
    beta_anchor = jnp.asarray(priors["beta"]["anchor_scale"], dtype=jnp.float64)
    z_beta = (jnp.asarray(theta_nm[6:12], dtype=jnp.float64) - jnp.log(beta_anchor)) / priors[
        "beta"
    ]["sd"]
    return {"contact_log": contact_log, "z_beta": z_beta}


def make_mcmc(
    sim_fn,
    y0: jax.Array,
    base_params: Dict[str, Any],
    targets: Dict[str, Any],
    priors: Dict[str, Any],
    theta_nm: jax.Array,
    num_warmup: int,
    num_samples: int,
    num_chains: int = 4,
    seed: int = 2026,
    chain_method: str = "parallel",
    target_accept_prob: float = 0.8,
    penalty_mode: str = "none",
    parameterization: str = "whiten",
    whiten_mean: Optional[jax.Array] = None,
    whiten_chol: Optional[jax.Array] = None,
) -> MCMC:
    init = nm_init_values(theta_nm, priors, parameterization=parameterization)
    kernel = NUTS(
        hcv_model,
        init_strategy=init_to_value(values=init),
        dense_mass=False,
        target_accept_prob=target_accept_prob,
    )
    mcmc = MCMC(
        kernel,
        num_warmup=num_warmup,
        num_samples=num_samples,
        num_chains=num_chains,
        chain_method=chain_method,
        progress_bar=True,
        jit_model_args=False,
    )
    mcmc.run(
        jax.random.PRNGKey(seed),
        sim_fn=sim_fn,
        y0=y0,
        base_params=base_params,
        targets=targets,
        priors=priors,
        obs_prev=targets["x_prev"],
        obs_pop=targets["prison_total"],
        penalty_mode=penalty_mode,
        parameterization=parameterization,
        whiten_mean=whiten_mean,
        whiten_chol=whiten_chol,
        extra_fields=("potential_energy", "accept_prob", "num_steps", "diverging"),
    )
    return mcmc
