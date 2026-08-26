"""JAX 模拟器与 legacy R/C++ 参考输出的一致性测试。"""

import numpy as np
import pandas as pd
import pytest

import jax
import jax.numpy as jnp

from src.simulator import (
    _age_shift,
    build_params_from_theta,
    equilibrium_metrics,
    j_summary,
)


def test_nm_summary_matches_legacy(sim_fn, y0, base_params, theta_nm):
    """在 legacy Nelder-Mead 点估计处，J 层摘要与 R 参考完全一致。"""
    th = jnp.asarray(theta_nm)
    pm = build_params_from_theta(th, base_params)
    y_T, y_T5, y_tgt = sim_fn(y0, pm)
    n_hat, p_sero, p_vir = j_summary(y_tgt)
    ref = pd.read_csv("tests/data/summary_nm.csv")
    assert np.max(np.abs(np.array(n_hat) - ref["N_hat"].values)) < 1e-8
    assert np.max(np.abs(np.array(p_sero) - ref["p_sero"].values)) < 1e-8
    assert np.array_equal(np.array(p_sero), np.array(p_vir))


def test_terminal_state_matches_legacy(sim_fn, y0, base_params, theta_nm):
    """T 与 T-5 的 384 隔室状态与 R 参考一致（浮点容差）。"""
    th = jnp.asarray(theta_nm)
    pm = build_params_from_theta(th, base_params)
    y_T, y_T5, _ = sim_fn(y0, pm)
    ref_T = pd.read_csv("tests/data/state_T.csv")
    ref_T5 = pd.read_csv("tests/data/state_T5.csv")
    assert np.max(np.abs(np.array(y_T).ravel() - ref_T.iloc[0, 1:].values)) < 1e-8
    assert np.max(np.abs(np.array(y_T5).ravel() - ref_T5.iloc[0, 1:].values)) < 1e-8


def test_equilibrium_pass_at_nm(sim_fn, y0, base_params, theta_nm):
    th = jnp.asarray(theta_nm)
    pm = build_params_from_theta(th, base_params)
    y_T, y_T5, _ = sim_fn(y0, pm)
    eq = equilibrium_metrics(y_T, y_T5)
    assert bool(np.asarray(eq["pass"]))
    assert float(eq["max_log_pop_ratio"]) < 0.01
    assert float(eq["max_prev_change"]) < 0.005


def test_simulator_jit_and_grad(sim_fn, y0, base_params, theta_nm):
    """JAX 兼容性：jit 与 grad 可正常计算。"""
    th = jnp.asarray(theta_nm)

    def summary_sum(t):
        pm = build_params_from_theta(t, base_params)
        _, _, y_tgt = sim_fn(y0, pm)
        n, p, _ = j_summary(y_tgt)
        return jnp.sum(n) + jnp.sum(p)

    val = summary_sum(th)
    grad = jax.jit(jax.grad(summary_sum))(th)
    assert np.isfinite(np.asarray(val))
    assert grad.shape == (12,)
    assert np.all(np.isfinite(np.asarray(grad)))


def test_age_shift_math():
    y = jnp.ones((2, 2, 2, 6))
    dt = jnp.asarray(1.0)
    out = _age_shift(y, dt)
    # 每个 (s,k,h) 行：i=0 减少 1/10；i=5 增加 1/10；总量守恒
    assert np.allclose(np.asarray(out[0, 0, 0, 0]), 0.9)
    assert np.allclose(np.asarray(out[0, 0, 0, 5]), 1.1)
    assert np.allclose(np.asarray(out.sum()), y.sum())
