# =============================================================================
# simulator.py — JAX 移植的 HCV PWID 4-层室 ODE 模拟器
#
# 对应 legacy C++ 实现（legacy/r_source/sim.cpp 与 setup.R）：
#   - 384 个隔室：4 层（D/J/F/X）x 4 肝病期（NC/CC/DC/HCC）
#     x 4 HCV 状态（u/a/c/t）x 6 年龄组；
#   - 连续时间 ODE，RK4 固定步长（dt = 1/365 年），步进后做 10 岁年龄带
#     推进（1/10 每年）并截断负值 —— 与 sim.cpp 的算子分裂完全一致；
#   - 频率依赖感染压：只有活跃吸毒层（D、F）参与；
#   - 全程 JAX 兼容（jit / grad / vmap），默认 float64（x64）。
# =============================================================================
from __future__ import annotations

import json
from functools import partial
from pathlib import Path
from typing import Any, Dict

import jax
import jax.numpy as jnp

jax.config.update("jax_enable_x64", True)

ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------------
# 配置加载
# ---------------------------------------------------------------------------
def _read_json(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def _strip_metadata(d: Dict[str, Any]) -> Dict[str, Any]:
    """移除 JSON 中以下划线开头的说明性字段，避免传入 JIT。"""
    return {k: v for k, v in d.items() if not str(k).startswith("_")}


def load_model_parameters(path: str | Path | None = None) -> Dict[str, Any]:
    path = Path(path) if path else ROOT / "config" / "model_parameters.json"
    return _strip_metadata(_read_json(path))


def load_simulation(path: str | Path | None = None) -> Dict[str, Any]:
    path = Path(path) if path else ROOT / "config" / "simulation.json"
    return _strip_metadata(_read_json(path))


def load_targets(path: str | Path | None = None) -> Dict[str, Any]:
    path = Path(path) if path else ROOT / "config" / "targets.json"
    return _strip_metadata(_read_json(path))


def load_priors(path: str | Path | None = None) -> Dict[str, Any]:
    path = Path(path) if path else ROOT / "config" / "priors.json"
    return _strip_metadata(_read_json(path))


def expand_y0(sim_cfg: Dict[str, Any]) -> jax.Array:
    """把 simulation.json 中的稀疏 y0 展开为 384 维向量。

    顺序与 C++ 一致：idx(s,k,h,i) = s*96 + (k-1)*24 + h*6 + i。
    布局复刻 legacy setup.R 的实际行为（见 simulation.json 注释）：
    D_u 在 flat index 1-6，D_a 在 7-12，241 在 D_c age1（index 12）。
    """
    y0 = jnp.zeros(384, dtype=jnp.float64)
    cfg = sim_cfg["y0"]
    y0 = y0.at[0].set(cfg["D_u_stage_NC_age1"])
    y0 = y0.at[1 + jnp.arange(5)].set(jnp.asarray(cfg["D_u_stage_NC_age2_6"]))
    y0 = y0.at[6 + jnp.arange(6)].set(jnp.asarray(cfg["D_a_stage_NC_age1_6"]))
    y0 = y0.at[12].set(cfg["D_c_stage_NC_age1"])
    return y0


def build_derived_params(p: Dict[str, Any]) -> Dict[str, Any]:
    """预计算基因型加权进展率（对应 sim.cpp 中 ptc_* 的计算）。"""
    rho = p["rho"]
    ptc = {
        "ptc12": (rho * p["r3_NC_CC"] + (1.0 - rho)) * p["p_NC_CC"],
        "ptc23": (rho * p["r3_CC_DC"] + (1.0 - rho)) * p["p_CC_DC"],
        "ptc24": (rho * p["r3_CC_HCC"] + (1.0 - rho)) * p["p_CC_HCC"],
        "ptc34": (rho * p["r3_DC_HCC"] + (1.0 - rho)) * p["p_DC_HCC"],
    }
    out = dict(p)
    out.update(ptc)
    return out


def build_params_from_theta(theta: jax.Array, base_params: Dict[str, Any]) -> Dict[str, Any]:
    """theta[0:6] = log contact row scale；theta[6:12] = log beta inflow scale。
    与 legacy likelihood.R::build_params 一致（C_contact 按行缩放，beta 整体缩放）。"""
    p = build_derived_params(base_params)
    contact_scale = jnp.exp(theta[0:6])
    beta_scale = jnp.exp(theta[6:12])
    C = jnp.asarray(p["C_contact"], dtype=jnp.float64)
    C_scaled = C * contact_scale[:, None]
    p = dict(p)
    p["C_contact"] = C_scaled
    p["beta"] = jnp.asarray(p["beta"], dtype=jnp.float64) * beta_scale
    return p


# ---------------------------------------------------------------------------
# ODE 右侧（对应 sim.cpp rhs()）
# ---------------------------------------------------------------------------
def _force_of_infection(y4: jax.Array, p: Dict[str, Any]) -> jax.Array:
    """频率依赖感染压：lambda_i = q * sum_j C(i,j) * infectious_j / active_j。
    active/infectious 仅统计活跃层 D(s=0) 与 F(s=2)。"""
    C = p["C_contact"]
    active = y4[0].sum(axis=(0, 1)) + y4[2].sum(axis=(0, 1))                # (6,)
    infectious = y4[0, :, 1:3, :].sum(axis=(0, 1)) + y4[2, :, 1:3, :].sum(
        axis=(0, 1)
    )
    ratio = jnp.where(active > 0.0, infectious / active, 0.0)
    return p["q"] * (C @ ratio)


def _trans_mult(t: jax.Array, p: Dict[str, Any]) -> jax.Array:
    m_min, m_max = p["m_min"], p["m_max"]
    return m_min + (m_max - m_min) / (1.0 + jnp.exp(-(t - p["m_t0"]) / p["m_tau"]))


def _mu_eff(p: Dict[str, Any]) -> jax.Array:
    """mu_eff(k,h,i)，形状 (4,4,6)。对应 sim.cpp mu_eff lambda。"""
    base = jnp.asarray(p["mu"], dtype=jnp.float64) * p["omega"]              # (6,)
    base = jnp.broadcast_to(base, (4, 4, 6))
    dc = p["mu_DC"] + base
    hcc = p["mu_HCC"] + base
    k_is_3 = (jnp.arange(4)[:, None, None] == 2).astype(jnp.float64)
    k_is_4 = (jnp.arange(4)[:, None, None] == 3).astype(jnp.float64)
    h_is_3 = (jnp.arange(4)[None, :, None] == 3).astype(jnp.float64)
    mu = (
        k_is_3 * jnp.where(h_is_3 == 1.0, p["psi_DC"] * dc, dc)
        + k_is_4 * jnp.where(h_is_3 == 1.0, p["psi_HCC"] * hcc, hcc)
        + (1.0 - k_is_3 - k_is_4) * base
    )
    return mu


def _rhs(y4: jax.Array, t: jax.Array, p: Dict[str, Any]) -> jax.Array:
    foi = _force_of_infection(y4, p)
    gam = _trans_mult(t, p) * foi                                            # (6,)
    mu = _mu_eff(p)                                                          # (4,4,6)
    svr = jnp.asarray(
        [p["alpha_NC"], p["alpha_NC"], p["alpha_DC_pos"], p["alpha_HCC"]],
        dtype=jnp.float64,
    )                                                                        # (4,)
    l1 = jnp.asarray(p["lambda1"], dtype=jnp.float64)
    l2 = jnp.asarray(p["lambda2"], dtype=jnp.float64)
    l3 = jnp.asarray(p["lambda3"], dtype=jnp.float64)
    pi = p["pi_recid"]
    iota1, iota2 = p["iota1"], p["iota2"]
    kappa = p["kappa"]

    prog_out = jnp.asarray(
        [p["ptc12"], p["ptc23"] + p["ptc24"], p["ptc34"], 0.0],
        dtype=jnp.float64,
    )                                                                        # (4,)
    age_ok = (jnp.arange(6) + 1.0) >= p["tau_min_age"]                       # (6,)
    treat_s = jnp.asarray(p["tau_stratum"], dtype=jnp.float64)[:, None] * age_ok[None, :]
    tau = jnp.asarray(p["tau"], dtype=jnp.float64)                           # (4,)
    treat_rate = tau[:, None, None] * treat_s[None, :, :]                    # (k,s,i)

    def prog_c_in(y_s: jax.Array) -> jax.Array:
        """进入慢性池的进展流入（同层内上一肝病期慢性患者）。形状 (4,6)。"""
        return jnp.stack(
            [
                jnp.zeros(6, dtype=jnp.float64),
                p["ptc12"] * y_s[0, 2, :],
                p["ptc23"] * y_s[1, 2, :],
                p["ptc24"] * y_s[1, 2, :] + p["ptc34"] * y_s[2, 2, :],
            ]
        )

    # ---- D 层（s=0）：从未被捕，活跃 ----
    yD = y4[0]
    Du, Da, Dc, Dt = yD[:, 0, :], yD[:, 1, :], yD[:, 2, :], yD[:, 3, :]
    trD = treat_rate[:, 0, :]                                                # (4,6)
    dD = jnp.zeros_like(yD)
    dD = dD.at[:, 0, :].add(
        jnp.concatenate([p["beta"][None, :], jnp.zeros((3, 6), dtype=jnp.float64)])
        + (kappa / iota1) * Da
        + (svr[:, None] / iota2) * Dt
        - (gam[None, :] + l1[None, :] + mu[:, 0, :]) * Du
    )
    dD = dD.at[:, 1, :].add(
        gam[None, :] * Du - (1.0 / iota1 + l1[None, :] + mu[:, 1, :]) * Da
    )
    dD = dD.at[:, 2, :].add(
        ((1.0 - kappa) / iota1) * Da
        + ((1.0 - svr[:, None]) / iota2) * Dt
        + prog_c_in(yD)
        - (prog_out[:, None] + trD + l1[None, :] + mu[:, 2, :]) * Dc
    )
    dD = dD.at[:, 3, :].add(
        trD * Dc - (1.0 / iota2 + l1[None, :] + mu[:, 3, :]) * Dt
    )
    # ---- F 层（s=2）：曾经被捕，活跃 ----
    yF = y4[2]
    yJ = y4[1]
    Fu, Fa, Fc, Ft = yF[:, 0, :], yF[:, 1, :], yF[:, 2, :], yF[:, 3, :]
    Ju, Ja, Jc, Jt = yJ[:, 0, :], yJ[:, 1, :], yJ[:, 2, :], yJ[:, 3, :]
    trF = treat_rate[:, 2, :]
    dF = jnp.zeros_like(yF)
    dF = dF.at[:, 0, :].add(
        pi * l2[None, :] * Ju
        + (kappa / iota1) * Fa
        + (svr[:, None] / iota2) * Ft
        - (gam[None, :] + l3[None, :] + mu[:, 0, :]) * Fu
    )
    dF = dF.at[:, 1, :].add(
        pi * l2[None, :] * Ja
        + gam[None, :] * Fu
        - (1.0 / iota1 + l3[None, :] + mu[:, 1, :]) * Fa
    )
    dF = dF.at[:, 2, :].add(
        pi * l2[None, :] * Jc
        + ((1.0 - kappa) / iota1) * Fa
        + ((1.0 - svr[:, None]) / iota2) * Ft
        + prog_c_in(yF)
        - (prog_out[:, None] + trF + l3[None, :] + mu[:, 2, :]) * Fc
    )
    dF = dF.at[:, 3, :].add(
        pi * l2[None, :] * Jt + trF * Fc - (1.0 / iota2 + l3[None, :] + mu[:, 3, :]) * Ft
    )
    # ---- J 层（s=1）：在押，不注射 ----
    dJ = (
        l1[None, None, :] * yD
        + l3[None, None, :] * yF
        - (l2[None, None, :] + mu) * yJ
    )
    trJ = treat_rate[:, 1, :]
    dJ = dJ.at[:, 1, :].add(-(1.0 / iota1) * Ja)
    dJ = dJ.at[:, 2, :].add(
        ((1.0 - kappa) / iota1) * Ja
        + ((1.0 - svr[:, None]) / iota2) * Jt
        + prog_c_in(yJ)
        - (prog_out[:, None] + trJ) * Jc
    )
    dJ = dJ.at[:, 3, :].add(trJ * Jc - (1.0 / iota2) * Jt)
    # ---- X 层（s=3）：退出吸毒/出狱后不再被捕 ----
    yX = y4[3]
    Xa, Xc, Xt = yX[:, 1, :], yX[:, 2, :], yX[:, 3, :]
    trX = treat_rate[:, 3, :]
    dX = (1.0 - pi) * l2[None, None, :] * yJ - mu * yX
    dX = dX.at[:, 1, :].add(-(1.0 / iota1) * Xa)
    dX = dX.at[:, 2, :].add(
        ((1.0 - kappa) / iota1) * Xa
        + ((1.0 - svr[:, None]) / iota2) * Xt
        + prog_c_in(yX)
        - (prog_out[:, None] + trX) * Xc
    )
    dX = dX.at[:, 3, :].add(trX * Xc - (1.0 / iota2) * Xt)
    return jnp.stack([dD, dJ, dF, dX])


def _age_shift(y4: jax.Array, dt: jax.Array) -> jax.Array:
    """10 岁年龄带推进：y_change = y[i]/10 * dt 移到 i+1 带（60+ 开放）。"""
    shift = y4[..., :5] / 10.0 * dt
    loss = jnp.concatenate([shift, jnp.zeros_like(y4[..., 5:])], axis=-1)
    gain = jnp.concatenate([jnp.zeros_like(y4[..., :1]), shift], axis=-1)
    return y4 - loss + gain


def _step(y4: jax.Array, t: jax.Array, dt: jax.Array, p: Dict[str, Any]) -> jax.Array:
    """一个固定步长：RK4 求导 -> 截断负值 -> 年龄推进（与 sim.cpp 一致）。"""
    k1 = _rhs(y4, t, p)
    k2 = _rhs(y4 + 0.5 * dt * k1, t + 0.5 * dt, p)
    k3 = _rhs(y4 + 0.5 * dt * k2, t + 0.5 * dt, p)
    k4 = _rhs(y4 + dt * k3, t + dt, p)
    y_new = y4 + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    y_new = jnp.maximum(y_new, 0.0)
    return _age_shift(y_new, dt)


# ---------------------------------------------------------------------------
# 主入口
# ---------------------------------------------------------------------------
def row_for_time(t: float, t_start: float, dt: float) -> int:
    return int(round((t - t_start) / dt))


def make_sim_fn(sim_cfg: Dict[str, Any]):
    """构造 JIT 编译的模拟函数。

    返回 sim_fn(y0, p) -> (y_T, y_T5, y_target)，其中：
      y_T      终点 T 状态（(4,4,4,6)）
      y_T5     终点前 t_lag 年的状态
      y_target 校准目标时刻的状态
    n_steps/行号在构建时固定（静态），保证 lax.scan 可 JIT。
    """
    t_start = float(sim_cfg["t_start"])
    t_end = float(sim_cfg["t_end"])
    dt = float(sim_cfg["dt"])
    n_steps = int(round((t_end - t_start) / dt)) + 1
    row_t5 = row_for_time(t_end - float(sim_cfg.get("equilibrium_t_lag", 5.0)), t_start, dt)
    row_target = row_for_time(float(sim_cfg.get("target_time", t_end)), t_start, dt)
    t_start_a = jnp.asarray(t_start, dtype=jnp.float64)
    dt_a = jnp.asarray(dt, dtype=jnp.float64)

    def sim_fn(y0: jax.Array, p: Dict[str, Any]) -> tuple[jax.Array, jax.Array, jax.Array]:
        y4 = y0.reshape(4, 4, 4, 6)

        def body(carry, step):
            y4c, t, y_t5, y_tgt = carry
            y_new = _step(y4c, t, dt_a, p)
            y_t5 = jnp.where(step == row_t5, y_new, y_t5)
            y_tgt = jnp.where(step == row_target, y_new, y_tgt)
            return (y_new, t + dt_a, y_t5, y_tgt), None

        (y_T, _, y_T5, y_tgt), _ = jax.lax.scan(
            body, (y4, t_start_a, y4, y4), jnp.arange(1, n_steps)
        )
        return y_T, y_T5, y_tgt

    return jax.jit(sim_fn)


def j_summary(y4: jax.Array) -> tuple[jax.Array, jax.Array, jax.Array]:
    """J 层（在押）年龄组摘要：N_hat、p_sero、p_viremic（均为 (6,)）。
    legacy 中 sero 与 viremic 使用相同的 h∈{1,2,3} 集合，因此二者相等。"""
    yJ = y4[1]
    n_hat = yJ.sum(axis=(0, 1))
    sero = yJ[:, 1:4, :].sum(axis=(0, 1))
    p_sero = jnp.where(n_hat > 0.0, sero / n_hat, 0.0)
    return n_hat, p_sero, p_sero


def equilibrium_metrics(y_T: jax.Array, y_T5: jax.Array) -> Dict[str, jax.Array]:
    """对应 legacy equilibrium.R 的 T vs T-5 稳定性门。"""
    n_T, p_T, _ = j_summary(y_T)
    n_T5, p_T5, _ = j_summary(y_T5)
    log_ratio = jnp.abs(jnp.log(n_T / n_T5))
    prev_change = jnp.abs(p_T - p_T5)
    total_log_ratio = jnp.abs(jnp.log(y_T.sum() / y_T5.sum()))

    hcv_T = y_T[:, 1:4, :].sum()
    hcv_T5 = y_T5[:, 1:4, :].sum()
    dc_T = y_T[:, 2, :].sum()
    dc_T5 = y_T5[:, 2, :].sum()
    hcc_T = y_T[:, 3, :].sum()
    hcc_T5 = y_T5[:, 3, :].sum()
    state_log_ratio = jnp.maximum(
        jnp.abs(jnp.log(hcv_T / hcv_T5)),
        jnp.maximum(
            jnp.abs(jnp.log(dc_T / dc_T5)),
            jnp.abs(jnp.log(hcc_T / hcc_T5)),
        ),
    )
    comp = jnp.abs(jnp.log(jnp.maximum(y_T, 1e-6) / jnp.maximum(y_T5, 1e-6)))
    max_comp_log_ratio = jnp.max(comp)

    pass_ = (
        jnp.all(log_ratio <= 0.01)
        & jnp.all(prev_change <= 0.005)
        & (state_log_ratio <= 0.01)
        & (max_comp_log_ratio <= 0.02)
    )
    return {
        "pass": pass_,
        "max_log_pop_ratio": jnp.max(log_ratio),
        "max_prev_change": jnp.max(prev_change),
        "total_log_ratio": total_log_ratio,
        "state_log_ratio": state_log_ratio,
        "max_comp_log_ratio": max_comp_log_ratio,
    }


def target_summary(y_T: jax.Array, y_T5: jax.Array, y_target: jax.Array):
    """一次模拟的完整摘要：目标时刻患病率/人口 + 平衡态指标。"""
    n_hat, p_sero, p_vir = j_summary(y_target)
    eq = equilibrium_metrics(y_T, y_T5)
    return {
        "N_hat": n_hat,
        "p_sero": p_sero,
        "p_viremic": p_vir,
        **eq,
    }
