"""配置 JSON 与初始条件的正确性测试。"""

import numpy as np


def test_model_parameters_shape(base_params):
    assert len(base_params["mu"]) == 6
    assert len(base_params["beta"]) == 6
    assert len(base_params["lambda1"]) == 6
    assert np.asarray(base_params["C_contact"]).shape == (6, 6)
    assert base_params["q"] > 0
    assert base_params["omega"] > 1


def test_targets_consistency(targets):
    x = np.asarray(targets["x_prev"])
    n = np.asarray(targets["prison_total"])
    p = np.asarray(targets["prev_supplied"])
    assert len(x) == len(n) == len(p) == 6
    # x 与 n/p 的重建误差应在四舍五入范围内
    assert np.all(np.abs(x / n - p) <= 0.5 / n)
    assert np.all(x <= n)


def test_y0_expansion_matches_legacy_layout(y0):
    """复刻 legacy setup.R 的实际布局（初始值相对名义隔室后移一格）。"""
    y = np.asarray(y0)
    assert y.shape == (384,)
    # flat index 0 = 0；1..5 = D_u age2..6；6..11 = D_a age1..6；12 = D_c age1
    assert y[0] == 0.0
    assert np.allclose(y[1:6], [119.0, 512.0, 595.0, 480.0, 381.0])
    assert np.allclose(y[6:12], [303.0, 26.0, 95.0, 164.0, 169.0, 183.0])
    assert y[12] == 241.0
    assert np.count_nonzero(y) == 12


def test_simulation_cfg(sim_cfg):
    assert sim_cfg["dt"] == 1.0 / 365.0
    assert sim_cfg["t_end"] == 150.0
    assert sim_cfg["target_time"] == 45.0
