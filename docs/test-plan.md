# 测试计划（Test Plan）

## 1. 目标

每次改动验证：

1. 配置 JSON 正确加载、初始条件与 legacy 布局一致；
2. JAX 模拟器与 legacy R/C++ 输出数值等价（容差 1e-8）；
3. 似然/先验与 legacy 标定输出一致；
4. 推断工具（摘要/比较/PPC）行为正确；
5. 端到端 NumPyro NUTS 与 PPC 流水线可运行（短时域冒烟）。

## 2. 测试层级

### 单元测试 —— `tests/unit/`

| 文件 | 覆盖内容 |
|---|---|
| `test_config.py` | 参数形状、目标一致性、y0 展开布局、模拟设置 |
| `test_simulator.py` | NM 点估计处与 R 参考的摘要/终态一致性、平衡态门、JIT/grad、年龄推进数学 |
| `test_likelihood.py` | 统计 NLL 与 legacy solutions.csv 对照、先验对数密度、边界因子 |
| `test_inference.py` | theta 长表、摘要、与 legacy 比较（KS）、PPC 摘要 |

### 集成测试 —— `tests/integration/`

| 文件 | 覆盖内容 |
|---|---|
| `test_nuts_smoke.py` | 短时域（5 年、7 天步长）NUTS 运行、样本有限、arviz 诊断；Predictive PPC 冒烟 |

## 3. 回归基准（reference data）

`tests/data/` 由 `legacy/validation/generate_reference.R` 用 R/C++ 模型生成：

- `summary_nm.csv`：NM 点估计处 t=45 的 J 层摘要；
- `state_T.csv` / `state_T5.csv`：t=150 / t=145 全 384 隔室状态；
- `traj_subset.csv`：每 1000 步的轨迹子集。

单元测试以 1e-8 绝对容差对照（实测偏差 ~1e-12）。

## 4. CI

`.github/workflows/ci.yml` 在 push/PR 时运行：

1. Python 3.13 + `requirements-dev.txt`；
2. `pytest tests/unit`；
3. `pytest tests/integration`（短时域冒烟，避免重计算）；
4. `compileall` 语法检查。

完整推断（500+ 样本 × 4 链）不在 CI 中运行，按需通过
`scripts/run_nuts.py` 或 `docker compose run infer` 执行。

## 5. 手动验证清单

- [ ] `pytest tests/unit -v` 全绿；
- [ ] `pytest tests/integration -v` 全绿；
- [ ] 完整 NUTS 运行：R-hat ∈ [0.995, 1.005]、ESS > 1000、发散数报告；
- [ ] PPC：观测值落入 95% CrI，平衡态通过率报告；
- [ ] 与 legacy 后验比较表（`compare_legacy.csv`）生成。
