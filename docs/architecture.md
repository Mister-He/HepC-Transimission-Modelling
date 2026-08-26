# 架构文档（Architecture）

## 1. 目录结构

```text
.
├── README.md / README.zh-CN.md   项目说明（中/英）
├── AGENTS.md                     面向 agent 的工作说明
├── memo.md                       全流程实验记录
├── config/                       固定参数 JSON（模型/模拟/目标/先验）
├── legacy/                       R/C++ 时代结果与源码的只读存档
├── docs/                         中文文档（需求/架构/测试计划/设计理念）
├── src/
│   ├── simulator.py              JAX 模拟器（sim.cpp 移植）
│   ├── model.py                  NumPyro 概率模型（先验+似然+惩罚）
│   ├── inference.py              诊断/PPC/与 legacy 比较的工具函数
│   └── __init__.py
├── scripts/
│   └── run_nuts.py               NUTS 推断流水线入口
├── tests/
│   ├── unit/                     单元测试（含与 R 参考的数值一致性）
│   └── integration/              集成测试（短时域 NUTS/PPC 冒烟）
├── results/numpyro_nuts/         本次推断输出（后验、诊断、PPC、图）
├── .github/workflows/ci.yml      CI
├── Dockerfile / docker-compose.yml
└── requirements*.txt
```

## 2. 数据流

```text
config/*.json
   │
   ▼
src/simulator.py ──> (y_T, y_T5, y_target)  ──> j_summary / equilibrium_metrics
   │  theta = [log contact scale, log beta scale]（12 维）
   ▼
src/model.py（NumPyro）
   先验（legacy）──> theta ──> 模拟 ──> p_hat / N_hat ──> 似然（Binomial + LogNormal）
                                              └──> 平衡态门（可选惩罚）
   ▼
scripts/run_nuts.py（NUTS，白化参数化，NM 点估计初值）
   ▼
收敛诊断（R-hat / ESS / 发散）→ 后验摘要 → PPC → 与 legacy 后验比较
   ▼
results/numpyro_nuts/*.csv / *.png
```

## 3. 模块说明

### 3.1 `src/simulator.py`

- `load_model_parameters / load_simulation / load_targets / load_priors`：
  读取 `config/*.json`（自动剔除 `_` 开头的说明字段）；
- `expand_y0`：把 JSON 中的稀疏初始条件展开为 384 维向量。**注意**：legacy
  `setup.R` 的 `idx()` 返回输出列索引（+2）却直接用于 `y0` 向量，导致初始值
  相对名义隔室后移一格；为数值等价，本模块复刻该布局（见
  `config/simulation.json` 注释）；
- `build_derived_params`：预计算基因型加权进展率（`ptc_*`）；
- `build_params_from_theta`：12 参数 → 缩放后的 `C_contact`（按行）与 `beta`；
- `_rhs / _step / _age_shift`：ODE 右侧、RK4 步、年龄推进；
- `make_sim_fn`：构造 JIT 函数 `sim_fn(y0, p) -> (y_T, y_T5, y_target)`，
  `lax.scan` 推进全部时间步，仅保留需要的三个状态（内存友好）；
- `j_summary`：J 层年龄组 `N_hat / p_sero`；
- `equilibrium_metrics`：T vs T-5 稳定性指标与通过判定。

### 3.2 `src/model.py`

- `hcv_model`：NumPyro 模型。支持两种参数化：
  - `"theta"`：直接在 12 维 log 参数空间采样（先验原生采样）；
  - `"whiten"`（默认）：在单位尺度 `u ~ N(0, I)` 空间采样，
    `theta = mu + chol(L) @ u`，其中 `mu`、`chol` 由 legacy 后验均值/协方差
    给出；先验以 `numpyro.factor` 计入。白化显著改善 NUTS 几何（避免步长
    塌缩），样本再变换回 theta 空间报告；
- 似然：`Binomial`（患病率）+ `LogNormal`（人口），与 legacy 一致；
- 平衡态惩罚：`penalty_mode ∈ {legacy, smooth, none}`。实测 legacy 后验
  100% 通过平衡态门，故默认 `none`（后处理过滤）；`legacy` 保留 1e6 硬惩罚
  原样；`smooth` 为平滑二次版本；
- `make_mcmc`：NUTS（对角质量矩阵）+ `init_to_value`（NM 点估计初值）。

### 3.3 `scripts/run_nuts.py`

流水线：

1. 读取配置与 legacy 点估计/后验；若用 `whiten` 参数化则计算白化变换；
2. `MCMC(NUTS(...))` 多链采样；
3. arviz 收敛诊断（R-hat / ESS bulk & tail / 发散数 / 接受率）；
4. `Predictive` 后验预测（目标级 + 数据级），按平衡态过滤；
5. 与 `legacy/posterior_samples_mcmc.csv` 逐参数比较（中位数/区间/KS）；
6. 输出 CSV 与图到 `results/numpyro_nuts/`。

## 4. 与 legacy 的关键差异（已量化）

| 项目 | legacy | 本仓库 | 影响 |
|---|---|---|---|
| 实现 | R + Rcpp/C++ | JAX（x64） | 数值一致（~1e-12） |
| 采样器 | 自适应 Metropolis | NUTS（白化） | 更高效、无梯度信息丢失 |
| 平衡态惩罚 | 1e6 硬惩罚 | 默认无（后处理过滤） | 后验质量区域 100% 通过门，等价 |
| 模拟步长 | 1/365 年 | 默认 1/365；推断主运行 7 天 | NM 摘要偏差 ≤ 0.27%（N）、0.0003（p） |
| 参数配置 | R 列表硬编码 | JSON | 可查看/修改 |

## 5. 部署与运行

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements-dev.txt

pytest tests/unit -v                 # 单元测试（含与 R 参考一致性）
pytest tests/integration -v          # 集成冒烟

python scripts/run_nuts.py \
  --num-warmup 500 --num-samples 1000 --num-chains 4 \
  --dt-days 7 --chain-method parallel --out-dir results/numpyro_nuts
```
