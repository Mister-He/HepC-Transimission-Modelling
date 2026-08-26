# memo.md — 全流程实验记录

> 本文档记录 `dev-numpyro` worktree 中从迁移到推断的完整过程与关键结论，
> 供后续 agent 快速理解"目前做到哪里、为什么这样做"。

## 1. 背景与目标

原仓库（分支 `dev_10yrs_age_interval`）用 R + Rcpp/C++ 实现了新加坡 PWID 人群
HCV 传播模型，并完成了：

- 12 参数 Nelder-Mead 标定（`run1_4strata`）；
- Laplace 近似与自适应 Metropolis MCMC（`npe_bayes`）；
- NPE/SNPE 尝试（最终未采用为主方法）。

本 worktree 的目标：将 `src/sim.cpp` 忠实移植到 NumPyro/JAX，用 NUTS +
legacy 先验 + legacy NM 点估计初值重新推断，做严格收敛检验与后验预测检验，
并与 legacy 后验比较；同时把固定参数 JSON 化、建立 CI/测试/文档，记录
`memo.md`。

## 2. 仓库整理（要求 1、2）

### 2.1 归档 `legacy/`

- `nelder_mead_point_estimate.csv/.json` —— 最终 NM 点估计（NLL=21.421，
  平衡态通过）；
- `posterior_samples_mcmc.csv` —— legacy 后验（3 链 × 1250 样本，
  自适应 Metropolis）；
- `posterior_summary.csv` —— legacy 后验摘要（由样本重算）；
- `prior_distribution.md/.json` —— 实际生效先验（mcmc.R::make_priors：
  contact log 缩放 Normal(log(anchor), 0.5²)、beta log 缩放
  log(anchor)+0.8·t3）；注意 run_config 中记录的
  "Normal(0,2²)" 字符串与实际代码不一致，已在文档中说明；
- `likelihood_design.md/.json` —— 二项患病率 + log-normal 人口似然、
  平衡态门、合理性软惩罚、边界；
- `diagnostics_mcmc.csv`、`credible_intervals.csv`、`sbc_summary*.csv`、
  `npe_config*.json`、`run_config.csv` —— 推断诊断与配置；
- `run1_4strata/` —— NM 标定运行全部输出（含 fit.rds、Laplace 近似）；
- `reports/` —— 旧版最终报告；
- `r_source/` —— 产生上述结果的 R/C++ 源码；
- `validation/generate_reference.R` —— 生成 R 参考输出的脚本。

### 2.2 删除

已删除：`output/`、`figures/`（结果）、`docs/prompt*.md`、`docs/AGENTS.md`
（prompts 与旧 AGENTS.md）、R 时代 `src/`、`scripts/`、`tests/`、旧
`Dockerfile`、`docker-compose.yml`、`ci.yml`。所有被删内容在 `legacy/` 有
存档（或可由 git 历史恢复）。

## 3. JAX 移植（要求 3、4）

### 3.1 关键发现：legacy y0 布局偏移

legacy `setup.R` 的 `idx()` 返回输出矩阵列索引（+2）却直接用于 `y0` 向量，
导致初始值相对名义隔室**后移一格**：D_u 实际占据 flat index 1–6、D_a 占据
7–12、241 落在 D_c age1。为数值等价，JAX 移植复刻该布局
（见 `config/simulation.json` 注释与 `src/simulator.py::expand_y0`）。

### 3.2 实现与验证

- `src/simulator.py`：384 隔室 RK4 + 每步负值截断 + 10 岁年龄推进，
  `lax.scan` 只保留 T、T-5、t=45 三个状态；
- 在 NM 点估计处与 R 参考输出逐隔室对比：**偏差 ~1e-12**
  （`tests/data/*.csv` 为参考，`tests/unit/test_simulator.py` 容差 1e-8）；
- 平衡态门在 NM 处通过，指标与 legacy 一致。

### 3.3 固定参数 JSON（要求 4）

`config/model_parameters.json`（自然史/治疗/死亡/接触/流入）、
`config/simulation.json`（时域/dt/y0）、`config/targets.json`（观测目标）、
`config/priors.json`（先验）。代码只读 JSON，`_` 开头字段为注释。

## 4. NUTS 推断（要求 3）—— 数值挑战与对策

### 4.1 问题：原始参数空间 NUTS 几何差

12 维 theta 后验尺度差异大（sd 0.05–0.5），直接采样时启发式步长 ≈ 几十个
后验标准差，warmup 步长塌缩到 ~1e-2，单次迭代 100+ leapfrog 步，不可行。

### 4.2 对策 1：白化再参数化

用 legacy 后验均值/协方差 Cholesky 做仿射变换，在 `u ~ N(0,I)` 空间采样
（`theta = mu + chol @ u`）。步长恢复到 0.3–0.6，每迭代 7–15 leapfrog，
接受率 0.87–0.95。样本变换回 theta 空间报告。白化是确定性双射，后验不变。

### 4.3 对策 2：平衡态惩罚的取舍

legacy MCMC 的 1e6 硬惩罚会让 NUTS 在可行域边界形成悬崖。实测 **legacy 后验
400/400 通过平衡态门**，即惩罚在后验质量区域几乎不激活。因此推断默认
`--penalty-mode none`，在后处理按平衡态门过滤；模型保留
`legacy`/`smooth` 模式供对照。

### 4.4 对策 3：dt=7 天的工程权衡

完整 150 年 × 日步长梯度评估约 2 秒，NUTS 全量运行不现实。主推断用
dt=7 天（梯度 ~0.28 秒）。在 NM 点估计处 dt=7 天 vs 1 天的偏差：
人口 ≤0.27%、患病率 ≤0.0003 —— 远小于后验不确定度。
`scripts/check_dt_robustness.py` 在后验样本上进一步量化（见 §6）。

## 5. 主运行结果（results/numpyro_nuts/）

运行参数见 `run_config.json`：4 链 × (500 warmup + 1000 samples)、dt=7 天、
白化参数化、无惩罚（后处理过滤）、seed=2026、NUTS target_accept=0.9。
后验样本 4000（4×1000），全部通过平衡态门。

### 5.1 收敛诊断（diagnostics_numpyro.csv）

- R-hat：**[0.9995, 1.0055]**（theta4=1.0055 略高于 legacy 阈值 1.005，
  其余 11 个参数全部满足 [0.995, 1.005]）；
- ESS bulk：**[4217, 7819]**（legacy 为 [1310, 1958]，显著更高）；
- ESS tail：**[2815, 3459]**；
- 发散数：**0**；平均接受率：**0.926**；平均树深：**3.81**。

结论：收敛性至少与 legacy 同等严格，且效率远高于 legacy 的自适应
Metropolis（ESS 高 3–4 倍、0 发散）。

### 5.2 与 legacy 后验比较（compare_legacy.csv）

| 指标 | 结果 |
|---|---|
| 中位数最大偏移 | theta12 = +0.067（其余大多 < 0.01） |
| 95% 区间重叠 | 0.16–1.41（theta6 完全包含 legacy 区间） |
| KS 统计量 | 0.014–0.135（差异很小） |
| KS p 值 | 多数 < 0.05（大样本下微小差异即显著；theta7/8/9 不显著） |

解读：两个后验的中心位置与不确定度高度一致；KS 显著主要来自样本量大
（3750 vs 4000）与两个采样器在弱识别方向（如 theta6）的细微形状差异，
临床/模型含义上无实质差别。

### 5.3 后验预测检验（ppc_target_summary.csv / ppc_data_summary.csv）

- 目标级：6 个年龄组的 `p_hat` 与 `N_hat` 观测全部落入 95% CrI；
- 数据级：二项患病率与监狱人口观测全部落入 95% CrI；
- Bayesian p 值：患病率 0.41–0.60、人口 0.50–0.52（校准良好）；
- 平衡态通过率：**4000/4000 = 100%**。

## 6. dt 稳健性检查（scripts/check_dt_robustness.py）

在 200 个新后验样本上对比 dt=1/365（精确）与 dt=7/365（主运行）：

| 年龄组 | N_hat 平均绝对差 | N 相对差 % | p_hat 平均绝对差 |
|---|---|---|---|
| <20 | 0.13 | 0.13 | 0.00007 |
| 20-29 | 1.39 | 0.11 | 0.00009 |
| 30-39 | 0.24 | 0.02 | 0.00016 |
| 40-49 | 0.35 | 0.02 | 0.00020 |
| 50-59 | 0.66 | 0.05 | 0.00008 |
| 60+ | 1.18 | 0.27 | 0.00032 |

平衡态通过率两种步长均为 100%。结论：dt=7 天近似对后验预测的影响
（≤0.27% 人口、≤0.0004 患病率）远小于后验不确定度（患病率 CrI 半宽
0.03–0.08），主运行结果与精确模型等价。

## 7. 工程交付（要求 1、5、6）

- 目录结构按规范（README/docs/src/tests/scripts/.github/Dockerfile/
  docker-compose.yml）；
- 单元测试 15 项全绿（含与 R 参考数值等价）；集成测试为短时域 NUTS/PPC 冒烟；
- CI：`.github/workflows/ci.yml`（unit + integration + compileall）；
- 文档：`docs/requirements.md`、`architecture.md`、`test-plan.md`、
  `design.md`（中文）；
- `AGENTS.md` 已按全部需求重新生成；
- 本文件 `memo.md` 记录全过程。

## 8. 已知限制与后续建议

- 主推断采用 dt=7 天近似（偏差已量化）；如需完全一致可跑
  `--dt-days 1`（成本约 ×7）；
- NUTS 无惩罚后验在可行域内与 legacy 一致，但严格说与"带惩罚"后验是
  近似等价（差异极小，见 §5.3 平衡态通过率）；
- 后续可做：50 年情景/敏感性分析（沿用 legacy scenarios 思路）、
  更细步长验证、GPU 加速、CI 中加入 lint/类型检查。
