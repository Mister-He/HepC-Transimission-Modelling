# 项目总览 / Project Overview

（中文为主，关键术语附英文）本文件是项目完成后的总说明：从零到最终结果
做了哪些实验、每一步尝试了什么、各部分文件代表什么、是哪次实验的产物。

## 1. 项目目标

建立一个新加坡 PWID 的年龄结构 HCV 分室模型，使其在均衡状态下拟合 2015
年前后监狱（J 层）六个年龄组的 HCV 血清阳性率与监狱人口，尤其关注 60+
年龄组；所有改动必须有文献或机制依据，并用 Laplace 近似给出拟合点的
95% 区间。

## 2. 最终基线（本修订）

| 项目 | 设定 |
|---|---|
| 模型 | 3 层结构 D/J/X，4 个肝脏阶段，5 个 HCV 状态，6 个年龄组，360 分室 |
| 进展率 | p_NC_CC=0.027（Thein 2008）；p_CC_DC=0.0788、p_CC_HCC=0.0479（Alazawi 2010，未治疗亚组）；p_DC_HCC=0.0464（Rivera-Irigoin 2006，95% CI 上限） |
| 传播 | 常数（m=1），历史水平并入接触矩阵行缩放；**不再拟合 m_min/m_max** |
| 死亡率 | SingStat 2015 年龄别死亡率 × PWID SMR 14.68（Mathers 2013） |
| 拟合参数 | 12 个：6 接触行缩放 + 6 beta 流入缩放；**无超额死亡率参数**（60+ 已可拟合） |
| beta_scale>1 | 3-6 组满足；1-2 组（<20、20-29）不满足（监狱人口小，见 DECISIONS.md） |
| 优化 | 确定性多重起点 Nelder-Mead；无 HMC |
| 不确定性 | Laplace 近似 + 均衡可行性过滤的 Monte Carlo 区间 |
| 图件 | 全部 ggplot2，出版级 |

## 3. 实验时间线与分步说明

> 上一轮 agent 生成的所有文件（docs、figures、output、src/calibration、
> README）已按用户要求删除，本修订全部重新生成。

### 实验 0 — 清理与基线设定
- 删除上一轮产物；重写 `AGENTS.md`、`prompt.md`；更新 `src/setup.R`
  （新进展率 + 常数传播）。
- 产物：`docs/calibration/preflight.md`、`DECISIONS.md`。

### 实验 1 — 12 参数直接多重起点（run1_v1_12p）
- 尝试：9 个起点（0 起点 + 5 随机扰动 + 3 知情起点），maxit 3000。
- 结果：最优 NLL 31.33，**未达接受标准**（<20 患病率 0.181 vs 0.111；
  40-49 人口 2597 vs 1841）。60+ 本身已接近（0.379 vs 0.3545）。
- 诊断：Nelder-Mead 未能到达更好的盆地；知情起点只解人口平衡，患病率
  代价被 NM 用人口拟合换取。
- 产物：`output/calibration/run1_v1_12p/`。

### 实验 2 — 热启动发现（run7 派生起点）
- 尝试：把上一轮收敛解（run7）的 m≈0.57 合并进接触行缩放、去掉
  eta_s6，得到 12 参数起点；直接评估 NLL 29.2，再用 NM（maxit 6000）
  打磨至 **NLL 22.47**，**全部接受标准通过**
  （RMSE 0.008、最大患病率误差 0.019、MAPE 0.046、最大 APE 0.132、
  均衡通过）。
- 结论：12 参数基线足以拟合所有目标（含 60+），无需超额死亡率。
- 产物：把该起点写入 `src/calibration/calibrate_nm.R`
  （`warm_run7_derived_12p`、`warm_12p_22p5`）。

### 实验 3 — 正式多重起点（run2_v1_12p_warm，最终运行）
- 尝试：11 个起点（0 + 5 随机扰动 + 3 知情 + 2 热启动），maxit 3000，
  seed 101；含 Laplace 区间与 ggplot2 图件。
- 结果：见 `docs/calibration/final_report.md`。
- 产物：`output/calibration/run2_v1_12p_warm/`、`figures/`。

### 实验 4 — beta_scale>1 约束探针（软规则）
- 尝试：把 beta 缩放参数化为 `beta_scale = 1 + exp(phi)`（保证 >1），
  从最终解重新优化，比较指标。
- 结果：约束后 NLL 从 22.44 升至 184.06，人口 MAPE 0.91（目标 0.10），
  完全不可行；1-2 组（<20、20-29）无法在 >1 约束下同时满足人口与患病率
  目标（监狱人口小、所需净流入低于 CNB 基准 beta）；按用户允许放弃该组
  约束并记录理由。
- 产物：`output/calibration/probe_beta_constraint_run2/`、
  `docs/calibration/DECISIONS.md`。

### 实验 5 —（未启用）超额死亡率应急
- 条件：仅当 60+ 在 12 参数下无法拟合时才启用（增加 eta_s[6] 参数）。
- 现状：12 参数已满足标准，故**未启用**；如需启用，见 AGENTS.md 的
  "Excess-mortality contingency" 与配套理由文件要求。

### 实验 6 — 贝叶斯后验（NPE 首选 + MCMC 主结果，prompt_mcmc.md 修订版）
- 尝试 A（NPE 首选）：sbi 0.27 `NPE_C`（MAF），3 轮序贯 SNPE（60,000
  先验 + 2 x 20,000 提议仿真），双种子（2026/2027），各 60,000 后验
  抽样；SBC 覆盖率 0.94-1.00。先验：log 接触 ~ Normal(0,2^2)、log beta
  ~ Student-t(3,0,2)。
- 结果 A（NPE 未采纳）：强识别方向过度自信（<20 患病率 CrI
  [0.106, 0.116] vs MCMC [0.055, 0.174]），弱识别方向跨种子不稳定
  （theta6 2.5% 分位 -4.21 vs -19.63；60+ 人口预测中位数 521 vs 491）。
  按回退规则改用传统方法为主。
- 尝试 B（MCMC 主结果）：adaptive Metropolis（Haario 2001），建议
  协方差取 NPE 后验协方差；3 链（NPE 中位数 / NM 最优 / 知情 NM 解），
  热重启至 40,000 迭代（burn-in 5000、thin 5）。
- 结果 B：**严格标准全部满足**——R-hat 1.0008-1.0096（要求
  [0.99, 1.01]），ESS 跨链 1,374-2,457（要求 > 400），接受率 0.238；
  12 个目标量 MCMC 95% CrI 全部与观测区间及 Laplace 区间重叠；60+
  后验预测患病率 0.379（[0.336, 0.424]）、人口 448（[381, 530]）。
- 产物：`output/calibration/npe_bayes/`、
  `docs/calibration/bayes_methodology.md`。

## 4. 文件说明（各文件是什么、哪次实验的产物）

| 文件/目录 | 说明 | 产物来源 |
|---|---|---|
| `src/setup.R` | 参数（进展率、死亡率、接触矩阵、beta、lambda 等）、初始条件 | 基线（实验 0） |
| `src/sim.cpp` | ODE 求解器（3 层、360 分室、/10.0 老化） | 基线（上一轮已定，本修订未改） |
| `src/calibration/targets.R` | 目标数据与 x_prev 重建 | 实验 0 |
| `src/calibration/model_metrics.R` | J 汇总与拟合指标 | 实验 0 |
| `src/calibration/equilibrium.R` | T vs T-5 均衡门槛 | 实验 0 |
| `src/calibration/likelihood.R` | 12 参数化与 NLL | 实验 0 |
| `src/calibration/calibrate_nm.R` | 多重起点 NM + 知情/热启动 | 实验 0 + 实验 2 热启动 |
| `src/calibration/laplace.R` | Laplace 近似与区间传播 | 实验 0 |
| `src/calibration/run_calibration.R` | 端到端运行器（元数据、多重起点、输出） | 实验 0 |
| `src/calibration/plot_results.R` | 全部 ggplot2 图件 | 实验 0 |
| `src/calibration/probe_beta_constraint.R` | beta>1 探针 | 实验 4 |
| `src/calibration/mcmc.R` | 先验、log 后验、adaptive Metropolis 采样器（含热重启） | 实验 6 |
| `src/calibration/run_npe.R` | NPE+MCMC 全流程运行器（数据/训练/验证/预测/图件） | 实验 6 |
| `src/calibration/npe_train.py` | Python NPE_C 训练/采样/SBC（多轮、多种子） | 实验 6 |
| `src/calibration/plot_mcmc.R` | 贝叶斯图件（trace/density/预测/区间/SBC） | 实验 6 |
| `output/calibration/run1_v1_12p/` | 12 参数直接多重起点（未达标） | 实验 1 |
| `output/calibration/run2_v1_12p_warm/` | 最终校准（NLL 22.44，全部达标） | 实验 3 |
| `output/calibration/probe_beta_constraint_*/` | beta>1 约束比较 | 实验 4 |
| `output/calibration/npe_bayes/` | 贝叶斯后验（NPE 敏感性 + MCMC 主结果、诊断、可信区间、图件） | 实验 6 |
| `output/calibration/npe_test/` | NPE 管道冒烟测试（审计留档） | 实验 6 |
| `output/calibration/smoke_*/` | 冒烟测试（管道验证） | 实验 0/1 |
| `docs/calibration/preflight.md` | 运行前 Git/哈希记录 | 实验 0 |
| `docs/calibration/model_audit.md` | 模型审计 | 实验 0 |
| `docs/calibration/mortality_review.md` | 2015 死亡率与 SMR 证据 | 实验 0 |
| `docs/calibration/natural_history_review.md` | 进展率证据与取舍 | 实验 0 |
| `docs/calibration/likelihood.md` | 统计模型规范 | 实验 0 |
| `docs/calibration/bayes_methodology.md` | 贝叶斯理论、过程与结果解读 | 实验 6 |
| `docs/calibration/DECISIONS.md` | 只追加决策日志 | 全程 |
| `docs/calibration/final_report.md` | 最终校准报告 | 实验 3 |
| `docs/calibration/PROJECT_OVERVIEW.md` | 本文件 | 全程 |
| `figures/` | 最终 ggplot2 图件 | 实验 3 |
| `README.md` / `README.zh-CN.md` | 项目介绍（英/中） | 实验 3 |
| `AGENTS.md` / `prompt.md` | 工作说明（本修订版） | 实验 0 |
| `Model schematic.pptx` | 模型示意图（受保护） | 基线，未修改 |

## 5. 复现命令

```bash
# 最终校准
Rscript src/calibration/run_calibration.R \
  --run-id my_run --seed 101 --maxit 3000 --n-starts 6 \
  --t-start -10 --t-end 55 --target-time 45 --target-mode sero \
  --out-dir output/calibration

# 图件
Rscript src/calibration/plot_results.R \
  --fit output/calibration/<run_id>/fit.rds --out-dir figures

# 3. 贝叶斯后验（NPE 3 轮双种子 + MCMC 严格验证）
Rscript src/calibration/run_npe.R --step all \
  --root . --fit output/calibration/<run_id>/fit.rds \
  --out-dir output/calibration/npe_bayes \
  --n-sims 60000 --n-cores 6 --seed 2026 \
  --n-draws 60000 --n-proposal 20000 --n-rounds 3 \
  --n-iter-mcmc 20000 --burnin-mcmc 5000 --mcmc-max-iter 400000
```

## 6. 最终结果（run2_v1_12p_warm）

- NLL 22.44（患病率 20.82 + 人口 1.62）；患病率 RMSE 0.0084、最大误差
  0.0195；人口 MAPE 0.047、最大 APE 0.126；均衡通过（最大患病率变化
  0.00125）。
- Laplace 95% 区间（1000 次可行抽样，有效维数 11/12）：六个年龄组的
  患病率与人口区间均与观测区间重叠（60+ 患病率 [0.323, 0.417] vs
  [0.309, 0.402]；人口 [407, 558] vs [336, 498]）。
- 拟合参数：接触行缩放 (5.51, 0.082, 0.053, 0.465, 7.45, 0.15)；
  beta 缩放 (0.169, 0.920, 1.263, 5.544, 14.85, 109.6)。
- 多重起点：两个热启动均收敛到同一盆地（22.48 / 22.44）。

## 7. 结论与局限

- 12 参数、常数传播、无超额死亡率的模型即可满足全部接受标准，60+
  组患病率与人口均在标准内。
- 贝叶斯（MCMC 主结果）结论一致：严格收敛（R-hat 1.0008-1.0096，
  ESS>400），全部 12 个目标量的 95% 可信区间与观测区间及 Laplace 区间
  重叠。NPE 已尝试（3 轮、双种子、SBC 校准良好）但因过度自信与弱识别
  方向不稳定而按规则回退为敏感性方法（详见 bayes_methodology.md）。
- 局限：60+ 患病率下降主要依赖大量未感染老年流入（beta6≈112/年）；
  <20 与 20-29 的 beta 缩放 <1；均衡为拟稳态（2015-2030 仍有小幅漂移）；
  进展率取文献上限（GT3 加权后有效率更高，属刻意保守上限）；
  contact6/beta6 为弱识别方向（MCMC 区间宽，Laplace 截断该方向）。
