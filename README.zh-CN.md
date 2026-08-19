# 新加坡注射吸毒人群（PWID）HCV 传播模型

一个针对新加坡注射吸毒人群（PWID）的年龄结构丙型肝炎（HCV）分室模型，
以监狱（羁押）PWID 群体的年龄别 HCV 血清阳性率和人口规模为校准目标。

## 1. 项目结构

```text
.
├── src/
│   ├── setup.R                # 参数、初始条件、辅助函数
│   ├── sim.cpp                # 编译的 ODE 求解器（4 层结构模型）
│   └── calibration/
│       ├── targets.R          # 校准目标数据与二项计数重建
│       ├── model_metrics.R    # J 层汇总提取与拟合指标
│       ├── equilibrium.R      # T 与 T-5 稳定性门槛
│       ├── likelihood.R       # 12 参数化与 NLL 目标函数
│       ├── calibrate_nm.R     # 多重起点 Nelder-Mead + 知情/热启动
│       ├── laplace.R          # Laplace 近似（Hessian、MC 区间）
│       ├── run_calibration.R  # 端到端校准主程序
│       ├── run_analysis.R     # 50 年敏感性投影（后验 CrI）
│       ├── scenarios.csv      # 敏感性场景清单（单一数据源）
│       ├── plot_results.R     # ggplot2 出版级图件（所有图）
│       ├── probe_beta_constraint.R  # beta_scale>1 可行性探针
│       └── （分析探针）        # 变体/敏感性探针
├── docs/calibration/
│   ├── preflight.md           # 运行前记录（Git SHA、哈希）
│   ├── model_audit.md         # 模型审计与失败诊断
│   ├── mortality_review.md    # 2015 年死亡率与 PWID SMR 证据
│   ├── natural_history_review.md  # 肝脏疾病进展率证据综述
│   ├── likelihood.md          # 统计模型规范
│   ├── DECISIONS.md           # 只追加的决策日志
│   ├── final_report.md        # 最终校准报告
│   └── PROJECT_OVERVIEW.md    # 全流程实验说明（试过什么、文件来源）
├── output/calibration/        # 每次运行的输出（run_config、solutions、predictions 等）
├── figures/                   # 最终 ggplot2 图件
├── README.md / README.zh-CN.md
├── AGENTS.md / prompt.md      # 面向 agent 的工作说明
└── Model schematic.pptx       # 模型示意图（受保护，不可修改）
```

## 2. 模型结构

```text
4 层结构： D（从未被捕的活跃 PWID）、J（目前在押）、
           F（曾被逮捕的活跃 PWID）、X（前 PWID，不再注射）
4 个肝脏阶段：NC、CC、DC、HCC
4 个 HCV 状态：u（易感/治愈后）、a（急性）、c（慢性）、t（治疗中）
6 个年龄组：<20、20-29、30-39、40-49、50-59、60+
共 384 个分室
```

- **传播**仅发生在活跃 PWID 中（D 与 F 注射，J 与 X 不注射）：频率依赖
  的感染力 `gamma_i = q * sum_j C_contact(i,j) * infectious_active_j /
  active_j`，
  传播乘子为**常数**（m = 1；历史传播水平由拟合的接触矩阵行缩放吸收）。
- **老化**：10 岁组距，每年 `y/10` 进入下一组（60+ 为开放组）。
- **逮捕/释放**：D→J（首次逮捕）速率 `lambda1[i]`；J→F 概率 `pi_recid`、
  J→X 概率 `1-pi_recid`，释放速率 `lambda2[i]`；F→J（再次逮捕）速率
  `lambda3[i]`。
- **死亡**：`mu[i] * omega`（omega = 14.68，PWID SMR）；失代偿性肝硬化
  附加 `mu_DC = 0.130/年`，HCC 附加 `mu_HCC = 0.430/年`。
- **进展率**（原始、基于文献）：`p_NC_CC = 0.027`（Thein 2008）；
  `p_CC_DC = 0.039`、`p_CC_HCC = 0.014`、`p_DC_HCC = 0.014`
  （Lim et al. 2018）。基因 3 型相对风险（Kanwal et al. 2014）保留。

## 3. 数据来源

- **目标数据**：2014-2016 樟宜监狱普筛的年龄别抗-HCV 血清阳性率与监狱
  人口（来源论文：Park et al.，medRxiv 10.1101/2025.10.24.25338708）。
- **背景死亡率**：新加坡统计局（SingStat）2015 年年龄别死亡率，映射到
  六个年龄组（见 `mortality_review.md`）。
- **PWID 超额死亡**：Mathers et al. 2013 汇总 SMR 14.68（PMID 23554523）。
- **进展率**：Thein et al. 2008（NC→CC）；Lim et al. 2018
  （CC→DC、CC→HCC、DC→HCC）。Alazawi et al. 2010 与 Rivera-Irigoin et al.
  2006 用于早期更高速率变体（见 PROJECT_OVERVIEW.md）。

## 4. 校准方法

- **12 个拟合参数**：6 个接触矩阵行缩放、6 个 beta 流入缩放
  （均为 `exp(theta)`）。不使用传播乘子参数，也不使用超额死亡率参数
  （60+ 组在 12 参数下已满足接受标准）。
- 目标函数：患病率二项 NLL + 人口对数正态 NLL（sigma_pop = 0.10）；
  均衡门槛（T 与 T-5）以罚项形式加入。
- 优化器：确定性多重起点 Nelder-Mead（`optim`，method = "Nelder-Mead"）。
  不使用 HMC。
- 不确定性：Laplace 近似（有限差分 Hessian、特征值截断 1e-4 的广义逆、
  1000 次按均衡可行性过滤的 Monte Carlo 抽样），逐年龄组与观测区间比较
  95% 区间。
- 接受标准：患病率 RMSE <= 0.02、最大患病率误差 <= 0.03、人口 MAPE
  <= 0.10、最大 APE <= 0.20、均衡通过。
- `beta_scale > 1` 为软性目标：3-6 组满足；<20 与 20-29 组因监狱人口
  规模小、所需净流入低于 CNB 基准 `beta` 而保持在 1 以下
  （见 DECISIONS.md）。

### 4.1 贝叶斯推断（MCMC）

在 NM 拟合之上进行贝叶斯后验采样。按修订版 `prompt_mcmc.md`，**先尝试
NPE/SNPE**（sbi `NPE_C`，3 轮序贯、双种子、SBC 校准），但未采纳：其在
强识别方向过度自信、弱识别方向跨种子不稳定。**传统 MCMC（adaptive
Metropolis-Hastings）为最终主方法**，收敛从严：R-hat ∈ [0.995, 1.005]
（实际 0.9997-1.0041）、合并 ESS > 1000（实际 1,310-1,958）、3 链 x
30,000 迭代、burn-in 5,000、thin 20。先验：log 接触缩放 ~ Normal(0, 2^2)；
log beta 缩放 ~
Student-t(3, 0, 2)。逐年龄组后验预测 95% 可信区间与 Laplace 及观测区间
重叠。详见 `prompt_mcmc.md`。

## 5. 如何复现

```bash
# 1. 完整校准（默认参数即可复现最终结果）
Rscript src/calibration/run_calibration.R \
  --run-id my_run --seed 101 --maxit 3000 --n-starts 6 \
  --t-start 0 --t-end 150 --target-time 45 --target-mode sero \
  --out-dir output/calibration

# 2. 出版级图件（ggplot2）
Rscript src/calibration/plot_results.R \
  --fit output/calibration/<run_id>/fit.rds --out-dir figures

# 3. 贝叶斯后验（NPE 3 轮 + MCMC 严格验证）
Rscript src/calibration/run_npe.R --step all \
  --root . --fit output/calibration/<run_id>/fit.rds \
  --out-dir output/calibration/npe_bayes \
  --n-sims 60000 --n-cores 6 --seed 2026 \
  --n-draws 60000 --n-proposal 20000 --n-rounds 3 \
  --n-iter-mcmc 20000 --burnin-mcmc 5000 --mcmc-max-iter 400000
```

每次运行输出：`run_config.csv`、`targets.csv`、`initial_values.csv`、
`optimization_history.csv`、`solutions.csv`、`predictions.csv`、
`residuals.csv`、`equilibrium.csv`、`laplace_intervals.csv`、
`laplace_diagnostics.csv`、`sessionInfo.txt`、`fit.rds` 及 ggplot2 图件。

### 5.1 敏感性分析（2017 起 50 年投影）

场景策略统一存放在一个 CSV（`src/calibration/scenarios.csv`）：各肝脏
阶段的治疗速率、GT3 比例、SVR/RBV 模式、可治疗层（D/J/F/X 标志）与
最小可治疗年龄组。每行一个策略；输出为 300 条后验抽样下的中位数与 95%
可信区间：

```bash
Rscript src/calibration/run_analysis.R \
  --root . --fit output/calibration/run1_4strata/fit.rds \
  --posterior output/calibration/npe_bayes/posterior_samples_mcmc.csv \
  --out-dir output/analysis --n-draws 300 --n-cores 4
```

输出：`scenario_summary.csv`、`scenario_key_years.csv`、
`fig_hcv_trajectories.png`、`fig_dc_hcc_trajectories.png`。

## 6. 结果

最终运行：`output/calibration/run1_4strata/`。

12 参数、无超额死亡率，满足全部接受标准。贝叶斯输出见
`output/calibration/npe_bayes/`；50 年敏感性分析见
`output/analysis/`。详见 `docs/calibration/final_report.md` 和
`docs/calibration/analysis_report.md`。
