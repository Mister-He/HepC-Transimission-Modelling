# legacy/ — 旧版（R/C++）标定与贝叶斯推断结果存档

本目录是分支 `dev_10yrs_age_interval` 上完成的最终标定/推断结果的只读存档，
供 NumPyro 迁移版（本仓库主代码）作为对照基线。

## 内容

| 路径 | 内容 |
|---|---|
| `nelder_mead_point_estimate.csv` / `.json` | 最终 Nelder-Mead 点估计（12 个 log 参数 + 反变换尺度） |
| `posterior_samples_mcmc.csv` | legacy 后验样本（3 链 × 1250，自适应 Metropolis） |
| `prior_distribution.md` / `.json` | 实际生效的先验分布 |
| `likelihood_design.md` / `.json` | 似然函数与惩罚项设计 |
| `diagnostics_mcmc.csv` | legacy MCMC 收敛诊断（R-hat/ESS） |
| `credible_intervals.csv` | NPE/MCMC/Laplace 与观测的目标摘要区间对比 |
| `sbc_summary*.csv`、`npe_config*.json` | NPE（SNPE）尝试的记录（未采用为主方法） |
| `run_config.csv` | NPE 运行配置 |
| `run1_4strata/` | 最终 NM 标定运行的全部输出（含 `fit.rds`、Laplace 近似） |
| `reports/` | 旧版最终报告（`final_report.md`、`analysis_report.md`） |
| `r_source/` | 产生上述结果的 R/C++ 源码（sim.cpp、setup.R、calibration/*.R、scripts/*） |

## 最终结果速览

- NM 点估计 NLL = 21.421，平衡态通过，prevalence RMSE ≈ 0.0051，
  population MAPE ≈ 0.036。
- MCMC 后验：R-hat ∈ [0.9997, 1.0041]，pooled ESS ∈ [1310, 1958]。
- 拟合时刻 `TARGET_TIME = 45`（模型年，对应 2015），`t_end = 150`，`dt = 1/365`。

注意：本目录为存档，不参与新的构建/测试流程。迁移版模型与推断代码位于
`src/`、`config/`、`scripts/`。
