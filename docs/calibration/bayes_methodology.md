# Bayesian inference for the HCV PWID model — methodology and results

Date: 2026-08-09. Status: final (output `output/calibration/npe_bayes/`).
This document explains the theory, the experimental process, and the
interpretation of the Bayesian posterior analysis performed on top of the
Nelder-Mead calibration (`run2_v1_12p_warm`, NLL 22.44).

## 1. 目标 / Goal

在冻结的模型结构与 12 个拟合参数上，用贝叶斯方法获得参数后验样本与 12
个目标量的 95% 可信区间，同时保留 Laplace 95% 区间。按 `prompt_mcmc.md`
（修订版）优先使用 NPE/SNPE 等现代方法，并用传统 MCMC 作为独立验证；
收敛标准从严：R-hat 必须在 [0.99, 1.01]，ESS 必须 > 400。

## 2. 理论 / Theory

### 2.1 Bayesian formulation

后验 `p(theta | x_obs) ∝ p(x_obs | theta) p(theta)`，参数为 12 个 log 尺度
参数。观测摘要 `x = (p_1..p_6, N_1..N_6)`（监狱血清阳性率与人口）由 ODE
模型确定性给出；观测不确定性来自二项/对数正态观测模型（与校准阶段一致：
Binomial seroprevalence + log-Normal population, sigma_pop = 0.10）。
均衡约束（T vs T-5）在后验预测阶段作为有效性过滤器（并进入 MCMC 的
log 后验罚项）。

### 2.2 Priors（依据模型特点与文献）

| 参数块 | 先验 | 理由 |
|---|---|---|
| log contact scales theta[1:6] | Normal(0, 2^2) | 基础接触矩阵是模型特有的标定猜测；弱信息先验，95% 尺度先验区间约 [0.02, 55] |
| log beta scales theta[7:12] | Student-t(3, 0, 2) | 基础 beta 锚定 CNB 官方新吸毒者数据但分年龄值为占位；重尾 t 先验以 1 为中位并容忍 60+ 大额流入（~110x） |

先验不取自 NM 拟合（避免双重使用数据）；NM 解只用于链起点与建议协方差
调参。变换：`z = theta/2`（先验化为 N(0,1)^6 x t(3)^6）；
`x_t = [logit(p), log(N)]` 用训练集均值/方差标准化并裁剪到 [-8, 8]。

### 2.3 NPE / SNPE（首选方法，已尝试）

- Neural Posterior Estimation：用条件归一化流（MAF）直接学习
  `p(z | x_t)`，训练数据为 `(theta~prior, x = simulator(theta))`。
- 采用 sbi 0.27 的 `NPE_C`（SNPE-C），**3 轮序贯细化**：第 1 轮从先验
  采样（60,000 条），第 2/3 轮从上一轮后验建议分布采样（各 20,000 条，
  R 端并行仿真），`force_first_round_loss=True` 保持密度比正确。
- 校验：**SBC**（simulation-based calibration，298 个留出先验点的后验
  秩统计量）、**双种子稳定性**（seed 2026/2027 独立训练）、以及与
  **MCMC 参考后验**的分位数对比。

### 2.4 Traditional MCMC（验证 / 最终主方法）

- 采样器：adaptive Metropolis-Hastings（Haario et al. 2001），建议协方差
  初值取 NPE 后验协方差（比 Laplace 协方差混合更快），`c_d = 2.38^2/d`。
- 多起点：3 条链分别从 NPE 后验中位数、NM 最优（warm_12p_22p5）与另一
  NM 解（population_informed_b6_c1）出发；支持热重启续跑，直到严格标准
  满足。
- 严格标准：R-hat ∈ [0.99, 1.01]（逐参数，另报 mpsrf）；ESS > 400
  （逐参数，跨链求和并分链报告）。

### 2.5 为何最终以 MCMC 为主（NPE 回退理由）

见第 4 节：NPE 在强识别方向给出过度自信（过窄）的区间，在弱识别方向
（contact6/beta5 等）跨种子严重不稳定；按 `prompt_mcmc.md` 第 1 条的
回退规则（"实在不好用再选择传统方法"），MCMC 后验（满足严格收敛标准）
被采纳为最终贝叶斯主结果，NPE 完整保留为敏感性分析。

## 3. 实验过程 / Process

| 阶段 | 设置 |
|---|---|
| 训练数据 | 60,000 先验抽样（59,682 条有效）；SBC 留出集 300（298 有效） |
| NPE | sbi 0.27 NPE_C (MAF)，3 轮 x 20,000 提议/轮，双种子（2026/2027），60,000 后验抽样/种子 |
| MCMC 验证 | adaptive Metropolis，3 链，块 20,000 迭代热重启；总 40,000 迭代、burn-in 5,000、thin 5 |
| 后验预测 | NPE 30,000 抽、MCMC 7,000 抽，均衡可行过滤 |
| 计算 | R 并行仿真 6 核；Python (venv /tmp/bayes-venv) torch 2.13 CPU |

所有配置与哈希见 `output/calibration/npe_bayes/run_config.csv`、
`npe_config.json`、`sessionInfo.txt`。

## 4. 结果 / Results

### 4.1 MCMC 验证（最终主方法）— 严格标准全部满足

| 指标 | 值 | 标准 |
|---|---:|---|
| R-hat（12 参数） | 1.0008 – 1.0096 | [0.99, 1.01] ✓ |
| mpsrf | 1.007 | 报告 |
| ESS（跨链求和） | 1,374 – 2,457 | > 400 ✓ |
| 每链 ESS | 332 – 855 | 报告 |
| 接受率 | 0.238 | 接近最优 0.234 |
| 总迭代 | 40,000（2 块热重启） | 记录 |

后验中位数与 NM 最优最大差 0.51（theta6，弱识别方向），其余参数一致。

### 4.2 MCMC 后验预测可信区间（主结果，与 Laplace、观测对比）

| 年龄组 | p MCMC 中位数 | p 95% CrI | p Laplace 95% CI | p 观测 95% CI | N MCMC 中位数 | N 95% CrI | N Laplace 95% CI | N 观测 95% CI |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| <20 | 0.103 | [0.055, 0.172] | [0.062, 0.187] | [0.060, 0.184] | 99.3 | [81.8, 120.9] | [81.3, 120.5] | [81.4, 120.4] |
| 20-29 | 0.173 | [0.153, 0.195] | [0.157, 0.199] | [0.153, 0.195] | 1255.2 | [1035.3, 1528.4] | [1031.1, 1536.9] | [1022.6, 1513.4] |
| 30-39 | 0.269 | [0.246, 0.292] | [0.251, 0.296] | [0.246, 0.292] | 1464.1 | [1196.3, 1778.1] | [1231.9, 1832.1] | [1205.9, 1784.6] |
| 40-49 | 0.432 | [0.409, 0.454] | [0.410, 0.455] | [0.408, 0.453] | 1923.1 | [1571.7, 2372.3] | [1584.6, 2326.0] | [1513.3, 2240.0] |
| 50-59 | 0.476 | [0.452, 0.500] | [0.457, 0.508] | [0.458, 0.506] | 1398.4 | [1191.6, 1628.9] | [1249.7, 1683.0] | [1338.2, 1980.5] |
| 60+ | 0.379 | [0.336, 0.424] | [0.323, 0.417] | [0.309, 0.402] | 448.5 | [380.9, 530.2] | [406.6, 558.0] | [336.2, 497.6] |

全部 12 个目标量的 MCMC 95% CrI 均与观测区间和 Laplace 区间重叠。

### 4.3 NPE 结果与回退依据（敏感性）

- SBC：各参数 95% 带覆盖率 0.94-1.00，90% 带 0.91-0.99，均值秩
  0.46-0.56 —— 平均校准良好。
- 跨种子稳定性：theta6（contact 行 6）2.5% 分位数 seed1 -4.21 vs seed2
  -19.63；后验预测 60+ 人口中位数 seed1 521 vs seed2 491、区间 [488,
  588] vs [310, 576] —— 弱识别方向不稳定。
- 与 MCMC 参考对比（过度自信）：<20 患病率 NPE 95% CrI [0.106, 0.116]
  远窄于 MCMC [0.055, 0.174]（甚至窄于 n=99 二项似然本身的 ±0.06 量级）；
  20-29 人口 [1232, 1263] vs [1035, 1529]。60+ 人口中位数 521 vs MCMC
  448。
- 结论：NPE（3 轮、双种子）不满足可靠性要求，按回退规则采用 MCMC 为
  主；NPE 结果保留在 `posterior_samples_npe*.csv`、
  `sbc_summary*.csv`、`fig_density.png` 中作为敏感性。

### 4.4 解读要点

- **60+ 组**：MCMC 后验预测患病率中位数 0.379（CrI [0.336, 0.424]）、
  人口 448（[381, 530]），与观测区间重叠；与 Laplace 结论一致——12
  参数模型无需超额死亡率即可容纳 60+ 目标。
- **beta_scale 后验**：<20 与 20-29 组后验中位数仍 < 1（0.169, 0.916），
  3-6 组 > 1；60+ 流入缩放中位数 ~106（后验 CrI 宽，弱识别）。
- **弱识别方向**：contact6/beta5/beta6 方向似然平坦，后验宽且由先验
  主导；Laplace 在该方向截断（有效维 11/12），MCMC 保留完整边际
  不确定性，NPE 在该方向不稳定——三类方法在此方向的行为差异已记录。
- **与 Laplace 的关系**：MCMC 是精确后验（采样误差内），Laplace 是
  高斯近似；二者在已识别方向上高度一致（60+ 区间重叠），差异集中在
  弱识别方向，报告两者并存、互为对照。

## 5. 复现 / Reproduction

```bash
# 数据生成 -> NPE 3 轮双种子 -> MCMC 严格验证 -> 预测 -> 图件
Rscript src/calibration/run_npe.R --step all \
  --root . --fit output/calibration/run2_v1_12p_warm/fit.rds \
  --out-dir output/calibration/npe_bayes \
  --n-sims 60000 --n-cores 6 --seed 2026 \
  --n-draws 60000 --n-proposal 20000 --n-rounds 3 \
  --n-iter-mcmc 20000 --burnin-mcmc 5000 --mcmc-max-iter 400000
```
