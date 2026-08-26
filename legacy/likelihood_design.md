# Legacy 似然函数设计（12 参数校准）

来源：`legacy/r_source/calibration/likelihood.R`、`targets.R`、
`model_metrics.R`、`equilibrium.R`。拟合目标是新加坡监狱（J 层）6 个年龄组
的 HCV 血清学患病率与监狱人口数。

## 1. 目标数据（targets.R）

| 年龄组 | 监狱人口 n | 患病率(供给) | x（二项计数，四舍五入） |
|---|---|---|---|
| <20 | 99 | 0.1118421 | 11 |
| 20-29 | 1244 | 0.1731044 | 215 |
| 30-39 | 1467 | 0.2684954 | 394 |
| 40-49 | 1841 | 0.4301165 | 792 |
| 50-59 | 1628 | 0.4821029 | 785 |
| 60+ | 409 | 0.3544304 | 145 |

## 2. 模型输出摘要（model_metrics.R）

- 模拟 150 年（`t = 0..150`，`dt = 1/365`，RK4 + 每步年龄推进 1/10/年）。
- 校准对比时刻 `TARGET_TIME = 45`（对应 2015 年）。
- `N_hat[i]` = J 层（在押）年龄组 i 总人数；`p_hat[i]` = J 层 HCV
  seroprevalence = (acute + chronic + on-treatment)/N_hat（`target_mode = "sero"`）。

## 3. 统计似然

```text
NLL(theta) = nll_prev + nll_pop

nll_prev = -sum_i Binomial(x_i | n_i, p_hat_i, log)
nll_pop  = 0.5 * sum_i ((log N_hat_i - log N_obs_i) / sigma_pop)^2
          （等价于 log N_hat ~ Normal(log N_obs, sigma_pop^2), sigma_pop = 0.10）
```

患病率在进入二项似然前被截断到 `[1e-10, 1-1e-10]`。

## 4. 惩罚项

### 4.1 平衡态门（equilibrium.R）

模拟终点 T=150 与 T-5=145 的比较，要求全 384 个隔室与 HCV/DC/HCC 总量稳定：

- `max_log_pop_ratio`（各年龄组人口对数比）<= 0.01
- `max_prev_change`（患病率绝对变化）<= 0.005
- `state_log_ratio`（HCV/DC/HCC 总量对数比）<= 0.01
- `max_comp_log_ratio`（384 隔室对数比，下限 1e-6）<= 0.02
- `total_log_ratio`（总人口对数比）<= 0.01（仅检查，不进入 pass 判定）

NM 目标函数：不通过时加入
`EQ_PENALTY_FACTOR = 1e6` × 各归一化越界量之和（log_pop/0.01 + prev/0.005 +
total/0.01 + state/0.01 + comp/0.02）。

MCMC `log_posterior`：不通过时只扣除
`1e6 * (max(log_pop_ratio/0.01,0) + max(prev_change/0.005,0) +
max(total_log_ratio/0.01,0))`。

### 4.2 合理性软惩罚（仅 NM 目标函数，系数 1e-2）

- contact 矩阵对角元素应单峰，峰值在 30-39 岁组；
- 50+/60+ 对角元素应小于 30-39 峰值的 1/2；
- beta 随年龄单调不减，且偏好 > 1（`0.01 * sum(max(1 - beta_scale,0)^2)`）。

### 4.3 参数边界（硬边界）

- `contact_scale = exp(theta[1:6]) ∈ [0.01, 1000]`
- `beta_scale = exp(theta[7:12]) ∈ [0.01, 1e6]`

越界时 NM 返回 `1e12 + 1e6 * barrier`；MCMC `log_posterior` 返回 `-Inf`。

## 5. 后验

- 方法：3 条链 × 30,000 次自适应 Metropolis（Haario et al. 2001），
  burn-in 5,000，thin 20，每条链保留 1,250 样本；
- 收敛标准：R-hat ∈ [0.995, 1.005]（实测 0.9997–1.0041），pooled ESS > 1000
  （实测 1,310–1,958）；
- 后验预测：每条样本重新模拟并过滤平衡态可行样本，得到 95% CrI。

机器可读版本见 `likelihood_design.json`。
