# Legacy 先验分布（12 参数）

最终 MCMC 后验（`posterior_samples_mcmc.csv`）使用的先验来自
`legacy/r_source/calibration/mcmc.R` 的 `make_priors()`。12 个拟合参数
`theta[1:12]` 定义在 log 尺度：

| 参数 | 先验 | 锚点（anchor） | 尺度 |
|---|---|---|---|
| `theta1..6`（log contact row scale） | Normal(log(anchor), 0.5^2) | `[6.3713, 0.0857, 0.0589, 0.5584, 9.5232, 0.2558]` | sd = 0.5 |
| `theta7..12`（log beta inflow scale） | log(anchor) + 0.8 * Student-t(3) | `[0.1670, 0.9443, 1.2297, 5.2385, 16.9018, 110.5976]` | scale = 0.8, df = 3 |

即：

```text
theta[1:6]  ~ Normal(log(contact_anchor), 0.5^2)
z[1:6]      ~ Student-t(df = 3)
theta[7:12] = log(beta_anchor) + 0.8 * z
```

## 说明与口径差异

- `run_config.csv` 中记录的先验字符串为
  `"log contact ~ Normal(0,2^2); log beta ~ Student-t(3,0,2)"`，这是早期
  文档化描述；**实际生效的先验**（`run_npe.R` 调用 `make_priors()` 默认参数）
  是上面以文献锚点为均值的 Normal / Student-t 混合先验。
- 设计意图：contact 锚点是模型自带的"合理猜测"接触矩阵，采用较弱先验；
  beta 锚点来自 CNB 官方新吸毒者数据，但逐年龄组是占位值，因此使用重尾
  Student-t 容忍数据要求的较大偏离（例如 60+ 组 inflow scale ~ 110）。

机器可读版本见 `prior_distribution.json`。
