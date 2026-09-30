# Singapore PWID HCV transmission model

这是一个包含监狱—社区流动和 PWID 传播动力学的 HCV 模型。当前代码路径
已经合并到 `scripts/`；`src/` 不再是运行时入口。

## 当前模型

- 状态为 `u`、`a`、`c`、`t`、`r`；`r` 表示 recovered。
- `a -> r` 为自发清除，`t -> r` 为 SVR12，`r -> a` 为再感染。
- `r` 保持 anti-HCV reactive，因此 seroprevalence 使用
  `(r+a+c+t)/(u+r+a+c+t)`。
- 每次 simulator 调用运行 150 年，最后五年检查 equilibrium；第 150
  年作为 2015 snapshot。
- 2016--2025 的新 drug-abuser inflow 来自
  `data/singstat_new_drug_abusers.csv`，最新年份数值延续到 2030。

## 校准

Target 全部写在 `data/calibration_targets.json`，当前只包含六个年龄组的
seroprevalence 和 population count；年龄 target 不乘 `0.74`，也不再使用
47% halfway-house/release-cohort target。

```bash
Rscript scripts/run_calibration.R \
  --run-id r_YYYYMMDD_agepop --maxit 1000 --n-starts 5 --seed 20260929
```

NPE 是主要的不确定性方法，不使用 Laplace：

```bash
Rscript scripts/run_npe.R --step all \
  --fit output/calibration/r_YYYYMMDD_agepop/fit.rds \
  --out-dir output/calibration/npe_agepop \
  --n-sims 20000 --n-cores 4 --n-rounds 2 \
  --python /path/to/python
```

Python 需要 `torch` 和 `sbi`。SBC 与 posterior-spread gate 未通过时，不应把
NPE output 当作最终 posterior；可按 `problems.md` 记录后再考虑 MCMC-like
fallback。

## Scenario 和测试

```bash
Rscript scripts/scenarios_2026.R \
  --fit output/calibration/r_YYYYMMDD_agepop/fit.rds
Rscript scripts/run_tests.R --scope all
```

Scenario 图、CSV 和可选的 Excel 文件写入 `output/scenarios_2026/`；校准图
写入 `figures/`。post-calibration analysis 的入口是
`scripts/post_calibration_analysis.R`，按要求不会在训练时自动运行。

详细说明见 [`docs/training.md`](docs/training.md)、
[`docs/evaluation.md`](docs/evaluation.md)、
[`docs/model_card.md`](docs/model_card.md) 和 [`answer.md`](answer.md)。
