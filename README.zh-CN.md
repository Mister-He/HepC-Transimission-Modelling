# Singapore PWID HCV transmission model

这是一个包含监狱—社区流动和 PWID 传播动力学的 HCV 模型。当前代码路径
已经合并到 `scripts/`；`src/` 不再是运行时入口。

## 当前模型

- 状态为 `u`、`a`、`c`、`t`、`r`；`r` 表示 recovered。
- `a -> r` 为自发清除，`t -> r` 为 SVR12，`r -> a` 为再感染。
- `r` 保持 anti-HCV reactive，因此 seroprevalence 使用
  `(r+a+c+t)/(u+r+a+c+t)`。
- 每次 simulator 调用运行 150 年，最后十年检查 equilibrium；第 150
  年作为 2015 snapshot。
- 2016--2025 的新 drug-abuser inflow 来自
  `data/singstat_new_drug_abusers.csv`，最新年份数值延续到 2030。

## 校准

Target 全部写在 `data/calibration_targets.json`，当前只包含六个年龄组的
seroprevalence 和 population count；年龄 target 不乘 `0.74`，也不再使用
47% halfway-house/release-cohort target。

```bash
Rscript scripts/run_calibration.R --run-id r_20261013_agepop \
  --maxit 3000 --n-starts 7 --sd-perturb 0.9 --seed 20261001 --cores 8
```

生产路径：Nelder–Mead 多起点拟合给出先验中心（`prior_centre()` 把最优解移到参数边界内至少 5 个先验 SD）→ 参考 MCMC（自适应随机游走 Metropolis，4 条链，每条 50,000 次评估，warmup 10,000，从第 5,000 次开始自适应，thin 10）给出生产后验 → scenario。单轮 NPE（50,000 个 LHS 先验模拟，5 个网络的集成）与 SBC 只做交叉验证，不进入生产路径：

```bash
Rscript scripts/run_mcmc_reference.R --cores 4   # 约 4.3 小时
Rscript scripts/run_npe.R --step all --n-sims 50000 --n-sbc 1000 --n-cores 14 --python /path/to/python
Rscript scripts/compare_npe_mcmc.R --policy-cores 14
```

Python 需要 `torch` 和 `sbi`。详见 `docs/training.md` 与 `output/calibration/README.md`。

## Scenario 和测试

```bash
Rscript scripts/scenarios_2026.R --root . --draws-per-seed 2000 --cores 14 --sample-seed 20261005
Rscript scripts/run_tests.R --scope all
```

Scenario 图、CSV 和可选的 Excel 文件写入 `output/scenarios_2026/`；校准图
写入 `figures/`。post-calibration analysis 的入口是
`scripts/post_calibration_analysis.R`，按要求不会在训练时自动运行。

详细说明见 [`docs/training.md`](docs/training.md)、
[`docs/evaluation.md`](docs/evaluation.md)、
[`docs/model_card.md`](docs/model_card.md) 和 [`answer.md`](answer.md)。
