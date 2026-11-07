# 敏感性分析：基因型 3 比例

作者要求（2026-11-04）。主分析假设 78% 的感染为基因型 3（`setup.R` `rho`，Zhao 2019）；这里改为 10%（`rho010/`）和 50%（`rho050/`），其余与主分析相同（J 从 2026 年开始，恒定流入，同一后验的 2,000 个抽样，`--sample-seed 20261005`）。先用 `scripts/gt3_params.R` 写出 `gt3_params.rds`，再运行：

```bash
Rscript scripts/gt3_params.R --root .
Rscript scripts/scenarios_2026.R --root . --draws-per-seed 2000 --cores 14 --sample-seed 20261005 \
  --param-file output/scenarios_2026_gt3/rho010/gt3_params.rds --out-dir output/scenarios_2026_gt3/rho010
```

- `gt3_params.csv/.rds`：所用的 `rho` 和随之改变的 SVR12 复合值。由 `scripts/model/genotype.R` 的 `rho_overrides()` 计算（`scripts/tornado.R` 也用它）：
  `alpha += (rho - 0.78) × (SVR_GT3 - SVR_other)`（NC/CC、DC 无 RBV、DC 有 RBV 与 HCC）。进展率的基因型加权在 `sim.cpp` 内由 `rho` 计算。
- **未重新校准。** 在 Nelder–Mead 最优点上，`rho` 0.10 使目标负对数似然增加 0.27（0.50 时 0.06），各年龄血清阳性率上升 0.4–2.2 个百分点，首次逮捕不变（2026-11-04 临时计算，未存档）。
- 其余文件与 `output/scenarios_2026/` 相同（不含图）。

| 2030 降幅（中位数） | 主分析 78% | 50% | 10% |
| --- | --- | --- | --- |
| J-only 发病 / 肝病死亡 | 51.5% / 22.1% | 51.4% / 21.2% | 51.2% / 19.7% |
| 8000/年 F25% 肝病死亡 | 44.0% | 44.1% | 44.2% |
| 无容量上限 肝病死亡 | 71.4% | 71.3% | 71.1% |

达标情景与主分析完全相同：5 个有限情景达到发病目标，只有无容量上限达到死亡目标。S1 Text 第 12 节、Table M。
