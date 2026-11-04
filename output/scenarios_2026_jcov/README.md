# J 治疗覆盖率网格（主分析，2026 年起）

`scripts/scenarios_2026.R` 的生成结果。与 `output/scenarios_2026/`（J 年启动概率 99.995%）使用同一后验、同一 2,000 个抽样（`--sample-seed 20261005`），只改 J 的年启动概率：

| 子目录 | J 年启动概率 | 命令参数 |
| --- | --- | --- |
| `j000/` | 0% | `--j-probability 0 --out-dir output/scenarios_2026_jcov/j000` |
| `j050/` | 50% | `--j-probability 0.5 --out-dir output/scenarios_2026_jcov/j050` |
| `j090/` | 90% | `--j-probability 0.9 --out-dir output/scenarios_2026_jcov/j090` |

每个子目录的文件与 `output/scenarios_2026/` 相同（不含图）。`jcov_heatmap.png/.svg`（`scripts/figures.py`）：17 个情景 × 4 个覆盖率的 2030 年发病（A）和死亡（B）降幅中位数，粗体为达标。

- 所有降幅随 J 覆盖率单调上升。
- 发病目标达标的有限情景数：0% 时 1 个（8000/年 F75%），50% 时 3 个，90% 与 99.995% 时 5 个。
- 任何覆盖率下都没有有限情景达到死亡目标（死亡 = DC/HCC 肝病死亡，2026-11-04 重跑）；无容量上限在所有覆盖率下都达到两个目标，0% 时死亡下降 65.4%，仅勉强达标（2026-11-04 ψ 更新后）。
- J 0% 的 J-only 即无治疗（降幅约 0）。
