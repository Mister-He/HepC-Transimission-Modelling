# 敏感性分析：J 治疗从 2019 年开始

`scripts/scenarios_2026.R --j-start-year 2019 --out-dir output/scenarios_2026_j2019`（2,000 个抽样，`--sample-seed 20261005`）的生成结果。与主分析（`output/scenarios_2026/`）的唯一差别：J 治疗（年启动概率 99.995%）从 2019 年 1 月 1 日开始；F/X 筛查仍从 2026 年开始。文件与主分析相同；图：`scenario_target_trajectories.png`（S1 Text Fig F）、`paper_fig2_trajectories.png`。

结果为 2026-11-04 重新校准后（死亡 = DC/HCC 肝病死亡）、ψ 更新后（只作用于已治愈者；ψ_DC 0.40，ψ_HCC 0.34）的重跑。

- J-only：发病下降 85.0%（达标），肝病死亡下降 40.8%。
- 所有 15 个有限情景达到发病目标；没有有限情景达到死亡目标（最高 8000/年 F25%，56.1%）；只有无容量上限同时达到两个目标（死亡 76.4%）。
