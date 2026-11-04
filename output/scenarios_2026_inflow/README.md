# 敏感性分析：观测的 2016–2025 年新吸毒者流入

`scripts/scenarios_2026.R --inflow observed --out-dir output/scenarios_2026_inflow`（2,000 个抽样，`--sample-seed 20261005`）的生成结果，由 2026-11-03 独立审查加入（`REVIEW_REPORT.md` R-02）。与主分析（`output/scenarios_2026/`）的唯一差别：2016–2025 年进入 D 的人数按 `data/singstat_new_drug_abusers.csv` 中各年龄组相对 2015 年的比例变化，2025 年的值保持到 2030 年；2015 年平衡态和校准不变。文件与主分析相同；图 `scenario_target_trajectories.png` 和 `paper_fig2_trajectories.png`（及 .svg）由 `python scripts/figures.py . --scen-dir output/scenarios_2026_inflow --only scenario_target_trajectories paper_fig2_trajectories` 绘制。

结果为 2026-11-04 重新校准后（肝病死亡为主结局）、ψ 更新后（只作用于已治愈者；ψ_DC 0.40，ψ_HCC 0.34）的重跑。

- J-only：发病下降 56.0%（95% UI 45.1%–61.9%；主分析 51.5%），肝病死亡下降 22.2%（主分析 22.1%）。
- 发病目标中位数达标的有限情景由 5 个增至 6 个（新增 2000/年 F75%，82.9%）；死亡目标仍只有无容量上限达标（71.5%）；8000/年 F25% 死亡下降 44.4%，双目标抽样 0%（主分析 0%）。
