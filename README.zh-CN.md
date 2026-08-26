# 新加坡 PWID 丙型肝炎传播模型 —— NumPyro/JAX 版

本仓库将 legacy R/C++ 版本的年龄结构 HCV 传播模型迁移到 **NumPyro / JAX**，
并使用 **NUTS** 进行贝叶斯推断。完整实验记录见 [memo.md](memo.md)，设计理念
见 [docs/design.md](docs/design.md)，架构说明见 [docs/architecture.md](docs/architecture.md)。

## 模型

- 4 层（D/J/F/X）× 4 肝病期（NC/CC/DC/HCC）× 4 HCV 状态（u/a/c/t）
  × 6 年龄组 = **384 个隔室**；
- 连续时间 ODE，固定步长 RK4（`dt = 1/365` 年），10 岁年龄带推进，
  模拟至 t = 150（目标时刻 t = 45，对应 2015）；
- 活跃吸毒者（D、F）间频率依赖感染压；
- **12 个拟合 log 参数**：6 个 contact 矩阵行缩放 + 6 个 beta 流入缩放，
  校准至新加坡监狱（J 层）HCV 血清学患病率与在押人口。

## 目录

```text
config/    固定参数 JSON（模型/模拟/目标/先验）
src/       JAX 模拟器（sim.cpp 移植）+ NumPyro 模型 + 推断工具
scripts/   run_nuts.py —— NUTS 推断流水线
tests/     单元测试（含与 R 参考数值等价）+ 集成冒烟
legacy/    R/C++ 时代结果与源码的只读存档
docs/      中文文档（需求/架构/测试计划/设计理念）
results/   推断输出（后验、诊断、PPC、图）
```

## 快速开始

```bash
python3.13 -m venv .venv && source .venv/bin/activate
pip install -r requirements-dev.txt

pytest tests/unit -v          # 含与 R 参考的数值一致性测试
pytest tests/integration -v   # 短时域 NUTS/PPC 冒烟

python scripts/run_nuts.py --num-warmup 500 --num-samples 1000 \
  --num-chains 4 --dt-days 7 --out-dir results/numpyro_nuts
```

并行链需设置 `XLA_FLAGS="--xla_force_host_platform_device_count=4"`。

## 核心结果（详见 memo.md 与 results/numpyro_nuts/）

- JAX 移植在 legacy NM 最优处与 R/C++ 参考偏差 ~1e-12；
- NUTS（白化参数化、NM 点估计初值）得到与 legacy 同标准（R-hat/ESS）收敛
  的后验；
- 后验预测覆盖观测目标；与 legacy 后验的逐参数比较见
  `results/numpyro_nuts/compare_legacy.csv`。

## 来源

由 legacy 分支 `dev_10yrs_age_interval` 迁移而来；旧结果与源码存档于
`legacy/`。
