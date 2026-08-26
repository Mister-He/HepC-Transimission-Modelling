# AGENTS.md

本文件面向后续在此仓库工作的 agent（包括 AI agent），依据项目全部需求生成。

## 1. 项目概述

新加坡 PWID 人群丙型肝炎（HCV）传播的年龄结构室模型，已从 R/C++ 迁移到
**NumPyro / JAX**。12 个 log 参数（6 个 contact 矩阵行缩放 + 6 个 beta 流入
缩放）校准到新加坡监狱（J 层）6 个年龄组的 HCV 血清学患病率与在押人口。
推断使用 **NUTS**，初值取自 legacy Nelder-Mead 点估计，先验与 legacy 一致。

## 2. 仓库结构（规范结构）

```text
README.md / README.zh-CN.md   项目说明
AGENTS.md                     本文件
memo.md                       全流程实验记录（先读这里！）
config/                       固定参数 JSON（模型/模拟/目标/先验）
legacy/                       R/C++ 时代结果与源码存档（只读）
docs/                         中文文档（requirements/architecture/test-plan/design）
src/                          Python 源码
  simulator.py                JAX 模拟器（sim.cpp 忠实移植）
  model.py                    NumPyro 模型（先验/似然/惩罚/白化参数化）
  inference.py                诊断/PPC/与 legacy 比较工具
scripts/run_nuts.py           NUTS 推断流水线
tests/unit + tests/integration pytest 测试
results/numpyro_nuts/         本次推断输出
.github/workflows/ci.yml      CI
Dockerfile / docker-compose.yml
requirements*.txt / pyproject.toml
```

## 3. 工作流与常用命令

```bash
# 环境（Python >= 3.11，推荐 3.13）
python3.13 -m venv .venv && source .venv/bin/activate
pip install -r requirements-dev.txt

# 测试
pytest tests/unit -v
pytest tests/integration -v

# 完整推断（重计算，按需运行）
python scripts/run_nuts.py --num-warmup 500 --num-samples 1000 \
  --num-chains 4 --dt-days 7 --chain-method parallel \
  --out-dir results/numpyro_nuts

# 并行链需先设置：XLA_FLAGS="--xla_force_host_platform_device_count=4"
```

## 4. 关键约定与注意点

1. **数值等价性是硬约束**：JAX 模拟器复刻了 legacy `sim.cpp` 的每一步
   （RK4 固定步长、每步负值截断、10 岁年龄带推进），甚至复刻了 legacy
   `setup.R` 的 `y0` 布局偏移。改模拟器必须跑
   `tests/unit/test_simulator.py`（对照 `tests/data/` 的 R 参考，容差 1e-8）。
2. **legacy/ 是只读存档**：记录旧结果与源码，不参与构建/测试；不要删除或
   修改其中的记录。
3. **参数一律走 config/*.json**：不要在新代码里硬编码模型参数；修改参数
   即编辑 JSON（`_` 开头字段为注释，加载时剔除）。
4. **推断默认参数**：白化参数化（`--parameterization whiten`）、无平衡态
   惩罚（`--penalty-mode none`，后处理按平衡态门过滤）、主运行 `--dt-days 7`
   （与日步长的偏差已量化：人口 ≤0.27%、患病率 ≤0.0003）。
5. **收敛标准**（与 legacy 同等严格）：R-hat ∈ [0.995, 1.005]、
   pooled ESS > 1000、记录发散数与接受率。
6. **文档与记录**：任何新的实验/决策必须更新 `memo.md`；新增设计决策写入
   `docs/design.md`；代码改动保持单元测试绿色。

## 5. 当前状态（截至 memo.md）

- JAX 移植完成并经 R 参考验证（~1e-12）；
- NUTS 全量推断已运行并输出到 `results/numpyro_nuts/`；
- 收敛诊断、PPC、与 legacy 后验比较结果见该目录与 `memo.md`；
- 后续工作可围绕：场景/敏感性分析、更细步长验证、CI 扩充等。
