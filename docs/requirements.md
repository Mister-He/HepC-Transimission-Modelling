# 需求文档（Requirements）

## 1. 目的与范围

本仓库实现新加坡 PWID（注射吸毒人群）中丙型肝炎（HCV）传播的年龄结构
室模型，并迁移至 **NumPyro / JAX** 技术栈：

1. 忠实移植 legacy R/C++ 模型（`legacy/r_source/sim.cpp`）为 JAX 实现；
2. 以 legacy 的 Nelder-Mead 点估计为初值、legacy 先验为先验，用 **NUTS**
   重新做贝叶斯推断；
3. 做同样严格的收敛性检验、后验预测检验（PPC），并与 legacy 后验比较；
4. 所有固定参数写入 JSON 配置，便于导入、查看与修改；
5. 建立 CI/CD、单元/集成测试与中文设计文档；
6. 全流程实验记录在 `memo.md`。

## 2. 模型需求

### 2.1 模型结构

- 4 个层（strata）：`D`（从未被捕、活跃）、`J`（在押）、`F`（曾经被捕、
  活跃）、`X`（退出吸毒/不再被捕）；
- 4 个肝病期：NC、CC、DC、HCC；
- 4 个 HCV 状态：`u`（易感/治愈后）、`a`（急性）、`c`（慢性）、`t`（治疗中）；
- 6 个年龄组：<20、20-29、30-39、40-49、50-59、60+；
- 共 `4 × 4 × 4 × 6 = 384` 个隔室；
- 平铺索引：`idx(s,k,h,i) = s*96 + (k-1)*24 + h*6 + i`。

### 2.2 传播与流动

- 仅活跃吸毒层（D、F）参与注射与感染；J、X 不注射；
- 频率依赖感染压：
  `gamma_i = transMult(t) * q * sum_j C_contact(i,j) * infectious_j / active_j`；
- 逮捕/释放：`D -> J`（首捕 `lambda1`）；`J -> F`（概率 `pi_recid`）或
  `J -> X`（概率 `1 - pi_recid`），释放率 `lambda2`；`F -> J`（再捕 `lambda3`）；
- 新易感者按年龄别 `beta[i]` 进入 `D_{u,NC,i}`；
- 治疗：`tau[k]`（各肝病期年治疗启动率）× 层/年龄资格
  （`tau_stratum`、`tau_min_age`）；
- 死亡：背景死亡率 `mu[i] * omega`（PWID SMR），DC/HCC 附加
  `mu_DC`/`mu_HCC`，SVR 后乘 `psi_DC`/`psi_HCC`。

### 2.3 数值方案（与 legacy 完全一致）

- 连续时间 ODE，固定步长 RK4，`dt = 1/365` 年；
- 每步 RK4 后：负值截断为 0，然后做 10 岁年龄带推进
  （`y_change = y[i]/10 * dt`，60+ 开放组）；
- 模拟时域 `t = 0..150`，校准目标时刻 `t = 45`（对应 2015 年）；
- 平衡态门：比较 `T=150` 与 `T-5=145` 的隔室/汇总稳定性。

### 2.4 移植约束

- **JAX 兼容**：全部计算可用 `jit` / `grad` / `vmap`；
- **数值等价**：在 legacy 点估计处与 R/C++ 参考输出的偏差 < 1e-8
  （实测 ~1e-12，见 `tests/data/` 与 `tests/unit/test_simulator.py`）。

## 3. 校准与推断需求

- **拟合参数**：12 个 log 参数 —— 6 个 contact 矩阵行缩放
  `theta[1:6]`、6 个 beta 流入缩放 `theta[7:12]`；
- **目标数据**：J 层 6 个年龄组的 HCV 血清学患病率（二项计数）与在押人口；
- **似然**：患病率 `Binomial(x_i | n_i, p_hat_i)` + 人口
  `log N_hat ~ Normal(log N_obs, sigma_pop^2)`（`sigma_pop = 0.10`）；
- **先验**（来自 legacy）：contact log 缩放
  `Normal(log(anchor), 0.5^2)`；beta log 缩放
  `log(anchor) + 0.8 * Student-t(3)`；
- **推断**：NUTS（多链），初值为 legacy NM 点估计；
- **收敛标准**（与 legacy 同等严格）：R-hat ∈ [0.995, 1.005]、
  pooled ESS > 1000、无/极少发散、接受率合理；
- **后验预测检验**：目标级（p_hat、N_hat）与数据级（二项计数、人口）PPC，
  观测落在 95% CrI 内，并按平衡态门过滤；
- **与 legacy 比较**：逐参数中位数/95% 区间/重叠/2 样本 KS 检验。

## 4. 工程需求

- 固定参数 JSON 化：`config/model_parameters.json`、
  `config/simulation.json`、`config/targets.json`、`config/priors.json`；
- 目录结构按项目规范（README、docs、src、tests、scripts、.github、Dockerfile、
  docker-compose.yml）；
- CI：单元测试 + 集成测试（短时域 NUTS/PPC 冒烟）自动运行；
- 文档：需求、架构、测试计划、设计理念（中文）；
- `AGENTS.md` 依据全部需求重新生成；
- `memo.md` 记录全流程实验。

## 5. 非功能需求

- 可复现性：固定随机种子、记录依赖版本与运行配置；
- 可维护性：模块划分清晰（simulator / model / inference / config）；
- 性能：完整 NUTS 运行可在数小时内在单机 CPU 完成（并行链 + 白化参数化 +
  7 天步长近似，偏差已量化）；
- 兼容性：Python ≥ 3.11，x64 精度（`jax_enable_x64 = True`）。
