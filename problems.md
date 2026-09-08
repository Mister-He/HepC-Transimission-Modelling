# preprocess.R 数据预处理问题清单

本文档记录对 data/sps.csv、src/preprocess.R 和当前
output/cleaned_data.rds 的审计结果。

当前数据审计结果：

- 原始数据：44,352 条记录，26,519 个唯一 IDENTIFER；
- 清洗后数据：43,023 条记录；减少的 1,329 条主要来自 incarceration 合并；
- ADM_DT、YEAR_OF_BIRTH 和 IDENTIFER 没有缺失；
- PERM_REL_DATE 缺失 5,725 条；
- 没有发现 PERM_REL_DATE < ADM_DT 的记录；
- 2 个 ID 的缺失释放日期不是该 ID 的最后一条记录，当前代码会删除其后续记录；
- 有 2,418 条释放日期晚于 2017-12-31，最晚为 2018-09-18。

## 1. 高优先级问题

### 1.1 观察截止日期没有统一处理

相关代码：[src/preprocess.R:97-107](src/preprocess.R#L97-L107)

当前代码只把 end == NA 当作删失：

    event = ifelse(is.na(end), 0, 1)

但是数据中有 2,418 条 end > 2017-12-31。如果研究观察窗口截止于
2017-12-31，这些记录在截止日时仍然是 right-censored，不应该被标记为
已观察到的 release event。

当前处理会造成：

- 使用 2018 年释放信息预测 2017 年以前的 release rate；
- t_injail 超过研究观察窗口；
- release/rearrest likelihood 使用了 cutoff 之后的信息；
- rearrest_dat 对最后一条记录的处理与 release 数据不一致。

建议保留原始日期，并建立观察窗口内的日期：

    cutoff <- as.Date("2017-12-31")

    data <- data %>%
      mutate(
        end_raw = PERM_REL_DATE,
        end_obs = case_when(
          is.na(end_raw) ~ cutoff,
          end_raw > cutoff ~ cutoff,
          TRUE ~ end_raw
        ),
        event = as.integer(!is.na(end_raw) & end_raw <= cutoff),
        t_injail = as.numeric(end_obs - start) + 1
      )

这里假设 PERM_REL_DATE 缺失代表截至 cutoff 尚未释放。如果缺失实际表示
“释放日期未知”，则不能自动填成 cutoff，而应作为未知记录单独处理。

### 1.2 快照筛选排除了所有缺失释放日期的在押人员

相关代码：[src/preprocess.R:140-143](src/preprocess.R#L140-L143)

当前代码为：

    filter(start <= time_end & end >= time_end)

在 R 中，NA >= time_end 的结果是 NA，而 dplyr::filter() 会丢弃该行。
因此，所有 end == NA 的人员都没有被计入快照。

四个日期的记录数如下：

| 快照日期 | 当前代码 | 被排除的缺失 end | 包含删失记录后 |
|---|---:|---:|---:|
| 2014-12-31 | 7,834 | 1,821 | 9,655 |
| 2015-12-31 | 6,688 | 3,090 | 9,778 |
| 2016-12-31 | 5,096 | 4,398 | 9,494 |
| 2017-12-31 | 2,418 | 5,725 | 8,143 |

因此 7,834 -> 6,688 -> 5,096 -> 2,418 的明显下降主要是筛选逻辑造成的，
不是 cleaned_data.rds 在逐年删除人员。

建议使用：

    filter(start <= time_end, is.na(end_raw) | end_raw >= time_end)

或者在已创建 end_obs 后使用：

    filter(start <= time_end, end_obs >= time_end)

### 1.3 缺失释放日期导致的记录截断会静默丢失数据

相关代码：[src/preprocess.R:122-135](src/preprocess.R#L122-L135)

当前代码没有“把删失记录移到最后”，而是保留缺失日期之前的记录，并删除
后续记录：

    if(any(is.na(temp$end))){
        loc = which(is.na(temp$end))
        temp = temp[c(1:loc), ]
    }

当前数据中有 2 个 ID 触发了该问题，共有 2 条后续记录被删除。

建议先生成异常表并停止处理：

    bad_censor <- data %>%
      arrange(IDENTIFER, start) %>%
      group_by(IDENTIFER) %>%
      summarise(
        bad_censor_order = any(is.na(end) & seq_along(end) < n()),
        .groups = "drop"
      ) %>%
      filter(bad_censor_order)

    if (nrow(bad_censor) > 0) {
      stop("A missing release date is not the last record for some individuals.")
    }

不能用自动截断替代对异常记录的调查。

## 2. incarceration 合并的问题

### 2.1 max(PERM_REL_DATE) 遇到缺失值可能错误地产生删失

相关代码：[src/preprocess.R:56-72](src/preprocess.R#L56-L72)

    PERM_REL_DATE = max(PERM_REL_DATE)

R 的 max() 默认 na.rm = FALSE。如果同一个合并组中同时存在有效释放日期和
缺失释放日期，结果通常会是 NA，整个合并后的 incarceration 会被当作删失。

更重要的是，缺失释放日期与后续入监记录之间不能简单判定为同一段或不同段，
应先检查数据逻辑，再决定如何合并。

建议：

1. 合并前检查同一 ID 的记录是否存在缺失 end 后又有后续记录；
2. 对异常 ID 单独核查；
3. 对正常的合并组使用明确的 end_obs 规则；
4. 不要依赖 max() 自动解决缺失日期。

### 2.2 合并时使用 first() 可能丢失 HCV 信息

当前代码使用：

    Screening_result = first(Screening_result)
    HepC_database = first(HepC_database)
    Date_diagnosed = first(Date_diagnosed)

如果第一条记录是 Unknown，而同一合并组后面的记录是 Positive，后面的阳性
信息会被丢弃。

建议使用显式聚合规则，例如：

    Screening_result = case_when(
      any(Screening_result == "Positive") ~ "Positive",
      any(Screening_result == "Negative") ~ "Negative",
      TRUE ~ "Not Known"
    )

Date_diagnosed 应保留最早的有效诊断日期；PRE EXISTING 应作为独立状态保留，
而不是仅依赖第一条记录。

### 2.3 合并规则需要明确日期边界

当前逻辑是：

    gap = as.numeric(ADM_DT - lag(PERM_REL_DATE))
    new_incarceration = ifelse(is.na(gap) | gap > 0, 1, 0)

这会把 ADM_DT == previous PERM_REL_DATE 的记录合并，把下一天入监的记录视为
新的 incarceration。这个规则本身可以合理，但必须明确采用的是闭区间日期还是
连续时间区间，并与 t_injail = end - start + 1 保持一致。

建议在合并前后记录：

- 原始记录数；
- 合并后记录数；
- 合并组数量；
- 每个 ID 是否仍有重叠区间；
- 是否存在非法日期顺序。

## 3. HCV 字段处理的问题

### 3.1 不能依赖字符串前导空格

相关代码：[src/preprocess.R:22-31](src/preprocess.R#L22-L31)

    HepC_database == " Hep C Positive"

该写法依赖原始数据中恰好有一个前导空格。建议先统一清理：

    HepC_database = trimws(HepC_database)

### 3.2 No Hep C Record 不应自动当作 Negative

原始 Hep C Database 包含：

    NIL              40,046
    Hep C Positive    4,210
    No Hep C Record      96

当前代码把除 Hep C Positive 之外的所有值都变成 Negative，因此
No Hep C Record 也会被当作阴性。

除非数据字典明确说明 No Hep C Record 等价于检测阴性，否则建议：

    HepC_database = case_when(
      HepC_database == "Hep C Positive" ~ "Positive",
      HepC_database == "NIL" ~ "Negative",
      HepC_database == "No Hep C Record" ~ "Not Known",
      TRUE ~ "Not Known"
    )

### 3.3 screening 空值应先标准化

screening 字段有 33,062 条空格值 " "，而不是标准的 NA。当前
TRUE ~ "Not Known" 恰好能处理它们，但这种处理不够显式。

建议：

    Screening_result <- trimws(Screening_result)
    Screening_result[Screening_result == ""] <- NA_character_

然后再进行 case_when()。

### 3.4 诊断日期的 90 天规则是强假设

相关代码：[src/preprocess.R:75-95](src/preprocess.R#L75-L95)

当前规则把入监前后 90 天内的诊断都视为入监时已感染：

    abs(as.numeric(ADM_DT - as.Date(Date_diagnosed, format = "%m/%d/%Y"))) <= 90

原始数据中只有 99 条可解析诊断日期，其中 9 条早于入监日期。该规则需要在
分析报告中明确说明，因为它：

- 同时把入监前和入监后 90 天纳入；
- 对入监后超过 90 天才诊断的人员标记为 Negative；
- 没有使用实际筛查日期或释放日期判断诊断是否发生在 incarceration 内。

建议将诊断日期先解析成独立的 Date_diagnosed_date 字段，并分别标记入监前、
入监期间、释放后和无诊断日期，再根据研究问题选择 90 天规则。

## 4. 快照和 target 构造的问题

### 4.1 注释说“2014 年以后”，但代码没有过滤

相关代码：[src/preprocess.R:33-38](src/preprocess.R#L33-L38)

当前数据中有 27,863 条记录的入监日期早于 2014 年，只有 15,993 条记录是
2014 年及以后入监。

这不一定是错误：如果需要完整历史来识别既往入监和再逮捕，早期记录应该保留。
但当前注释与实际代码不一致，必须明确以下两种分析对象之一：

1. 保留全部历史，用于构造完整的 incarceration/re-arrest history；
2. 只分析 2014 年以后入监者，并明确过滤 ADM_DT >= 2014-01-01。

建议不要只修改注释，而要在代码中显式定义 cohort。

### 4.2 使用 slice(n()) 选择人员记录不够明确

当前 prevalence target 是：

    filter(Screening_result != "Not Known") %>%
    group_by(IDENTIFER) %>%
    slice(n())

这会先删除 Unknown 记录，再选择最后一条记录。因此，如果某人的最后一次
incarceration 没有 screening result，代码会回退到更早的一条已知记录，而不是
严格选择最后一次或筛查窗口内的记录。

建议：

- 先定义明确的 screening/index record；
- 再按日期、screening 时间或目标窗口选择一条记录；
- 对同一天多条记录规定 tie-breaking 规则；
- 使用 n_distinct(IDENTIFER) 检查最终是否一人一条。

### 4.3 年龄计算忽略生日

相关代码：[src/preprocess.R:109-120](src/preprocess.R#L109-L120)

当前年龄为：

    adm_age = year(start) - YEAR_OF_BIRTH

这不是精确周岁。生日尚未到时会高估 1 岁，生日已过时才应加 1。虽然大部分
人员不会因此跨越十岁边界，但 19/20、29/30、39/40 等边界人员可能被分到
错误年龄组。

如果只有出生年份而没有出生日期，应明确这是 calendar-year age，而不是精确年龄；
如果有出生日期，则应使用完整日期计算年龄。

### 4.4 人口 target 应明确是记录数还是人数

当前代码最终使用：

    group_by(c_age_grp) %>% reframe(total = n())

在当前四个快照中，清洗后的 active record 数与 unique ID 数基本一致，但代码
本身没有保证这一点。如果未来仍有重叠 incarceration，n() 会把一个人计数多次。

建议使用：

    group_by(c_age_grp) %>%
    summarise(total = n_distinct(IDENTIFER), .groups = "drop")

并在此之前检查每个快照中是否存在一个 ID 多条 active interval。

## 5. release/rearrest 数据构造的问题

### 5.1 release rate 使用了 cutoff 之后的 release event

release_dat 直接使用：

    t_injail = data$t_injail
    event = data$event

由于当前 event 只检查 is.na(end)，2018 年释放的 2,418 条记录被当作真实
release event。应先按 cutoff 生成 end_obs 和 event，再进行 Lexis splitting。

### 5.2 rearrest 最后一条记录的 censor 逻辑不完整

相关代码：[src/preprocess.R:294-315](src/preprocess.R#L294-L315)

当前条件是：

    if(temp$end < as.Date("2017-12-31") & !is.na(temp$end))

问题包括：

- end == 2017-12-31 被排除，因为使用了 <；
- end > 2017-12-31 没有被统一视为 cutoff censoring；
- end == NA 的最后一次入监不应自动产生释放后的观察区间；
- 若前一条记录是异常的缺失 end，start[i] - end[i-1] 可能产生错误或缺失的
  release-to-rearrest duration。

建议使用 end_obs 和 event，并只对 cutoff 前真实释放的最后一条记录建立释放后的
right-censored interval。

### 5.3 t_injail + 1 是建模假设，不是单纯清洗

相关代码：[src/preprocess.R:97-104](src/preprocess.R#L97-L104)

    t_injail = 1 + as.numeric(end - start)

这会把入监日和释放日都视为完整暴露日，并将同日入监/释放记录赋予 1 天暴露。
这可以作为离散日暴露的定义，但 release Stan model 使用的是连续时间 hazard，
因此每条记录额外加一天可能带来系统性偏差。

建议在以下两种定义中选择一种并保持一致：

- 连续时间模型：使用 end - start；
- 离散日模型：使用 end - start + 1，并在方法中明确说明。

## 6. 工作流和可复现性问题

### 6.1 脚本多次清空环境

相关代码：[src/preprocess.R:4-6](src/preprocess.R#L4-L6) 和
[src/preprocess.R:180-182](src/preprocess.R#L180-L182)

    rm(list = ls())
    gc()

这会删除调用者环境中的对象，使 preprocess.R 不适合作为模块 source，也让前面
阶段创建的对象不能安全复用。

建议把脚本拆成函数，或者使用独立环境；不要在分析脚本中全局清空环境。

### 6.2 相对路径依赖当前工作目录

例如：

    read.csv("data/sps.csv")
    saveRDS(data, "output/cleaned_data.rds")
    cmdstan_model("release_rate.stan")

从非项目根目录运行脚本时可能读写失败。建议统一定义 ROOT，并使用
file.path(ROOT, ...)。

### 6.3 缺少预处理审计日志

当前脚本没有系统记录：

- 原始行数与清洗后行数；
- 合并了多少条 incarceration；
- 删除了多少异常记录；
- 每个字段的 missing/unknown 数量；
- 每个 cutoff 快照的 record 数与 unique-person 数；
- 被截断、右删失和 cutoff 后释放的记录数。

建议在保存 RDS 前写出 preprocess_audit.csv，并在异常数量不为零时发出 warning
或停止运行。

## 7. 四个快照人口下降的最终判断

当前 cleaned_data.rds 的快照记录数为：

    2014-12-31: 7,834
    2015-12-31: 6,688
    2016-12-31: 5,096
    2017-12-31: 2,418

主要原因是快照条件排除了 end 缺失记录，而缺失记录会随时间累积：

    2014: 1,821 条
    2015: 3,090 条
    2016: 4,398 条
    2017: 5,725 条

把 right-censored records 纳入后，快照计数为：

    9,655 -> 9,778 -> 9,494 -> 8,143

因此，原先的剧烈下降主要是程序性低估。修正删失后仍然存在的下降不能直接
解释为真实监狱人口下降，因为 2017 年入监记录只有约 3,226 条，而 2014–2016
年约为 4,400 条，可能还反映了数据覆盖、入监量或记录定义变化。

## 8. 建议的修复顺序

1. 统一 2017-12-31 cutoff，并正确处理 end 缺失和 end 超过 cutoff；
2. 修正快照条件，避免排除 right-censored active records；
3. 删除静默截断逻辑，改为异常检测和人工核查；
4. 修正 event、t_injail、release rate 和 rearrest rate 的时间窗口；
5. 修正 HCV 数据的 trim、Unknown/Negative 区分和合并聚合规则；
6. 明确定义 2014 年前历史记录、screening/index record 和 cohort；
7. 检查年龄边界以及 record count 与 unique-person count；
8. 增加审计日志、固定项目根路径，并移除 rm(list = ls())。

在这些问题修复前，当前 cleaned_data.rds 可以用于初步探索，但不建议直接把
当前快照人口、release event 或 rearrest exposure 解释为严格的 2014–2017
观察数据。
