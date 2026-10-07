#!/usr/bin/env python3
"""Build docs/main.docx and docs/Supplementary Information.docx from generated results.

Usage: python3 scripts/build_paper.py [ROOT]
Needs: output/scenarios_2026/*.csv and paper_fig*.png (scripts/scenarios_2026.R, then
scripts/paper_figures.R), docs/references.bib, docs/plos.csl, pandoc and Rscript+readxl.
Every number in the text is read from the result files; edit prose here, not in Word.
"""
import csv, json, sys, os, subprocess, tempfile
ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.getcwd())
TMP = tempfile.mkdtemp()
_byseed = os.path.join(TMP, "byseed.csv")
subprocess.run(["Rscript", "-e", f'write.csv(readxl::read_xlsx("{ROOT}/output/scenarios_2026/scenario_results.xlsx",sheet="BySeed"),"{_byseed}",row.names=FALSE)'], check=True)
S = os.path.join(ROOT, "output/scenarios_2026")
def rd(p): return list(csv.DictReader(open(os.path.join(ROOT, p))))
summ = {r["scenario_id"]: r for r in rd("output/scenarios_2026/scenario_2030_summary.csv")}
traj = rd("output/scenarios_2026/scenario_trajectories.csv")
nt = rd("output/scenarios_2026/no_treatment_trajectory.csv")
ann = rd("output/scenarios_2026/scenario_annual_services_costs.csv")
treat = rd("output/scenarios_2026/scenario_treatment_summary.csv")
cost_ev = rd("output/scenarios_2026/cost_evidence.csv")
F = float

def ui(m, lo, hi, scale=1, d=1, pre="", suf=""):
    f = lambda v: f"{pre}{F(v)*scale:,.{d}f}{suf}"
    return f"{f(m)} (95% UI: {f(lo)}–{f(hi)})"
def col(row, c, **k): return ui(row[c], row[c + "_lo"], row[c + "_hi"], **k)
def pct(row, c, d=1): return col(row, c, scale=100, d=d, suf="%")
def red_from_ratio(row, c, d=1):  # reduction = 1 - ratio; swap bounds
    return ui(1 - F(row[c]), 1 - F(row[c + "_hi"]), 1 - F(row[c + "_lo"]), scale=100, d=d, suf="%")
def ratio_pct(row, c, d=1): return col(row, c, scale=100, d=d, suf="%")
def tr(sid, yr): return next(r for r in traj if r["scenario_id"] == sid and int(F(r["year"])) == yr)
def ntr(yr): return next(r for r in nt if int(F(r["year"])) == yr)
def an(sid, yr, st): return next(r for r in ann if r["scenario_id"] == sid and int(F(r["year"])) == yr and r["stratum"] == st)
def money(row, c): return col(row, c, scale=1e-6, d=1, pre="S$", suf=" million")
def prob(sid): return f"{100*F(summ[sid]['joint_pass_probability']):.0f}%"
def red(sid, c): return pct(summ[sid], c)
def med(sid, c, scale=100, d=1): return f"{F(summ[sid][c])*scale:.{d}f}"

j = summ["J_only"]
nt30, nt25, nt19 = ntr(2030), ntr(2025), ntr(2019)
j19, j25, j26, j30 = tr("J_only", 2019), tr("J_only", 2025), tr("J_only", 2026), tr("J_only", 2030)
cF, cX = an("J_only", 2026, "F"), an("J_only", 2026, "X")
xshare = F(cX["chronic_stock_start"]) / (F(cX["chronic_stock_start"]) + F(cF["chronic_stock_start"]))
fx_death_share_2026 = F(j26["fx_mortality"]) / F(j26["mortality"])
fx_death_share_2030 = F(j30["fx_mortality"]) / F(j30["mortality"])
f_pool = F(cF["screening_eligible"])
# quota use over five years
def fx_tot(sid, st, c):
    return sum(F(r[c]) for r in ann if r["scenario_id"] == sid and r["stratum"] == st)
fq_8025 = fx_tot("S8000_F25", "F", "actual_screened") / fx_tot("S8000_F25", "F", "screening_quota")
fq_8075 = fx_tot("S8000_F75", "F", "actual_screened") / fx_tot("S8000_F75", "F", "screening_quota")
finite = [k for k, v in summ.items() if v["scenario_type"] == "finite"]
n_inc = sum(1 for k in finite if summ[k]["who_incidence_pass"] == "TRUE")
n_mort = sum(1 for k in finite if summ[k]["who_mortality_pass"] == "TRUE")
n_joint = sum(1 for k in finite if summ[k]["who_joint_pass"] == "TRUE")
ndraw = int(F(j["n_posterior_draws"]))
def inc_per_pp(sid):
    d = (F(summ[sid]["mortality_reduction"]) - F(j["mortality_reduction"])) * 100
    return F(summ[sid]["incremental_discounted_cost_vs_j_only_sgd"]) / 1e6 / d

def label(sid):
    r = summ[sid]
    if sid == "J_only": return "Prison treatment only (J only)"
    if sid == "SInf": return "J + no screening cap"
    return f"J + {int(F(r['screening_total'])):,}/year, F {round(100*F(r['f_share']))}%"

table1_ids = ["J_only", "S600_F50", "S2000_F25", "S2000_F75", "S5000_F25", "S5000_F50", "S5000_F75",
              "S8000_F25", "S8000_F50", "S8000_F75", "SInf"]
def t1row(sid):
    r = summ[sid]
    return (f"| {label(sid)} | {r['incidence_reduction_text']} | {r['mortality_reduction_text']} | "
            f"{prob(sid)} | {r['fx_starts_text'] if sid!='J_only' else '0'} | {col(r,'total_discounted_cost_2026_2030_sgd',scale=1e-6,d=1)} |")

main = f"""---
title: "Prison treatment and community screening to reduce hepatitis C infections and deaths among people who inject drugs in Singapore: A modelling study"
---

**Short title:** Hepatitis C treatment after prison release in Singapore

Yichen He^1^\\*

^1^ [Affiliation to be added]

\\* Corresponding author: [email to be added]

# Abstract

## Background

Hepatitis C virus (HCV) infection in Singapore is concentrated among people who inject drugs (PWID), who often cycle between the community and prison. Singapore expanded treatment in prison from 2019. We asked whether prison treatment alone could meet the World Health Organization (WHO) 2030 targets, and how community screening should be split after release.

## Methods and findings

We built an age structured HCV model with four groups: active PWID never imprisoned (D), people in prison (J), and released people who still inject (F) or have stopped (X). We calibrated it to 2015 prison data with Bayesian methods. We assumed treatment in prison from 2019. From 2026 we compared 15 community screening scenarios of 600 to 8,000 tests a year, with 25%, 50% or 75% for F and the rest for X. Results are medians with 95% uncertainty intervals (UI) across {ndraw:,} posterior draws. Without treatment, new infections in 2030 were {ratio_pct(nt30, 'incidence_ratio_2015')} of the 2015 level. With prison treatment alone, new infections fell by {red('J_only', 'incidence_reduction')} by 2030, but HCV deaths fell by only {red('J_only', 'mortality_reduction')}. In 2026, {100*xshare:.0f}% of chronic infections among released people were in X. Only {n_joint} of 15 scenarios met both WHO targets. At 8,000 tests a year, the incidence reduction was {med('S8000_F25','incidence_reduction')}% for every split, while the mortality reduction fell from {med('S8000_F25','mortality_reduction')}% to {med('S8000_F75','mortality_reduction')}% as the F share rose from 25% to 75%. The model covers only PWID and assumes ideal linkage to treatment.

## Conclusions

In this model, prison treatment sharply reduces new HCV infections, but the fall in HCV deaths lags. Most remaining deaths occur among released people who no longer inject. Community screening is needed to reach the mortality target. The F group is small, so screening beyond its size adds no benefit. Most community capacity should go to X.

# Author summary

## Why was this study done?

- Singapore aims to eliminate hepatitis C as a public health threat. Most infections are among people who inject drugs, and many of them pass through prison.
- Treatment in prison expanded in 2019. It is not known whether this alone can meet the WHO 2030 targets for new infections and deaths.
- Community services have limited screening capacity. It is not clear how this capacity should be shared between people who still inject after release and people who have stopped.

## What did the researchers do and find?

- We built a model of hepatitis C transmission and liver disease among people with a history of injecting drugs in Singapore. We calibrated it to prison data and carried uncertainty in all model parameters through to the results.
- Prison treatment alone reduced new infections by about three quarters by 2030. It reduced hepatitis C deaths by less than half.
- By 2026, nearly all chronic infections and deaths among released people were in those who had stopped injecting. Shifting screening towards people who still inject did not improve infections once that small group was covered. It did reduce the fall in deaths.

## What do these findings mean?

- Prison treatment is a strong tool for cutting transmission. Its effect on deaths is slow because deaths come from liver disease that has built up over many years.
- Reaching the mortality target needs large scale community screening, mainly among people who have stopped injecting.
- These results come from a model with ideal links from testing to treatment. Real programmes will need data on uptake, referral and cost before they set targets.

# Introduction

Hepatitis C virus (HCV) infection causes cirrhosis, liver cancer and death. Direct acting antivirals (DAAs) cure more than 95% of people treated [@pawlotsky2020]. The World Health Organization (WHO) has set targets for 2030. These are an 80% reduction in new HCV infections and a 65% reduction in HCV deaths compared with 2015 [@who2022strategy; @who2023validation]. Many countries are not on track [@dore2021].

In Singapore, HCV is concentrated among people who inject drugs (PWID) [@chaillon2021; @cda_hcv]. An earlier study estimated that about 11,000 people in Singapore have ever injected drugs [@teo2019]. Many of them move between the community, halfway houses and prison [@chaillon2021]. Prison is a key point of contact with health services. Singapore expanded DAA treatment in prison in 2019. In that year 681 people were treated in one prison cohort, with a cure rate above 99% [@prison2021]. Sofosbuvir and velpatasvir were later subsidised for wider use [@tan2023; @ace_daa]. Community programmes have also started. One programme plans to screen 600 high risk people a year [@moh2025ett; @moh2024ett]. A national HCV programme now coordinates this work [@cda_programme].

Treating people who inject drugs can prevent onward transmission [@martin2013; @lim2018]. Deaths are a different matter. They come from cirrhosis and liver cancer, which develop over decades [@thein2008]. Cure lowers but does not remove the risk of liver cancer in people with cirrhosis [@morgan2013; @vandermeer2012]. A programme can therefore cut infections quickly and still miss the mortality target. A previous Singapore model estimated the treatment needed for elimination [@chaillon2021]. It did not separate released people who still inject from those who have stopped. These two groups play different roles. The first drives transmission. The second may hold most of the burden of advanced liver disease.

We used a calibrated model to answer three questions. First, how much does prison treatment from 2019 reduce new infections and deaths by 2030? Second, how much community screening is needed to meet both WHO targets? Third, how should a fixed screening capacity be shared between released people who still inject and those who have stopped?

# Methods

## Study design and reporting

This is a modelling study. It used published and aggregate data only, so ethics approval was not required. The analysis was not prospectively registered. We report it following the CHEERS 2022 checklist where it applies [@cheers2022] (S1 Text, Section 7).

## Model structure

We built a deterministic compartmental model of HCV transmission and liver disease (Fig 1). The population includes everyone with a history of injecting drugs in Singapore. It is split into four groups. Active PWID who have never been imprisoned form group D. People in prison form group J. Released people who still inject form group F. Released people who have stopped injecting form group X. People enter D through new drug use. They move to J on first arrest. On release they move to F or X. People in F can be arrested again. The model has six age groups and four liver stages. The stages are no cirrhosis, compensated cirrhosis, decompensated cirrhosis and hepatocellular carcinoma. Each person is susceptible, acutely infected, chronically infected, on treatment, or recovered.

Transmission occurs only among active injectors in D and F. People in J and X are not exposed. Recovery does not give immunity, so recovered people in D and F can be reinfected. Liver disease progresses at rates taken from the literature [@thein2008; @lim2018; @kanwal2014]. Cure stops progression before cirrhosis and reduces it after cirrhosis [@morgan2013; @vandermeer2012]. Treatment is 12 weeks of sofosbuvir and velpatasvir. Cure rates depend on genotype and liver stage [@pawlotsky2020; @xiao2025]. People with HCV have higher death rates from decompensated cirrhosis and liver cancer [@lim2018; @mathers2013]. S1 Text gives all equations, parameter values and sources.

## Data and calibration

We fitted the model to 2015 data from the Singapore Prison Service. These were seroprevalence and prison population counts in six age groups. We estimated 12 parameters. Six scaled the contact rates by age. Six scaled the rate of entry into injecting by age. Each model run started in model year 0 and ran for 150 years to reach equilibrium. The end of that period represents 2015. From 2016, entry into injecting followed the official annual number of new drug users by age [@singstat_drug]. The 2025 value was held constant to 2030.

We used neural posterior estimation to obtain the posterior distribution of the 12 parameters [@greenberg2019; @papamakarios2017; @sbi2020]. We trained the estimator on 3,000 simulations. We trained it twice with different random seeds. We accepted the result only if simulation based calibration showed good coverage [@talts2018] and the posterior was not overly narrow. Both seeds passed these checks. A separate Markov chain Monte Carlo run gave nearly the same posterior. S1 Text, Section 2 gives details.

## Treatment and screening scenarios

All scenarios share the same history up to 2025. There was no treatment from 2016 to 2018. From 2019, people with chronic infection in prison started treatment at a high rate. This rate corresponds to an annual probability of 0.99995 when there are no competing events. We call this the prison only scenario (J only). It is a model assumption, not an observed treatment volume. We also ran a counterfactual with no treatment at any time.

From 2026, we added community screening for F and X. Each year a fixed number of screening tests was available. We tested 600, 1,000, 2,000, 5,000 and 8,000 tests a year. The 600 level matches the target of the existing community programme [@moh2025ett]. Each total was split between F and X in three ways. F received 25%, 50% or 75% and X received the rest. Unused tests in one group were not passed to the other. This gave 15 scenarios. We also ran an ideal reference in which everyone eligible in F and X was screened every year.

Screening worked as follows. At the start of each year, we counted the people in F and X who were not already on treatment. The number screened was the smaller of the quota and this eligible count. The share of those screened who had chronic infection set the maximum number of treatment starts that year. Identified people then started treatment at a rate equal to an annual probability of 0.99. The model ensured that starts never exceeded the maximum. No one in D was treated in any scenario.

## Outcomes

The two main outcomes were the reductions in 2030 compared with 2015 in annual new infections and annual HCV deaths. New infections include reinfections. HCV deaths are excess deaths from decompensated cirrhosis and liver cancer. These are event counts per year in the modelled population. They are not rates per person. A scenario met a target if its median reduction was at least 80% for infections or 65% for deaths. We also report the share of posterior draws in which both targets were met. Secondary outcomes were screening tests, treatment starts and costs from 2026 to 2030.

Costs covered antibody tests, RNA tests, one consultation per treatment start and a full course of sofosbuvir and velpatasvir. The drug cost was S$40,405 per course [@wong2019; @sgd_rate]. Test and consultation prices came from public price lists [@nhg_dsc; @parkway_rna]. These are model prices. They are not Singapore payer costs. Costs were discounted at 3% a year to 2026. Treatment starts were not discounted. S1 Text, Section 3 lists all cost inputs.

## Uncertainty analysis

We drew 1,000 parameter sets at random from each seed's posterior, giving {ndraw:,} draws. For each draw we ran the full analysis. This covered the 2015 equilibrium, the shared history to 2025 and all scenarios to 2030. All draws passed the equilibrium check. We report the median and the 2.5th and 97.5th percentiles across draws as the 95% uncertainty interval (UI). The UI therefore reflects uncertainty in the calibrated parameters and differences between the two training seeds. It does not reflect uncertainty in uptake, linkage to care or prices. We used R with C++ code through Rcpp [@rcpp2011] and Python with the sbi package [@sbi2020].

# Results

## Calibration

The model reproduced the age specific seroprevalence and population targets. All 95% predictive intervals overlapped the target intervals (S1 Text, Section 2). In 2015, the model estimated 20,915 people with a history of injecting drugs. Of these, 1,513 were in F and 14,680 were in X.

## Prison treatment reduced new infections but deaths lagged

Without any treatment, new infections rose after 2015 as entry into injecting increased. They remained high to 2030 (Fig 2). In 2030 annual new infections were {ratio_pct(nt30, 'incidence_ratio_2015')} of the 2015 level. HCV deaths were {ratio_pct(nt30, 'mortality_ratio_2015')} of the 2015 level.

Prison treatment from 2019 changed the course of new infections. By 2025 they were {ratio_pct(j25, 'incidence_ratio_2015')} of the 2015 level. By 2030 they had fallen by {red('J_only', 'incidence_reduction')}. This was close to the 80% target. Deaths fell more slowly. By 2025 they were still {ratio_pct(j25, 'mortality_ratio_2015')} of the 2015 level. By 2030 the reduction was {red('J_only', 'mortality_reduction')}, well short of the 65% target.

The reason lies in where infection and disease sit after release. In 2026, F held {col(cF, 'chronic_stock_start', d=0)} people with chronic infection. X held {col(cX, 'chronic_stock_start', d=0)}. X thus held {100*xshare:.0f}% of chronic infections among released people. In the J only scenario, {100*fx_death_share_2026:.0f}% of HCV deaths in 2026 occurred in F and X. By 2030 this share was {100*fx_death_share_2030:.0f}%. People in X do not inject and do not transmit. Prison treatment does not reach them after release.

## Community screening needed to meet both targets

Adding community screening increased both reductions (Table 1). The incidence target was met at the median in {'all 15' if n_inc==15 else str(n_inc)+' of 15'} scenarios. The mortality target was met in only {n_mort}. Only {n_joint} scenarios met both. These were 5,000 tests a year with 25% to F, and 8,000 tests a year with 25% or 50% to F. At 5,000 tests a year with 25% to F, the mortality reduction was {red('S5000_F25', 'mortality_reduction')}. Both targets were met in only {prob('S5000_F25')} of draws, so this result is fragile. At 8,000 tests a year with 25% to F, both targets were met in {prob('S8000_F25')} of draws. The ideal reference met both targets with a mortality reduction of {red('SInf', 'mortality_reduction')}.

At 600 tests a year, the planned size of the current community programme, the mortality reduction was at most {max(F(summ[k]['mortality_reduction']) for k in ['S600_F25','S600_F50','S600_F75'])*100:.1f}%. This was a gain of under 5 percentage points over prison treatment alone.

## How screening is split between F and X matters

The split between F and X changed the results in a consistent way (Fig 3). More screening for F gave a larger fall in new infections, up to a ceiling. More screening for X gave a larger fall in deaths.

The F group is small. In 2026 about {round(f_pool,-1):,.0f} people in F were eligible for screening each year. Once the F quota exceeded this number, extra tests could not be used (S1 Text, Fig S5). At 8,000 tests a year with 25% to F, {100*fq_8025:.0f}% of the F quota was used. With 75% to F, only {100*fq_8075:.0f}% was used. In all three splits at this capacity, the incidence reduction was {med('S8000_F25','incidence_reduction')}%. The extra F tests were idle, while the X tests they replaced were lost. The mortality reduction fell from {red('S8000_F25','mortality_reduction')} with 25% to F to {red('S8000_F75','mortality_reduction')} with 75% to F.

At lower capacities there was a real trade off. At 2,000 tests a year, moving from 25% to 75% for F raised the incidence reduction from {med('S2000_F25','incidence_reduction')}% to {med('S2000_F75','incidence_reduction')}%. It lowered the mortality reduction from {med('S2000_F25','mortality_reduction')}% to {med('S2000_F75','mortality_reduction')}%. A simple rule follows. F should receive enough tests to cover its eligible group each year, about 1,550. All remaining capacity should go to X.

## Treatment volume and cost

Treatment starts in F and X over 2026 to 2030 ranged from {summ['S600_F75']['fx_starts_text']} at 600 tests a year with 75% to F, to {summ['S8000_F25']['fx_starts_text']} at 8,000 tests a year with 25% to F (Table 1). Starts were highest in 2026 and fell each year as the stock of chronic infection was treated. The discounted cost of the three scenarios that met both targets ranged from S${F(summ['S8000_F50']['total_discounted_cost_2026_2030_sgd'])/1e6:.1f} million to S${F(summ['S8000_F25']['total_discounted_cost_2026_2030_sgd'])/1e6:.1f} million over five years. Drugs made up 95% to 96% of this cost. Compared with prison treatment alone, each extra percentage point of mortality reduction cost about S${(inc_per_pp('S5000_F25')+inc_per_pp('S8000_F25'))/2:.1f} million in the two scenarios with 25% to F at 5,000 and 8,000 tests a year. The two training seeds gave nearly the same results (S1 Text, Table S10).

**Table 1. HCV outcomes in 2030 and services from 2026 to 2030 in selected scenarios.**

| Scenario | Incidence reduction | Mortality reduction | Draws meeting both targets | F and X treatment starts | Discounted cost, 2026 to 2030 (S$ million) |
| --- | --- | --- | --- | --- | --- |
{chr(10).join(t1row(s) for s in table1_ids)}

Reductions compare 2030 with 2015. Values are medians (95% UI) across {ndraw:,} posterior draws. Targets are at least 80% for incidence and at least 65% for mortality. Costs are in 2026 Singapore dollars, discounted at 3% a year, and include prison treatment costs from 2026. Treatment starts count courses, not people. S1 Text, Table S8 gives all 17 scenarios.

# Discussion

In this model, prison treatment from 2019 cut new HCV infections among people who inject drugs in Singapore by about three quarters by 2030. It did much less for deaths. By 2026, nearly all chronic infections and HCV deaths among released people were in those who had stopped injecting. Community screening was needed to meet the mortality target. The split of screening between released people who still inject and those who have stopped mattered as much as the total amount.

Prison treatment works on transmission because most active injectors pass through prison. Curing them before release removes infection from the injecting network. This matches earlier work on treatment as prevention [@martin2013; @lim2018]. Deaths follow a different path. They arise from cirrhosis and liver cancer that built up over many years. Most of this burden sits in people who stopped injecting long ago. They are rarely in prison and do not transmit. Cure also does not remove all risk in people who already have cirrhosis [@morgan2013; @vandermeer2012]. So the mortality curve bends late. This lag is a feature of HCV natural history. It is not a failure of prison treatment.

Our findings on the F and X split have direct use for planning. A programme that targets people who still inject may seem the best choice for elimination. Our results show this holds only up to a point. The F group is small and quickly covered. After that, extra tests for F add nothing. The same tests given to X would prevent deaths. Planners should therefore set the F quota by the size of the eligible F group. They should direct the rest to people who have stopped injecting. These people may be harder to reach. They may not use drug services and may not see themselves as at risk. Outreach through primary care, addiction services and records of past imprisonment may be needed.

The scale required is large. The scenarios that met both targets needed 5,000 to 8,000 tests a year. That is 8 to 13 times the planned capacity of the current community programme [@moh2025ett]. They also assumed that every person identified starts treatment. Real linkage is weaker. In a Singapore trial among halfway house residents, only 12.6% of people with antibodies started treatment [@hsiang2022]. Care navigation and simpler models of care can raise uptake [@clink2022; @eckhardt2022]. Even then, real programmes will likely need more screening than our model suggests.

Our study has several strengths. We propagated parameter uncertainty through every scenario using {ndraw:,} posterior draws from two independently trained estimators. We modelled screening capacity explicitly, so treatment depends on how many people are found. We separated released people by injecting status, which earlier work did not do [@chaillon2021].

The study also has limitations. First, the model covers only people with a history of injecting drugs. It cannot show elimination in the whole population [@who2023validation]. Second, we calibrated the model to one year of prison data. We did not fit it to observed treatment numbers after 2019. The prison treatment assumption may therefore overstate or understate what happened. Third, we assumed perfect linkage from screening to treatment and no loss to follow up. Fourth, people in X cannot be reinfected or return to injecting. Some do relapse, which would add transmission. Fifth, entry into injecting uses counts of new drug users from law enforcement as a proxy. Sixth, costs use list prices and an older drug price. Actual payer costs are not public. Seventh, our uncertainty intervals cover model parameters only.

In conclusion, prison treatment has done much to reduce HCV transmission among people who inject drugs in Singapore. It will not, by itself, meet the WHO mortality target by 2030. Large scale community screening is needed. Most of it should reach people who have stopped injecting, after enough is set aside to cover those who still inject.

# Supporting information

**S1 Text. Supplementary Information.** Model equations and parameters, calibration methods and diagnostics, scenario definitions, cost inputs, full scenario results with 95% uncertainty intervals, results by training seed, cost sensitivity and the CHEERS 2022 checklist.

# Acknowledgments

[To be added.]

# Data availability

All model code, inputs and generated results are in the project repository [repository link to be added]. The analysis uses aggregate data only.

# Funding

[To be added.]

# Competing interests

[To be added.]

# Author contributions

[To be added.]

# Figure legends

**Fig 1. Model structure.** People enter group D when they start injecting. They move to prison (J) on arrest. On release they move to F if they keep injecting or to X if they stop. People in F can be arrested again. Transmission occurs only in D and F. Treatment in prison starts in 2019. Community screening for F and X starts in 2026.

![Fig 1. Model structure.](OUTPUT_DIR/scenarios_2026/paper_fig1_schematic.png){{width=100%}}

**Fig 2. Annual new HCV infections and HCV deaths relative to 2015.** (A) New infections, including reinfections. (B) HCV deaths from decompensated cirrhosis and liver cancer. Lines show medians and shaded bands show 95% UIs across {ndraw:,} posterior draws. Dashed horizontal lines show the WHO targets of 20% and 35% of the 2015 level. Dotted vertical lines mark 2019 and 2026. In panel A the two lines with 8,000 tests a year overlap.

![Fig 2. Annual new HCV infections and HCV deaths relative to 2015.](OUTPUT_DIR/scenarios_2026/paper_fig2_trajectories.png){{width=100%}}

**Fig 3. Reductions in 2030 by screening capacity and share allocated to F.** (A) Incidence. (B) Mortality. Points show medians and bars show 95% UIs. Dashed lines show the 80% and 65% targets. Grey lines show the prison only scenario and the ideal reference.

![Fig 3. Reductions in 2030 by screening capacity and share allocated to F.](OUTPUT_DIR/scenarios_2026/paper_fig3_allocation.png){{width=100%}}

# References
"""

# ---------------- Supplementary Information ----------------
calib = rd("output/calibration/r_20260929_agepop/predictions.csv")
ovl = rd("figures/calibration_interval_overlap.csv")
npe = rd("output/calibration/npe_agepop_recal/npe_evaluation.csv")
npe2 = rd("output/calibration/npe_agepop_recal/npe_evaluation_seed2.csv")
mcmc = rd("output/calibration/mcmc_recal/mcmc_parameter_summary.csv")
inflow = rd("data/singstat_new_drug_abusers.csv")
cfg = json.load(open(os.path.join(ROOT, "data/scenario_policy_config.json")))


byseed = list(csv.DictReader(open(_byseed)))
def quant(xs, p):
    xs = sorted(xs); i = (len(xs) - 1) * p; lo = int(i); hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (i - lo)
post = {}
for k, fn in (("1", "posterior_samples_npe.csv"), ("2", "posterior_samples_npe_seed2.csv")):
    post[k] = rd("output/calibration/npe_agepop_recal/" + fn)
theta_names = ([f"Contact scale, age {a}" for a in ["<20", "20–29", "30–39", "40–49", "50–59", "60+"]] +
               [f"Entry scale, age {a}" for a in ["<20", "20–29", "30–39", "40–49", "50–59", "60+"]])
def qtxt(xs, d=3): return f"{quant(xs,.5):.{d}f} ({quant(xs,.025):.{d}f} to {quant(xs,.975):.{d}f})"
post_rows = []
for i in range(12):
    t = f"theta{i+1}"
    a = [F(r[t]) for r in post["1"]]; b = [F(r[t]) for r in post["2"]]
    m = next(r for r in mcmc if r["parameter"] == t)
    post_rows.append(f"| θ{i+1} | {theta_names[i]} | {qtxt(a)} | {qtxt(b)} | {qtxt(a+b)} | {F(m['median']):.3f} ({F(m['q025']):.3f} to {F(m['q975']):.3f}) |")
sbc_rows = [f"| θ{i+1} | {F(npe[i]['coverage_95']):.3f} | {F(npe[i]['coverage_90']):.3f} | {F(npe[i]['posterior_sd']):.3f} | {F(npe2[i]['coverage_95']):.3f} | {F(npe2[i]['coverage_90']):.3f} | {F(npe2[i]['posterior_sd']):.3f} |" for i in range(12)]
seros = [o for o in ovl if o["measure"] == "seroprevalence"]; pops = [o for o in ovl if o["measure"] == "population"]
cal_rows = []
for c, s_, p_ in zip(calib, seros, pops):
    cal_rows.append(f"| {c['age_group']} | {100*F(c['seroprevalence_target']):.1f}% | {100*F(c['seroprevalence_hat']):.1f}% | {100*F(s_['NPE_lower']):.1f}% to {100*F(s_['NPE_upper']):.1f}% | {F(c['population_target']):.0f} | {F(c['population_hat']):.0f} | {F(p_['NPE_lower']):.0f} to {F(p_['NPE_upper']):.0f} |")
ages = ["<20", "20-29", "30-39", "40-49", "50-59", "60+"]
years = sorted({r["year"] for r in inflow})
inflow_rows = [f"| {y} | " + " | ".join(next(r["new_drug_abusers"] for r in inflow if r["year"] == y and r["age_group"] == a) for a in ages) + " |" for y in years]
def s6row(sid):
    r = summ[sid]
    return (f"| {label(sid)} | {r['incidence_reduction_text']} | {r['mortality_reduction_text']} | {prob(sid)} | "
            f"{pct(r,'fx_incidence_reduction')} | {pct(r,'fx_mortality_reduction')} |")
def s7row(sid):
    r = summ[sid]
    return (f"| {label(sid)} | {col(r,'fx_actual_screened',d=0)} | {col(r,'fx_rna_tests',d=0)} | {r['fx_starts_text']} | {r['j_starts_text']} | "
            f"{money(r,'total_discounted_cost_2026_2030_sgd')} | {money(r,'incremental_discounted_cost_vs_j_only_sgd')} |")
def s9row(sid):
    r = summ[sid]
    return "| " + label(sid) + " | " + " | ".join(col(r, f"cost_x{m}_sgd", scale=1e-6, d=1, pre="S$") for m in ("0_5", "1", "1_5")) + " |"
order = list(summ.keys())
def seedrow(r):
    f3 = lambda c: ui(r[c], r[c+"_lo"], r[c+"_hi"], scale=100, d=1, suf="%")
    return f"| {label(r['scenario_id'])} | {r['seed']} | {f3('incidence_reduction')} | {f3('mortality_reduction')} | {100*F(r['who_joint_pass_prob']):.0f}% |"
seed_ids = ["J_only", "S5000_F25", "S8000_F25", "S8000_F50", "S8000_F75", "SInf"]
cost_rows = [f"| {c['item']} | S${F(c['base_cost_sgd']):,.2f} | {c['unit']} | {c['source_currency']} {c['source_unit_price']} ({c['source_price_year']}) |" for c in cost_ev]

si = f"""---
title: "Supplementary Information"
subtitle: "Prison treatment and community screening to reduce hepatitis C infections and deaths among people who inject drugs in Singapore: A modelling study"
---

This document is S1 Text. It gives the technical details behind the main text. Section 1 describes the model. Section 2 covers data and calibration. Section 3 defines the scenarios and costs. Section 4 reports full results. Section 5 describes the uncertainty analysis and software. Section 6 lists further limitations. Section 7 is the CHEERS 2022 checklist.

# 1. Model structure

## 1.1 Population groups and states

The model is a deterministic system of ordinary differential equations. It has 480 compartments. These come from 4 population groups, 4 liver stages, 5 HCV states and 6 age groups.

- **Population groups.** D: active PWID never imprisoned. J: people in prison. F: people released from prison who still inject. X: people released from prison who have stopped injecting.
- **Liver stages.** Stage 1: no cirrhosis. Stage 2: compensated cirrhosis. Stage 3: decompensated cirrhosis. Stage 4: hepatocellular carcinoma.
- **HCV states.** u: susceptible. a: acute infection. c: chronic infection. t: on treatment. r: recovered after spontaneous clearance or cure.
- **Age groups.** <20, 20–29, 30–39, 40–49, 50–59 and 60+ years.

People enter D in the susceptible state and stage 1 at an age specific rate $\\beta_i(t)$. People in D move to J at the first arrest rate $\\lambda^{{(1)}}_i$. People in J are released at rate $\\lambda^{{(2)}}_i$. A share $\\pi = 0.59$ return to injecting and enter F. The rest enter X. People in F return to J at the re-arrest rate $\\lambda^{{(3)}}_i$. Arrest and release rates are fixed inputs based on Singapore Prison Service records. People in X stay in X.

## 1.2 Transmission

Transmission occurs only in D and F. The force of infection on age group $i$ is

$$\\gamma_i(t) = q\\; m(t) \\sum_j C_{{ij}} \\frac{{I_j(t)}}{{N_j(t)}}$$

where $q$ is the per event transmission probability [@boelen2014], $C_{{ij}}$ is the calibrated contact matrix, $I_j$ is the number of people in states a and c in D and F of age $j$, and $N_j$ is the number of people in D and F of age $j$. The time multiplier $m(t)$ was fixed at 1. People on treatment are not infectious. Susceptible and recovered people in D and F are infected at rate $\\gamma_i(t)$. People in J and X are not exposed.

## 1.3 Natural history and treatment

Acute infection lasts on average $\\iota_1 = 180$ days [@fasano2024; @tanaka2000]. A share $\\kappa = 0.26$ clears infection spontaneously [@micallef2006]. The rest become chronically infected. Chronic infection progresses through the four liver stages. Progression rates are weighted by genotype, with 78% genotype 3 [@zhao2019]. Genotype 3 has higher progression rates [@kanwal2014]. The effective rate for transition $A \\to B$ is

$$\\tilde p_{{A\\to B}} = \\rho\\, r_3^{{A\\to B}} p^{{A\\to B}} + (1-\\rho)\\, p^{{A\\to B}}.$$

Treatment is sofosbuvir and velpatasvir for $\\iota_2 = 84$ days. People on treatment are cured with probability $\\bar\\alpha_k$ and return to state c otherwise. In the scenarios, decompensated cirrhosis used the cure rate without ribavirin, to match the costed regimen. After cure, people in stage 1 do not progress. People in stage 2 progress at reduced rates [@vandermeer2012; @morgan2013]. People in stage 3 keep the same risk of liver cancer. Treatment lowers excess mortality in stages 3 and 4 by factors $\\psi^{{DC}}$ and $\\psi^{{HCC}}$ [@premkumar2024; @dang2020]. Background mortality varies by age and is multiplied by a standardised mortality ratio for PWID [@mathers2013]. Excess mortality applies in stages 3 and 4 [@lim2018].

The state transitions within each group, stage and age are

$$a \\xrightarrow{{\\kappa/\\iota_1}} r,\\quad a \\xrightarrow{{(1-\\kappa)/\\iota_1}} c,\\quad c \\xrightarrow{{\\tau}} t,\\quad t \\xrightarrow{{\\bar\\alpha/\\iota_2}} r,\\quad t \\xrightarrow{{(1-\\bar\\alpha)/\\iota_2}} c,\\quad u, r \\xrightarrow{{\\gamma_i(t)}} a$$

As an example, the equation for chronic infection in group F, stage 1, age $i$ is

$$\\frac{{dF_{{c,1,i}}}}{{dt}} = \\pi\\lambda^{{(2)}}_i J_{{c,1,i}} + \\frac{{1-\\kappa}}{{\\iota_1}} F_{{a,1,i}} + \\frac{{1-\\bar\\alpha_1}}{{\\iota_2}} F_{{t,1,i}} - \\left(\\tau_{{F}} + \\lambda^{{(3)}}_i + \\omega\\mu_i + \\tilde p_{{1\\to2}}\\right) F_{{c,1,i}}$$

The other 479 equations follow the same pattern. The full set is in the project file `docs/equations.tex` and is implemented in `scripts/model/sim.cpp`.

## 1.4 Parameters

**Table S1. Natural history and treatment parameters.**

| Parameter | Description | Value | Source |
| --- | --- | --- | --- |
| $q$ | Transmission probability per sharing event | 0.0057 | [@boelen2014] |
| $\\kappa$ | Probability of spontaneous clearance | 0.26 | [@micallef2006] |
| $\\iota_1$ | Mean duration of acute infection | 180 days | [@fasano2024] |
| $\\iota_2$ | Duration of treatment | 84 days | [@pawlotsky2020] |
| $\\rho$ | Proportion with genotype 3 | 0.78 | [@zhao2019] |
| $\\bar\\alpha_1$ | Cure rate, stages 1 and 2 | 0.979 | [@pawlotsky2020] |
| $\\bar\\alpha_3^-$ | Cure rate, stage 3, no ribavirin | 0.549 | [@pawlotsky2020; @xiao2025] |
| $\\bar\\alpha_3^+$ | Cure rate, stage 3, with ribavirin | 0.882 | [@pawlotsky2020; @xiao2025] |
| $p^{{NC\\to CC}}$ | Progression, no cirrhosis to compensated cirrhosis (per year) | 0.027 | [@thein2008] |
| $p^{{CC\\to DC}}$ | Progression, compensated to decompensated cirrhosis | 0.039 | [@lim2018] |
| $p^{{CC\\to HCC}}$ | Progression, compensated cirrhosis to liver cancer | 0.014 | [@lim2018] |
| $p^{{DC\\to HCC}}$ | Progression, decompensated cirrhosis to liver cancer | 0.014 | [@lim2018] |
| $r_3^{{NC\\to CC}}, r_3^{{CC\\to DC}}$ | Relative risk for genotype 3, fibrosis | 1.36 | [@kanwal2014] |
| $r_3^{{CC\\to HCC}}, r_3^{{DC\\to HCC}}$ | Relative risk for genotype 3, liver cancer | 1.93 | [@kanwal2014] |
| $\\phi^{{CC\\to DC}}$ | Relative progression after cure, compensated to decompensated | 0.07 | [@vandermeer2012] |
| $\\phi^{{CC\\to HCC}}$ | Relative progression after cure, compensated cirrhosis to liver cancer | 0.23 | [@morgan2013] |
| $\\phi^{{DC\\to HCC}}$ | Relative progression after cure, decompensated cirrhosis to liver cancer | 1.00 | Assumed |
| $\\omega$ | Standardised mortality ratio for PWID | 14.68 | [@mathers2013] |
| $\\mu^{{DC}}$ | Excess mortality, decompensated cirrhosis (per year) | 0.13 | [@lim2018] |
| $\\mu^{{HCC}}$ | Excess mortality, liver cancer (per year) | 0.43 | [@lim2018] |
| $\\psi^{{DC}}$ | Relative excess mortality after treatment, stage 3 | 0.45 | [@premkumar2024] |
| $\\psi^{{HCC}}$ | Relative excess mortality after treatment, stage 4 | 0.37 | [@dang2020] |
| $\\pi$ | Share returning to injecting after release | 0.59 | Based on prison records |

Genotype specific cure rates were 0.976 (genotype 3) and 0.988 (other) in stages 1 and 2. In stage 3 they were 0.45 and 0.90 without ribavirin and 0.87 and 0.92 with ribavirin [@pawlotsky2020; @xiao2025]. Liver cancer used the stage 3 cure rate because no specific trial data exist.

## 1.5 Main assumptions

- There is no lasting immunity after clearance or cure. Reinfection follows the same force of infection as first infection [@pawlotsky2020].
- Liver disease does not progress during acute infection.
- Fibrosis does not regress after cure.
- People in X do not relapse to injecting and cannot be reinfected.
- People on treatment are not infectious.
- The genotype mix is fixed over time.

# 2. Data and calibration

## 2.1 Targets

The calibration targets were 2015 Singapore Prison Service data on HCV seroprevalence and the number of people in prison with a history of drug use, in six age groups. Seroprevalence in the model is the share of people in states a, c, t and r among all people in J. The prevalence terms used a binomial likelihood. The population terms used a log normal likelihood with standard deviation 0.10 on the log scale.

## 2.2 Entry into injecting

From 2016 to 2025, entry into D followed the annual number of new drug users by age published by the Singapore Department of Statistics [@singstat_drug] (Table S2). Each age specific series was multiplied by a calibrated scale. The 2025 value was held constant to 2030. These counts come from law enforcement records and are a proxy for entry into injecting.

**Table S2. New drug users by age group, Singapore, 2015 to 2025.**

| Year | <20 | 20–29 | 30–39 | 40–49 | 50–59 | 60+ |
| --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(inflow_rows)}

## 2.3 Calibrated parameters and priors

We estimated 12 parameters on the log scale. Parameters θ1 to θ6 multiply the contact rates of the six age groups. Parameters θ7 to θ12 multiply the entry rates of the six age groups. A deterministic fit with seven Nelder Mead starts gave the centre of the prior. The best fit had a prevalence root mean square error of 0.0087 and a population mean absolute percentage error of 3.8%. Every model evaluation ran from model year 0 to 150. The last five years had to be at equilibrium. The prior for each parameter was normal on the log scale, centred on the best fit, with standard deviation 0.20.

## 2.4 Neural posterior estimation

We used neural posterior estimation (NPE-C) [@greenberg2019] with a masked autoregressive flow [@papamakarios2017], implemented in the sbi package [@sbi2020]. The network had 128 hidden features and 10 transforms. We trained it for up to 500 epochs on 3,000 simulations from the prior. Simulated data included binomial noise on prevalence and log normal noise on population. We trained the network twice with different seeds and drew 5,000 posterior samples from each.

We checked each network with simulation based calibration (SBC) on 300 test simulations [@talts2018]. We accepted a network if coverage of the 95% and 90% intervals was close to nominal for every parameter and if no posterior standard deviation was below 0.05. Both seeds passed (Table S3).

**Table S3. Simulation based calibration and posterior spread by training seed.**

| Parameter | Seed 1, 95% coverage | Seed 1, 90% coverage | Seed 1, posterior SD | Seed 2, 95% coverage | Seed 2, 90% coverage | Seed 2, posterior SD |
| --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(sbc_rows)}

## 2.5 Posterior distribution

Table S4 gives the posterior of each parameter by seed and pooled. As a check we also ran a Markov chain Monte Carlo sampler with four chains of 120,000 iterations each, after 20,000 burn in iterations. All parameters had an effective sample size above 10,000 and R-hat below 1.001. MCMC medians differed from pooled NPE medians by at most 0.02 on the log scale.

**Table S4. Posterior median (95% interval) of the calibrated log scale parameters.**

| Parameter | Description | NPE seed 1 | NPE seed 2 | NPE pooled | MCMC |
| --- | --- | --- | --- | --- | --- |
{chr(10).join(post_rows)}

## 2.6 Model fit

Table S5 compares the targets with the best deterministic fit and the NPE 95% predictive interval. All 12 NPE and MCMC predictive intervals overlapped the target intervals. Fig S1 shows the fit.

**Table S5. Calibration targets and model fit, 2015.**

| Age group | Seroprevalence target | Best fit | NPE 95% interval | Population target | Best fit | NPE 95% interval |
| --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(cal_rows)}

![Fig S1. Calibration of age specific seroprevalence and prison population in 2015. NPE and MCMC predictive intervals are compared with target intervals.](FIGURES_DIR/calibration_agepop.png){{width=100%}}

In 2015 the model estimated 1,386 (95% interval: 1,299–1,483) people in D, 3,340 (3,128–3,571) in J, 1,513 (1,407–1,631) in F and 14,680 (13,511–15,989) in X. The total was 20,915 (19,360–22,643) people with a history of injecting drugs. Of these, 4,241 (3,919–4,613) were viraemic.

# 3. Scenarios and costs

## 3.1 Shared history

All scenarios start from the same 2015 equilibrium. There was no treatment from 2016 to 2018. From 2019 the chronic infection treatment rate in J was $\\tau_J = -\\ln(1 - 0.99995)$ per year with no annual cap. People already on treatment completed their course. People who failed treatment or were reinfected could start again. No one in D was treated in any scenario. The no treatment counterfactual kept $\\tau = 0$ in all groups up to 2030.

## 3.2 Community screening quota

From 2026, each year began with a screening plan for F and X. For group $s$ we computed $N_s = u + a + c + r$, which excludes people on treatment. The number screened was $S_s = \\min(Q_s, N_s)$, where $Q_s$ is the annual quota. The number of antibody tests was $S_s$. The number of RNA tests was $S_s (a + c + r)/N_s$. The number of people with chronic infection identified was $K_s = S_s\\, c/N_s$. $K_s$ was the maximum number of treatment starts in group $s$ that year. Identified people started treatment at rate $\\tau_{{FX}} = -\\ln(1 - 0.99)$ per year. The ODE solver tracked cumulative starts and switched off the $c \\to t$ flow in group $s$ once starts reached $K_s$. Unused quota was not moved between F and X or carried to the next year. In the ideal reference, $Q_s$ was infinite, so everyone eligible was screened once a year.

This is an aggregate approximation. The model does not track individual screening history. The same person can be counted as screened in more than one year. Acute infections were not treated.

**Table S6. Scenario definitions.**

| Scenario | Annual tests | F quota per year | X quota per year |
| --- | --- | --- | --- |
| J only | 0 | 0 | 0 |
""" + "\n".join(f"| {label(k)} | {int(F(summ[k]['screening_total'])):,} | {F(summ[k]['f_screening_quota']):,.0f} | {F(summ[k]['x_screening_quota']):,.0f} |" for k in order if summ[k]['scenario_type']=='finite') + f"""
| J + no screening cap | Unlimited | All eligible | All eligible |

## 3.3 Outcomes

New infections are the annual flow from states u and r to state a in D and F, evaluated at the end of each calendar year. HCV deaths are the annual excess deaths in stages 3 and 4. Both are event counts per year. The reduction in 2030 is one minus the ratio of the 2030 value to the 2015 value. We also report these reductions for F and X combined.

## 3.4 Costs

Costs were counted for F, J and X from 2026 to 2030 (Table S7). In F and X, each person screened had one antibody test. Each expected antibody positive person had one RNA test. Each treatment start had one consultation and one full drug course. In J, each antibody positive person entering prison had one antibody test. Each treatment start in J had one RNA test, one consultation and one drug course. Costs were counted in the year they occurred and discounted as $C_y / 1.03^{{(y-2026)}}$.

The drug price was USD 2,577 per week for 12 weeks [@wong2019]. This is USD 30,924 per course. We converted it at the 2025 average rate of 1.3066 Singapore dollars per US dollar [@sgd_rate]. Prices were not inflated. These are model prices and not Singapore payer costs. Follow up visits, SVR12 tests, fibrosis assessment and programme staff were not included.

**Table S7. Unit costs.**

| Item | Unit cost | Unit | Source price |
| --- | --- | --- | --- |
{chr(10).join(cost_rows)}

Sources: drug price [@wong2019], antibody test and consultation [@nhg_dsc], RNA test [@parkway_rna].

# 4. Full results

## 4.1 Targets in all scenarios

**Table S8. Reductions in 2030 compared with 2015, all 17 scenarios.**

| Scenario | Incidence reduction | Mortality reduction | Draws meeting both targets | F and X incidence reduction | F and X mortality reduction |
| --- | --- | --- | --- | --- | --- |
{chr(10).join(s6row(k) for k in order)}

Values are medians (95% UI) across {ndraw:,} posterior draws.

![Fig S2. Trajectories of new infections and HCV deaths relative to 2015 (left) and reductions in 2030 (right) in all 17 scenarios. Medians across posterior draws.](OUTPUT_DIR/scenarios_2026/scenario_target_trajectories.png){{width=100%}}

![Fig S3. Attainment of the 2030 targets by scenario. Cells show median reductions. Joint status requires both targets.](OUTPUT_DIR/scenarios_2026/scenario_target_heatmap.png){{width=90%}}

## 4.2 Services and costs

**Table S9. Services and discounted costs, 2026 to 2030.**

| Scenario | F and X people screened | F and X RNA tests | F and X treatment starts | J treatment starts | Total discounted cost | Incremental cost vs J only |
| --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(s7row(k) for k in order)}

Screening counts may include the same person in more than one year. Treatment starts count courses.

![Fig S4. Annual treatment starts and screening derived treatment quotas in F and X, by capacity and allocation. Medians across posterior draws.](OUTPUT_DIR/scenarios_2026/scenario_treatment_starts.png){{width=100%}}

![Fig S5. Share of the five year screening quota used in F and X. The F group is small, so a large F share leaves part of its quota unused.](OUTPUT_DIR/scenarios_2026/paper_fig4_f_pool.png){{width=90%}}

![Fig S6. Cost composition of the prison only scenario and the scenarios that met both targets. Shares use the median cost of each item.](OUTPUT_DIR/scenarios_2026/scenario_cost_composition_pies.png){{width=100%}}

![Fig S7. Discounted costs and treatment starts by group, 2026 to 2030. Segments show the median for each group. Labels show the median total.](OUTPUT_DIR/scenarios_2026/scenario_cumulative_costs.png){{width=100%}}

## 4.3 Results by training seed

**Table S10. Selected results by NPE training seed.**

| Scenario | Seed | Incidence reduction | Mortality reduction | Draws meeting both targets |
| --- | --- | --- | --- | --- |
{chr(10).join(seedrow(r) for sid in seed_ids for r in byseed if r['scenario_id']==sid)}

Each seed contributed 1,000 draws.

## 4.4 Cost sensitivity

Table S11 scales all unit costs by 0.5, 1 and 1.5. Epidemiological results do not change.

**Table S11. Total discounted cost (S$ million) under unit cost multipliers.**

| Scenario | 0.5 times | 1 times | 1.5 times |
| --- | --- | --- | --- |
{chr(10).join(s9row(k) for k in order)}

# 5. Uncertainty analysis and software

We sampled 1,000 parameter sets at random without replacement from each seed's 5,000 posterior samples, using sampling seed 20261005. For each set we ran the full pipeline. This included the 150 year equilibrium, the equilibrium check, the shared history from 2016 to 2025, the no treatment counterfactual and all 17 scenarios to 2030. All {ndraw:,} draws passed the equilibrium check. Derived quantities such as reductions, costs and incremental costs were computed within each draw before summarising. We report the median and the 2.5th and 97.5th percentiles. A scenario met a target if its median reduction met the threshold. We also report the share of draws meeting both targets.

The model was solved in C++ through Rcpp [@rcpp2011] with a daily time step. NPE used Python, PyTorch and sbi [@sbi2020]. The full analysis of {ndraw:,} draws took about 18 minutes on 14 processor cores. The project test suite (62 checks) passed. Running the scenario code at the pooled posterior median reproduced the earlier single point results exactly.

# 6. Further limitations

- **Population scope.** The model includes only people with a history of injecting drugs. It does not include other routes of infection or the general population. National elimination needs national data [@who2023validation; @cda_hcv].
- **Prison treatment history.** The 2019 prison treatment rate is an assumption. The model was not fitted to observed treatment numbers. In a Singapore prison cohort, 681 people were treated in 2019 [@prison2021]. The model gave {col(an('J_only',2026,'J'),'treatment_starts',d=0)} prison treatment starts in 2026 because most prevalent cases had been treated by then.
- **Linkage to care.** Screening, confirmation and treatment start were assumed to be complete. Real linkage is lower [@hsiang2022; @clink2022].
- **Group X.** People in X cannot relapse or be reinfected. Relapse would add transmission and reinfection.
- **Event counts.** Outcomes are annual event counts in the modelled population. They are not rates per person.
- **Costs.** Prices are public list prices and an older drug price. Many programme costs are not included.
- **Uncertainty.** The UIs reflect calibrated parameters only. Fixed natural history parameters, uptake and price were not varied.

# 7. CHEERS 2022 checklist

**Table S12. CHEERS 2022 items and where they are reported.**

| No. | Item | Location |
| --- | --- | --- |
| 1 | Title | Title page |
| 2 | Abstract | Abstract |
| 3 | Background and objectives | Introduction |
| 4 | Health economic analysis plan | Methods, Study design (no prospective plan) |
| 5 | Study population | Methods, Model structure |
| 6 | Setting and location | Introduction, Methods |
| 7 | Comparators | Methods, Treatment and screening scenarios |
| 8 | Perspective | Methods, Outcomes (health service service costs only) |
| 9 | Time horizon | Methods, Outcomes (2026 to 2030) |
| 10 | Discount rate | Methods, Outcomes; S1 Text 3.4 |
| 11 | Selection of outcomes | Methods, Outcomes |
| 12 | Measurement of outcomes | Methods, Outcomes; S1 Text 3.3 |
| 13 | Valuation of outcomes | Not applicable (no utility weights) |
| 14 | Measurement and valuation of resources and costs | S1 Text 3.4 |
| 15 | Currency, price date and conversion | S1 Text 3.4 |
| 16 | Rationale and description of model | Methods; S1 Text 1 |
| 17 | Analytics and assumptions | S1 Text 1.5, 3 |
| 18 | Characterising heterogeneity | Results by group and age structure in model |
| 19 | Characterising distributional effects | Not assessed |
| 20 | Characterising uncertainty | Methods, Uncertainty analysis; S1 Text 4.4, 5 |
| 21 | Approach to engagement with patients and others affected by the study | No engagement |
| 22 | Study parameters | S1 Text Tables S1 to S7 |
| 23 | Summary of main results | Results; Table 1 |
| 24 | Effect of uncertainty | Results; S1 Text 4.3, 4.4 |
| 25 | Effect of engagement with patients and others affected by the study | Not applicable |
| 26 | Study findings, limitations, generalisability and current knowledge | Discussion |
| 27 | Source of funding | Funding statement |
| 28 | Conflicts of interest | Competing interests statement |

# References
"""
si = si.replace("health service service costs", "health service costs")
DOCS = os.path.join(ROOT, "docs")
for name, text, out in (("main.md", main, "main.docx"), ("supplementary.md", si, "Supplementary Information.docx")):
    text = text.replace("OUTPUT_DIR/", os.path.join(ROOT, "output") + "/").replace("FIGURES_DIR/", os.path.join(ROOT, "figures") + "/")
    md = os.path.join(TMP, name); open(md, "w").write(text)
    subprocess.run(["pandoc", md, "--citeproc", "--bibliography", os.path.join(DOCS, "references.bib"),
                    "--csl", os.path.join(DOCS, "plos.csl"), "-o", os.path.join(DOCS, out)], check=True)
print("Wrote docs/main.docx and docs/Supplementary Information.docx")
