#!/usr/bin/env python3
# =============================================================================
# build_paper.py — builds docs/main.docx and docs/Supplementary Information.docx from
# generated results; every number in the text is read from the result files.
#
# Usage: python3 scripts/build_paper.py [ROOT]. Needs output/scenarios_2026/, output/calibration/,
# docs/references.bib, docs/plos.csl and pandoc. Edit prose here, not in Word.
# =============================================================================
import csv, json, re, sys, os, subprocess, tempfile
ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.getcwd())
TMP = tempfile.mkdtemp()
S = os.path.join(ROOT, "output/scenarios_2026")
def rd(p): return list(csv.DictReader(open(os.path.join(ROOT, p))))
SETUP = open(os.path.join(ROOT, "scripts/model/setup.R")).read()
def setup_val(name):  # numeric value(s) of a `name = ...` entry in setup.R
    v = re.search(rf"^\s*{name}\s*=\s*(c\([^)]*\)|[0-9.eE+-]+)", SETUP, re.M).group(1)
    return [float(x) for x in re.findall(r"[0-9.]+(?:[eE][+-]?[0-9]+)?", v)]
PI, HR, VREF, PSI_DC = (setup_val(n)[0] for n in ("pi_recid", "hr_viraemic", "viraemic_share_ref", "psi_DC"))
F_FIRST = setup_val("f_first")
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
def ratio_pct(row, c, d=1): return col(row, c, scale=100, d=d, suf="%")
def tr(sid, yr): return next(r for r in traj if r["scenario_id"] == sid and int(F(r["year"])) == yr)
def ntr(yr): return next(r for r in nt if int(F(r["year"])) == yr)
def an(sid, yr, st): return next(r for r in ann if r["scenario_id"] == sid and int(F(r["year"])) == yr and r["stratum"] == st)
def money(row, c): return col(row, c, scale=1e-6, d=1, pre="S$", suf=" million")
def prob(sid): return f"{100*F(summ[sid]['joint_pass_probability']):.0f}%"
def red(sid, c): return pct(summ[sid], c)
def med(sid, c, scale=100, d=1): return f"{F(summ[sid][c])*scale:.{d}f}"

j = summ["J_only"]
j_inc = 100 * F(j["incidence_reduction"])
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
def scen_phrase(ids):
    by = {}
    for k in ids:
        by.setdefault(int(F(summ[k]["screening_total"])), []).append(round(100 * F(summ[k]["f_share"])))
    parts = []
    for tot in sorted(by):
        sh = [f"{x}%" for x in sorted(by[tot])]
        sh_txt = sh[0] if len(sh) == 1 else ", ".join(sh[:-1]) + " or " + sh[-1]
        parts.append(f"{tot:,} tests a year with {sh_txt} to F")
    return " and ".join(parts)
mort_finite = [k for k in finite if summ[k]["who_mortality_pass"] == "TRUE"]
j_mort = 100 * F(j["mortality_reduction"])
sinf_inc = 100 * F(summ["SInf"]["incidence_reduction"])
max_inc = max(100 * F(summ[k]["incidence_reduction"]) for k in finite)
d_inc_share = 1 - F(summ["SInf"]["fx_incidence"]) / F(summ["SInf"]["incidence"])  # D share of 2030 infections
# The narrative below is written for these results; stop the build if they change.
assert n_joint == 0 and n_inc == 0 and sinf_inc < 80 and mort_finite, "results changed: revise the narrative"
assert abs(j_inc - j_mort) < 10, "prison treatment no longer cuts infections and deaths alike: revise the narrative"
def mprob(sid): return f"{100*F(summ[sid]['who_mortality_pass_prob']):.0f}%"
pass_caps = sorted({int(F(summ[k]["screening_total"])) for k in mort_finite})
cap_txt = (f"{pass_caps[0]:,} tests a year" if len(pass_caps) == 1 else f"{pass_caps[0]:,} to {pass_caps[-1]:,} tests a year")
mult_txt = (f"about {round(pass_caps[0] / 600)} times" if len(pass_caps) == 1 else f"{round(pass_caps[0] / 600)} to {round(pass_caps[-1] / 600)} times")
inc8 = [F(summ[f"S8000_F{x}"]["incidence_reduction"]) * 100 for x in (25, 50, 75)]
inc8_txt = f"{min(inc8):.1f}%" if round(min(inc8), 1) == round(max(inc8), 1) else f"{min(inc8):.1f}% to {max(inc8):.1f}%"
prior_dom = [x for x in rd("output/calibration/mcmc_reference/mcmc_parameter_summary.csv") if F(x["sd"]) / 0.20 > 0.8]
n_dom_words = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six"}.get(len(prior_dom), str(len(prior_dom)))
def drug_share(sid):
    rr = [x for x in ann if x["scenario_id"] == sid]
    return sum(F(x["discounted_cost_daa_sgd"]) for x in rr) / sum(F(x["discounted_cost_total_sgd"]) for x in rr)
ds = [drug_share(k) for k in mort_finite]
ds_txt = (f"{100*min(ds):.0f}%" if round(100*min(ds)) == round(100*max(ds)) else f"{100*min(ds):.0f}% to {100*max(ds):.0f}%")
pass_costs = [F(summ[k]["total_discounted_cost_2026_2030_sgd"]) / 1e6 for k in mort_finite]
strata = {r["measure"]: r for r in rd("output/calibration/mcmc_reference/model_prediction_2015_strata.csv")}
CMP = "output/calibration/npe_vs_mcmc"
acc_sum = {r["run"]: r for r in rd(CMP + "/accuracy_summary.csv")}
pol_cmp = rd(CMP + "/policy_comparison.csv")
npe_l = acc_sum["npe_lhs"]
npe_w = F(npe_l["mean_width95_ratio"])
npe_pol = max(F(r["error"]) for r in pol_cmp if r["run"] == "npe_lhs" and r["error_unit"] == "relative")
npe_pp = max(F(r["error"]) for r in pol_cmp if r["run"] == "npe_lhs" and r["error_unit"] != "relative")
npe_ss = json.load(open(os.path.join(ROOT, "output/calibration/npe_lhs/simulation_summary.json")))
npe_cfg = json.load(open(os.path.join(ROOT, "output/calibration/npe_lhs/npe_config.json")))
npe_sbc = [r for r in rd("output/calibration/npe_lhs/sbc_summary.csv")]
mc_sbc = {r["parameter"]: r for r in rd(CMP + "/mcmc_sbc/mcmc_sbc_summary.csv")}["all"]
ref_cfg = json.load(open(os.path.join(ROOT, "output/calibration/mcmc_reference/mcmc_run_config.json")))
def rng(k, d=2): v = [F(r[k]) for r in npe_sbc]; return f"{min(v):.{d}f} to {max(v):.{d}f}"
sbc_uniform = all(r["rank_shape"] == "uniform" for r in npe_sbc)
w_txt = (f"on average {100*(1-npe_w):.0f}% narrower than" if npe_w < 0.995 else
         f"on average {100*(npe_w-1):.0f}% wider than" if npe_w > 1.005 else "about as wide as")
def stx(m, full=True):
    r = strata[m]
    return (f"{F(r['median']):,.0f} (95% interval: {F(r['lower_95']):,.0f}–{F(r['upper_95']):,.0f})" if full
            else f"{F(r['median']):,.0f} ({F(r['lower_95']):,.0f}–{F(r['upper_95']):,.0f})")
assert summ["S5000_F25"]["who_mortality_pass"] != "TRUE", "S5000_F25 now meets the mortality target: revise the narrative"
s5_txt = (f"At 5,000 tests a year with 25% to F, the median mortality reduction was {red('S5000_F25', 'mortality_reduction')}, "
          f"just below the target. The mortality target was met in {mprob('S5000_F25')} of draws, so this scenario sits on the threshold.")
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
            f"{mprob(sid)} | {r['fx_starts_text'] if sid!='J_only' else '0'} | {col(r,'total_discounted_cost_2026_2030_sgd',scale=1e-6,d=1)} |")

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

We built an age structured HCV model with four groups: active PWID never imprisoned (D), people in prison (J), and released people who still inject (F) or have stopped (X). We calibrated it to 2015 prison data with Bayesian methods. We assumed treatment in prison from 2019. From 2026 we compared 15 community screening scenarios of 600 to 8,000 tests a year, with 25%, 50% or 75% for F and the rest for X. Results are medians with 95% uncertainty intervals (UI) across {ndraw:,} posterior draws. Without treatment, new infections in 2030 were {ratio_pct(nt30, 'incidence_ratio_2015')} of the 2015 level. With prison treatment alone, new infections fell by {red('J_only', 'incidence_reduction')} by 2030 and HCV deaths by {red('J_only', 'mortality_reduction')}. In 2026, {100*xshare:.0f}% of chronic infections among released people were in X. No scenario met the incidence target: even screening everyone eligible after release cut new infections by only {sinf_inc:.0f}%, because {100*d_inc_share:.0f}% of new infections in 2030 occurred among people who had never been imprisoned. The mortality target was met at the median only at {scen_phrase(mort_finite)}. At 8,000 tests a year, the incidence reduction was {inc8_txt} for every split, while the mortality reduction fell from {med('S8000_F25','mortality_reduction')}% to {med('S8000_F75','mortality_reduction')}% as the F share rose from 25% to 75%. The model covers only PWID and assumes ideal linkage to treatment.

## Conclusions

In this model, prison treatment alone roughly halves new HCV infections and HCV deaths by 2030, short of both targets. Most remaining deaths occur among released people who no longer inject, so large scale community screening is needed to reach the mortality target. The F group is small, so screening beyond its size adds no benefit, and most community capacity should go to X. The incidence target cannot be met through prison and post release services, because most new infections occur in people who have never been imprisoned.

# Author summary

## Why was this study done?

- Singapore aims to eliminate hepatitis C as a public health threat. Most infections are among people who inject drugs, and many of them pass through prison.
- Treatment in prison expanded in 2019. It is not known whether this alone can meet the WHO 2030 targets for new infections and deaths.
- Community services have limited screening capacity. It is not clear how this capacity should be shared between people who still inject after release and people who have stopped.

## What did the researchers do and find?

- We built a model of hepatitis C transmission and liver disease among people with a history of injecting drugs in Singapore. We calibrated it to prison data and carried uncertainty in all model parameters through to the results.
- Prison treatment alone reduced new infections by {j_inc:.0f}% and hepatitis C deaths by {j_mort:.0f}% by 2030, short of the WHO targets of 80% and 65%.
- By 2026, nearly all chronic infections and deaths among released people were in those who had stopped injecting. Shifting screening towards people who still inject did not improve infections once that small group was covered. It did reduce the fall in deaths.
- No strategy met the target for new infections. Most new infections occurred among people who inject but have never been imprisoned, whom prison and post release services do not reach.

## What do these findings mean?

- Prison treatment reaches only people who are imprisoned. Cutting new infections further needs testing and prevention for people who inject before they are ever imprisoned.
- Reaching the mortality target needs large scale community screening, mainly among people who have stopped injecting.
- These results come from a model with ideal links from testing to treatment. Real programmes will need data on uptake, referral and cost before they set targets.

# Introduction

Hepatitis C virus (HCV) infection causes cirrhosis, liver cancer and death. Direct acting antivirals (DAAs) cure more than 95% of people treated [@pawlotsky2020]. The World Health Organization (WHO) has set targets for 2030. These are an 80% reduction in new HCV infections and a 65% reduction in HCV deaths compared with 2015 [@who2022strategy; @who2023validation]. Many countries are not on track [@dore2021].

In Singapore, HCV is concentrated among people who inject drugs (PWID) [@chaillon2021; @cda_hcv]. An earlier study estimated that about 11,000 people in Singapore have ever injected drugs [@teo2019]. Many of them move between the community, halfway houses and prison [@chaillon2021]. Prison is a key point of contact with health services. Singapore expanded DAA treatment in prison in 2019. In that year 681 people were treated in one prison cohort, with a cure rate above 99% [@prison2021]. Sofosbuvir and velpatasvir were later subsidised for wider use [@tan2023; @ace_daa]. Community programmes have also started. One programme plans to screen 600 high risk people a year [@moh2025ett; @moh2024ett]. A national HCV programme now coordinates this work [@cda_programme].

Treating people who inject drugs can prevent onward transmission [@martin2013; @lim2018]. Deaths are a different matter. They come from cirrhosis and liver cancer, which develop over decades [@thein2008]. Cure lowers but does not remove the risk of liver cancer in people with cirrhosis [@morgan2013; @vandermeer2012]. A programme can therefore cut infections quickly and still miss the mortality target. A previous Singapore model estimated the treatment needed for elimination [@chaillon2021]. It did not separate released people who still inject from those who have stopped. These two groups play different roles. The first drives transmission. The second may hold most of the burden of advanced liver disease.

We used a calibrated model to answer three questions. First, how much does prison treatment from 2019 reduce new infections and deaths by 2030? Second, how much community screening is needed to meet the WHO targets? Third, how should a fixed screening capacity be shared between released people who still inject and those who have stopped?

# Methods

## Study design and reporting

This is a modelling study. It used published and aggregate data only, so ethics approval was not required. The analysis was not prospectively registered. We report it following the CHEERS 2022 checklist where it applies [@cheers2022] (S1 Text, Section 7).

## Model structure

We built a deterministic compartmental model of HCV transmission and liver disease (Fig 1). The population includes everyone with a history of injecting drugs in Singapore. It is split into four groups. Active PWID who have never been imprisoned form group D. People in prison form group J. Released people who still inject form group F. Released people who have stopped injecting form group X. People enter D through new drug use. They move to J on first arrest. On release they move to F or X. People in F can be arrested again. The model has six age groups and four liver stages. The stages are no cirrhosis, compensated cirrhosis, decompensated cirrhosis and hepatocellular carcinoma. Each person is susceptible, acutely infected, chronically infected, on treatment, or recovered.

Transmission occurs only among active injectors in D and F. People in J and X are not exposed. Recovery does not give immunity, so recovered people in D and F can be reinfected. Liver disease progresses at rates taken from the literature [@thein2008; @lim2018; @kanwal2014]. Cure stops progression before cirrhosis and reduces it after cirrhosis [@morgan2013; @vandermeer2012]. Treatment is 12 weeks of sofosbuvir and velpatasvir. Cure rates depend on genotype and liver stage [@pawlotsky2020; @xiao2025]. People with HCV have higher death rates from decompensated cirrhosis and liver cancer [@lim2018; @mathers2013]. People with chronic infection also die {HR:.2f} times as fast from other causes [@lee2012]. Release and re-arrest rates and the yearly number of first arrests came from Singapore Prison Service records, and the first arrest rate was set so that the model reproduces that number. S1 Text gives all equations, parameter values and sources.

## Data and calibration

We fitted the model to 2015 data from the Singapore Prison Service. These were seroprevalence and prison population counts in six age groups. We estimated 12 parameters. Six scaled the contact rates by age. Six set the rate of entry into injecting by age, as an excess above the minimum entry that supplies the observed first arrests. Each model run started in model year 0 and ran for 150 years to reach equilibrium. The end of that period represents 2015. From 2016, entry into injecting followed the official annual number of new drug users by age [@singstat_drug]. The 2025 value was held constant to 2030.

A multi start Nelder Mead fit gave the best fitting parameter values. We used them as the centre of a normal prior with standard deviation 0.20 on the log scale. Any best fit value closer than five prior standard deviations to a parameter bound was moved inward to that distance, so that the prior lies within the bounds. The likelihood was binomial for the prevalence counts and log normal for the population counts and for the 2015 entry into injecting, compared with the official count of new drug users. Parameter sets that did not reach equilibrium were excluded. We sampled the posterior with an adaptive random walk Metropolis algorithm [@haario2001]. We ran four chains from independent prior draws, with 50,000 model runs each. All parameters reached rank normalised R-hat below 1.01 [@vehtari2021]. We also trained single round neural posterior estimation [@greenberg2019; @sbi2020] on {npe_ss['n_attempted']:,} prior simulations drawn by Latin hypercube sampling and compared it with MCMC. Simulation based calibration on {int(npe_sbc[0]['n_sbc']):,} simulated datasets showed well calibrated neural posteriors (95% coverage {rng('coverage_95')}). Policy results agreed within {100*npe_pol:.1f}%. The neural posterior's 95% intervals were {w_txt} those from MCMC. We used MCMC, which samples the exact posterior and whose error falls with longer runs. S1 Text, Section 2 gives details.

## Treatment and screening scenarios

All scenarios share the same history up to 2025. There was no treatment from 2016 to 2018. From 2019, people with chronic infection in prison started treatment at a high rate. This rate corresponds to an annual probability of 0.99995 when there are no competing events. We call this the prison only scenario (J only). It is a model assumption, not an observed treatment volume. We also ran a counterfactual with no treatment at any time.

From 2026, we added community screening for F and X. Each year a fixed number of screening tests was available. We tested 600, 1,000, 2,000, 5,000 and 8,000 tests a year. The 600 level matches the target of the existing community programme [@moh2025ett]. Each total was split between F and X in three ways. F received 25%, 50% or 75% and X received the rest. Unused tests in one group were not passed to the other. This gave 15 scenarios. We also ran an ideal reference in which everyone eligible in F and X was screened every year.

Screening worked as follows. At the start of each year, we counted the people in F and X who were not already on treatment. The number screened was the smaller of the quota and this eligible count. The share of those screened who had chronic infection set the maximum number of treatment starts that year. Identified people then started treatment at a rate equal to an annual probability of 0.99. The model ensured that starts never exceeded the maximum. No one in D was treated in any scenario.

## Outcomes

The two main outcomes were the reductions in 2030 compared with 2015 in annual new infections and annual HCV deaths. New infections include reinfections. HCV deaths are excess deaths from decompensated cirrhosis and liver cancer. These are event counts per year in the modelled population. They are not rates per person. A scenario met a target if its median reduction was at least 80% for infections or 65% for deaths. We also report the share of posterior draws in which each target, and both together, were met. Secondary outcomes were screening tests, treatment starts and costs from 2026 to 2030.

Costs covered antibody tests, RNA tests, one consultation per treatment start and a full course of sofosbuvir and velpatasvir. The drug cost was S$40,405 per course [@wong2019; @sgd_rate]. Test and consultation prices came from public price lists [@nhg_dsc; @parkway_rna]. These are model prices. They are not Singapore payer costs. Costs were discounted at 3% a year to 2026. Treatment starts were not discounted. S1 Text, Section 3 lists all cost inputs.

## Uncertainty analysis

We drew {ndraw:,} parameter sets at random from the MCMC posterior. For each draw we ran the full analysis. This covered the 2015 equilibrium, the shared history to 2025 and all scenarios to 2030. All draws passed the equilibrium check. We report the median and the 2.5th and 97.5th percentiles across draws as the 95% uncertainty interval (UI). The UI therefore reflects uncertainty in the calibrated parameters. It does not reflect uncertainty in uptake, linkage to care or prices. We used R with C++ code through Rcpp [@rcpp2011].

# Results

## Calibration

The model reproduced the age specific seroprevalence and population targets. All 95% predictive intervals overlapped the target intervals (S1 Text, Section 2). In 2015, the model estimated {stx('population_ever_PWID')} people with a history of injecting drugs. Of these, {stx('population_F', False)} were in F and {stx('population_X', False)} were in X.

## Prison treatment alone met neither target

Without any treatment, new infections rose after 2015 as entry into injecting increased. They remained high to 2030 (Fig 2). In 2030 annual new infections were {ratio_pct(nt30, 'incidence_ratio_2015')} of the 2015 level. HCV deaths were {ratio_pct(nt30, 'mortality_ratio_2015')} of the 2015 level.

Prison treatment from 2019 changed the course of new infections. By 2025 they were {ratio_pct(j25, 'incidence_ratio_2015')} of the 2015 level. By 2030 they had fallen by {red('J_only', 'incidence_reduction')}. {"This was close to" if j_inc >= 70 else "This was well short of"} the 80% target. Deaths fell at a similar pace. By 2025 they were {ratio_pct(j25, 'mortality_ratio_2015')} of the 2015 level. By 2030 the reduction was {red('J_only', 'mortality_reduction')}, well short of the 65% target.

The reason lies in where infection and disease sit after release. In 2026, F held {col(cF, 'chronic_stock_start', d=0)} people with chronic infection. X held {col(cX, 'chronic_stock_start', d=0)}. X thus held {100*xshare:.0f}% of chronic infections among released people. In the J only scenario, {100*fx_death_share_2026:.0f}% of HCV deaths in 2026 occurred in F and X. By 2030 this share was {100*fx_death_share_2030:.0f}%. People in X do not inject and do not transmit. Prison treatment does not reach them after release.

## Community screening can meet the mortality target but not the incidence target

Adding community screening increased both reductions (Table 1). No scenario met the incidence target. The largest reduction in a finite scenario was {max_inc:.1f}%, and screening everyone eligible in F and X every year gave {red('SInf', 'incidence_reduction')}. In that ideal reference, {100*d_inc_share:.0f}% of new infections in 2030 occurred in D, among people who inject but have never been imprisoned. Neither prison treatment nor screening after release reaches them, and no draw met the incidence target in any scenario.

The mortality target was met at the median in {n_mort} of 15 scenarios: {scen_phrase(mort_finite)}. At 8,000 tests a year with 25% to F, the mortality target was met in {mprob('S8000_F25')} of draws. With 50% to F, the share was {mprob('S8000_F50')}. {s5_txt} The ideal reference reduced mortality by {red('SInf', 'mortality_reduction')}.

At 600 tests a year, the planned size of the current community programme, the mortality reduction was at most {max(F(summ[k]['mortality_reduction']) for k in ['S600_F25','S600_F50','S600_F75'])*100:.1f}%. This was a gain of under 5 percentage points over prison treatment alone.

## How screening is split between F and X matters

The split between F and X changed the results in a consistent way (Fig 3). More screening for F gave a larger fall in new infections, up to a ceiling. More screening for X gave a larger fall in deaths.

The F group is small. In 2026 about {round(f_pool,-1):,.0f} people in F were eligible for screening each year. Once the F quota exceeded this number, extra tests could not be used (S1 Text, Fig S5). At 8,000 tests a year with 25% to F, {100*fq_8025:.0f}% of the F quota was used. With 75% to F, only {100*fq_8075:.0f}% was used. In all three splits at this capacity, the incidence reduction was {inc8_txt}. The extra F tests were idle, while the X tests they replaced were lost. The mortality reduction fell from {red('S8000_F25','mortality_reduction')} with 25% to F to {red('S8000_F75','mortality_reduction')} with 75% to F.

At lower capacities there was a real trade off. At 2,000 tests a year, moving from 25% to 75% for F raised the incidence reduction from {med('S2000_F25','incidence_reduction')}% to {med('S2000_F75','incidence_reduction')}%. It lowered the mortality reduction from {med('S2000_F25','mortality_reduction')}% to {med('S2000_F75','mortality_reduction')}%. A simple rule follows. F should receive enough tests to cover its eligible group each year, about 1,550. All remaining capacity should go to X.

## Treatment volume and cost

Treatment starts in F and X over 2026 to 2030 ranged from {summ['S600_F75']['fx_starts_text']} at 600 tests a year with 75% to F, to {summ['S8000_F25']['fx_starts_text']} at 8,000 tests a year with 25% to F (Table 1). Starts were highest in 2026 and fell each year as the stock of chronic infection was treated. {f"The discounted cost of the finite scenario that met the mortality target was S${pass_costs[0]:.1f} million" if len(mort_finite) == 1 else f"The discounted cost of the {len(mort_finite)} finite scenarios that met the mortality target ranged from S${min(pass_costs):.1f} million to S${max(pass_costs):.1f} million"} over five years. Drugs made up {ds_txt} of this cost. Compared with prison treatment alone, each extra percentage point of mortality reduction cost about S${(inc_per_pp('S5000_F25')+inc_per_pp('S8000_F25'))/2:.1f} million in the two scenarios with 25% to F at 5,000 and 8,000 tests a year. Single round neural posterior estimation gave policy results within {100*npe_pol:.1f}% of these values (S1 Text, Table S4).

**Table 1. HCV outcomes in 2030 and services from 2026 to 2030 in selected scenarios.**

| Scenario | Incidence reduction | Mortality reduction | Draws meeting the mortality target | F and X treatment starts | Discounted cost, 2026 to 2030 (S$ million) |
| --- | --- | --- | --- | --- | --- |
{chr(10).join(t1row(s) for s in table1_ids)}

Reductions compare 2030 with 2015. Values are medians (95% UI) across {ndraw:,} posterior draws. Targets are at least 80% for incidence and at least 65% for mortality. No draw met the incidence target in any scenario. Costs are in 2026 Singapore dollars, discounted at 3% a year, and include prison treatment costs from 2026. Treatment starts count courses, not people. S1 Text, Table S9 gives all 17 scenarios.

# Discussion

In this model, prison treatment from 2019 cut new HCV infections among people who inject drugs in Singapore by {j_inc:.0f}% and HCV deaths by {j_mort:.0f}% by 2030. By 2026, nearly all chronic infections and HCV deaths among released people were in those who had stopped injecting. Community screening was needed to meet the mortality target. The split of screening between released people who still inject and those who have stopped mattered as much as the total amount. No strategy we modelled met the incidence target, because most new infections occurred in people who had not yet been imprisoned.

Prison treatment works on transmission among people who pass through prison. Curing them before release removes infection from that part of the injecting network, which matches earlier work on treatment as prevention [@martin2013; @lim2018]. But people inject for years before their first arrest, and in this model they acquired {100*d_inc_share:.0f}% of new infections in 2030 even when everyone eligible after release was screened. Reducing these infections needs testing and prevention in the community before imprisonment. Deaths follow a different path. They arise from cirrhosis and liver cancer that built up over many years. Most of this burden sits in people who stopped injecting long ago. They are rarely in prison and do not transmit. Cure also does not remove all risk in people who already have cirrhosis [@morgan2013; @vandermeer2012].

Our findings on the F and X split have direct use for planning. A programme that targets people who still inject may seem the best choice for elimination. Our results show this holds only up to a point. The F group is small and quickly covered. After that, extra tests for F add nothing. The same tests given to X would prevent deaths. Planners should therefore set the F quota by the size of the eligible F group. They should direct the rest to people who have stopped injecting. These people may be harder to reach. They may not use drug services and may not see themselves as at risk. Outreach through primary care, addiction services and records of past imprisonment may be needed.

The scale required is large. The scenarios that met the mortality target needed {cap_txt}. That is {mult_txt} the planned capacity of the current community programme [@moh2025ett]. Even 5,000 tests a year with 25% to F met the mortality target in only {mprob('S5000_F25')} of draws. All scenarios also assumed that every person identified starts treatment. Real linkage is weaker. In a Singapore trial among halfway house residents, only 12.6% of people with antibodies started treatment [@hsiang2022]. Care navigation and simpler models of care can raise uptake [@clink2022; @eckhardt2022]. Even then, real programmes will likely need more screening than our model suggests.

Our study has several strengths. We propagated parameter uncertainty through every scenario using {ndraw:,} draws from a well converged MCMC posterior. We checked the posterior against an independent method, neural posterior estimation. We modelled screening capacity explicitly, so treatment depends on how many people are found. We separated released people by injecting status, which earlier work did not do [@chaillon2021].

The study also has limitations. First, the model covers only people with a history of injecting drugs. It cannot show elimination in the whole population [@who2023validation]. Second, we calibrated the model to one year of prison data. We did not fit it to observed treatment numbers after 2019. The prison treatment assumption may therefore overstate or understate what happened. Third, we assumed perfect linkage from screening to treatment and no loss to follow up. Fourth, people in X cannot be reinfected or return to injecting. Some do relapse, which would add transmission. Fifth, entry into injecting uses counts of new drug users from law enforcement as a proxy, and the same counts anchor the number of people who inject without ever being imprisoned. That number sets the incidence result that prison and post release services cannot reach. Sixth, costs use list prices and an older drug price. Actual payer costs are not public. Seventh, our uncertainty intervals cover model parameters only. The prior was centred on a fit to the same data, and for {n_dom_words} parameters the data narrowed the prior little.

In conclusion, prison treatment has done much to reduce HCV among people who inject drugs in Singapore, but by itself it will not meet either WHO target by 2030. Large scale community screening after release can meet the mortality target, if most of it reaches people who have stopped injecting after enough is set aside for those who still inject. Meeting the incidence target needs testing and prevention that reach people who inject before they are ever imprisoned.

# Supporting information

**S1 Text. Supplementary Information.** Model equations and parameters, calibration methods and diagnostics, the comparison with neural posterior estimation, scenario definitions, cost inputs, full scenario results with 95% uncertainty intervals, cost sensitivity and the CHEERS 2022 checklist.

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
calib = rd("output/calibration/r_20261014_agepop/predictions.csv")
hr1 = rd("output/calibration/hr_sensitivity/hr_1/predictions.csv")
hr_sens = [rd(f"output/calibration/hr_sensitivity/hr_{h}/predictions.csv") for h in ("1.11", "1.23", "1.77")]
nm_cfg = json.load(open(os.path.join(ROOT, "output/calibration/r_20261014_agepop/run_config.json")))
nm_sol = rd("output/calibration/r_20261014_agepop/solutions.csv")
nm_best = min(nm_sol, key=lambda r: F(r["objective"]))
ovl = rd("figures/calibration_interval_overlap.csv")
REF = "output/calibration/mcmc_reference"
ref_diag = [r for r in rd(REF + "/mcmc_diagnostics.csv") if r["subset"] == "all_retained"]
ref_cs = rd(REF + "/mcmc_chain_summary.csv")
for r in ref_cs:
    r["acceptance"] = r.get("acceptance") or r.get("acceptance_forward_budget")
inflow = rd("data/singstat_new_drug_abusers.csv")
cfg = json.load(open(os.path.join(ROOT, "data/scenario_policy_config.json")))


def quant(xs, p):
    xs = sorted(xs); i = (len(xs) - 1) * p; lo = int(i); hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (i - lo)
ref_draws = rd(REF + "/mcmc_draws.csv")
theta_names = ([f"Contact scale, age {a}" for a in ["<20", "20–29", "30–39", "40–49", "50–59", "60+"]] +
               [f"Excess entry scale, age {a}" for a in ["<20", "20–29", "30–39", "40–49", "50–59", "60+"]])
def qtxt(xs, d=3): return f"{quant(xs,.5):.{d}f} ({quant(xs,.025):.{d}f} to {quant(xs,.975):.{d}f})"
post_rows = []; sd_ratio = []
for i in range(12):
    t = f"theta{i+1}"
    a = [F(r[t]) for r in ref_draws]
    m = sum(a) / len(a); sd = (sum((v - m) ** 2 for v in a) / (len(a) - 1)) ** 0.5
    post_rows.append(f"| θ{i+1} | {theta_names[i]} | {F(nm_best[t]):.3f} | {ref_cfg['prior_centre'][i]:.3f} | {qtxt(a)} | {sd:.3f} | {sd/0.20:.2f} |")
    sd_ratio.append(sd / 0.20)
moved = [f"θ{i+1}" for i in range(12) if abs(ref_cfg['prior_centre'][i] - F(nm_best[f"theta{i+1}"])) > 1e-9]
prior_dom = [f"θ{i+1}" for i in range(12) if sd_ratio[i] > 0.8]
assert moved and prior_dom
pd_txt = ", ".join(prior_dom[:-1]) + " and " + prior_dom[-1] if len(prior_dom) > 1 else prior_dom[0]
diag_rows = [f"| θ{i+1} | {F(d['rhat']):.4f} | {F(d['ess_bulk']):,.0f} | {F(d['ess_tail']):,.0f} | {F(d['mcse_mean'])/F(d['sd']):.3f} |" for i, d in enumerate(ref_diag)]
def bench_row(k):
    a = acc_sum[k]
    pc = [r for r in pol_cmp if r["run"] == k and r["error_unit"] == "relative"]
    pp = [r for r in pol_cmp if r["run"] == k and r["error_unit"] != "relative"]
    name = f"Single round NPE, {npe_cfg['n_ensemble']} network ensemble, Latin hypercube sampling"
    return (f"| {name} | {int(F(a['forward_runs'])):,} ({int(F(a['valid_simulations'])):,}) | {F(a['max_abs_smd']):.2f} | "
            f"{F(a['max_abs_endpoint']):.2f} | {F(a['mean_width95_ratio']):.3f} ({F(a['min_width95_ratio']):.2f} to {F(a['max_width95_ratio']):.2f}) | "
            f"{100*F(a['share_pass_all']):.0f}% | {100*max(F(r['error']) for r in pc):.1f}% | {max(F(r['error']) for r in pp):.1f} |")
bench_rows = [bench_row("npe_lhs")]
seros = [o for o in ovl if o["measure"] == "seroprevalence"]; pops = [o for o in ovl if o["measure"] == "population"]
n_ovl = sum(o["MCMC_overlap"] == "TRUE" for o in ovl)
miss_txt = "" if n_ovl == 12 else " (not: " + ", ".join(f"{o['measure']} {o['age_group']}" for o in ovl if o["MCMC_overlap"] != "TRUE") + ")"
cal_rows = []
for c, s_, p_ in zip(calib, seros, pops):
    cal_rows.append(f"| {c['age_group']} | {100*F(c['seroprevalence_target']):.1f}% | {100*F(c['seroprevalence_hat']):.1f}% | {100*F(s_['MCMC_lower']):.1f}% to {100*F(s_['MCMC_upper']):.1f}% | {F(c['population_target']):.0f} | {F(c['population_hat']):.0f} | {F(p_['MCMC_lower']):.0f} to {F(p_['MCMC_upper']):.0f} |")
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

People enter D in the susceptible state and stage 1 at an age specific rate $\\beta_i(t)$. People in D move to J at the first arrest rate $\\lambda^{{(1)}}_i$. People in J are released at rate $\\lambda^{{(2)}}_i$. A share $\\pi = {PI:.2f}$ return to injecting and enter F. The rest enter X. People in F return to J at the re-arrest rate $\\lambda^{{(3)}}_i$. People in X stay in X.

Release and re-arrest rates and $\\pi$ were fitted with Stan to Singapore Prison Service records of everyone held for drug consumption. The modelled population is the people with a known HCV screening result, who were assumed to be the PWID; their 2015 prison counts are the calibration targets. The re-arrest model includes the 2015 prison steady state, with the same ageing (0.1 per year) and background mortality as the model. It gives the yearly first arrest flow $f_i$ ({', '.join(f'{x:g}' for x in F_FIRST)} people by age). The first arrest rate is derived from it: in the HCV free steady state of D, $D_i = (\\beta_i + 0.1 D_{{i-1}} - f_i)/(\\omega\\mu_i + a_i)$, with ageing $a_i = 0.1$ ($a_6 = 0$), and $\\lambda^{{(1)}}_i = f_i / D_i$. Entry is parameterised as the minimum that supplies $f_i$ plus a positive excess (Section 2.3), so every parameter set can supply the first arrests.

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

$$\\frac{{dF_{{c,1,i}}}}{{dt}} = \\pi\\lambda^{{(2)}}_i J_{{c,1,i}} + \\frac{{1-\\kappa}}{{\\iota_1}} F_{{a,1,i}} + \\frac{{1-\\bar\\alpha_1}}{{\\iota_2}} F_{{t,1,i}} - \\left(\\tau_{{F}} + \\lambda^{{(3)}}_i + h\\,\\omega_0\\mu_i + \\tilde p_{{1\\to2}}\\right) F_{{c,1,i}}$$

Background mortality is $\\omega_0\\mu_i$, where $\\mu_i$ is the Singapore resident death rate. Chronic infection (state c) multiplies it by $h = {HR:.2f}$. This is the multivariable adjusted hazard ratio for extrahepatic death with detectable HCV RNA, against anti-HCV negative people, in the REVEAL-HCV cohort (95% CI 1.23 to 1.77; people with undetectable RNA had 0.90) [@lee2012]. We used the extrahepatic ratio because liver deaths are explicit in stages 3 and 4, and we applied it to state c only because the cohort measured chronic carriers. $\\omega_0 = \\omega/(1 + (h-1)v)$ keeps the average PWID rate at the pooled standardised mortality ratio $\\omega$ [@mathers2013], with $v = {VREF:.3f}$ the chronic viraemic share of PWID at the 2015 best fit (fixed point; held constant so that treating one person does not change the death rate of others). These extra deaths are not counted as HCV-related deaths.

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
| $\\psi^{{DC}}$ | Relative excess mortality after treatment, stage 3 | {PSI_DC:.2f} | [@premkumar2024] |
| $h$ | Hazard ratio for non-liver death, chronic infection | {HR:.2f} | [@lee2012] |
| $v$ | Chronic viraemic share of PWID, 2015 (fixed point) | {VREF:.3f} | Model fit |
| $\\psi^{{HCC}}$ | Relative excess mortality after treatment, stage 4 | 0.37 | [@dang2020] |
| $\\pi$ | Share returning to injecting after release | {PI:.2f} | Fitted to prison records |
| $f_i$ | First arrests per year, by age | {', '.join(f'{x:g}' for x in F_FIRST)} | Fitted to prison records |

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

We estimated 12 parameters on the log scale. Parameters θ1 to θ6 multiply the contact rates of the six age groups. Parameters θ7 to θ12 set the entry rates of the six age groups. Entry into age group $i$ is the minimum that supplies its first arrests, $\\max(0, f_i - 0.1 D_{{i-1}})$, plus an excess $e^{{\\theta_{{6+i}}}}\\beta_i$, so every parameter set can supply the first arrest flow. Where first arrests use up the entry, the excess sets the size of the never imprisoned group D; at age 60+, where people ageing in supply the first arrests, it is the whole entry. A deterministic Nelder Mead fit with {len(nm_sol)} starts gave the centre of the prior. The starts were an earlier anchor, {len(nm_sol)-2} random perturbations of it (standard deviation {nm_cfg['sd_perturb']} on the log scale) and a population informed start. Each start ran for up to {nm_cfg['maxit']:,} iterations and was then restarted from its own solution, up to {nm_cfg['nm_restarts']} times, while the objective still fell by more than 0.001. The best fit had an objective of {F(nm_best['objective']):.2f}, a prevalence root mean square error of {F(nm_best['prevalence_rmse']):.4f} and a population mean absolute percentage error of {100*F(nm_best['population_mape']):.1f}%. Every model evaluation ran from model year 0 to 150. The last ten years had to be at equilibrium. The prior for each parameter was normal on the log scale with standard deviation 0.20. Its centre was the best fit, except that any value closer than five prior standard deviations (1.0 on the log scale) to a parameter bound was moved inward to that distance. This changed only {', '.join(moved)}, which the best fit placed at its lower bound, where the likelihood is flat (Table S5). Without this step about half of the prior mass for that parameter would lie outside the bounds. Contact scales were bounded to 0.01 to 1,000 and entry scales to 0.01 to 10⁶.

## 2.4 Posterior sampling

The posterior combined the prior, a binomial likelihood for the six prevalence counts, a log normal likelihood (standard deviation 0.10) for the six population counts and a log normal likelihood (standard deviation 1.0) comparing the 2015 entry into D with the official count of new drug users [@singstat_drug]. Because first arrests are fixed by the prison records, the prison data alone do not determine how many people inject without ever being imprisoned; the entry term anchors that number to the official count, allowing roughly a sevenfold difference either way. Parameter sets outside the bounds or whose simulation failed the equilibrium check had zero posterior density. We sampled it with adaptive random walk Metropolis [@haario2001]. Four chains started from independent prior draws. The initial proposal covariance was 0.01 times the prior variance. Adaptation began after {ref_cfg['adapt_start']:,} proposals, was refreshed every {ref_cfg['adapt_every']} proposals and was frozen after warmup. A 5% share of proposals always used the initial covariance. Each chain ran {ref_cfg['max_evals_per_chain']:,} model evaluations, of which {ref_cfg['n_warmup']:,} were warmup, giving {sum(int(F(r['retained_total'])) for r in ref_cs):,} retained iterations. Every {ref_cfg['thin_storage']}th iteration was stored ({len(ref_draws):,} draws), and diagnostics used all retained iterations. Acceptance was {min(F(r['acceptance']) for r in ref_cs):.2f} to {max(F(r['acceptance']) for r in ref_cs):.2f}. All parameters met the convergence checks (Table S3) [@vehtari2021].

**Table S3. Convergence diagnostics of the MCMC posterior.**

| Parameter | Rank normalised R-hat | Bulk ESS | Tail ESS | MCSE of mean / posterior SD |
| --- | --- | --- | --- | --- |
{chr(10).join(diag_rows)}

## 2.5 Comparison with neural posterior estimation

We checked the MCMC posterior against single round neural posterior estimation (NPE) [@greenberg2019], implemented in the sbi package [@sbi2020]. Both methods used the same prior, likelihood, solver and equilibrium restriction. We drew {npe_ss['n_attempted']:,} parameter sets from the prior by Latin hypercube sampling. The design had one point in each of {npe_ss['n_attempted']:,} equal probability strata for every parameter, and each column was mapped through the normal inverse distribution function of that parameter's prior. Of these, {npe_ss['n_valid']:,} simulations were valid. Exclusions were {npe_ss['n_excluded_bounds']:,} draws outside the bounds, {npe_ss['n_excluded_equilibrium']:,} equilibrium failures and {npe_ss['n_excluded_nonfinite']:,} non-finite simulations; excluded draws are dropped, not replaced. The estimator was an equal weight ensemble of {npe_cfg['n_ensemble']} masked autoregressive flows, each with 50 hidden features and 5 transforms [@papamakarios2017], trained in one round with the standard loss from different random seeds. A single network gave a small seed dependent bias in one parameter, which the ensemble removed.

We checked the neural posterior by simulation based calibration [@talts2018] on {int(npe_sbc[0]['n_sbc']):,} datasets simulated from independent prior draws. For each parameter, 95% interval coverage was {rng('coverage_95', 3)} and 90% coverage {rng('coverage_90', 3)}. {"No" if sbc_uniform else "Some"} rank histogram{"" if sbc_uniform else "s"} showed systematic skew (U shape, central peak, shift or Bonferroni corrected chi-square departure). The predeclared gate required 95% coverage of at least 0.85 and 90% coverage of at least 0.75 for every parameter. As a smaller direct check, we reran the MCMC on {mc_sbc['n_datasets']} of these simulated datasets, with one chain of 30,000 evaluations each. Pooled 95% coverage was {F(mc_sbc['coverage_95']):.2f} and 90% coverage {F(mc_sbc['coverage_90']):.2f} over {12*int(mc_sbc['n_datasets'])} ranks; this guards only against gross miscalibration.

Table S4 compares the neural posterior with the MCMC posterior. The full comparison is in the project repository (`output/calibration/npe_vs_mcmc/report.md`).

**Table S4. Single round NPE compared with the MCMC posterior.**

| Method | Prior simulations (valid) | Worst absolute mean difference (posterior SD) | Worst 95% endpoint difference (posterior SD) | 95% interval width ratio, mean (range) | Parameters passing all checks | Worst policy median error | Worst joint pass probability error (percentage points) |
| --- | --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(bench_rows)}

Checks: absolute mean difference at most 0.1 SD, 95% interval endpoints within 0.2 SD and width ratio between 0.9 and 1.1. Policy errors compare 1,000 neural posterior states with the {ndraw:,} MCMC states used for the scenarios.

## 2.6 Posterior distribution and model fit

Table S5 gives the posterior of each parameter. For {pd_txt} the posterior standard deviation was more than 0.8 times the prior standard deviation, so these parameters remain mostly informed by the prior.

**Table S5. Nelder Mead best fit, prior centre and posterior median (95% interval) of the log scale parameters.**

| Parameter | Description | Nelder Mead best fit | Prior centre | Posterior median (95% interval) | Posterior SD | Posterior SD / prior SD |
| --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(post_rows)}

Table S6 compares the targets with the best deterministic fit and the 95% posterior predictive interval. {"All 12" if n_ovl == 12 else f"{n_ovl} of 12"} predictive intervals overlapped the target intervals{miss_txt}. Fig S1 shows the fit. At age 60+, the best fit gave a seroprevalence of {100*F(calib[-1]['seroprevalence_hat']):.1f}% against a target of {100*F(calib[-1]['seroprevalence_target']):.1f}%. Without the excess mortality ($h = 1$) the 60+ fit was {100*F(hr1[-1]['seroprevalence_hat']):.1f}%, and over $h$ from 1.11 to 1.77 it was {100*min(F(x[-1]['seroprevalence_hat']) for x in hr_sens):.1f}% to {100*max(F(x[-1]['seroprevalence_hat']) for x in hr_sens):.1f}% (3 starts each).

**Table S6. Calibration targets and model fit, 2015.**

| Age group | Seroprevalence target | Best fit | Posterior 95% interval | Population target | Best fit | Posterior 95% interval |
| --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(cal_rows)}

![Fig S1. Calibration of age specific seroprevalence and prison population in 2015. Target intervals are compared with 95% posterior predictive intervals from the MCMC posterior and, for contrast, single round neural posterior estimation.](FIGURES_DIR/calibration_agepop.png){{width=100%}}

In 2015 the model estimated {stx('population_D')} people in D, {stx('population_J', False)} in J, {stx('population_F', False)} in F and {stx('population_X', False)} in X. The total was {stx('population_ever_PWID', False)} people with a history of injecting drugs. Of these, {stx('active_hcv_ever_PWID', False)} were viraemic.

# 3. Scenarios and costs

## 3.1 Shared history

All scenarios start from the same 2015 equilibrium. There was no treatment from 2016 to 2018. From 2019 the chronic infection treatment rate in J was $\\tau_J = -\\ln(1 - 0.99995)$ per year with no annual cap. People already on treatment completed their course. People who failed treatment or were reinfected could start again. No one in D was treated in any scenario. The no treatment counterfactual kept $\\tau = 0$ in all groups up to 2030.

## 3.2 Community screening quota

From 2026, each year began with a screening plan for F and X. For group $s$ we computed $N_s = u + a + c + r$, which excludes people on treatment. The number screened was $S_s = \\min(Q_s, N_s)$, where $Q_s$ is the annual quota. The number of antibody tests was $S_s$. The number of RNA tests was $S_s (a + c + r)/N_s$. The number of people with chronic infection identified was $K_s = S_s\\, c/N_s$. $K_s$ was the maximum number of treatment starts in group $s$ that year. Identified people started treatment at rate $\\tau_{{FX}} = -\\ln(1 - 0.99)$ per year. The ODE solver tracked cumulative starts and switched off the $c \\to t$ flow in group $s$ once starts reached $K_s$. Unused quota was not moved between F and X or carried to the next year. In the ideal reference, $Q_s$ was infinite, so everyone eligible was screened once a year.

This is an aggregate approximation. The model does not track individual screening history. The same person can be counted as screened in more than one year. Acute infections were not treated.

**Table S7. Scenario definitions.**

| Scenario | Annual tests | F quota per year | X quota per year |
| --- | --- | --- | --- |
| J only | 0 | 0 | 0 |
""" + "\n".join(f"| {label(k)} | {int(F(summ[k]['screening_total'])):,} | {F(summ[k]['f_screening_quota']):,.0f} | {F(summ[k]['x_screening_quota']):,.0f} |" for k in order if summ[k]['scenario_type']=='finite') + f"""
| J + no screening cap | Unlimited | All eligible | All eligible |

## 3.3 Outcomes

New infections are the annual flow from states u and r to state a in D and F, evaluated at the end of each calendar year. HCV deaths are the annual excess deaths in stages 3 and 4. Both are event counts per year. The reduction in 2030 is one minus the ratio of the 2030 value to the 2015 value. We also report these reductions for F and X combined.

## 3.4 Costs

Costs were counted for F, J and X from 2026 to 2030 (Table S8). In F and X, each person screened had one antibody test. Each expected antibody positive person had one RNA test. Each treatment start had one consultation and one full drug course. In J, each antibody positive person entering prison had one antibody test. Each treatment start in J had one RNA test, one consultation and one drug course. Costs were counted in the year they occurred and discounted as $C_y / 1.03^{{(y-2026)}}$.

The drug price was USD 2,577 per week for 12 weeks [@wong2019]. This is USD 30,924 per course. We converted it at the 2025 average rate of 1.3066 Singapore dollars per US dollar [@sgd_rate]. Prices were not inflated. These are model prices and not Singapore payer costs. Follow up visits, SVR12 tests, fibrosis assessment and programme staff were not included.

**Table S8. Unit costs.**

| Item | Unit cost | Unit | Source price |
| --- | --- | --- | --- |
{chr(10).join(cost_rows)}

Sources: drug price [@wong2019], antibody test and consultation [@nhg_dsc], RNA test [@parkway_rna].

# 4. Full results

## 4.1 Targets in all scenarios

**Table S9. Reductions in 2030 compared with 2015, all 17 scenarios.**

| Scenario | Incidence reduction | Mortality reduction | Draws meeting both targets | F and X incidence reduction | F and X mortality reduction |
| --- | --- | --- | --- | --- | --- |
{chr(10).join(s6row(k) for k in order)}

Values are medians (95% UI) across {ndraw:,} posterior draws.

![Fig S2. Trajectories of new infections and HCV deaths relative to 2015 (left) and reductions in 2030 (right) in all 17 scenarios. Medians across posterior draws.](OUTPUT_DIR/scenarios_2026/scenario_target_trajectories.png){{width=100%}}

![Fig S3. Attainment of the 2030 targets by scenario. Cells show median reductions. Joint status requires both targets.](OUTPUT_DIR/scenarios_2026/scenario_target_heatmap.png){{width=90%}}

## 4.2 Services and costs

**Table S10. Services and discounted costs, 2026 to 2030.**

| Scenario | F and X people screened | F and X RNA tests | F and X treatment starts | J treatment starts | Total discounted cost | Incremental cost vs J only |
| --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(s7row(k) for k in order)}

Screening counts may include the same person in more than one year. Treatment starts count courses.

![Fig S4. Annual treatment starts and screening derived treatment quotas in F and X, by capacity and allocation. Medians across posterior draws.](OUTPUT_DIR/scenarios_2026/scenario_treatment_starts.png){{width=100%}}

![Fig S5. Share of the five year screening quota used in F and X. The F group is small, so a large F share leaves part of its quota unused.](OUTPUT_DIR/scenarios_2026/paper_fig4_f_pool.png){{width=90%}}

![Fig S6. Cost composition of the prison only scenario and the scenarios that met both targets. Shares use the median cost of each item.](OUTPUT_DIR/scenarios_2026/scenario_cost_composition_pies.png){{width=100%}}

![Fig S7. Discounted costs and treatment starts by group, 2026 to 2030. Segments show the median for each group. Labels show the median total.](OUTPUT_DIR/scenarios_2026/scenario_cumulative_costs.png){{width=100%}}

## 4.3 Cost sensitivity

Table S11 scales all unit costs by 0.5, 1 and 1.5. Epidemiological results do not change.

**Table S11. Total discounted cost (S$ million) under unit cost multipliers.**

| Scenario | 0.5 times | 1 times | 1.5 times |
| --- | --- | --- | --- |
{chr(10).join(s9row(k) for k in order)}

# 5. Uncertainty analysis and software

We sampled {ndraw:,} parameter sets at random without replacement from the {len(ref_draws):,} MCMC posterior draws, using sampling seed 20261005. For each set we ran the full pipeline. This included the 150 year equilibrium, the equilibrium check, the shared history from 2016 to 2025, the no treatment counterfactual and all 17 scenarios to 2030. All {ndraw:,} draws passed the equilibrium check. Derived quantities such as reductions, costs and incremental costs were computed within each draw before summarising. We report the median and the 2.5th and 97.5th percentiles. A scenario met a target if its median reduction met the threshold. We also report the share of draws meeting both targets.

The model was solved in C++ through Rcpp [@rcpp2011] with a daily time step. MCMC ran in R. Neural posterior estimation used Python, PyTorch and sbi [@sbi2020]. The MCMC posterior took about 3.8 hours on four processor cores. The scenario analysis of {ndraw:,} draws took about 17 minutes on 14 cores. The project test suite (63 checks) passed.

# 6. Further limitations

- **Population scope.** The model includes only people with a history of injecting drugs. It does not include other routes of infection or the general population. National elimination needs national data [@who2023validation; @cda_hcv].
- **Prison treatment history.** The 2019 prison treatment rate is an assumption. The model was not fitted to observed treatment numbers. In a Singapore prison cohort, 681 people were treated in 2019 [@prison2021]. The model gave {col(an('J_only',2026,'J'),'treatment_starts',d=0)} prison treatment starts in 2026 because most prevalent cases had been treated by then.
- **Linkage to care.** Screening, confirmation and treatment start were assumed to be complete. Real linkage is lower [@hsiang2022; @clink2022].
- **Group X.** People in X cannot relapse or be reinfected. Relapse would add transmission and reinfection.
- **Event counts.** Outcomes are annual event counts in the modelled population. They are not rates per person.
- **Costs.** Prices are public list prices and an older drug price. Many programme costs are not included.
- **Excess mortality.** The hazard ratio for non-liver deaths comes from a general population cohort in Taiwan, not from PWID [@lee2012]. A cohort in Naples found 1.11 (0.82 to 1.49) [@piselli2021]. In a PWID cohort, all-cause mortality was lower with undetectable RNA (0.54), but the comparison group was self selected [@cepeda2022].
- **Prison rates.** The first arrest flow comes from a steady state balance of the 2015 prison population. The same prison counts are calibration targets, so the population targets carry little extra information.
- **Uncertainty.** The UIs reflect calibrated parameters only. Fixed natural history parameters, uptake and price were not varied.
- **Prior.** The prior was centred on a deterministic fit to the same 12 observations. For {pd_txt} the posterior standard deviation was {100*min(r for r in sd_ratio if r > 0.8):.0f}% to {100*max(sd_ratio):.0f}% of the prior standard deviation, so the 0.20 prior width partly sets their uncertainty.

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
| 20 | Characterising uncertainty | Methods, Uncertainty analysis; S1 Text 2.4, 2.5, 4.3, 5 |
| 21 | Approach to engagement with patients and others affected by the study | No engagement |
| 22 | Study parameters | S1 Text Tables S1 to S8 |
| 23 | Summary of main results | Results; Table 1 |
| 24 | Effect of uncertainty | Results; S1 Text 2.5, 4.3 |
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
