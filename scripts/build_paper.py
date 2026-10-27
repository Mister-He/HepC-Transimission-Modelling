#!/usr/bin/env python3
# =============================================================================
# build_paper.py — builds docs/main.docx and docs/Supplementary Information.docx from
# generated results; every number in the text is read from the result files.
#
# Usage: python3 scripts/build_paper.py [ROOT]. Needs output/scenarios_2026/, output/calibration/,
# docs/references.bib, docs/plos.csl and pandoc. Edit prose here, not in Word.
# =============================================================================
import csv, json, math, re, sys, os, subprocess, tempfile
ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.getcwd())
TMP = tempfile.mkdtemp()
S = os.path.join(ROOT, "output/scenarios_2026")
def rd(p): return list(csv.DictReader(open(os.path.join(ROOT, p))))
SETUP = open(os.path.join(ROOT, "scripts/model/setup.R")).read()
def setup_val(name):  # numeric value(s) of a `name = ...` entry in setup.R
    v = re.search(rf"^\s*{name}\s*=\s*(c\([^)]*\)|[0-9.eE+-]+)", SETUP, re.M).group(1)
    return [float(x) for x in re.findall(r"[0-9.]+(?:[eE][+-]?[0-9]+)?", v)]
PI, SMR_V, SMR_NV, PSI_DC = (setup_val(n)[0] for n in ("pi_recid", "smr_viraemic", "smr_nonviraemic", "psi_DC"))
LAMBDA3 = setup_val("lambda3")
TG = json.load(open(os.path.join(ROOT, "data/calibration_targets.json")))
summ = {r["scenario_id"]: r for r in rd("output/scenarios_2026/scenario_2030_summary.csv")}
traj = rd("output/scenarios_2026/scenario_trajectories.csv")
nt = rd("output/scenarios_2026/no_treatment_trajectory.csv")
ann = rd("output/scenarios_2026/scenario_annual_services_costs.csv")
treat = rd("output/scenarios_2026/scenario_treatment_summary.csv")
cost_ev = rd("output/scenarios_2026/cost_evidence.csv")
F = float
WORDS = 'zero one two three four five six seven eight nine ten eleven twelve thirteen'.split()  # small counts in prose

def ui(m, lo, hi, scale=1, d=1, pre="", suf=""):
    f = lambda v: f"{pre}{F(v)*scale:,.{d}f}{suf}"
    return f"{f(m)} (95% UI: {f(lo)}–{f(hi)})"
def col(row, c, **k): return ui(row[c], row[c + "_lo"], row[c + "_hi"], **k)
def pct(row, c, d=1): return col(row, c, scale=100, d=d, suf="%")
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
finite = [k for k, v in summ.items() if v["scenario_type"] == "finite"]
n_inc = sum(1 for k in finite if summ[k]["who_incidence_pass"] == "TRUE")
n_mort = sum(1 for k in finite if summ[k]["who_mortality_pass"] == "TRUE")
ndraw = int(F(j["n_posterior_draws"]))
j_mort = 100 * F(j["mortality_reduction"])
d_inc_share = 1 - F(summ["SInf"]["fx_incidence"]) / F(summ["SInf"]["incidence"])  # D share of 2030 infections
# The narrative below is written for these results; stop the build if they change.
mort_finite = [k for k in finite if summ[k]["who_mortality_pass"] == "TRUE"]
assert n_inc == 15 and j_inc >= 80 and j_mort < 65, "results changed: revise the narrative"
assert set(mort_finite) == {"S5000_F25", "S8000_F25", "S8000_F50"}, "mortality-passing scenarios changed: revise the narrative"
def scen_phrase(ids):  # "5,000 tests a year with 25% to F and 8,000 with 25% or 50% to F" style list
    caps = sorted({int(F(summ[k]["screening_total"])) for k in ids})
    parts = [f"{c:,}" + (" tests a year" if i == 0 else "") + " with " +
             " or ".join(f"{x}%" for x in sorted(round(100 * F(summ[k]["f_share"])) for k in ids if int(F(summ[k]["screening_total"])) == c)) +
             " to F" for i, c in enumerate(caps)]
    return " and ".join(parts)
assert summ["SInf"]["who_mortality_pass"] == "TRUE" and j_inc - j_mort > 20, "results changed: revise the narrative"
strata = {r["measure"]: r for r in rd("output/calibration/mcmc_reference/model_prediction_2015_strata.csv")}
def mprob(sid): return f"{100*F(summ[sid]['who_mortality_pass_prob']):.0f}%"
def iprob(sid): return f"{100*F(summ[sid]['who_incidence_pass_prob']):.0f}%"
best_mort = max(finite, key=lambda k: F(summ[k]["mortality_reduction"]))
assert best_mort == "S8000_F25", "best finite scenario changed: revise the narrative"
ref_cfg = json.load(open(os.path.join(ROOT, "output/calibration/mcmc_reference/mcmc_run_config.json")))
PRIOR_SDS = ref_cfg["prior_sd"]
c_post = next(r for r in rd("output/calibration/mcmc_reference/mcmc_parameter_summary.csv") if r["parameter"] == "theta13")
c_med, c_lo, c_hi = (math.exp(F(c_post[k])) for k in ("median", "q025", "q975"))
C_TXT = f"{c_med:.2f} (95% interval {c_lo:.2f} to {c_hi:.2f})"
def trunc_sd(m, sd, lo, hi):  # SD of Normal(m, sd) truncated to [lo, hi]
    phi = lambda z: math.exp(-z * z / 2) / math.sqrt(2 * math.pi); Phi = lambda z: 0.5 * (1 + math.erf(z / math.sqrt(2)))
    a, b = (lo - m) / sd, (hi - m) / sd; Z = Phi(b) - Phi(a)
    return sd * math.sqrt(1 + (a * phi(a) - b * phi(b)) / Z - ((phi(a) - phi(b)) / Z) ** 2)
PRIOR_TSD = [trunc_sd(m, sd, lo, hi) for m, sd, lo, hi in zip(ref_cfg["prior_mean"], PRIOR_SDS, ref_cfg["theta_lo"], ref_cfg["theta_hi"])]
c_sd_ratio = F(c_post["sd"]) / PRIOR_TSD[12]
def drug_share(sid):
    rr = [x for x in ann if x["scenario_id"] == sid]
    return sum(F(x["discounted_cost_daa_sgd"]) for x in rr) / sum(F(x["discounted_cost_total_sgd"]) for x in rr)
ds_txt = f"{100*drug_share('S5000_F25'):.0f}%"
CMP = "output/calibration/npe_vs_mcmc"
acc_sum = {r["run"]: r for r in rd(CMP + "/accuracy_summary.csv")}
pol_cmp = rd(CMP + "/policy_comparison.csv")
npe_l = acc_sum["npe_lhs"]
npe_w = F(npe_l["mean_width95_ratio"])
npe_pol_med = max(F(r["error"]) for r in pol_cmp if r["run"] == "npe_lhs" and r["error_unit"] == "relative"
                  and "cost" not in r["outcome"])
npe_sup = dict(l.split() for l in open(os.path.join(ROOT, CMP, "npe_support.txt")))
npe_sup = [int(npe_sup["npe_draws_checked"]), int(npe_sup["npe_draws_equilibrium_pass"])]
npe_jp = F(next(r for r in pol_cmp if r["scenario"] == "S5000_F25" and r["outcome"] == "joint_pass_probability")["error"])
npe_ss = json.load(open(os.path.join(ROOT, "output/calibration/npe_lhs/simulation_summary.json")))
npe_cfg = json.load(open(os.path.join(ROOT, "output/calibration/npe_lhs/npe_config.json")))
npe_sbc = [r for r in rd("output/calibration/npe_lhs/sbc_summary.csv")]
mc_sbc = {r["parameter"]: r for r in rd(CMP + "/mcmc_sbc/mcmc_sbc_summary.csv")}["all"]
def rng(k, d=2): v = [F(r[k]) for r in npe_sbc]; return f"{min(v):.{d}f} to {max(v):.{d}f}"
sbc_uniform = all(r["rank_shape"] == "uniform" for r in npe_sbc)
w_txt = (f"on average {100*(1-npe_w):.0f}% narrower than" if npe_w < 0.995 else
         f"on average {100*(npe_w-1):.0f}% wider than" if npe_w > 1.005 else "about as wide as")
def stx(m, full=True):
    r = strata[m]
    return (f"{F(r['median']):,.0f} (95% interval: {F(r['lower_95']):,.0f}–{F(r['upper_95']):,.0f})" if full
            else f"{F(r['median']):,.0f} ({F(r['lower_95']):,.0f}–{F(r['upper_95']):,.0f})")
assert max(F(summ[k]['mortality_reduction']) for k in ['S600_F25','S600_F50','S600_F75']) - F(j['mortality_reduction']) < 0.05, "600/yr gain is no longer under 5 points"
fq_8075 = fx_tot("S8000_F75", "F", "actual_screened") / fx_tot("S8000_F75", "F", "screening_quota")
f_chronic_share = F(cF["chronic_stock_start"]) / f_pool  # chronic share of people eligible for screening in F
x_chronic_share = F(cX["chronic_stock_start"]) / F(cX["screening_eligible"])
J_OBS_2015 = sum(int(r["n_in_prison"]) for r in rd("output/prison_2015_counts.csv"))  # screened PWID in prison, 2015 (preprocess.R); not a target
nm_sol = rd("output/calibration/r_20261027_agepop/solutions.csv")
nm_best = min(nm_sol, key=lambda r: F(r["objective"]))
max_abs_err = F(nm_best["prevalence_max_abs_err"])
REF_DIAG = [r for r in rd("output/calibration/mcmc_reference/mcmc_diagnostics.csv") if r["subset"] == "all_retained"]
# Claims the text makes; stop the build if any stops holding.
assert 50 <= j_mort < 60, "J-only mortality reduction is no longer 'about half'"
robust = [k for k in mort_finite if F(summ[k]["joint_pass_probability"]) >= 0.95]
assert robust == ["S8000_F25"], "robust passing set changed: revise the narrative"
assert 0.9 <= F(summ["S8000_F50"]["joint_pass_probability"]) < 0.95 and 0.7 <= F(summ["S5000_F25"]["joint_pass_probability"]) < 0.9, "pass robustness changed: revise the narrative"
assert min(mort_finite, key=lambda k: F(summ[k]["total_discounted_cost_2026_2030_sgd"])) == "S5000_F25", "cheapest median pass changed"
assert summ["S5000_F50"]["who_mortality_pass"] != "TRUE" and F(summ["S5000_F50"]["joint_pass_probability"]) < 0.2, "5,000/yr F50% changed: revise the narrative"
assert summ["S8000_F75"]["who_mortality_pass"] != "TRUE" and summ["S5000_F75"]["who_mortality_pass"] != "TRUE", "75% F scenarios changed: revise the narrative"
assert fq_8075 > 0.999, "F quota is no longer used in full at 8,000/yr F75%: revise the narrative"
for cap in (600, 1000, 2000, 5000, 8000):
    ids = [f"S{cap}_F{x}" for x in (25, 50, 75)]
    m = [F(summ[k]["mortality_reduction"]) for k in ids]; i = [F(summ[k]["incidence_reduction"]) for k in ids]
    assert m == sorted(m, reverse=True) and i == sorted(i), f"F/X ordering changed at {cap}/yr"
assert x_chronic_share > f_chronic_share, "X no longer has the higher chronic share"
assert xshare > 0.8 and fx_death_share_2026 > 0.9, "X/FX shares changed: revise 'most'/'nearly all'"
assert max(F(d["rhat"]) for d in REF_DIAG) < 1.01, "MCMC R-hat is no longer below 1.01"
j_model_2015 = F(strata["population_J"]["median"])
assert j_model_2015 > J_OBS_2015, "modelled prison population is no longer larger than observed"
assert (SMR_V, SMR_NV) == (9.52, 4.73), "SMRs changed: update the confidence intervals in the text"
SCEN_SEED = re.search(r"sampling seed (\d+)", open(os.path.join(S, "calibration_source.txt")).read()).group(1)
ovl = rd("figures/calibration_interval_overlap.csv")
n_ovl = sum(o["MCMC_overlap"] == "TRUE" for o in ovl)
ent = [math.exp(F(r["median"])) for r in rd("output/calibration/mcmc_reference/mcmc_parameter_summary.csv") if r["parameter"] in [f"theta{i}" for i in range(7, 12)]]
entry_lo, entry_hi = min(ent), max(ent)
assert len({r["n_valid"] for r in pol_cmp}) == 1, "policy comparison rows use different NPE states"
pol_n = f"{int(F(pol_cmp[0]['n_valid'])):,} valid neural posterior states (of {int(F(pol_cmp[0]['n_requested'])):,})"
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

We built an age structured HCV model with four groups: active PWID never imprisoned (D), people in prison (J), and released people who still inject (F) or have stopped (X). We calibrated it with Bayesian methods to age specific HCV seroprevalence in prison in 2015, the official 2015 counts of new drug users by age (taken as first arrests) and a published estimate of the number of people who inject drugs. Entry into injecting and the first arrest rate were calibrated and held at their 2015 values in every year. We assumed treatment in prison from 2019. From 2026 we compared 15 community screening scenarios of 600 to 8,000 tests a year, with 25%, 50% or 75% for F and the rest for X. Results are medians with 95% uncertainty intervals (UI) across {ndraw:,} posterior draws. Without treatment, new infections in 2030 were {pct(nt30, 'incidence_ratio_2015')} of the 2015 level. With prison treatment alone, new infections fell by {red('J_only', 'incidence_reduction')} by 2030, meeting the incidence target in {iprob('J_only')} of draws, but HCV deaths fell by only {red('J_only', 'mortality_reduction')}. In 2026, {100*xshare:.0f}% of chronic infections among released people were in X, and {100*fx_death_share_2026:.0f}% of HCV deaths occurred in F and X. Every screening scenario met the incidence target. The mortality target was met at the median only at {scen_phrase(mort_finite)}. Only 8,000 tests a year with 25% to F met both targets in at least 95% of draws ({prob('S8000_F25')}); 8,000 with 50% did so in {prob('S8000_F50')} and 5,000 with 25% in {prob('S5000_F25')}. At every capacity, giving more tests to X gave larger falls in deaths. The model covers only PWID, assumes ideal linkage to treatment and was not fitted to the size of the prison population.

## Conclusions

In this model, prison treatment alone cuts new HCV infections by more than four fifths by 2030 but cuts HCV deaths by only about half. Most remaining deaths occur among released people, and most chronic infection after release is in people who no longer inject. Community screening met both targets at the median from 5,000 tests a year, with three quarters of tests for released people who have stopped injecting; it did so reliably only at 8,000 tests a year. That is {WORDS[round(5000/600)]} to {WORDS[round(8000/600)]} times the planned capacity of the current programme.

# Author summary

## Why was this study done?

- Singapore aims to eliminate hepatitis C as a public health threat. Most infections are among people who inject drugs, and many of them pass through prison.
- Treatment in prison expanded in 2019. It is not known whether this alone can meet the WHO 2030 targets for new infections and deaths.
- Community services have limited screening capacity. It is not clear how this capacity should be shared between people who still inject after release and people who have stopped.

## What did the researchers do and find?

- We built a model of hepatitis C transmission and liver disease among people with a history of injecting drugs in Singapore. We calibrated it to prison seroprevalence, counts of new drug users and an estimate of the number of people who inject drugs, and carried uncertainty in the calibrated parameters through to the results.
- Prison treatment alone reduced new infections by {j_inc:.0f}% and hepatitis C deaths by {j_mort:.0f}% by 2030. It met the WHO target of 80% for infections but not the target of 65% for deaths.
- By 2026, nearly all HCV deaths occurred among released people, and most of their chronic infections were in those who had stopped injecting. Community screening met both targets at the median with 5,000 tests a year when three quarters went to people who had stopped injecting, and with 8,000 tests a year when at least half did.

## What do these findings mean?

- Prison treatment works well against transmission. Deaths come from liver disease that built up over decades and fall more slowly.
- Reaching the mortality target reliably needs community screening about {WORDS[round(8000/600)]} times larger than current plans, aimed mainly at people who have stopped injecting.
- These results come from a model with ideal links from testing to treatment, calibrated to one year of data. Real programmes will need data on uptake, referral and cost before they set targets.

# Introduction

Hepatitis C virus (HCV) infection causes cirrhosis, liver cancer and death. Direct acting antivirals (DAAs) cure more than 95% of people treated [@pawlotsky2020]. The World Health Organization (WHO) has set targets for 2030. These are an 80% reduction in new HCV infections and a 65% reduction in HCV deaths compared with 2015 [@who2022strategy; @who2023validation]. Many countries are not on track [@dore2021].

In Singapore, HCV is concentrated among people who inject drugs (PWID) [@chaillon2021; @cda_hcv]. A network scale-up survey in 2017 estimated that about 11,000 people in Singapore inject drugs (95% CI 6,500 to 17,000) [@teo2019]. Many of them move between the community, halfway houses and prison [@chaillon2021]. Prison is a key point of contact with health services. Singapore expanded DAA treatment in prison in 2019. In that year 681 people were treated in one prison cohort, with a cure rate above 99% [@prison2021]. Sofosbuvir and velpatasvir were later subsidised for wider use [@tan2023; @ace_daa]. Community programmes have also started. One programme plans to screen 600 high risk people a year [@moh2025ett; @moh2024ett]. A national HCV programme now coordinates this work [@cda_programme].

Treating people who inject drugs can prevent onward transmission [@martin2013; @lim2018]. Deaths are a different matter. They come from cirrhosis and liver cancer, which develop over decades [@thein2008]. Cure lowers but does not remove the risk of liver cancer in people with cirrhosis [@morgan2013; @vandermeer2012]. A programme can therefore cut infections quickly and still miss the mortality target. A previous Singapore model estimated the treatment needed for elimination [@chaillon2021]. It did not separate released people who still inject from those who have stopped. These two groups play different roles. The first drives transmission. The second may hold most of the burden of advanced liver disease.

We used a calibrated model to answer three questions. First, how much does prison treatment from 2019 reduce new infections and deaths by 2030? Second, how much community screening is needed to meet the WHO targets? Third, how should a fixed screening capacity be shared between released people who still inject and those who have stopped?

# Methods

## Study design and reporting

This is a modelling study. It used published and aggregate data only, so ethics approval was not required. The analysis was not prospectively registered. We report it following the CHEERS 2022 checklist where it applies [@cheers2022] (S1 Text, Section 7).

## Model structure

We built a deterministic compartmental model of HCV transmission and liver disease (Fig 1). The population includes everyone with a history of injecting drugs in Singapore. It is split into four groups. Active PWID who have never been imprisoned form group D. People in prison form group J. Released people who still inject form group F. Released people who have stopped injecting form group X. People enter D through new drug use. They move to J on first arrest. On release they move to F or X. People in F can be arrested again. The model has six age groups and four liver stages. The stages are no cirrhosis, compensated cirrhosis, decompensated cirrhosis and hepatocellular carcinoma. Each person is susceptible, acutely infected, chronically infected, on treatment, or recovered.

Transmission occurs only among active injectors in D and F. People in J and X are not exposed. Recovery does not give immunity, so recovered people in D and F can be reinfected. Liver disease progresses at rates taken from the literature [@thein2008; @lim2018; @kanwal2014]. Cure stops progression before cirrhosis and reduces it after cirrhosis [@morgan2013; @vandermeer2012]. Treatment is 12 weeks of sofosbuvir and velpatasvir. Cure rates depend on genotype and liver stage [@pawlotsky2020; @xiao2025]. People with HCV have higher death rates from decompensated cirrhosis and liver cancer [@lim2018]. Other deaths occur at the Singapore resident death rate multiplied by standardised mortality ratios for people who use drugs: {SMR_V:.2f} for people with viraemic infection and {SMR_NV:.2f} for everyone else [@vansanten2014]. The viraemic share therefore changes over time, and cure lowers a person's death rate. Release and re-arrest rates came from Singapore Prison Service records. These records do not observe people who have never been arrested, so the first arrest rate was a calibrated multiple $c$ of the re-arrest rate, and entry into injecting was a calibrated multiple of the official counts of new drug users. S1 Text gives all equations, parameter values and sources.

## Data and calibration

We fitted the model to three sets of 2015 data. The first was HCV seroprevalence among people in prison in six age groups, from the Singapore Prison Service, with a beta likelihood centred on the model value and a fixed precision of {TG['beta_precision']}. The second was the official number of new drug users by age in 2015 [@singstat_drug], which we compared with the modelled number of first arrests from D, with a negative binomial likelihood (size {TG['first_arrests_nb_size']}). The third was the network scale-up estimate of {TG['community_pwid']:,} people who inject drugs in Singapore in 2017 [@teo2019], which we compared with the modelled number of active injectors outside prison (D and F), with a log normal likelihood whose standard deviation matched the published 95% interval. The size of the prison population was not fitted. We estimated 13 parameters. Six scaled the contact rates by age. Six scaled the entry into injecting by age, as multiples of at least 1 of the official counts of new drug users. The last, $c$, multiplied the re-arrest rate to give the first arrest rate, and was bounded below at 0.5. Each model run started in model year 0 and ran for 150 years to reach equilibrium. The end of that period represents 2015, and with constant entry and no treatment before 2019 it also stands for the 2017 survey year. Entry was held at its 2015 value from 2016 to 2030, so the projections show the effect of treatment and screening, not of changes in drug use.

The prior did not depend on the data. Each log contact scale had a normal prior with mean 0 and standard deviation 1.5. Each log entry scale had a half normal prior with standard deviation 1, so entry was at least the official count. The log of $c$ had a normal prior with mean 0 and standard deviation 0.5, centred on equal first and re-arrest rates. All priors were truncated at the bounds. A multi start Nelder Mead search found the posterior mode. Parameter sets that did not reach equilibrium were excluded. We sampled the posterior with an adaptive random walk Metropolis algorithm [@haario2001]. We ran {WORDS[ref_cfg['n_chains']]} chains from independent prior draws, with {ref_cfg['max_evals_per_chain']:,} model runs each. All parameters reached rank normalised R-hat below 1.01 [@vehtari2021]. We also trained single round neural posterior estimation [@greenberg2019; @sbi2020] on {npe_ss['n_attempted']:,} prior simulations drawn by Latin hypercube sampling and compared it with MCMC. S1 Text, Section 2 gives details.

## Treatment and screening scenarios

All scenarios share the same history up to 2025. There was no treatment from 2016 to 2018. From 2019, people with chronic infection in prison started treatment at a high rate. This rate corresponds to an annual probability of 0.99995 when there are no competing events. We call this the prison only scenario (J only). It is a model assumption, not an observed treatment volume. We also ran a counterfactual with no treatment at any time.

From 2026, we added community screening for F and X. Each year a fixed number of screening tests was available. We tested 600, 1,000, 2,000, 5,000 and 8,000 tests a year. The 600 level matches the target of the existing community programme [@moh2025ett]. Each total was split between F and X in three ways. F received 25%, 50% or 75% and X received the rest. Unused tests in one group were not passed to the other. This gave 15 scenarios. We also ran an ideal reference in which everyone eligible in F and X was screened every year.

Screening worked as follows. At the start of each year, we counted the people in F and X who were not already on treatment. The number screened was the smaller of the quota and this eligible count. The share of those screened who had chronic infection set the maximum number of treatment starts that year. Identified people then started treatment at a rate equal to an annual probability of 0.99. The model ensured that starts never exceeded the maximum. No one in D was treated in any scenario.

## Outcomes

The two main outcomes were the reductions in 2030 compared with 2015 in annual new infections and annual HCV deaths. New infections include reinfections. HCV deaths are excess deaths from decompensated cirrhosis and liver cancer, plus the excess non-liver deaths of people with viraemic infection over the non-viraemic death rate. These are event counts per year in the modelled population. They are not rates per person. A scenario met a target if its median reduction was at least 80% for infections or 65% for deaths. We also report the share of posterior draws in which each target, and both together, were met. Secondary outcomes were screening tests, treatment starts and costs from 2026 to 2030.

Costs covered antibody tests, RNA tests, one consultation per treatment start and a full course of sofosbuvir and velpatasvir. The drug cost was S$40,405 per course [@wong2019; @sgd_rate]. Test and consultation prices came from public price lists [@nhg_dsc; @parkway_rna]. These are model prices. They are not Singapore payer costs. Costs were discounted at 3% a year to 2026. Treatment starts were not discounted. S1 Text, Section 3 lists all cost inputs.

## Uncertainty analysis

We drew {ndraw:,} parameter sets at random from the MCMC posterior. For each draw we ran the full analysis. This covered the 2015 equilibrium, the shared history to 2025 and all scenarios to 2030. All draws passed the equilibrium check. We report the median and the 2.5th and 97.5th percentiles across draws as the 95% uncertainty interval (UI). The UI therefore reflects uncertainty in the calibrated parameters. It does not reflect uncertainty in uptake, linkage to care or prices. We used R with C++ code through Rcpp [@rcpp2011].

# Results

## Calibration

The model reproduced the age specific seroprevalence targets; the best fit was within {100*max_abs_err:.1f} percentage points of every target. All {n_ovl} 95% posterior predictive intervals overlapped the target intervals (S1 Text, Section 2). Modelled first arrests matched the counts of new drug users at ages under 30 but were higher at older ages, because people who start injecting young age before their first arrest and entry could not fall below the official counts. The model had {stx('population_current_PWID')} active injectors outside prison in 2015, below the survey estimate of {TG['community_pwid']:,} but within its interval. The first arrest multiplier was {C_TXT}, well above its prior median of about 1.05, but its posterior was only modestly narrower than its prior (standard deviation {100*c_sd_ratio:.0f}% of the prior's). Entry was close to the official counts: the median scale was {entry_lo:.2f} to {entry_hi:.2f} below age 60. Because the prison population size was not fitted, the modelled prison population in 2015 was {stx('population_J', False)}, compared with {J_OBS_2015:,} people with a known screening result in prison records. In 2015, the model estimated {stx('population_ever_PWID')} people with a history of injecting drugs. Of these, {stx('population_D', False)} were in D, {stx('population_F', False)} in F and {stx('population_X', False)} in X.

## Prison treatment alone met the incidence target but not the mortality target

Without any treatment, new infections in 2030 were {pct(nt30, 'incidence_ratio_2015')} of the 2015 level and HCV deaths {pct(nt30, 'mortality_ratio_2015')} (Fig 2).

Prison treatment from 2019 changed the course of new infections. By 2025 they were {pct(j25, 'incidence_ratio_2015')} of the 2015 level. By 2030 they had fallen by {red('J_only', 'incidence_reduction')}, which met the 80% target at the median and in {iprob('J_only')} of draws. Deaths fell more slowly. By 2025 they were {pct(j25, 'mortality_ratio_2015')} of the 2015 level. By 2030 the reduction was {red('J_only', 'mortality_reduction')}, well short of the 65% target.

The reason lies in where infection and disease sit after release. In 2026, F held {col(cF, 'chronic_stock_start', d=0)} people with chronic infection. X held {col(cX, 'chronic_stock_start', d=0)}. X thus held {100*xshare:.0f}% of chronic infections among released people. In the J only scenario, {100*fx_death_share_2026:.0f}% of HCV deaths in 2026 occurred in F and X. By 2030 this share was {100*fx_death_share_2030:.0f}%. People in X do not inject and do not transmit. Prison treatment does not reach them after release.

## Community screening can meet both targets

Adding community screening increased both reductions (Table 1). Every scenario met the incidence target at the median, including 600 tests a year ({red('S600_F25', 'incidence_reduction')} with 25% to F; target met in {iprob('S600_F25')} of draws). Screening everyone eligible in F and X every year gave {red('SInf', 'incidence_reduction')}. In that ideal reference, {100*d_inc_share:.0f}% of the remaining new infections in 2030 occurred in D, among people who inject but have never been imprisoned.

The mortality target was met at the median in {n_mort} of 15 scenarios: {scen_phrase(mort_finite)}. Only 8,000 tests a year with 25% to F met both targets in at least 95% of draws ({prob('S8000_F25')}). With 50% to F at that capacity both were met in {prob('S8000_F50')} of draws, and at 5,000 tests a year with 25% to F ({red('S5000_F25', 'mortality_reduction')}) in {prob('S5000_F25')}. At 5,000 tests a year with 50% to F the median reduction was {red('S5000_F50', 'mortality_reduction')}, below the target, and both targets were met in only {prob('S5000_F50')} of draws. With 75% to F the target was missed at both capacities ({med('S5000_F75','mortality_reduction')}% and {med('S8000_F75','mortality_reduction')}%). Screening everyone eligible in F and X every year reduced mortality by {red('SInf', 'mortality_reduction')}. At 600 tests a year, the planned size of the current community programme, the mortality reduction was at most {max(F(summ[k]['mortality_reduction']) for k in ['S600_F25','S600_F50','S600_F75'])*100:.1f}%, a gain of under 5 percentage points over prison treatment alone.

## How screening is split between F and X matters

The split between F and X changed the results in a consistent way (Fig 3). More screening for F gave a larger fall in new infections. More screening for X gave a larger fall in deaths, at every capacity.

The two groups differ in where chronic infection is. In 2026 about {round(f_pool,-1):,.0f} people in F were eligible for screening, of whom {100*f_chronic_share:.0f}% had chronic infection. In X, {100*x_chronic_share:.0f}% of those eligible did. A test given to X therefore found more people to treat. At 8,000 tests a year, the mortality reduction fell from {red('S8000_F25','mortality_reduction')} with 25% to F to {red('S8000_F75','mortality_reduction')} with 75% to F, while the incidence reduction rose from {med('S8000_F25','incidence_reduction')}% to {med('S8000_F75','incidence_reduction')}%. The F quota was used in full in every scenario, so the difference comes from what each test finds, not from unused tests. At 2,000 tests a year, moving from 25% to 75% for F raised the incidence reduction from {med('S2000_F25','incidence_reduction')}% to {med('S2000_F75','incidence_reduction')}% and lowered the mortality reduction from {med('S2000_F25','mortality_reduction')}% to {med('S2000_F75','mortality_reduction')}%. Because every scenario already met the incidence target, the mortality target is the binding one, and it favours X: at 5,000 and 8,000 tests a year, only the splits with 25% or 50% to F met it.

## Treatment volume and cost

Treatment starts in F and X over 2026 to 2030 ranged from {summ['S600_F75']['fx_starts_text']} at 600 tests a year with 75% to F, to {summ['S8000_F25']['fx_starts_text']} at 8,000 tests a year with 25% to F (Table 1). Starts were highest in 2026 and fell each year as the stock of chronic infection was treated. The cheapest scenario that met both targets at the median, 5,000 tests a year with 25% to F, cost {money(summ['S5000_F25'], 'total_discounted_cost_2026_2030_sgd')} over five years; drugs made up {ds_txt} of it. The only scenario that met both targets in at least 95% of draws, 8,000 tests a year with 25% to F, cost {money(summ['S8000_F25'], 'total_discounted_cost_2026_2030_sgd')}. Compared with prison treatment alone, each extra percentage point of mortality reduction cost about S${inc_per_pp('S5000_F25'):.1f} million at 5,000 tests a year and S${inc_per_pp('S8000_F25'):.1f} million at 8,000. Single round neural posterior estimation gave median reductions within {100*npe_pol_med:.1f}% of the MCMC results (S1 Text, Table S4).

**Table 1. HCV outcomes in 2030 and services from 2026 to 2030 in selected scenarios.**

| Scenario | Incidence reduction | Mortality reduction | Draws meeting the mortality target | F and X treatment starts | Discounted cost, 2026 to 2030 (S$ million) |
| --- | --- | --- | --- | --- | --- |
{chr(10).join(t1row(s) for s in table1_ids)}

Reductions compare 2030 with 2015. Values are medians (95% UI) across {ndraw:,} posterior draws. Targets are at least 80% for incidence and at least 65% for mortality. Every scenario met the incidence target at the median; the mortality target was met at {scen_phrase(mort_finite)}. Costs are in 2026 Singapore dollars, discounted at 3% a year, and include prison treatment costs from 2026. Treatment starts count courses, not people. S1 Text, Table S9 gives all 17 scenarios.

# Discussion

In this model, prison treatment from 2019 cut new HCV infections among people who inject drugs in Singapore by {j_inc:.0f}% and HCV deaths by {j_mort:.0f}% by 2030. This met the incidence target but not the mortality target. By 2026, nearly all HCV deaths occurred among released people, and most chronic infections among released people were in those who had stopped injecting. Community screening of 5,000 tests a year or more met the mortality target at the median, but only when most tests went to released people who had stopped injecting, and only 8,000 tests a year met it reliably.

Prison treatment works on transmission among people who pass through prison. Curing them before release removes infection from that part of the injecting network, which matches earlier work on treatment as prevention [@martin2013; @lim2018]. In this model people who inject are first arrested about {c_med:.1f} times as fast as released people are re-arrested, so few people who start injecting go long without being imprisoned, and prison reaches most of the network; even with screening of everyone eligible after release, {100*d_inc_share:.0f}% of the remaining infections in 2030 occurred in people never imprisoned. The first arrest rate rests on treating the official counts of new drug users as first arrests of people who inject. These counts also include people who do not inject, so the true first arrest flow of injectors is likely smaller. A slower first arrest would leave more infections out of reach of prison. Deaths follow a different path. They arise from cirrhosis and liver cancer that built up over many years. Most of this burden sits in people who stopped injecting long ago. They do not transmit and are not reached by prison treatment after release. Cure also does not remove all risk in people who already have cirrhosis [@morgan2013; @vandermeer2012].

Our findings on the F and X split have direct use for planning. A programme that targets people who still inject may seem the best choice for elimination. In this model it improved incidence, which was already on target, at the cost of deaths, which were not. People in X held most chronic infection, so tests given to them found more people to treat. Planners should therefore direct most community screening to people who have stopped injecting. These people may be harder to reach. They may not use drug services and may not see themselves as at risk. Outreach through primary care, addiction services and records of past imprisonment may be needed.

The scale required is large. The only robust scenario, 8,000 tests a year, is {round(8000/600)} times the planned capacity of the current community programme [@moh2025ett], and 5,000 tests a year, which met the targets at the median, is {round(5000/600)} times. All scenarios also assumed that every person identified starts treatment. Real linkage is weaker. In a Singapore trial among halfway house residents, only 12.6% of people with antibodies started treatment [@hsiang2022]. Care navigation and simpler models of care can raise uptake [@clink2022; @eckhardt2022]. Even then, real programmes will likely need more screening than our model suggests.

Our study has several strengths. We propagated parameter uncertainty through every scenario using {ndraw:,} draws from a well converged MCMC posterior. We modelled screening capacity explicitly, so treatment depends on how many people are found. We separated released people by injecting status, which earlier work did not do [@chaillon2021]. Mortality outside the liver depended on each person's infection status, using mortality ratios for people who use drugs [@vansanten2014], so the share of deaths from viraemic people changed as treatment progressed.

The study also has limitations. First, the model covers only people with a history of injecting drugs. It cannot show elimination in the whole population [@who2023validation]. Second, we calibrated the model to one year of data: prison seroprevalence, counts of new drug users and one small survey estimate of the number of people who inject drugs (199 respondents) [@teo2019]. We did not fit the size of the prison population, and the model's 2015 prison population ({stx('population_J', False)}) was much larger than the {J_OBS_2015:,} people with a known screening result in prison records. Absolute numbers of infections, deaths and treatments are therefore uncertain beyond what the UIs show; the reductions relative to 2015 are less affected. We did not fit observed treatment numbers after 2019. Third, we assumed perfect linkage from screening to treatment and no loss to follow up. Fourth, people in X cannot be reinfected or return to injecting. Some do relapse, which would add transmission. Fifth, the first arrest rate and entry into injecting are both informed by the 2015 counts of new drug users from law enforcement. These counts include people who use drugs without injecting, and entry was not allowed to fall below them, so modelled first arrests were too high at older ages. Entry was held constant to 2030; later changes in drug use are not modelled. Sixth, the mortality ratios come from an Amsterdam cohort grouped by antibody status rather than viraemia [@vansanten2014]. HCV deaths include the excess non-liver deaths of viraemic people, which cure removes in the model; if part of that excess reflects other risks of people who use drugs, the mortality benefit of treatment is overstated. Seventh, costs use list prices and an older drug price. Actual payer costs are not public. Eighth, our uncertainty intervals cover model parameters only.

In conclusion, prison treatment has done much to reduce HCV among people who inject drugs in Singapore. In this model it meets the WHO incidence target by 2030 but not the mortality target. Community screening after release can meet both targets from about 5,000 tests a year, and reliably at 8,000, if most tests go to people who have stopped injecting, but that is far larger than current plans.

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

**Fig 2. Annual new HCV infections and HCV deaths relative to 2015.** (A) New infections, including reinfections. (B) HCV deaths: liver deaths and the excess non-liver deaths of viraemic people. Lines show medians and shaded bands show 95% UIs across {ndraw:,} posterior draws. Dashed horizontal lines show the WHO targets of 20% and 35% of the 2015 level. Dotted vertical lines mark 2019 and 2026.

![Fig 2. Annual new HCV infections and HCV deaths relative to 2015.](OUTPUT_DIR/scenarios_2026/paper_fig2_trajectories.png){{width=100%}}

**Fig 3. Reductions in 2030 by screening capacity and share allocated to F.** (A) Incidence. (B) Mortality. Points show medians and bars show 95% UIs. Dashed lines show the 80% and 65% targets. Grey lines show the prison only scenario and the ideal reference.

![Fig 3. Reductions in 2030 by screening capacity and share allocated to F.](OUTPUT_DIR/scenarios_2026/paper_fig3_allocation.png){{width=100%}}

# References
"""

# ---------------- Supplementary Information ----------------
calib = rd("output/calibration/r_20261027_agepop/predictions.csv")
nm_cfg = json.load(open(os.path.join(ROOT, "output/calibration/r_20261027_agepop/run_config.json")))
nm_same = [r for r in nm_sol if F(r["objective"]) - F(nm_best["objective"]) < 1]
assert len(nm_same) == len(nm_sol) - 1, "Nelder-Mead start agreement changed: revise 'the other'"
REF = "output/calibration/mcmc_reference"
ref_diag = REF_DIAG
ref_cs = rd(REF + "/mcmc_chain_summary.csv")
for r in ref_cs:
    r["acceptance"] = r.get("acceptance") or r.get("acceptance_forward_budget")
inflow = rd("data/singstat_new_drug_abusers.csv")
cfg = json.load(open(os.path.join(ROOT, "data/scenario_policy_config.json")))
# Values the prose states in words; stop if the inputs behind them change.
assert (cfg["j_annual_start_probability"], cfg["fx_annual_start_probability"]) == (0.99995, 0.99), "start probabilities changed"
assert (cfg["sof_vel_usd_per_week"], cfg["sof_vel_weeks"], cfg["usd_to_sgd"]) == (2577, 12, 1.3066), "drug price inputs changed"
assert abs(F(next(c for c in cost_ev if "SOF" in c["item"])["base_cost_sgd"]) - 2577 * 12 * 1.3066) < 0.01, "cost_evidence.csv drug price differs from the config"
dv = dict(l.split() for l in open(os.path.join(S, "draw_validity.txt")))
assert dv["draws_valid"] == dv["draws_requested"], "not every scenario draw passed equilibrium: revise 'All draws passed'"


def quant(xs, p):
    xs = sorted(xs); i = (len(xs) - 1) * p; lo = int(i); hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (i - lo)
ref_draws = rd(REF + "/mcmc_draws.csv")
AGES = ["<20", "20–29", "30–39", "40–49", "50–59", "60+"]
theta_names = ([f"Contact scale, age {a}" for a in AGES] + [f"Entry scale, age {a}" for a in AGES] +
               ["First arrest multiplier $c$"])
prior_txt = ["N(0, 1.5)"] * 6 + ["N⁺(0, 1)"] * 6 + ["N(0, 0.5)"]
def qtxt(xs, d=3): return f"{quant(xs,.5):.{d}f} ({quant(xs,.025):.{d}f} to {quant(xs,.975):.{d}f})"
post_rows = []; sd_ratio = []
for i in range(13):
    t = f"theta{i+1}"
    a = [F(r[t]) for r in ref_draws]
    m = sum(a) / len(a); sd = (sum((v - m) ** 2 for v in a) / (len(a) - 1)) ** 0.5
    post_rows.append(f"| θ{i+1} | {theta_names[i]} | {F(nm_best[t]):.3f} | {prior_txt[i]} | {qtxt(a)} | {sd:.3f} | {sd/PRIOR_TSD[i]:.2f} |")
    sd_ratio.append(sd / PRIOR_TSD[i])
prior_dom = [f"θ{i+1}" for i in range(13) if sd_ratio[i] > 0.7]
assert prior_dom == ["θ5", "θ6", "θ13"], "prior-dominated set changed: revise the narrative"
pd_txt = ", ".join(prior_dom[:-1]) + " and " + prior_dom[-1] if len(prior_dom) > 1 else prior_dom[0]
diag_rows = [f"| θ{d['parameter'][5:]} | {F(d['rhat']):.4f} | {F(d['ess_bulk']):,.0f} | {F(d['ess_tail']):,.0f} | {F(d['mcse_mean'])/F(d['sd']):.3f} |" for d in ref_diag]
def bench_row(k):
    a = acc_sum[k]
    pc = [r for r in pol_cmp if r["run"] == k and r["error_unit"] == "relative"]
    pp = [r for r in pol_cmp if r["run"] == k and r["error_unit"] != "relative"]
    name = f"Single round NPE, {npe_cfg['n_ensemble']} network ensemble, Latin hypercube sampling"
    return (f"| {name} | {int(F(a['forward_runs'])):,} ({int(F(a['valid_simulations'])):,}) | {F(a['max_abs_smd']):.2f} | "
            f"{F(a['max_abs_endpoint']):.2f} | {F(a['mean_width95_ratio']):.3f} ({F(a['min_width95_ratio']):.2f} to {F(a['max_width95_ratio']):.2f}) | "
            f"{100*F(a['share_pass_all']):.0f}% | {100*max(F(r['error']) for r in pc):.1f}% | {max(F(r['error']) for r in pp):.1f} |")
bench_rows = [bench_row("npe_lhs")]
miss_txt = "" if n_ovl == len(ovl) else " (not: " + ", ".join(f"{o['measure']} {o['age_group']}" for o in ovl if o["MCMC_overlap"] != "TRUE") + ")"
MEAS = {"seroprevalence": ("Seroprevalence", "seroprevalence_hat", lambda v: f"{100*v:.1f}%"),
        "first_arrests": ("First arrests per year", "first_arrests_hat", lambda v: f"{v:,.0f}"),
        "community_pwid": ("Active injectors outside prison (D + F)", "community_pwid_hat", lambda v: f"{v:,.0f}")}
def cal_row(o):
    name, hat, f = MEAS[o["measure"]]
    c = calib[0] if o["measure"] == "community_pwid" else next(x for x in calib if x["age_group"] == o["age_group"])
    return (f"| {name} | {o['age_group']} | {f(F(o['target']))} ({f(F(o['target_lower']))} to {f(F(o['target_upper']))}) | "
            f"{f(F(c[hat]))} | {f(F(o['MCMC_lower']))} to {f(F(o['MCMC_upper']))} |")
cal_rows = [cal_row(o) for o in ovl]
ages = ["<20", "20-29", "30-39", "40-49", "50-59", "60+"]
years = ["2015"]
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

People enter D in the susceptible state and stage 1 at an age specific rate $s_i\\beta_i$, where $\\beta_i$ is the 2015 count of new drug users and $s_i \\ge 1$ a calibrated entry scale. People in D move to J at the first arrest rate $\\lambda^{{(1)}}_i = c\\,\\lambda^{{(3)}}_i$. People in J are released at rate $\\lambda^{{(2)}}_i$. A share $\\pi = {PI:.2f}$ return to injecting and enter F. The rest enter X. People in F return to J at the re-arrest rate $\\lambda^{{(3)}}_i$. People in X stay in X.

Release and re-arrest rates and $\\pi$ were fitted with Stan to Singapore Prison Service records of everyone held for drug consumption, as survival models of time to release and of time to re-arrest with a share $\\pi$ who ever return. The modelled population is the people with a known HCV screening result, who were assumed to be the PWID; their 2015 seroprevalence is the calibration target. The records cannot tell true first offenders apart and never observe people who have not been arrested, so they do not identify the first arrest rate. It was set to $\\lambda^{{(1)}}_i = c\\,\\lambda^{{(3)}}_i$, with one multiplier $c \\ge 0.5$ for all ages, calibrated (Section 2.3).

## 1.2 Transmission

Transmission occurs only in D and F. The force of infection on age group $i$ is

$$\\gamma_i(t) = q\\; m(t) \\sum_j C_{{ij}} \\frac{{I_j(t)}}{{N_j(t)}}$$

where $q$ is the per event transmission probability [@boelen2014], $C_{{ij}}$ is the calibrated contact matrix, $I_j$ is the number of people in states a and c in D and F of age $j$, and $N_j$ is the number of people in D and F of age $j$. The time multiplier $m(t)$ was fixed at 1. People on treatment are not infectious. Susceptible and recovered people in D and F are infected at rate $\\gamma_i(t)$. People in J and X are not exposed.

## 1.3 Natural history and treatment

Acute infection lasts on average $\\iota_1 = 180$ days [@fasano2024; @tanaka2000]. A share $\\kappa = 0.26$ clears infection spontaneously [@micallef2006]. The rest become chronically infected. Chronic infection progresses through the four liver stages. Progression rates are weighted by genotype, with 78% genotype 3 [@zhao2019]. Genotype 3 has higher progression rates [@kanwal2014]. The effective rate for transition $A \\to B$ is

$$\\tilde p_{{A\\to B}} = \\rho\\, r_3^{{A\\to B}} p^{{A\\to B}} + (1-\\rho)\\, p^{{A\\to B}}.$$

Treatment is sofosbuvir and velpatasvir for $\\iota_2 = 84$ days. People on treatment are cured with probability $\\bar\\alpha_k$ and return to state c otherwise. In the scenarios, decompensated cirrhosis used the cure rate without ribavirin, to match the costed regimen. After cure, people in stage 1 do not progress. People in stage 2 progress at reduced rates [@vandermeer2012; @morgan2013]. People in stage 3 keep the same risk of liver cancer. Treatment lowers excess mortality in stages 3 and 4 by factors $\\psi^{{DC}}$ and $\\psi^{{HCC}}$ [@premkumar2024; @dang2020]. Background mortality varies by age and is multiplied by a standardised mortality ratio for people who use drugs that depends on viraemic status [@vansanten2014]. Excess mortality applies in stages 3 and 4 [@lim2018].

The state transitions within each group, stage and age are

$$a \\xrightarrow{{\\kappa/\\iota_1}} r,\\quad a \\xrightarrow{{(1-\\kappa)/\\iota_1}} c,\\quad c \\xrightarrow{{\\tau}} t,\\quad t \\xrightarrow{{\\bar\\alpha/\\iota_2}} r,\\quad t \\xrightarrow{{(1-\\bar\\alpha)/\\iota_2}} c,\\quad u, r \\xrightarrow{{\\gamma_i(t)}} a$$

As an example, the equation for chronic infection in group F, stage 1, age $i$ is

$$\\frac{{dF_{{c,1,i}}}}{{dt}} = \\pi\\lambda^{{(2)}}_i J_{{c,1,i}} + \\frac{{1-\\kappa}}{{\\iota_1}} F_{{a,1,i}} + \\frac{{1-\\bar\\alpha_1}}{{\\iota_2}} F_{{t,1,i}} - \\left(\\tau_{{F}} + \\lambda^{{(3)}}_i + s_v\\mu_i + \\tilde p_{{1\\to2}}\\right) F_{{c,1,i}}$$

Background mortality is $s\\,\\mu_i$, where $\\mu_i$ is the Singapore resident death rate and $s$ is an all-cause standardised mortality ratio for people who use drugs against the general population, from the Amsterdam Cohort Studies (1985 to 2012) [@vansanten2014]. People in the viraemic states (a, c, t) have $s_v = {SMR_V:.2f}$ (95% CI 8.06 to 11.26; HCV monoinfected). People in states u and r have $s_n = {SMR_NV:.2f}$ (3.47 to 6.45; HCV and HIV uninfected). Recovered people take the uninfected ratio because extrahepatic mortality after loss of viraemia is close to that of uninfected people [@lee2012]. Liver deaths in stages 3 and 4 are added on top. The rate is set by each person's state, so the viraemic share of the population, and with it average mortality, changes over time and falls with treatment. The non-liver deaths of viraemic people above the non-viraemic rate, $(s_v - s_n)\\mu_i$, are counted as HCV deaths together with the liver deaths; the deaths at rate $s_n\\mu_i$ are not.

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
| $s_v$ | Standardised mortality ratio, viraemic (a, c, t) | {SMR_V:.2f} | [@vansanten2014] |
| $s_n$ | Standardised mortality ratio, non-viraemic (u, r) | {SMR_NV:.2f} | [@vansanten2014] |
| $\\mu^{{DC}}$ | Excess mortality, decompensated cirrhosis (per year) | 0.13 | [@lim2018] |
| $\\mu^{{HCC}}$ | Excess mortality, liver cancer (per year) | 0.43 | [@lim2018] |
| $\\psi^{{DC}}$ | Relative excess mortality after treatment, stage 3 | {PSI_DC:.2f} | [@premkumar2024] |
| $\\psi^{{HCC}}$ | Relative excess mortality after treatment, stage 4 | 0.37 | [@dang2020] |
| $\\pi$ | Share returning to injecting after release | {PI:.2f} | Fitted to prison records |
| $\\lambda^{{(3)}}_i$ | Re-arrest rate (per year), by age | {', '.join(f'{x:.3f}' for x in LAMBDA3)} | Fitted to prison records |
| $c$ | First arrest multiplier, $\\lambda^{{(1)}}_i = c\\,\\lambda^{{(3)}}_i$ | {C_TXT} | Calibrated (θ13) |

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

There were 13 calibration targets in three groups, all for 2015 and all independent.

- **Seroprevalence.** HCV seroprevalence among people in prison in six age groups, from Singapore Prison Service screening data. Seroprevalence in the model is the share of people in states a, c, t and r among all people in J. Each observed seroprevalence $y_i$ had the likelihood $y_i \\sim \\mathrm{{Beta}}\\big(p_i\\phi,\\;(1-p_i)\\phi\\big)$, where $p_i$ is the model value and $\\phi = {TG['beta_precision']}$ is a fixed precision (a standard deviation of about 0.045 at $p_i = 0.3$).
- **First arrests.** The official number of new drug users by age in 2015 (Table S2) [@singstat_drug], $n_i$, was compared with the modelled annual flow from D to J, $f_i = c\\,\\lambda^{{(3)}}_i D_i$, as $n_i \\sim \\mathrm{{NegBin}}(\\text{{mean}}\\ f_i,\\ \\text{{size}}\\ {TG['first_arrests_nb_size']})$. The size allows for overdispersion and for the count of 1 at age 60+.
- **Active injectors outside prison.** A network scale-up survey of 199 Singapore residents in {TG['community_pwid_year']} estimated {TG['community_pwid']:,} people who inject drugs (95% CI 6,500 to 17,000) [@teo2019]. It was compared with the modelled $D + F$ as $\\log y \\sim \\mathrm{{N}}(\\log(D+F),\\ {TG['community_pwid_log_sd']})$, with the standard deviation taken from the interval. We counted D and F because the survey asked about people living in Singapore who inject drugs; people in X have stopped injecting, and the authors note that injectors in prison are hard to reach through respondents' contacts. The model value is the 2015 equilibrium, which also stands for 2017 because entry is constant and there was no treatment before 2019.

## 2.2 Entry into injecting

Entry into D was the 2015 number of new drug users by age published by the Singapore Department of Statistics [@singstat_drug] (Table S2), times a calibrated entry scale of at least 1 for each age group. These counts come from law enforcement records. Almost everyone who starts injecting is eventually arrested, so the counts are a lower bound on entry; they also include people who do not inject. Entry set the equilibrium up to 2015 and was held constant from 2016 to 2030, so the projections isolate the effect of treatment and screening.

**Table S2. New drug users by age group, Singapore, 2015 (first arrest target and base entry into D).**

| Year | <20 | 20–29 | 30–39 | 40–49 | 50–59 | 60+ |
| --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(inflow_rows)}

## 2.3 Calibrated parameters and priors

We estimated 13 parameters on the log scale. Parameters θ1 to θ6 multiply the contact rates of the six age groups relative to the contact matrix in the model code. Parameters θ7 to θ12 are the log entry scales $\\log s_i$. Parameter θ13 is $\\log c$, the first arrest multiplier. The prior was fixed before fitting and did not use the data: θ1 to θ6 were each normal with mean 0 and standard deviation 1.5 (a 95% range of 0.05 to 19 times the contact rate); θ7 to θ12 were half normal with standard deviation 1 (median entry about 2.0 times the official count); and θ13 was normal with mean 0 and standard deviation 0.5, centred on equal first and re-arrest rates. The prior was truncated to the bounds: contact scales 0.01 to 1,000, entry scales 1 to 100 and $c$ 0.5 to 100. A deterministic Nelder Mead search with {len(nm_sol)} starts found the posterior mode. The starts were θ = 0 and {len(nm_sol)-1} random perturbations of it (standard deviation {nm_cfg['sd_perturb']} on the log scale). Each start ran for up to {nm_cfg['maxit']:,} iterations and was then restarted from its own solution, up to {nm_cfg['nm_restarts']} times, while the objective still fell by more than 0.001. {WORDS[len(nm_same)].capitalize()} of the {len(nm_sol)} starts reached the same mode (negative log posterior {F(nm_best['objective']):.2f} to {max(F(r['objective']) for r in nm_same):.2f}; best prevalence root mean square error {F(nm_best['prevalence_rmse']):.4f}); the other stopped in a local optimum ({max(F(r['objective']) for r in nm_sol):.1f}). At the mode the entry scales for ages 30 and over were at their lower bound of 1. Every model evaluation ran from model year 0 to 150. The last ten years had to be at equilibrium.

## 2.4 Posterior sampling

The posterior combined the prior and the three likelihood terms (Section 2.1). Parameter sets outside the bounds or whose simulation failed the equilibrium check had zero posterior density. We sampled it with adaptive random walk Metropolis [@haario2001]. {WORDS[ref_cfg['n_chains']].capitalize()} chains started from independent prior draws. The initial proposal covariance was {ref_cfg['init_scale']**2:g} times the prior variance. Adaptation began after {ref_cfg['adapt_start']:,} proposals, was refreshed every {ref_cfg['adapt_every']} proposals and was frozen after warmup. A {100*ref_cfg['safety_p']:.0f}% share of proposals always used the initial covariance. Each chain ran {ref_cfg['max_evals_per_chain']:,} model evaluations, of which {ref_cfg['n_warmup']:,} were warmup, giving {sum(int(F(r['retained_total'])) for r in ref_cs):,} retained iterations. Every {ref_cfg['thin_storage']}th iteration was stored ({len(ref_draws):,} draws), and diagnostics used all retained iterations. Acceptance was {min(F(r['acceptance']) for r in ref_cs):.2f} to {max(F(r['acceptance']) for r in ref_cs):.2f}. All parameters met the convergence checks (Table S3) [@vehtari2021].

**Table S3. Convergence diagnostics of the MCMC posterior.**

| Parameter | Rank normalised R-hat | Bulk ESS | Tail ESS | MCSE of mean / posterior SD |
| --- | --- | --- | --- | --- |
{chr(10).join(diag_rows)}

## 2.5 Comparison with neural posterior estimation

We checked the MCMC posterior against single round neural posterior estimation (NPE) [@greenberg2019], implemented in the sbi package [@sbi2020]. Both methods used the same prior, likelihood, solver and equilibrium restriction. We drew {npe_ss['n_attempted']:,} parameter sets from the truncated prior by Latin hypercube sampling. Of these, {npe_ss['n_valid']:,} simulations were valid. Exclusions were {npe_ss['n_excluded_bounds']:,} draws outside the bounds, {npe_ss['n_excluded_equilibrium']:,} equilibrium failures and {npe_ss['n_excluded_nonfinite']:,} non-finite simulations; excluded draws are dropped, not replaced. The simulated targets received the same observation noise as the likelihood: beta noise on seroprevalence, negative binomial noise on first arrests and log normal noise on the number of active injectors. The estimator was an equal weight ensemble of {npe_cfg['n_ensemble']} masked autoregressive flows, each with 50 hidden features and 5 transforms [@papamakarios2017], trained in one round with the standard loss from different random seeds. Posterior draws outside the bounds were rejected.

Simulation based calibration [@talts2018] on {int(npe_sbc[0]['n_sbc']):,} datasets simulated from independent prior draws gave 95% interval coverage of {rng('coverage_95', 3)} and 90% coverage of {rng('coverage_90', 3)} for each parameter. {"No rank histogram showed systematic skew." if sbc_uniform else "The rank histograms of " + ", ".join(f"θ{r['parameter'][5:]}" for r in npe_sbc if r["rank_shape"] != "uniform") + " were skewed (shifts of about 4% of the rank range, or excess mass in the tails), so the neural posterior failed the predeclared quality gate. It served only as a cross-check; every reported result uses the MCMC posterior."} As a smaller direct check, we reran the MCMC on {mc_sbc['n_datasets']} of these simulated datasets, with one chain of 45,000 evaluations each. Pooled 95% coverage was {F(mc_sbc['coverage_95']):.2f} and 90% coverage {F(mc_sbc['coverage_90']):.2f} over {13*int(mc_sbc['n_datasets'])} ranks; this guards only against gross miscalibration.

At the observed data, however, the neural posterior did not match the MCMC posterior (Table S4). Of {npe_sup[0]:,} random neural posterior draws, {npe_sup[0]-npe_sup[1]:,} ({100*(1-npe_sup[1]/npe_sup[0]):.0f}%) failed the equilibrium check, which the MCMC treats as impossible; they were dropped before the comparison. Its 95% intervals were {w_txt} the MCMC intervals, and {"no parameter" if F(npe_l['share_pass_all']) == 0 else f"{100*F(npe_l['share_pass_all']):.0f}% of parameters"} passed all three checks. Calibration on average across prior datasets does not guarantee accuracy at one dataset; with noisy summaries and a wide prior, the network spends its capacity over the whole prior range. The policy medians were less sensitive: they differed from the MCMC results by at most {100*npe_pol_med:.1f}%, although the share of draws meeting both targets at 5,000 tests a year with 25% to F differed by {npe_jp:.0f} percentage points. We used the MCMC posterior for all results. The full comparison is in the project repository (`output/calibration/npe_vs_mcmc/report.md`).

**Table S4. Single round NPE compared with the MCMC posterior.**

| Method | Prior simulations (valid) | Worst absolute mean difference (posterior SD) | Worst 95% endpoint difference (posterior SD) | 95% interval width ratio, mean (range) | Parameters passing all checks | Worst policy error | Worst joint pass probability error (percentage points) |
| --- | --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(bench_rows)}

Checks: absolute mean difference at most 0.1 SD, 95% interval endpoints within 0.2 SD and width ratio between 0.9 and 1.1. Policy errors compare {pol_n} with the {ndraw:,} MCMC states used for the scenarios; the worst policy error is for discounted cost.

## 2.6 Posterior distribution and model fit

Table S5 gives the posterior of each parameter. For {pd_txt} the posterior standard deviation was more than 0.7 times the prior standard deviation, so these parameters remain mostly informed by the prior. The ratio uses the standard deviation of the truncated prior. For $c$ (θ13) the posterior moved well away from the prior (median {c_med:.2f} against a prior median of about 1.05) even though its spread narrowed little.

**Table S5. Nelder Mead posterior mode, prior and posterior median (95% interval) of the log scale parameters.**

| Parameter | Description | Posterior mode | Prior | Posterior median (95% interval) | Posterior SD | Posterior SD / prior SD |
| --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(post_rows)}

Table S6 compares the targets with the posterior mode and the 95% posterior predictive interval (with observation noise) of each modelled target. {f"All {n_ovl}" if n_ovl == len(ovl) else f"{n_ovl} of {len(ovl)}"} intervals overlapped the target intervals{miss_txt}. Fig S1 shows the fit.

**Table S6. Calibration targets and model fit, 2015.**

| Target | Age group | Observed (95% range of the likelihood) | Posterior mode | Posterior predictive 95% interval |
| --- | --- | --- | --- | --- |
{chr(10).join(cal_rows)}

![Fig S1. Calibration targets, 2015: seroprevalence in prison by age, first arrests by age, and active injectors outside prison. Target values with the 95% range of the likelihood are compared with 95% posterior predictive intervals from the MCMC posterior and, for contrast, single round neural posterior estimation.](FIGURES_DIR/calibration_agepop.png){{width=80%}}

In 2015 the model estimated {stx('population_D')} people in D, {stx('population_J', False)} in J, {stx('population_F', False)} in F and {stx('population_X', False)} in X. The total was {stx('population_ever_PWID', False)} people with a history of injecting drugs. Of these, {stx('active_hcv_ever_PWID', False)} were viraemic. Only D + F ({stx('population_current_PWID', False)}) was fitted, to the survey estimate. The modelled prison population is much larger than the {J_OBS_2015:,} people with a known screening result in prison records.

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

New infections are the annual flow from states u and r to state a in D and F, evaluated at the end of each calendar year. HCV deaths are the annual excess deaths in stages 3 and 4 plus $(s_v - s_n)\\mu_i$ times the number of viraemic people (states a, c and t) in every group, stage and age. Both are event counts per year. The reduction in 2030 is one minus the ratio of the 2030 value to the 2015 value. We also report these reductions for F and X combined.

## 3.4 Costs

Costs were counted for F, J and X from 2026 to 2030 (Table S8). In F and X, each person screened had one antibody test. Each expected antibody positive person had one RNA test. Each treatment start had one consultation and one full drug course. In J the same rule applied to people entering prison: every entrant not on treatment had one antibody test and every expected antibody positive entrant one RNA test. Each treatment start in J had one consultation and one drug course. Costs were counted in the year they occurred and discounted as $C_y / 1.03^{{(y-2026)}}$.

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

![Fig S5. Share of the five year screening quota used in F and X. Every quota was used in full.](OUTPUT_DIR/scenarios_2026/paper_fig4_f_pool.png){{width=90%}}

![Fig S6. Cost composition of the prison only scenario and the scenarios that met both targets at the median. Shares use the median cost of each item.](OUTPUT_DIR/scenarios_2026/scenario_cost_composition_pies.png){{width=100%}}

![Fig S7. Discounted costs and treatment starts by group, 2026 to 2030. Segments show the median for each group. Labels show the median total.](OUTPUT_DIR/scenarios_2026/scenario_cumulative_costs.png){{width=100%}}

## 4.3 Cost sensitivity

Table S11 scales all unit costs by 0.5, 1 and 1.5. Epidemiological results do not change.

**Table S11. Total discounted cost (S$ million) under unit cost multipliers.**

| Scenario | 0.5 times | 1 times | 1.5 times |
| --- | --- | --- | --- |
{chr(10).join(s9row(k) for k in order)}

# 5. Uncertainty analysis and software

We sampled {ndraw:,} parameter sets at random without replacement from the {len(ref_draws):,} MCMC posterior draws, using sampling seed {SCEN_SEED}. For each set we ran the full pipeline. This included the 150 year equilibrium, the equilibrium check, the shared history from 2016 to 2025, the no treatment counterfactual and all 17 scenarios to 2030. All {ndraw:,} draws passed the equilibrium check. Derived quantities such as reductions, costs and incremental costs were computed within each draw before summarising. We report the median and the 2.5th and 97.5th percentiles. A scenario met a target if its median reduction met the threshold. We also report the share of draws meeting both targets.

The model was solved in C++ through Rcpp [@rcpp2011] with a daily time step. MCMC ran in R. Neural posterior estimation used Python, PyTorch and sbi [@sbi2020]. The MCMC posterior took about {F(ref_cfg['method_sec'])/3600:.1f} hours on {WORDS[ref_cfg['cores']]} processor cores. The project test suite passed.

# 6. Further limitations

- **Population scope.** The model includes only people with a history of injecting drugs. It does not include other routes of infection or the general population. National elimination needs national data [@who2023validation; @cda_hcv].
- **Prison treatment history.** The 2019 prison treatment rate is an assumption. The model was not fitted to observed treatment numbers. In a Singapore prison cohort, 681 people were treated in 2019 [@prison2021]. The model gave {col(an('J_only',2026,'J'),'treatment_starts',d=0)} prison treatment starts in 2026 because most prevalent cases had been treated by then.
- **Linkage to care.** Screening, confirmation and treatment start were assumed to be complete. Real linkage is lower [@hsiang2022; @clink2022].
- **Group X.** People in X cannot relapse or be reinfected. Relapse would add transmission and reinfection.
- **Event counts.** Outcomes are annual event counts in the modelled population. They are not rates per person.
- **Costs.** Prices are public list prices and an older drug price. Many programme costs are not included.
- **Mortality ratios.** The standardised mortality ratios come from the Amsterdam Cohort Studies of people who use drugs, 1985 to 2012, against the Dutch general population [@vansanten2014]. They group people by HCV antibody status, so the viraemic ratio (HCV monoinfected, {SMR_V:.2f}) includes people who cleared infection and liver deaths, which the model also counts separately. No cohort reports these ratios by viraemic status.
- **Population size.** The prison population size was not fitted. The modelled 2015 prison population ({stx('population_J', False)}) was much larger than the {J_OBS_2015:,} people with a known screening result in prison records, so absolute counts of infections, deaths, screening and treatment, and costs, are uncertain beyond the UIs.
- **First arrest rate and entry.** The prison records cannot identify the first arrest rate of people who have never been imprisoned. It was a multiple $c \\ge 0.5$ of the re-arrest rate ({C_TXT}), informed by the official counts of new drug users and the survey estimate of people who inject drugs. The counts include people who do not inject, and the survey had 199 respondents [@teo2019]. Entry could not fall below the counts, so modelled first arrests at ages 30 and over were higher than observed. A lower first arrest rate would leave more people in D, who are never reached by prison treatment.
- **Uncertainty.** The UIs reflect calibrated parameters only. Fixed natural history parameters, uptake and price were not varied.
- **Prior.** With 13 targets and 13 parameters, the prior matters. For {pd_txt} the posterior standard deviation was {100*min(r for r in sd_ratio if r > 0.7):.0f}% to {100*max(sd_ratio):.0f}% of the prior standard deviation.

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
| 8 | Perspective | Methods, Outcomes (health service costs only) |
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
DOCS = os.path.join(ROOT, "docs")
for name, text, out in (("main.md", main, "main.docx"), ("supplementary.md", si, "Supplementary Information.docx")):
    text = text.replace("OUTPUT_DIR/", os.path.join(ROOT, "output") + "/").replace("FIGURES_DIR/", os.path.join(ROOT, "figures") + "/")
    md = os.path.join(TMP, name); open(md, "w").write(text)
    subprocess.run(["pandoc", md, "--citeproc", "--bibliography", os.path.join(DOCS, "references.bib"),
                    "--csl", os.path.join(DOCS, "plos.csl"), "-o", os.path.join(DOCS, out)], check=True)
print("Wrote docs/main.docx and docs/Supplementary Information.docx")
