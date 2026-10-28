---
title: "Prison treatment and post-release screening for hepatitis C among people who inject drugs in Singapore: A modelling study"
---

**Short title:** Hepatitis C treatment after prison release in Singapore

Yichen He\*

[Affiliation: to be provided by the authors]

\* Corresponding author. E-mail: [to be provided by the authors]

# Abstract

## Background

Hepatitis C virus (HCV) infection in Singapore is concentrated among people who inject drugs (PWID), many of whom cycle through prison, where treatment expanded from 2019. We asked whether prison treatment alone can meet the World Health Organization 2030 targets (80% fewer new infections and 65% fewer HCV deaths than in 2015), and how limited community screening should be shared after release. [claim:C001] [evidence:E001,E002,E003]

## Methods and findings

We built an age-structured model of PWID never imprisoned (D), in prison (J), and released and still injecting (F) or stopped (X), with prison rates estimated from records of 20,720 people. We calibrated it by Markov chain Monte Carlo to 2015 prison seroprevalence, first arrests and a survey estimate of active injectors, and projected 15 policies screening 600 to 8,000 released people a year from 2026 (25%–75% of tests to F) over 2,000 posterior draws. [claim:C002] [evidence:E004,E005,E006,E007,E001,E008,E009,E010]

Prison treatment alone reduced new infections in 2030 by 84.4% (95% uncertainty interval [UI] 67.9%–90.4%) but HCV deaths by only 54.4% (95% UI 51.6%–57.0%). In 2026, 89% of chronic infections among released people were in X, and F and X accounted for 97% of HCV deaths. All screening policies met the incidence target; the mortality target was met only with 5,000 tests a year and 25% to F (66.7%, 95% UI 63.8%–69.1%; both targets in 87% of draws) or 8,000 tests with 25% or 50% to F. Moving tests from F to X increased the fall in deaths and reduced the fall in infections. The cheapest passing policy cost S$86.2 million (95% UI 70.9–106.4) over 2026–2030, mostly for drugs. Limitations include a modelled prison population (9,741) far above the 3,367 screened people recorded, weakly identified entry parameters, charge-based costs and parameter-only uncertainty. [claim:C003] [evidence:E008,E011,E012,E013,E014,E015,E016]

## Conclusions

In this model, prison treatment is likely to meet the 2030 target for new HCV infections among PWID in Singapore but not the target for deaths, which required screening several thousand released people a year, mostly people who have stopped injecting. [claim:C004] [evidence:E008]

# Author summary

## Why was this study done?

- Most hepatitis C infections in Singapore are among people who inject drugs, and many of them spend time in prison, where treatment expanded from 2019. [claim:C005] [evidence:E001,E017,E018]

- It is unclear whether prison treatment alone can meet the WHO 2030 elimination targets, which ask for large falls in both new infections and deaths. [claim:C006] [evidence:E003]

- Community screening capacity is limited, and it is not known whether it should go to released people who still inject or to those who have stopped.

## What did the researchers do and find?

- We built a model of hepatitis C among people who inject drugs that follows them into and out of prison, and fitted it to prison and national data from 2015. [claim:C007] [evidence:E006]

- In the model, prison treatment alone cut new infections by about 84% by 2030 but deaths by only about 54%, because most deaths occur among released people who have stopped injecting, whom prison treatment does not reach again. [claim:C008] [evidence:E008,E012]

- The death target was met only when 5,000 or more released people were screened each year, with most tests going to people who have stopped injecting. [claim:C009] [evidence:E008]

## What do these findings mean?

- Prison treatment can prevent most new infections, but it cannot by itself prevent most hepatitis C deaths.

- Finding and treating people who have stopped injecting, through community testing outside prison, is needed to reduce deaths.

- These results come from a model with important uncertainties, including the size of the community population and the costs, and should be checked against programme data.

# Abbreviations

CHEERS, Consolidated Health Economic Evaluation Reporting Standards; D, people who inject drugs and have never been imprisoned; F, released people who still inject; HCV, hepatitis C virus; J, people in prison; MCMC, Markov chain Monte Carlo; PWID, people who inject drugs; UI, uncertainty interval; WHO, World Health Organization; X, released people who have stopped injecting.

# Introduction

WHO set targets for 2030 to eliminate viral hepatitis as a public health threat, including large reductions in new HCV infections and in HCV deaths compared with 2015 [@who2022strategy; @who2023validation]. [claim:C010] [evidence:E003,E019]
Direct-acting antivirals cure most people treated [@pawlotsky2020], and treating people who inject drugs can also prevent onward transmission [@martin2013; @dore2021]. [claim:C011] [evidence:E020,E021,E022]

In Singapore, people who inject drugs are a key population for HCV [@winslow2007; @cda_hcv], and their number has been estimated at about 11,000 [@teo2019]. [claim:C012] [evidence:E006,E017,E023,E024]
Many of them pass through prison, and treatment in prison with direct-acting antivirals expanded from 2019 [@prison2021; @tan2023]. [claim:C013] [evidence:E001,E018,E025]

Prison offers a setting where people can be found and treated, but linkage to care after release is often poor [@clink2022; @eckhardt2022]. [claim:C014] [evidence:E026,E027]
Community testing programmes in Singapore have been described in policy statements [@moh2024ett; @moh2025ett; @cda_programme], and point-of-care testing can improve linkage to care [@hsiang2022]. [claim:C015] [evidence:E028,E029,E030,E031]
A previous model examined what is needed to eliminate HCV in Singapore [@chaillon2021]. Our model adds the distinction between released people who return to injecting and those who stop. [claim:C016] [evidence:E032]

This distinction matters for policy. People who still inject keep transmitting and are often re-arrested, so prison treatment reaches them again. People who have stopped do not transmit, but they may carry long-standing infections that progress to cirrhosis and liver cancer, and they are less likely to return to prison. [claim:C017] [evidence:E009,E033]

We used a calibrated transmission model to estimate whether prison treatment alone can meet the WHO 2030 targets for PWID in Singapore, and how a fixed annual number of community screening tests should be split between released people who still inject and those who have stopped. [claim:C018] [evidence:E001,E002]

# Methods

## Study design

This is a modelling study with a costing of the screening and treatment policies. We report it following the Consolidated Health Economic Evaluation Reporting Standards 2022 (CHEERS 2022) [@cheers2022] where they apply (S1 Checklist). [claim:C019] [evidence:E034]

[Protocol registration: to be confirmed by the authors. No registered protocol is recorded in the project repository.]

## Model structure

We built a deterministic compartmental model of PWID in Singapore (Fig 1). People enter the model when they start injecting and are placed in group D (injecting, never imprisoned). First arrest moves them to J (in prison). On release, a share returns to injecting (F) and the rest stop (X). People in F can be re-arrested; people in X are not. [claim:C020] [evidence:E033,E004]

Each group was divided into six age groups (<20, 20–29, 30–39, 40–49, 50–59 and 60 years or older), four liver disease stages (non-cirrhotic, compensated cirrhosis, decompensated cirrhosis and hepatocellular carcinoma) and five HCV states: susceptible, acute, chronic, on treatment and recovered. Recovered people, after spontaneous clearance or cure, remain antibody positive and can be reinfected. [claim:C021] [evidence:E033,E006]

Transmission occurs only among people who inject (D and F). The force of infection in each age group depends on an age contact matrix and on the share of infectious (acute or chronic) people in each age group of D and F. There is no transmission in prison or after people stop injecting. [claim:C022] [evidence:E033]

Natural history, treatment efficacy and mortality parameters came from the literature (S1 Text, Table A). Background mortality was the 2015 Singapore death rate by age, multiplied by a standardised mortality ratio for people who use drugs of 9.52 for viraemic people and 4.73 for others [@vansanten2014]. Decompensated cirrhosis and hepatocellular carcinoma add liver deaths, which are lower after cure. [claim:C023] [evidence:E004,E033,E035]

The model was solved with a fourth-order Runge–Kutta method with a daily step. Each run started from a small seeded population and ran for 150 model years, so that 2015 was a steady state; runs whose last ten years changed by more than set tolerances were excluded (S1 Text, Section 1). [claim:C024] [evidence:E033,E036,E006]

## Prison rates

We estimated release and re-arrest rates by age from Singapore Prison Service records of 20,720 people held for drug consumption with admissions from 2000 to 2017 (33,557 prison stays), using Bayesian survival models fitted in Stan (S1 Text, Section 3). Release rates were 0.33 to 0.64 per year by age. A cure-type model gave a 72.4% probability of eventual re-arrest (95% credible interval 71.4%–73.4%), and re-arrest rates of 0.27 to 0.55 per year among those re-arrested. [claim:C025] [evidence:E009,E037,E004,E010]

We assumed that people who are eventually re-arrested return to injecting (F) and that the others stop (X). The first-arrest rate from D was the re-arrest rate multiplied by a calibrated factor c. [claim:C026] [evidence:E004,E033]

## Entry and calibration targets

Entry into D by age was the number of new drug abusers recorded by the Central Narcotics Bureau in 2015 (257, 644, 292, 94, 21, 1 by age group; 1,309 in total) [@singstat_drug], multiplied by a calibrated scale of at least 1, and was held constant in every year. [claim:C027] [evidence:E005,E038,E039]

We calibrated the model to three 2015 targets (S1 Text, Table C): HCV antibody prevalence by age among people in prison at the end of 2015 with a result from the December 2014 to February 2016 screening exercise (3,367 people); the same Central Narcotics Bureau counts as the annual number of first arrests by age; and 11,000 active injectors in the community (95% confidence interval 6,500–17,000), a network scale-up estimate [@teo2019], compared with the modelled number in D and F. [claim:C028] [evidence:E006,E014,E005,E024]

## Calibration

We estimated 13 parameters on the log scale: six multipliers of the contact matrix (one per age group), six entry scales and the first-arrest factor c. Priors were fixed before fitting: normal with mean 0 and standard deviation 1.5 for the contact multipliers, half-normal with standard deviation 1 for the entry scales, and normal with mean 0 and standard deviation 0.5 for log c, each truncated to wide bounds (S1 Text, Table D). [claim:C029] [evidence:E038]

The likelihood combined a beta distribution for each seroprevalence (precision 100), a negative binomial distribution for each first-arrest count (size 50) and a log-normal distribution for the community estimate (log-scale standard deviation 0.2453, taken from its confidence interval). [claim:C030] [evidence:E038,E006]

A multistart Nelder–Mead search found the posterior mode. We then sampled the posterior with an adaptive random-walk Metropolis algorithm [@haario2001]: 8 chains started from independent prior draws, each with 60,000 model runs, of which 20,000 were warm-up. Convergence was checked with rank-normalised R-hat and effective sample sizes [@vehtari2021]. [claim:C031] [evidence:E040,E007,E041,E042]

We checked the sampler with simulation-based calibration [@talts2018] on simulated data sets, and compared the posterior with neural posterior estimation [@greenberg2019; @sbi2020] as a cross-check (S1 Text, Section 4). [claim:C032] [evidence:E043,E044,E045]

## Policy scenarios

All scenarios share the same history: no treatment until 2018 and treatment in prison from 2019, at a rate that gives a 99.995% annual probability of starting treatment for people with chronic infection in J. From 2026, we added community screening for F and X. We compared prison treatment only (J only), 15 finite policies (600, 1,000, 2,000, 5,000, 8,000 tests per year, with 25%, 50% or 75% given to F and the rest to X), and uncapped screening of everyone in F and X each year as a reference. [claim:C033] [evidence:E001,E002]

The smallest capacity, 600 tests per year, was set with reference to a public community testing programme [@moh2024ett]; the other capacities and shares are illustrative and are not announced policies. [claim:C034] [evidence:E001,E028]

At the start of each year, the tests given to F or X were applied to people not on treatment, in proportion to their HCV states. Each antibody-positive person had an RNA test, and the number of chronic infections found set that year's treatment capacity in the group. Within that capacity, people with chronic infection started treatment at a rate that gives a 99% annual probability. Sofosbuvir–velpatasvir without ribavirin was assumed, and people in D were never treated. [claim:C035] [evidence:E001,E002,E004]

## Outcomes

The primary outcomes were the reductions in 2030, compared with 2015, of annual new HCV infections among PWID and of annual HCV deaths. New infections were the annual number of infections in D and F. HCV deaths were liver deaths in decompensated cirrhosis and hepatocellular carcinoma, plus the excess non-liver deaths of viraemic people over the rate of non-viraemic people. [claim:C036] [evidence:E002]

A policy met a target if the median reduction was at least 80% for infections or 65% for deaths. We also report the share of posterior draws in which both targets were met. [claim:C037] [evidence:E002,E001]

## Costs

We costed antibody tests (S$25.31), RNA tests (S$307), one consultation per treatment start (S$25) and a 12-week sofosbuvir–velpatasvir course (S$40,405) in 2026 Singapore dollars, from published charges and a published price (S1 Text, Table H) [@wong2019; @nhg_dsc; @parkway_rna; @sgd_rate]. [claim:C038] [evidence:E016,E001,E046,E047,E048,E049]

Costs for 2026–2030 in all three groups (J, F and X) were discounted at 3% a year to 2026. Unit prices were varied together by a factor of 0.5 and 1.5 in a sensitivity analysis. These are charges, not costs to the health system. [claim:C039] [evidence:E001,E016]

## Uncertainty

For each of 2,000 posterior draws we ran every scenario and computed every outcome. We report medians and 95% uncertainty intervals (UI; 2.5th to 97.5th percentiles) across draws. These intervals reflect parameter uncertainty only. [claim:C040] [evidence:E008,E050]

## Software

The model was written in C++ and called from R through Rcpp [@rcpp2011]; rate models used Stan. Code and generated results are in the project repository. [claim:C041] [evidence:E051]

## Ethics statement

[Ethics statement: to be completed by the authors. Record-level Singapore Prison Service data were used to estimate prison rates and the 2015 seroprevalence targets; the approving ethics committee, the approval number and the consent arrangements are not recorded in the project repository.] [claim:C042] [evidence:E010]

# Results

## Calibration and the 2015 population

The posterior chains converged (largest R-hat 1.004; smallest bulk effective sample size 2,158). The model reproduced the age pattern of seroprevalence in prison, and the 95% posterior intervals of all 13 modelled targets overlapped the target ranges (Fig 2). [claim:C043] [evidence:E052,E053]

At ages 30 and over, modelled first arrests were higher than the recorded counts (for example 324 against 292 at ages 30–39 at the posterior mode), because entry could not fall below the recorded counts and the entry scales settled at their lower bound (posterior medians 1.09 to 1.15 below age 60). [claim:C044] [evidence:E054,E015]

The modelled number of active injectors in the community was 7,957 (7,062–9,298), below the survey estimate of 11,000 but within its range. The first-arrest factor c was 2.78 (95% credible interval 1.55–5.49), so first arrest was faster than re-arrest. [claim:C045] [evidence:E013,E015,E006]

In 2015 the model had 57,742 people who had ever injected (95% credible interval 52,350–66,850): 1,130 in D, 9,741 in J, 6,802 in F and 40,035 in X. The modelled prison population was not a calibration target and was much larger than the 3,367 people with a screening result recorded in prison in 2015. [claim:C046] [evidence:E013,E014]

In 2015 the model had 637 (547–755) new infections and 307 (261–368) HCV deaths per year. Without treatment, both stayed at their 2015 level until 2030. [claim:C047] [evidence:E012,E055]

## Prison treatment alone

With treatment in prison from 2019, new infections fell steadily (Fig 3). By 2025 they were 34% of the 2015 level, and by 2030 the reduction was 84.4% (95% UI 67.9%–90.4%), meeting the incidence target at the median and in 78% of draws. [claim:C048] [evidence:E012,E008]

HCV deaths fell less: by 54.4% (95% UI 51.6%–57.0%) in 2030. No posterior draw met the mortality target with prison treatment alone. [claim:C049] [evidence:E008]

The reason is where infections remain. Released people who still inject are often re-arrested and treated again in prison, but people who have stopped injecting are never re-arrested in the model. In 2026, under prison treatment only, 89% of the 4,176 people with chronic infection in F and X were in X, and F and X together accounted for 97% of all HCV deaths. [claim:C050] [evidence:E011,E012]

People in D, who had never been imprisoned, accounted for 12% of new infections in 2030 under prison treatment only; no policy in this study reached them. [claim:C051] [evidence:E008]

## Community screening policies

All 15 finite screening policies met the incidence target at the median, but only three met the mortality target (Fig 4 and Table 1). With 8,000 tests a year, the fall in deaths was 72.1% (95% UI 69.1%–74.5%) with 25% to F and 67.9% (95% UI 64.9%–70.5%) with 50% to F; both targets were met in 99.9% and 97% of draws. [claim:C052] [evidence:E008]

With 5,000 tests a year and 25% to F, deaths fell by 66.7% (95% UI 63.8%–69.1%), and both targets were met in 87% of draws. With 50% to F at the same capacity, deaths fell by 63.7% (95% UI 60.7%–66.2%) and both targets were met in only 16% of draws. [claim:C053] [evidence:E008]

At every capacity, giving a larger share of tests to F increased the fall in new infections and reduced the fall in deaths. All quotas were used in full in every year and group. [claim:C054] [evidence:E008,E011]

This trade-off follows from the two groups. At the start of 2026, 6,985 people in F were eligible for screening, of whom 6.5% had chronic infection; in X, 40,413 were eligible and 9.2% had chronic infection. A test given to X therefore found more chronic infections, whereas a test given to F found fewer infections but in people who could transmit. [claim:C055] [evidence:E011]

Uncapped screening of everyone in F and X each year reduced new infections by 99.9% (95% UI 99.5%–100.0%) and deaths by 88.7% (95% UI 87.8%–89.6%). Even then, 12% of the remaining infections in 2030 were in D. [claim:C056] [evidence:E008]

## Costs

Discounted costs for 2026–2030 were S$34.3 million (95% UI 25.5–48.9) with prison treatment only. The cheapest policy that met both targets at the median, 5,000 tests a year with 25% to F, cost S$86.2 million (95% UI 70.9–106.4), or S$51.4 million (95% UI 42.3–61.9) more than prison treatment only. [claim:C057] [evidence:E008]

Drugs made up 95% of the cost of this policy, so the total scales almost directly with the drug price. With 8,000 tests a year and 25% to F, the cost was S$109.4 million (95% UI 90.9–133.0); S1 Text, Section 8 gives all policies. [claim:C058] [evidence:E011,E008]

# Discussion

In this model of PWID in Singapore, prison treatment from 2019 was enough to reduce new HCV infections by more than four fifths by 2030, but it reduced HCV deaths by only about half. The gap arises because released people account for nearly all HCV deaths, and most of their chronic infections are among people who have stopped injecting, whom prison treatment does not reach again. Meeting the mortality target required community screening of several thousand released people a year, with most tests given to people who have stopped injecting. [claim:C059] [evidence:E008,E011]

The two WHO targets pull screening in different directions. Tests given to people who still inject prevent more infections, while tests given to people who have stopped prevent more deaths. Because prison treatment already achieves the incidence target at the median, the binding target in Singapore is mortality, and screening policy should be judged mainly by the deaths it prevents. [claim:C060] [evidence:E008]

This suggests that the people who matter most for HCV mortality are not those most visible to drug services. People who have stopped injecting may no longer be in contact with prison or drug services, and reaching them may need testing in primary care or other community settings [@hsiang2022]. [claim:C061] [evidence:E031]

Part of the modelled HCV mortality is excess non-liver death among viraemic people, based on mortality ratios from a cohort of people who use drugs [@vansanten2014]. Evidence that viraemic infection raises mortality from causes other than liver disease supports this [@lee2012; @cepeda2022], but the excess is not found in every cohort [@piselli2021]; with a smaller excess, deaths would depend more on cirrhosis and liver cancer. [claim:C062] [evidence:E035,E056,E057,E058]

We extend earlier modelling for Singapore [@chaillon2021] by separating released people by injecting status and by modelling screening capacity explicitly. This separation gives a reason why the mortality target may be harder to meet than the incidence target. [claim:C063] [evidence:E032]

## Strengths and limitations

The model was calibrated with a full Bayesian method that converged well, and every policy result carries posterior uncertainty from 2,000 draws. Prison flows were estimated from individual prison records rather than assumed, and screening capacity was modelled explicitly, so treatment depends on how many people are found. [claim:C064] [evidence:E052,E008,E009]

The study has important limitations. First, the modelled prison population (9,741) was about 3 times the number of screened people recorded in prison in 2015. The prison size was not fitted, and the model may keep people in prison or in contact with it for too long, which would favour prison treatment. [claim:C065] [evidence:E013,E014]

Second, the entry and first-arrest parameters were weakly identified. Central Narcotics Bureau counts include people who do not inject, the entry scales settled at their lower bound, and the posterior of c was only modestly narrower than its prior (standard deviation 75% of the prior's). The size of the community population that never enters prison is therefore uncertain. [claim:C066] [evidence:E015,E007,E005]

Third, we assumed that released people who are eventually re-arrested return to injecting and that the others stop permanently, and that there is no transmission in prison. Relapse after a long gap would put more infections back into F. [claim:C067] [evidence:E004,E033]

Fourth, screening was an aggregate approximation: tests were spread across each group in proportion to its HCV states, with perfect linkage within capacity and no loss to follow-up. Real programmes will find and treat fewer people for the same number of tests. [claim:C068] [evidence:E002]

Fifth, costs were published charges and a published drug price carried forward without inflation, not payer costs, and follow-up visits were excluded. Our uncertainty intervals cover model parameters only, not structural assumptions, prices or uptake. [claim:C069] [evidence:E016,E001]

Finally, the seroprevalence targets could not be reproduced exactly from the current processing of the prison records (largest difference 0.5 percentage points), and the neural posterior used as a cross-check failed its calibration check (S1 Text, Section 4). [claim:C070] [evidence:E006,E014,E059]

## Conclusions

In this model, treatment in prison is likely to meet the 2030 target for new HCV infections among people who inject drugs in Singapore, but not the target for HCV deaths. Reducing deaths requires finding and treating people who have stopped injecting after release, which needs community testing on a scale of thousands of people a year. [claim:C071] [evidence:E008]

# Data availability

[Data availability statement: to be completed by the authors. Model code, aggregate inputs and all generated results are in the project repository (link and DOI to be added). The record-level Singapore Prison Service data cannot be shared by the authors; the access conditions and the contact for access requests must be stated.]

# Funding

[Funding statement: to be completed by the authors, including grant numbers and the role of any funder, or the PLOS statement for unfunded work.]

# Competing interests

[Competing interests: to be declared by each author.]

# Author contributions

[CRediT author contributions: to be completed by the authors.]

# Use of artificial intelligence

[AI use statement: to be completed and approved by the authors. Draft: Claude (Anthropic) was used to write and revise analysis code and to draft and revise the text, figures and supporting information, with procedures from Scientific Agent Skills [@kassis2026]. The authors are to state that they checked all content, data, code, references and conclusions, and take responsibility for them.] [claim:C072] [evidence:E060]

# Acknowledgments

[Acknowledgments: to be completed by the authors.]

# References

::: {#refs}
:::

# Figure captions

**Fig 1. Model structure.** People who inject drugs enter D (never imprisoned) and move to J (prison) at first arrest. On release they move to F (still injecting) or X (stopped). People in F can be re-arrested. Transmission occurs only in D and F. Treatment is given in J from 2019 and community screening is offered to F and X from 2026. [claim:C073] [evidence:E033,E001]

**Fig 2. Calibration to the 2015 targets.** (A) HCV antibody prevalence in prison by age. (B) First arrests by age (Central Narcotics Bureau new drug abusers). (C) Active injectors in the community (D and F). Black: target with the 95% range of its likelihood; blue: 95% posterior interval of the model value from the MCMC posterior; orange: neural posterior estimate (cross-check only). [claim:C074] [evidence:E053]

**Fig 3. New infections and HCV deaths, 2015–2030.** Annual new HCV infections (A) and HCV deaths (B) relative to 2015, without treatment, with prison treatment only, and with 8,000 community tests a year from 2026 given 25% or 75% to F. Each value is the annual rate on 31 December of the labelled year, drawn at 1 January of the next year, so the falls start at the vertical lines marking prison treatment (2019) and community screening (2026). Lines are medians and bands are 95% uncertainty intervals across 2,000 posterior draws. Dashed lines mark the WHO targets. [claim:C075] [evidence:E008,E012]

**Fig 4. Reductions in 2030 by screening capacity and share of tests given to F.** (A) New infections. (B) HCV deaths. Points are medians and bars are 95% uncertainty intervals. Dashed lines mark the targets; the dotted line is prison treatment only and the solid line is uncapped screening. [claim:C076] [evidence:E008]

# Tables

**Table 1. Reductions in 2030 compared with 2015, probability of meeting both WHO targets, and discounted costs for 2026–2030.** [claim:C077] [evidence:E008]

| Policy | Fall in new infections, % (95% UI) | Fall in HCV deaths, % (95% UI) | Draws meeting both targets | Cost, S$ million (95% UI) | [claim:C078] [evidence:E008]
|---|---|---|---|---|
| J only | 84.4% (67.9%–90.4%) | 54.4% (51.6%–57.0%) | 0.0% | 34.3 (25.5–48.9) | [claim:C079] [evidence:E008]
| J + 600/year, 25% to F | 85.6% (70.0%–91.2%) | 56.1% (53.4%–58.7%) | 0.0% | 41.6 (32.1–56.8) | [claim:C080] [evidence:E008]
| J + 600/year, 50% to F | 86.7% (71.9%–91.9%) | 55.7% (52.9%–58.3%) | 0.0% | 39.4 (30.1–54.4) | [claim:C081] [evidence:E008]
| J + 600/year, 75% to F | 87.7% (73.7%–92.5%) | 55.3% (52.4%–57.9%) | 0.0% | 37.0 (28.1–51.9) | [claim:C082] [evidence:E008]
| J + 1,000/year, 25% to F | 86.3% (71.3%–91.7%) | 57.2% (54.5%–59.7%) | 0.0% | 46.2 (36.3–62.1) | [claim:C083] [evidence:E008]
| J + 1,000/year, 50% to F | 88.0% (74.3%–92.8%) | 56.5% (53.7%–59.1%) | 0.0% | 42.6 (33.1–57.8) | [claim:C084] [evidence:E008]
| J + 1,000/year, 75% to F | 89.6% (77.1%–93.7%) | 55.8% (52.9%–58.4%) | 0.0% | 38.7 (29.6–53.6) | [claim:C085] [evidence:E008]
| J + 2,000/year, 25% to F | 88.0% (74.3%–92.8%) | 59.8% (57.2%–62.3%) | 0.0% | 57.3 (46.2–74.1) | [claim:C086] [evidence:E008]
| J + 2,000/year, 50% to F | 91.0% (79.7%–94.6%) | 58.5% (55.7%–61.1%) | 0.0% | 50.2 (40.1–66.3) | [claim:C087] [evidence:E008]
| J + 2,000/year, 75% to F | 93.4% (84.4%–96.2%) | 57.0% (54.1%–59.7%) | 0.0% | 42.6 (33.4–57.4) | [claim:C088] [evidence:E008]
| J + 5,000/year, 25% to F | 92.3% (82.0%–95.5%) | 66.7% (63.8%–69.1%) | 87.2% | 86.2 (70.9–106.4) | [claim:C089] [evidence:E008]
| J + 5,000/year, 50% to F | 96.7% (91.3%–98.2%) | 63.7% (60.7%–66.2%) | 15.8% | 70.3 (57.7–87.5) | [claim:C090] [evidence:E008]
| J + 5,000/year, 75% to F | 98.8% (96.1%–99.4%) | 60.0% (57.0%–62.8%) | 0.0% | 52.4 (42.7–66.7) | [claim:C091] [evidence:E008]
| J + 8,000/year, 25% to F | 95.3% (88.2%–97.3%) | 72.1% (69.1%–74.5%) | 99.9% | 109.4 (90.9–133.0) | [claim:C092] [evidence:E008]
| J + 8,000/year, 50% to F | 99.0% (96.8%–99.6%) | 67.9% (64.9%–70.5%) | 97.3% | 87.4 (72.6–106.7) | [claim:C093] [evidence:E008]
| J + 8,000/year, 75% to F | 99.8% (99.2%–100.0%) | 62.4% (59.5%–65.2%) | 3.6% | 61.8 (51.1–75.8) | [claim:C094] [evidence:E008]
| J + uncapped screening | 99.9% (99.5%–100.0%) | 88.7% (87.8%–89.6%) | 100.0% | 209.6 (173.4–259.6) | [claim:C095] [evidence:E008]

UI, uncertainty interval across posterior draws. Costs are 2026 Singapore dollars for J, F and X, discounted at 3% a year to 2026. [claim:C096] [evidence:E001]

# Supporting information

**S1 Text. Technical appendix.** Model equations, parameters, prison rate estimation, calibration and its checks, the 2015 model population, scenario implementation, costs and additional results. [claim:C097] [evidence:E033]

**S1 Checklist. CHEERS 2022 checklist.** [claim:C098] [evidence:E034]
