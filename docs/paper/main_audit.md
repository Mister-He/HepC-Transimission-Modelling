---
title: "Whom to test, where and how intensively: prison treatment and post-release screening to eliminate hepatitis C among people who inject drugs in Singapore. A modelling study with multi-criteria decision analysis"
---

**Short title:** Hepatitis C elimination after prison release in Singapore

Yichen He\*

[Affiliation: to be provided by the authors]

\* Corresponding author. E-mail: [to be provided by the authors]

# Abstract

## Background

Hepatitis C virus (HCV) infection in Singapore is concentrated among people who inject drugs (PWID), many of whom cycle through prison. Eliminating HCV among PWID needs decisions on whom to test and treat, through which delivery point and how intensively. We asked whether treatment in prison, alone or with community screening of released people from 2026, can meet the World Health Organization 2030 elimination targets (80% fewer new infections and 65% fewer HCV deaths than in 2015), and which strategy is preferred when health benefits, the WHO targets and the programme budget are weighed together. [claim:C001] [evidence:E001,E002,E003]

## Methods and findings

We built an age-structured model of PWID never imprisoned (D), in prison (J), and released and still injecting (F) or stopped (X), with prison rates estimated from records of 20,720 people. We calibrated it by Markov chain Monte Carlo to 2015 prison seroprevalence, first arrests of new heroin and methamphetamine users and a survey estimate of active injectors. In the main analysis nothing changed from 2015 until 2026, when treatment in prison and 15 community screening policies (600 to 8,000 tests a year, 25%–75% to F) began; we varied prison treatment coverage and, in a sensitivity analysis, started prison treatment in 2019. We projected to 2050 because liver deaths lag infection by decades, counted infections and liver deaths averted against prison treatment only (the status quo) and against no treatment, and chose among strategies with Pareto efficiency and TOPSIS, as in a multi-criteria decision analysis. Results are over 2,000 posterior draws. [claim:C002] [evidence:E004,E005,E006,E007,E008,E001,E009,E010,E011,E012,E013]

Prison treatment alone from 2026 reduced new infections in 2030 by 51.5% (95% uncertainty interval [UI] 40.6%–57.3%) and HCV deaths by 22.1% (95% UI 21.2%–23.0%). Five of 15 screening policies met the incidence target, but none met the mortality target; the best, 8,000 tests a year with 25% to F, reduced deaths by 44.0% (95% UI 41.4%–45.9%). Only uncapped screening of all released people met both targets, reducing deaths by 71.4% (95% UI 71.3%–71.6%), but not if cure lowered the risk of liver death in decompensated cirrhosis less than assumed. Moving tests from F to X increased the fall in deaths and reduced the fall in infections. Lower prison treatment coverage reduced every effect. Had prison treatment started in 2019, it alone would have reduced new infections by 85.0% (95% UI 72.7%–89.7%), meeting the incidence target, but no capped policy would have met the mortality target; the best reduced deaths by 56.1% (95% UI 54.0%–57.7%). Over 2026–2050, 8,000 tests a year with 75% to F averted 546 (95% UI 454–648) liver deaths more than prison treatment alone for an extra programme budget of S$43.7 million (95% UI 14.3–59.0), and was ranked first by TOPSIS with equal weights in 99% of draws; uncapped screening averted the most liver deaths and was preferred when the WHO targets or liver deaths were weighted more. Limitations include treating new heroin and methamphetamine users as the people who inject, weakly identified entry parameters, liver disease inputs from other settings, to which the fall in deaths was sensitive, charge-based costs, and uncertainty intervals that cover the calibrated parameters only. [claim:C003] [evidence:E009,E014,E015,E016,E006,E017,E012,E018,E019]

## Conclusions

In this model, prison treatment and capped community screening could not meet the 2030 target for HCV deaths among PWID in Singapore, even with prison treatment from 2019; with treatment from 2026, the incidence target needs thousands of community tests a year in addition to prison treatment. Over a longer horizon, several thousand community tests a year after release, most of them for people who still inject, gave the best balance of infections and liver deaths averted against the programme budget, while screening every released person each year is needed if meeting the 2030 mortality target is the priority. [claim:C004] [evidence:E009,E015,E018,E020,E013]

# Author summary

## Why was this study done?

- Most hepatitis C infections in Singapore are among people who inject drugs, and many of them spend time in prison. [claim:C005] [evidence:E021,E022]

- It is unclear whether treatment in prison, with or without community screening after release, can meet the WHO 2030 elimination targets, which ask for large falls in both new infections and deaths. [claim:C006] [evidence:E003]

- Community screening capacity is limited, and it is not known whether it should go to released people who still inject or to those who have stopped, or how many tests a year are worth their cost.

## What did the researchers do and find?

- We built a model of hepatitis C among people who inject drugs that follows them into and out of prison, and fitted it to prison and national data from 2015. [claim:C007] [evidence:E007]

- If treatment starts in 2026, prison treatment alone cut new infections by about 52% and deaths by about 22% by 2030. With up to 8,000 community tests a year, several policies met the infection target but none met the death target. [claim:C008] [evidence:E009]

- If prison treatment had started in 2019, it alone would have met the infection target, but only testing every released person every year would have met the death target. [claim:C009] [evidence:E015]

- Looking to 2050, and weighing infections and liver deaths prevented, the WHO targets and the programme budget together, 8,000 community tests a year with most given to people who still inject was the preferred strategy; testing every released person each year prevented the most liver deaths but cost much more. [claim:C010] [evidence:E018,E020,E013,E023,E012]

## What do these findings mean?

- Time matters: the later treatment starts, the less any policy can achieve by 2030. Deaths from liver disease fall slowly after cure, because people with advanced liver disease remain at risk. [claim:C011] [evidence:E009,E015]

- Prison treatment prevents many new infections, but most hepatitis C deaths are among released people who have stopped injecting, who need community testing.

- Because liver disease takes decades to develop, most of the deaths prevented by treatment now occur after 2030, so strategies should be judged over a longer horizon as well as against the 2030 targets. [claim:C012] [evidence:E023,E012,E024]

- These results come from a model with important uncertainties, including the size of the community population and the costs, and should be checked against programme data.

# Abbreviations

CHEERS, Consolidated Health Economic Evaluation Reporting Standards; D, people who inject drugs and have never been imprisoned; DAA, direct-acting antiviral; F, released people who still inject; HCV, hepatitis C virus; J, people in prison; MCMC, Markov chain Monte Carlo; PWID, people who inject drugs; TOPSIS, Technique for Order of Preference by Similarity to Ideal Solution; UI, uncertainty interval; WHO, World Health Organization; X, released people who have stopped injecting.

# Introduction

WHO set targets for 2030 to eliminate viral hepatitis as a public health threat, including large reductions in new HCV infections and in HCV deaths compared with 2015 [@who2022strategy; @who2023validation]. [claim:C013] [evidence:E003,E025]
Direct-acting antivirals cure most people treated [@pawlotsky2020], and treating people who inject drugs can also prevent onward transmission [@martin2013; @dore2021]. [claim:C014] [evidence:E026,E027,E028]

In Singapore, people who inject drugs are a key population for HCV [@winslow2007; @cda_hcv], and their number has been estimated at about 11,000 [@teo2019]. [claim:C015] [evidence:E007,E021,E029,E030]
Many of them pass through prison, where treatment with direct-acting antivirals has been offered [@prison2021; @tan2023]. [claim:C016] [evidence:E022,E031]

Prison offers a setting where people can be found and treated, but linkage to care after release is often poor [@clink2022; @eckhardt2022]. [claim:C017] [evidence:E032,E033]
Community testing programmes in Singapore have been described in policy statements [@moh2024ett; @moh2025ett; @cda_programme], and point-of-care testing can improve linkage to care [@hsiang2022]. [claim:C018] [evidence:E034,E035,E036,E037]
A previous model examined what is needed to eliminate HCV in Singapore [@chaillon2021]. Our model adds the distinction between released people who return to injecting and those who stop. [claim:C019] [evidence:E038]

This distinction matters for policy. People who still inject keep transmitting and are often re-arrested, so prison treatment reaches them again. People who have stopped do not transmit, but they may carry long-standing infections that progress to cirrhosis and liver cancer, and they are less likely to return to prison. [claim:C020] [evidence:E010,E039]

Elimination policy for PWID therefore turns on three questions: whom to test and treat (people who have never been imprisoned, people in prison, or released people who still inject or have stopped), through which delivery point (prison or the community after release), and how intensively (how many community tests a year). These choices trade off several criteria at once: the two WHO targets, the infections and deaths averted, and the programme budget. Multi-criteria decision analysis makes such trade-offs explicit [@thokala2016]: strategies that are worse on every criterion than another can be set aside (Pareto efficiency), and the remainder ranked by their closeness to an ideal strategy (TOPSIS) [@hwang1981]. [claim:C021] [evidence:E040,E041]

The time horizon also matters. Chronic HCV infection progresses to cirrhosis over decades [@thein2008], so most of the liver deaths that treatment in 2026 prevents will occur after 2030, and a horizon that ends in 2030 understates the value of treating people with long-standing infection. [claim:C022] [evidence:E024]

We used a calibrated transmission model to estimate whether prison treatment and community screening starting in 2026 can meet the WHO 2030 targets for PWID in Singapore, how a fixed annual number of community screening tests should be split between released people who still inject and those who have stopped, how the answer depends on prison treatment coverage, and what an earlier start in prison would have achieved. We then projected the strategies to 2050 and selected among them with Pareto efficiency and TOPSIS, weighing infections and liver deaths averted, the WHO targets and the programme budget. [claim:C023] [evidence:E001,E002,E013]

# Methods

## Study design

This is a modelling study with a costing of the screening and treatment policies. We report it following the Consolidated Health Economic Evaluation Reporting Standards 2022 (CHEERS 2022) [@cheers2022] where they apply (S1 Checklist). [claim:C024] [evidence:E042]

[Protocol registration: to be confirmed by the authors. No registered protocol is recorded in the project repository.]

## Model structure

We built a deterministic compartmental model of PWID in Singapore (Fig 1). People enter the model when they start injecting and are placed in group D (injecting, never imprisoned). First arrest moves them to J (in prison). On release, a share returns to injecting (F) and the rest stop (X). People in F can be re-arrested; people in X are not. [claim:C025] [evidence:E039,E004]

Each group was divided into six age groups (<20, 20–29, 30–39, 40–49, 50–59 and 60 years or older), four liver disease stages (non-cirrhotic, compensated cirrhosis, decompensated cirrhosis and hepatocellular carcinoma) and five HCV states: susceptible, acute, chronic, on treatment and recovered. Recovered people, after spontaneous clearance or cure, remain antibody positive and can be reinfected. [claim:C026] [evidence:E039,E007]

Transmission occurs only among people who inject (D and F). The force of infection in each age group depends on an age contact matrix and on the share of infectious (acute or chronic) people in each age group of D and F. There is no transmission in prison or after people stop injecting. [claim:C027] [evidence:E039]

Natural history, treatment efficacy and mortality parameters came from the literature (S1 Text, Table A). Each cause of death had one route in the model. Deaths from causes other than liver disease occurred at the 2015 Singapore death rate by age [@singstat_m810141] (ages 60–74 for the oldest group) multiplied by 4.73, the standardised mortality ratio of people who use drugs without HCV or HIV infection [@vansanten2014], and for viraemic people also by 1.47, the hazard ratio of extrahepatic death with viraemic HCV infection [@lee2012]. Liver deaths occurred only in decompensated cirrhosis and hepatocellular carcinoma, at lower rates after cure. [claim:C028] [evidence:E004,E039,E043,E044,E045]

The model was solved with a fourth-order Runge–Kutta method with a daily step. Each run started from a small seeded population and ran for 150 model years, so that 2015 was a steady state; runs whose last ten years changed by more than set tolerances were excluded (S1 Text, Section 1). [claim:C029] [evidence:E039,E046,E007]

## Prison rates

We estimated release and re-arrest rates by age from Singapore Prison Service records of 20,720 people held for drug consumption with admissions from 2000 to 2017 (33,557 prison stays), using Bayesian survival models fitted in Stan (S1 Text, Section 3). Release rates were 0.33 to 0.64 per year by age. A cure-type model gave a 72.4% probability of eventual re-arrest (95% credible interval 71.4%–73.4%), and re-arrest rates of 0.27 to 0.55 per year among those re-arrested. [claim:C030] [evidence:E010,E047,E004,E011]

We assumed that people who are eventually re-arrested return to injecting (F) and that the others stop (X). The first-arrest rate from D was the re-arrest rate multiplied by a calibrated factor c. [claim:C031] [evidence:E004,E039]

## Entry and calibration targets

The Central Narcotics Bureau recorded 1,309 new drug abusers in 2015 (257, 644, 292, 94, 21, 1 by age group) [@singstat_drug]. Of the 1,309 new abusers that year, 113 mainly used heroin and 1,014 methamphetamine [@singstat_m891521], the drugs most often injected in Singapore; we took this share (86.1%) as the share of new abusers who inject, the same at every age, because the drug table has no age breakdown. [claim:C032] [evidence:E005,E006,E004,E048,E049]
Entry into D by age was these counts multiplied by 86.1% (1,127 a year in total) and by a calibrated scale of at least 1, and was held constant in every year. [claim:C033] [evidence:E004,E050,E006]

We calibrated the model to three 2015 targets (S1 Text, Table C): HCV antibody prevalence by age among people in prison at the end of 2015 with a result from the December 2014 to February 2016 screening exercise (3,367 people); first arrests by age, the same counts multiplied by 86.1% and rounded (221, 554, 251, 81, 18, 1); and 11,000 active injectors in the community (95% confidence interval 6,500–17,000), a network scale-up estimate [@teo2019], compared with the modelled number in D and F. [claim:C034] [evidence:E007,E051,E005,E006,E030]

## Calibration

We estimated 13 parameters on the log scale: six multipliers of the contact matrix (one per age group), six entry scales and the first-arrest factor c. Priors were fixed before fitting: normal with mean 0 and standard deviation 1.5 for the contact multipliers, half-normal with standard deviation 1 for the entry scales, and normal with mean 0 and standard deviation 0.5 for log c, each truncated to wide bounds (S1 Text, Table D). [claim:C035] [evidence:E050]

The likelihood combined a beta distribution for each seroprevalence (precision 100), a negative binomial distribution for each first-arrest count (size 50) and a log-normal distribution for the community estimate (log-scale standard deviation 0.2453, taken from its confidence interval). [claim:C036] [evidence:E050,E007]

A multistart Nelder–Mead search found the posterior mode. We then sampled the posterior with an adaptive random-walk Metropolis algorithm [@haario2001]: 8 chains started from independent prior draws, each with 60,000 model runs, of which 20,000 were warm-up. Convergence was checked with rank-normalised R-hat and effective sample sizes [@vehtari2021]. [claim:C037] [evidence:E052,E008,E053,E054]

We checked the sampler with simulation-based calibration [@talts2018] on simulated data sets (S1 Text, Section 4). [claim:C038] [evidence:E055]

To validate the model, we checked that the 95% posterior predictive interval of every calibration target overlapped the target's range, and compared the modelled 2015 prison population, which was not a calibration target, with the number of people held in prison for drug consumption at the end of 2015. The model code was checked by unit and integration tests, including a comparison with an independent R implementation of the model equations. No observed series of new infections or HCV deaths among PWID in Singapore was available for external validation. [claim:C039] [evidence:E056,E051,E057]

## Policy scenarios

In the main analysis, nobody was treated before 2026, so the model stayed at its 2015 equilibrium until then. From 2026, people with chronic infection in J were treated at a rate that gives a 99.995% annual probability of starting treatment, and community screening of F and X began. We compared no treatment, prison treatment only (J only), 15 finite policies (600, 1,000, 2,000, 5,000, 8,000 tests per year, with 25%, 50% or 75% given to F and the rest to X), and uncapped screening of everyone in F and X each year as a reference. [claim:C040] [evidence:E001,E002]

Two comparators were used. No treatment, in which nobody is treated in any group, is the counterfactual for the infections and deaths averted by each strategy. Prison treatment only (J only) is the status quo comparator: it represents continuing the treatment offered in prison without community screening after release, and all incremental costs and outcomes of community screening are measured against it. [claim:C041] [evidence:E001,E002,E012]

To show how much depends on prison treatment, we repeated all 17 policies with an annual probability of starting treatment in J of 0%, 50%, 90% and 99.995% (68 combinations). [claim:C042] [evidence:E001,E014,E002]

In a sensitivity analysis, prison treatment started in 2019 instead, at the main coverage, and community screening still started in 2026. [claim:C043] [evidence:E015,E002]

The smallest capacity, 600 tests per year, was set with reference to a public community testing programme [@moh2024ett]; the other capacities and shares are illustrative and are not announced policies. [claim:C044] [evidence:E001,E034]

At the start of each year, the tests given to F or X were applied to people not on treatment, in proportion to their HCV states. Each antibody-positive person had an RNA test, and the number of chronic infections found set that year's treatment capacity in the group. Within that capacity, people with chronic infection started treatment at a rate that gives a 99% annual probability. Sofosbuvir–velpatasvir without ribavirin was assumed, and people in D were never treated. [claim:C045] [evidence:E001,E002,E004]

The care cascade was ideal in the main analysis: every chronic infection found by community screening was linked to treatment within the year's capacity, every person who started treatment completed it, and cure followed at the efficacy reported in trials. In a sensitivity analysis, only 80% or 60% of the chronic infections found in F and X were linked to treatment; treatment completion was not varied separately. [claim:C046] [evidence:E002,E058,E004]

## Time horizon

The WHO targets were assessed in 2030. Because HCV-related cirrhosis and liver cancer develop over decades [@thein2008], we also projected every strategy to 2050 and report cumulative outcomes over 2026–2030 and 2026–2050; 2026–2050 was the main horizon for strategy selection. Entry of new injectors was held at its 2015 level to 2050, as before 2026, and the prison rates and model parameters were held constant. [claim:C047] [evidence:E001,E002,E023,E004,E024]

## Outcomes

The primary outcomes were the reductions in 2030, compared with 2015, of annual new HCV infections among PWID and of annual HCV deaths. New infections were the number of infections in D and F during the calendar year. HCV deaths were the liver deaths in decompensated cirrhosis and hepatocellular carcinoma during the year, the deaths from HCV-related cirrhosis and liver cancer that the WHO target counts. The 2015 values were the yearly numbers at the calibrated equilibrium. As a secondary outcome, we added the extrahepatic excess deaths of viraemic people (S1 Text, Section 10). [claim:C048] [evidence:E002,E003]

A policy met a target if the median reduction was at least 80% for infections or 65% for deaths. We also report the share of posterior draws in which both targets were met. [claim:C049] [evidence:E002,E001]

For each posterior draw we also summed new infections and liver deaths over 2026–2030 and 2026–2050, and subtracted them from the totals without treatment and with J only to give the infections and liver deaths averted. Averted outcomes were not discounted. Incident cases of decompensated cirrhosis and liver cancer were not counted separately. [claim:C050] [evidence:E002,E012]

## Costs

Costs are the programme budget only: testing, consultations and direct-acting antiviral (DAA) drugs. They exclude the costs of caring for people with decompensated cirrhosis or liver cancer and any savings from the disease that treatment averts, so they are not the net costs of a strategy to the health system. [claim:C051] [evidence:E002,E001]

We costed antibody tests (S$25.31), RNA tests (S$307), one consultation per treatment start (S$25) and a 12-week sofosbuvir–velpatasvir course (S$40,405) in 2026 Singapore dollars, from published charges and a published price (S1 Text, Table H) [@wong2019; @nhg_dsc; @parkway_rna; @sgd_rate]. [claim:C052] [evidence:E059,E001,E060,E061,E062,E063]

Programme budgets for 2026–2030 and 2026–2050 in all three groups (J, F and X) were discounted at 3% a year to 2026. Unit prices were varied together by a factor of 0.5 and 1.5 in a sensitivity analysis. These are charges, not costs to the health system. The cost per infection or liver death averted was the incremental programme budget against J only divided by the outcomes averted against J only, both discounted at the same rate; it was not computed when fewer than one outcome was averted. [claim:C053] [evidence:E001,E059,E002,E012]

## Strategy selection

We compared strategies on five criteria at each horizon: three benefits (infections averted and liver deaths averted against no treatment, and the probability of meeting both WHO targets in 2030) and two costs (the discounted programme budget and the number of treatment courses). We used two strategy sets: the 17 policies of the main analysis, and the 68 combinations of policy and prison treatment coverage. Each criterion was rescaled to 0–1 across the set, with the costs reverse-coded so that 1 was best. [claim:C054] [evidence:E018,E020,E013,E012]

We applied three decision rules: (a) the strategy with the largest health benefit whose programme budget was within a cap; (b) the strategy with the lowest cost per outcome averted against no treatment; and (c) TOPSIS among the Pareto-efficient strategies, ranking each by its relative closeness to the ideal (best on every criterion) and away from the anti-ideal (worst on every criterion) [@hwang1981]. We also drew the cost-effectiveness frontier for each benefit, flagging strategies that were dominated (another cost no more and gained no less) or extendedly dominated (a mix of two other strategies would do better) [@cantor1994]. [claim:C055] [evidence:E018,E020,E013,E041,E064]

The main analysis used posterior medians and equal weights. To carry the uncertainty through, we repeated Pareto efficiency, TOPSIS and the budget rule within each of the 2,000 posterior draws and report the probability that each strategy was Pareto-efficient and that it was ranked first. We varied the weights by doubling or removing each criterion in turn, by using only the benefits or only the costs, and by drawing 1,000 weight vectors uniformly from the simplex. [claim:C056] [evidence:E018,E020,E013,E019]

## Uncertainty

For each of 2,000 posterior draws we ran every scenario and computed every outcome. We report medians and 95% uncertainty intervals (UI; 2.5th to 97.5th percentiles) across draws, and the share of draws meeting a target. These describe uncertainty in the 13 calibrated parameters given the fixed inputs, not the probability that a target will be met. They do not include uncertainty in screening uptake, linkage to care, treatment completion, drug prices or model structure. [claim:C057] [evidence:E009,E065]

To show how much the fixed inputs matter, we held the calibrated parameters at the posterior mode and set each fixed epidemiological input in turn to the low and high end of its literature range, without refitting (S1 Text, Section 13). In further sensitivity analyses, entry of new injectors from 2016 to 2025 followed the observed yearly counts of new drug abusers instead of the 2015 level (S1 Text, Section 11), and the share of infections with genotype 3 was 10% or 50% instead of 78% (S1 Text, Section 12). These analyses were run to 2030; the coverage and linkage analyses were run to 2050 (S1 Text, Section 17). [claim:C058] [evidence:E066,E067,E068,E069,E002,E058]

## Software

The model was written in C++ and called from R through Rcpp [@rcpp2011]; rate models used Stan. Code and generated results are in the project repository. [claim:C059] [evidence:E070]

## Ethics statement

[Ethics statement: to be completed by the authors. Record-level Singapore Prison Service data were used to estimate prison rates and the 2015 seroprevalence targets; the approving ethics committee, the approval number and the consent arrangements are not recorded in the project repository.] [claim:C060] [evidence:E011]

# Results

## Calibration and the 2015 population

The posterior chains converged (largest R-hat 1.005; smallest bulk effective sample size 1,641). The model reproduced the age pattern of seroprevalence in prison, and the 95% posterior intervals of all 13 modelled targets overlapped the target ranges (Fig 2). [claim:C061] [evidence:E071,E056]

At ages 30 and over, modelled first arrests were higher than the targets (for example 279 against 251 at ages 30–39 at the posterior mode), because entry could not fall below the scaled counts and the entry scales settled near their lower bound (posterior medians 1.08 to 1.15 below age 60). [claim:C062] [evidence:E072,E016]

The modelled number of active injectors in the community was 7,332 (6,513–8,610), below the survey estimate of 11,000 but within its range. The first-arrest factor c was 2.71 (95% credible interval 1.54–5.23), so first arrest was faster than re-arrest. [claim:C063] [evidence:E057,E016,E007]

In 2015 the model had 59,599 people who had ever injected (95% credible interval 53,968–69,269): 1,000 in D, 8,768 in J, 6,303 in F and 43,510 in X. The modelled prison population was not a calibration target; it was close to the 8,725 people held in prison for drug consumption at the end of 2015, the population whose release and re-arrest rates the model uses. [claim:C064] [evidence:E057,E051]

In 2015 the model had 529 (454–628) new infections and 157 (130–191) HCV deaths per year. Without treatment, both stayed at their 2015 level until 2030 (within 0.09% until 2025). [claim:C065] [evidence:E073,E074]

## Prison treatment alone

With treatment in prison from 2026, 3,189 people started treatment in J in 2026, as the infections accumulated since 2015 were treated, and fewer in later years. New infections fell by 51.5% (95% UI 40.6%–57.3%) by 2030 and HCV deaths by 22.1% (95% UI 21.2%–23.0%) (Fig 3). No posterior draw met either target. [claim:C066] [evidence:E075,E073,E009]

In 2026, 73% of the 8,027 people with chronic infection in F and X were in X, and F and X together accounted for 90% of all HCV deaths. People in X are never re-arrested in the model, so prison treatment does not reach them. [claim:C067] [evidence:E075,E073,E039]

People in D, who had never been imprisoned, accounted for 15% of new infections in 2030 under prison treatment only; no policy in this study reached them. [claim:C068] [evidence:E009]

## Community screening policies

Five of the 15 finite screening policies met the incidence target at the median (Fig 4 and Table 1): 8,000 tests a year with any share to F, and 5,000 tests with 50% or 75% to F. None met the mortality target. The largest fall in deaths among finite policies was 44.0% (95% UI 41.4%–45.9%), with 8,000 tests and 25% to F, and no finite policy met both targets in any posterior draw. [claim:C069] [evidence:E009]

At every capacity, giving a larger share of tests to F increased the fall in new infections and reduced the fall in deaths. All quotas were used in full. [claim:C070] [evidence:E009,E075]

At the start of 2026, 6,318 people in F were eligible for screening, of whom 33.8% had chronic infection; in X, 43,683 were eligible and 13.5% had chronic infection. A test given to F therefore found more chronic infections, in people who could transmit, but people in F are re-arrested and can be treated in prison, whereas people in X can be reached only by community screening; this is a likely reason why tests given to X prevented more deaths. [claim:C071] [evidence:E075,E010,E039]

Uncapped screening of everyone in F and X each year reduced new infections by 99.8% (95% UI 98.7%–99.9%) and deaths by 71.4% (95% UI 71.3%–71.6%), and was the only policy that met both targets. Even then, 12% of the remaining infections in 2030 were in D. In the one-way analysis of the fixed inputs, the only input that stopped it meeting the mortality target was the relative risk of liver death after cure in decompensated cirrhosis: at the less favourable end of its range, deaths fell by 59.4% (S1 Text, Section 13). [claim:C072] [evidence:E009,E076,E017]

## Prison treatment coverage

Every reduction rose with the coverage of prison treatment (Fig 5 and Table 2). Without prison treatment, deaths fell by 29.6% (95% UI 26.3%–32.0%) with 8,000 tests and 25% to F, and only 8,000 tests with 75% to F met the incidence target, reducing new infections by 88.0% (95% UI 82.8%–90.7%). At 50% coverage, three finite policies met the incidence target, and at 90% the same five as at full coverage. [claim:C073] [evidence:E014,E009]

At no coverage did a finite policy meet the mortality target. Uncapped screening met it at every coverage, though only just when there was no prison treatment, reducing deaths by 65.4% (95% UI 65.1%–65.7%). [claim:C074] [evidence:E014,E009]

## Sensitivity analysis: prison treatment from 2019

Had prison treatment started in 2019, prison treatment alone would have reduced new infections by 85.0% (95% UI 72.7%–89.7%) and deaths by 40.8% (95% UI 39.5%–42.1%) by 2030, meeting the incidence target at the median. All 15 finite policies would then have met the incidence target, but none would have met the mortality target: the largest fall in deaths was 56.1% (95% UI 54.0%–57.7%), with 8,000 tests and 25% to F. Only uncapped screening met both targets, reducing deaths by 76.4% (95% UI 76.1%–76.7%) (S1 Text, Section 9). [claim:C075] [evidence:E015]

## Programme budget, 2026–2030

Discounted programme budgets for 2026–2030 were S$211.6 million (95% UI 180.3–252.4) with prison treatment only, mostly for drugs and mostly in 2026. The cheapest policy that met the incidence target, 5,000 tests a year with 75% to F, cost S$244.3 million (95% UI 208.5–288.3), or S$32.5 million (95% UI 27.7–37.8) more than prison treatment only. [claim:C076] [evidence:E009,E075]

Uncapped screening, the only policy that met both targets, cost S$484.1 million (95% UI 408.9–583.6). Drugs made up 98% of the cost of 5,000 tests a year with 25% to F, so totals scale almost directly with the drug price; S1 Text, Section 8 gives all policies. [claim:C077] [evidence:E075,E009]

## Long-term outcomes, 2026–2050

Liver deaths kept falling after 2030 under every strategy (Fig 6). The median yearly number of liver deaths fell to 35% of the 2015 level in 2028 with uncapped screening, 2035 with 8,000 tests and 25% to F, 2041 with 8,000 tests and 75% to F, and 2048 with prison treatment only; new infections under prison treatment only fell to 20% of the 2015 level in 2036. [claim:C078] [evidence:E023,E012]

Against no treatment, prison treatment only averted 9,892 (95% UI 8,279–11,843) new infections and 1,651 (95% UI 1,397–1,983) liver deaths over 2026–2050; 93% of these liver deaths averted occurred after 2030. Against prison treatment only, uncapped screening averted a further 1,433 (95% UI 1,176–1,788) liver deaths, 8,000 tests with 25% to F 947 (95% UI 782–1,138) and with 75% to F 546 (95% UI 454–648) (S1 Text, Section 15). [claim:C079] [evidence:E023,E012]

At every capacity, giving 25% of tests to F averted more liver deaths over 2026–2050 than giving 75%, while giving 75% averted more infections and needed fewer treatment courses, because infections prevented in people who still inject do not need treatment later. [claim:C080] [evidence:E023,E012]

## Cost per outcome averted

Over 2026–2050, the extra programme budget against prison treatment only was S$43.7 million (95% UI 14.3–59.0) with 8,000 tests and 75% to F, S$140.7 million (95% UI 111.9–170.8) with 25% to F and S$274.1 million (95% UI 224.9–338.8) with uncapped screening. The cost per liver death averted was S$117,000 (95% UI S$41,000–S$148,000), S$214,000 (95% UI S$186,000–S$225,000) and S$256,000 (95% UI S$232,000–S$269,000), respectively. [claim:C081] [evidence:E023,E012]

With 600, 1,000 or 2,000 tests a year and 75% to F, the programme budget over 2026–2050 was about the same as with prison treatment only (95% UI of the difference including zero), because the infections they prevented did not need treatment later; these policies averted more infections and liver deaths at no clear extra cost. [claim:C082] [evidence:E023,E012]

## Strategy selection

Over 2026–2050, the cost-effectiveness frontier for liver deaths averted ran, in order of cost, through 1,000/year, 75% to F; 2,000/year, 75% to F; 5,000/year, 75% to F; 8,000/year, 75% to F; 8,000/year, 50% to F; uncapped screening (Fig 7). Prison treatment only was dominated: 1,000/year, 75% to F cost no more and averted more infections and liver deaths. Every policy giving 25% of tests to F was dominated or extendedly dominated; for infections averted the frontier ran through the 75%-to-F policies from 1,000 tests and uncapped screening. By rule (b), the lowest cost per outcome averted against no treatment was with 5,000/year, 75% to F for both infections and liver deaths. [claim:C083] [evidence:E018,E020,E013]

Eight of the 17 strategies were Pareto-efficient at the median. With equal weights, TOPSIS ranked J + 8,000/year, 75% to F first (closeness 0.565; first in 99.2% of posterior draws), followed by J + 5,000/year, 75% to F (0.559) and J + uncapped screening (0.551) (Fig 8 and Table 3). Over 2026–2030 alone, uncapped screening, the only strategy that met both targets, ranked first. [claim:C084] [evidence:E018,E020,E013]

The ranking depended on the weights. Doubling the weight on liver deaths averted or on meeting both targets made uncapped screening first; doubling the weight on the budget or on treatment courses made 5,000/year, 75% to F first. Across 1,000 random weight vectors, the strategy ranked first most often was J + uncapped screening (48%), then J + 5,000/year, 75% to F (18%), J + 8,000/year, 75% to F (17%) and J + 2,000/year, 75% to F (15%). [claim:C085] [evidence:E018,E020,E013,E019]

With a cap on the programme budget over 2026–2050, the strategy averting the most liver deaths was 5,000/year, 75% to F with a cap of S$300 million; 5,000/year, 50% to F with a cap of S$350 million; 8,000/year, 50% to F with a cap of S$400 million; 8,000/year, 25% to F with a cap of S$450–500 million; uncapped screening with a cap of S$550 million (rule a; S1 Text, Section 16). [claim:C086] [evidence:E018,E020,E013,E077]

When prison treatment coverage was also a choice, TOPSIS ranked uncapped screening with full prison coverage first in 100% of draws, and the top four were uncapped screening at each prison coverage. For every screening policy with more than one Pareto-efficient coverage level, TOPSIS closeness rose with prison coverage. [claim:C087] [evidence:E018,E020,E013]

## Linkage to care

If only 80% or 60% of the chronic infections found in F and X were linked to treatment, every finite policy had smaller effects. With 60% linkage, 8,000 tests and 25% to F reduced liver deaths in 2030 by 36.6% (95% UI 34.6%–38.1%) (main analysis 44.0%) and averted 723 (95% UI 598–862) liver deaths over 2026–2050 against prison treatment only. Uncapped screening still met both targets at 80% and 60% linkage (fall in liver deaths 71.0% (95% UI 70.8%–71.1%) and 68.9% (95% UI 68.7%–69.1%)), because people missed in one year were screened again the next (S1 Text, Section 17). [claim:C088] [evidence:E058,E009]

# Discussion

In this model of PWID in Singapore, starting HCV treatment in 2026 left too little time to meet the WHO 2030 mortality target with prison treatment and capped community screening. Prison treatment from 2026 reduced new infections by about half and deaths by about a fifth. Up to 8,000 community tests a year could meet the incidence target, but only screening of every released person each year met the mortality target, and that depended on how much cure lowers the risk of liver death in decompensated cirrhosis. [claim:C089] [evidence:E009]

Timing mattered most for new infections. With the same prison treatment starting in 2019, prison treatment alone met the incidence target, but even then no capped policy met the mortality target. Each year without treatment adds to the pool of chronic infections that progress to cirrhosis and liver cancer. Cure lowers but does not remove the risk of liver death in people who already have decompensated cirrhosis or liver cancer, and people with compensated cirrhosis can still progress, so liver deaths fall much more slowly than infections after treatment begins. [claim:C090] [evidence:E015,E009]

The two WHO targets pull screening in different directions. Tests given to people who still inject prevent more infections, while tests given to people who have stopped prevent more deaths. Because the mortality target was harder to meet in every analysis, a policy aimed at the 2030 targets should be judged mainly by the deaths it prevents. [claim:C091] [evidence:E009,E015]

Over a longer horizon the choice looks different. Most of the liver deaths that treatment averts occur after 2030, and preventing infections in people who still inject also saves the treatment they would later need. Over 2026–2050, policies giving most community tests to F lay on the cost-effectiveness frontier, the smallest of them at no clear extra budget over prison treatment only, and 8,000 tests a year with 75% to F was ranked first with equal weights in almost every posterior draw. Giving more tests to X averted more liver deaths but cost more, because it found and treated more long-standing infections; uncapped screening averted the most liver deaths and was the only strategy to meet the 2030 targets, but at the highest budget. [claim:C092] [evidence:E023,E012,E018,E020,E013]

In terms of whom, where and how intensively, the analysis gives a consistent picture. Whom: people released from prison, both those who still inject (to prevent infections) and those who have stopped (to prevent liver deaths); people who have never been imprisoned were not reached by any strategy. Where: prison treatment is the base on which community screening builds, and higher prison coverage was preferred in every comparison. How intensively: several thousand community tests a year, well above the 600 a year of the smallest policy, which was set with reference to a public programme, and screening every released person if the 2030 mortality target is to be met. Which strategy is preferred depends on the weight placed on meeting the 2030 targets and on averting liver deaths relative to the budget, and the decision analysis makes that dependence explicit rather than resolving it. [claim:C093] [evidence:E009,E018,E020,E013,E014]

This suggests that the people who matter most for HCV mortality are not those most visible to drug services. People who have stopped injecting may no longer be in contact with prison or drug services, and reaching them may need testing in primary care or other community settings [@hsiang2022]. [claim:C094] [evidence:E037]

We counted only liver deaths, as the WHO target does. Viraemic infection may also raise mortality from other causes [@lee2012; @cepeda2022], although the excess is not found in every cohort [@piselli2021]. Cure removes this excess at once, so counting it as well gave larger falls in deaths (prison treatment only 25.8%, uncapped screening 80.6%; S1 Text, Section 10); still, no finite policy reached 65%. [claim:C095] [evidence:E009,E045,E078,E079]

We extend earlier modelling for Singapore [@chaillon2021] by separating released people by injecting status and by modelling screening capacity explicitly. This separation gives a reason why the mortality target may be harder to meet than the incidence target. [claim:C096] [evidence:E038]

## Strengths and limitations

The model was calibrated with a full Bayesian method that converged well, and every policy result carries posterior uncertainty from 2,000 draws. Prison flows were estimated from individual prison records rather than assumed, and screening capacity was modelled explicitly, so treatment depends on how many people are found. [claim:C097] [evidence:E071,E009,E010]

The study has important limitations. First, the seroprevalence targets came from the 3,367 people in prison at the end of 2015 who had a screening result, out of the 8,725 held for drug consumption. The modelled prison population (8,768) matched the number held, although it was not fitted, but the people screened may not represent everyone in prison. [claim:C098] [evidence:E057,E051]

Second, the main analysis assumes that nobody was treated before 2026. Treatment in prison has in fact been offered [@prison2021; @tan2023], but its start and coverage are not recorded in the project data. The sensitivity analysis with full coverage from 2019 shows how much an earlier start changes the results; the real history probably lies between the two analyses. [claim:C099] [evidence:E015,E001,E022,E031]

Third, we did not observe who injects. We treated new drug abusers whose main drug was heroin or methamphetamine (86.1% of 1,309 in 2015) as the people who start injecting, at the same share in every age group because the drug table has no age breakdown. The records give the main drug, not the route of use: users of these drugs who smoke rather than inject inflate entry, and people who inject other drugs are missed. Entry could rise above this share through the calibrated entry scales but could not fall below it, and it was held at its 2015 level in every year; a fall in new injectors after 2015 would lower infections without treatment, and a rise would raise them. [claim:C100] [evidence:E005,E006,E004,E050,E039]

Fourth, the entry and first-arrest parameters were weakly identified. The entry scales settled near their lower bound, and the posterior of c was only modestly narrower than its prior (standard deviation 72% of the prior's), as was the contact multiplier for ages 50–59. The first-arrest rate was assumed to have the same age pattern as the re-arrest rate, because the records contain only people who had already been arrested. [claim:C101] [evidence:E016,E008,E004]

Fifth, we assumed that released people who are eventually re-arrested return to injecting and that the others stop permanently and are never re-arrested, and that there is no transmission in prison. Relapse after a long gap would put more infections back into F. [claim:C102] [evidence:E004,E039]

Sixth, screening was an aggregate approximation: tests were spread across each group in proportion to its HCV states, with perfect linkage within capacity and no loss to follow-up, and people who had never been imprisoned were never screened or treated. Real programmes will find and treat fewer people for the same number of tests. The linkage sensitivity analysis shows that lower linkage reduces the effect of every finite policy; treatment completion and reinfection after cure in the community were not varied separately. [claim:C103] [evidence:E002,E058]

Seventh, the model assumed that the epidemic was in a steady state in 2015, with constant prison rates, mixing that depends only on age, and the same mortality ratios for people who have stopped injecting as for those who still inject. Costs were published charges and a published drug price carried forward without inflation, not payer costs, and follow-up visits were excluded. Our uncertainty intervals cover the calibrated parameters only, not the fixed inputs, structural assumptions, prices or uptake, so they are too narrow. The one-way analysis shows that the fall in deaths depends most on the liver disease inputs, such as the risk of liver death after cure in decompensated cirrhosis (with uncapped screening, 59.4% to 78.9% across its range), which come from cohorts outside Singapore. S1 Text, Section 14 (Table O) lists every assumption with its basis and its likely effect on the results. [claim:C104] [evidence:E039,E046,E004,E059,E001,E017]

Eighth, the projection to 2050 holds entry of new injectors, prison rates, prices and treatment efficacy at their current values for 25 years, so the long-horizon results are scenarios rather than forecasts. The programme budget excludes the cost of caring for advanced liver disease and so understates the value of strategies that avert liver deaths. The decision analysis depends on the criteria and weights chosen: infections averted and treatment courses are correlated, the benefit criteria were not discounted, and TOPSIS ranks can change when strategies are added to or removed from the set, as the different first choices in the two strategy sets show. We did not value outcomes in quality-adjusted life years, by design. [claim:C105] [evidence:E023,E012,E018,E020,E013,E019]

Finally, the seroprevalence targets could not be reproduced exactly from the current processing of the prison records (largest difference 0.5 percentage points). [claim:C106] [evidence:E007,E051]

## Conclusions

In this model, treatment for PWID in Singapore that starts in 2026 can meet the 2030 target for new HCV infections only with prison treatment and thousands of community tests a year, and cannot meet the target for HCV deaths without screening every released person; even then the result depends on liver disease inputs from other settings. Reducing deaths requires finding and treating people who have stopped injecting after release, and every year of delay makes the target harder to reach. Judged over 2026–2050 on infections and liver deaths averted, the WHO targets and the programme budget together, prison treatment combined with several thousand community tests a year, most of them for released people who still inject, was the preferred strategy unless meeting the 2030 mortality target is given priority. [claim:C107] [evidence:E009,E015,E018,E020,E013]

# Data availability

[Data availability statement: to be completed by the authors. Model code, aggregate inputs and all generated results are in the project repository (link and DOI to be added). The record-level Singapore Prison Service data cannot be shared by the authors; the access conditions and the contact for access requests must be stated.]

# Funding

[Funding statement: to be completed by the authors, including grant numbers and the role of any funder, or the PLOS statement for unfunded work.]

# Competing interests

[Competing interests: to be declared by each author.]

# Author contributions

[CRediT author contributions: to be completed by the authors.]

# Use of artificial intelligence

[AI use statement: to be completed and approved by the authors. Draft: Claude (Anthropic) was used to write and revise analysis code and to draft and revise the text, figures and supporting information, with procedures from Scientific Agent Skills [@kassis2026]. The authors are to state that they checked all content, data, code, references and conclusions, and take responsibility for them.] [claim:C108] [evidence:E080]

# Acknowledgments

[Acknowledgments: to be completed by the authors.]

# References

::: {#refs}
:::

# Figure captions

**Fig 1. Model structure.** People who inject drugs enter D (never imprisoned) and move to J (prison) at first arrest. On release they move to F (still injecting) or X (stopped). People in F can be re-arrested. Transmission occurs only in D and F. In the main analysis, treatment is given in J and community screening is offered to F and X from 2026. [claim:C109] [evidence:E039,E001]

**Fig 2. Calibration to the 2015 targets.** (A) HCV antibody prevalence in prison by age. (B) First arrests by age (Central Narcotics Bureau new drug abusers, scaled to heroin and methamphetamine users). (C) Active injectors in the community (D and F). Black: target with the 95% range of its likelihood; blue: 95% posterior interval of the model value from the MCMC posterior. [claim:C110] [evidence:E056]

**Fig 3. New infections and HCV deaths, 2015–2030.** Annual new HCV infections (A) and HCV deaths (B) relative to 2015, without treatment, with prison treatment only, and with 8,000 community tests a year given 25% or 75% to F, all from 2026. Yearly counts are drawn at mid-year and joined by a monotone curve; the vertical line marks 2026. Lines are medians and bands are 95% uncertainty intervals across 2,000 posterior draws. Dashed lines mark the WHO targets. [claim:C111] [evidence:E009,E073]

**Fig 4. Reductions in 2030 by screening capacity and share of tests given to F.** (A) New infections. (B) HCV deaths. Points are medians and bars are 95% uncertainty intervals. Dashed lines mark the targets; the dotted line is prison treatment only and the solid line is uncapped screening. [claim:C112] [evidence:E009]

**Fig 5. Reductions in 2030 by prison treatment coverage.** Median reductions in new infections (A) and HCV deaths (B) compared with 2015, for each policy (rows) and annual probability of starting treatment in prison from 2026 (columns). Bold values meet the target. [claim:C113] [evidence:E014,E009]

**Fig 6. New infections and liver deaths, 2015–2050.** Annual new HCV infections (A) and liver deaths from decompensated cirrhosis and hepatocellular carcinoma (B) relative to 2015, for no treatment, prison treatment only, 8,000 community tests a year with 25% or 75% to F, and uncapped screening. Lines are medians and bands are 95% uncertainty intervals across 2,000 posterior draws. Dashed lines mark the WHO 2030 targets; dotted lines mark 2026 and 2030. Entry of new injectors is held at its 2015 level throughout. [claim:C114] [evidence:E023,E012]

**Fig 7. Cost-effectiveness frontiers, 2026–2050.** Discounted programme budget (testing, consultations and drugs; 2026 Singapore dollars) against new infections averted (A) and liver deaths averted (B) compared with no treatment, for the 17 policies (posterior medians). Numbered points are on the frontier; squares are extendedly dominated and crosses dominated. Costs exclude care for advanced liver disease. [claim:C115] [evidence:E018,E020,E013]

**Fig 8. TOPSIS ranking of the Pareto-efficient strategies, 2026–2050.** Relative closeness to the ideal strategy with equal weights on infections averted, liver deaths averted, the probability of meeting both WHO targets, the programme budget and treatment courses. (A) The 17 policies at full prison treatment coverage. (B) The 68 combinations of policy and prison treatment coverage (at most 15 shown). Labels give the share of posterior draws in which each strategy ranked first and was Pareto-efficient. [claim:C116] [evidence:E018,E020,E013]

# Tables

**Table 1. Reductions in 2030 compared with 2015, probability of meeting both WHO targets, and discounted costs for 2026–2030 (main analysis).** [claim:C117] [evidence:E009]

| Policy | Fall in new infections, % (95% UI) | Fall in HCV deaths, % (95% UI) | Draws meeting both targets | Cost, S$ million (95% UI) | [claim:C118] [evidence:E009]
|---|---|---|---|---|
| J only | 51.5% (40.6%–57.3%) | 22.1% (21.2%–23.0%) | 0.0% | 211.6 (180.3–252.4) | [claim:C119] [evidence:E009]
| J + 600/year, 25% to F | 54.9% (44.0%–60.4%) | 24.2% (23.2%–25.1%) | 0.0% | 223.6 (190.8–264.8) | [claim:C120] [evidence:E009]
| J + 600/year, 50% to F | 58.2% (47.4%–63.6%) | 23.8% (22.9%–24.8%) | 0.0% | 221.1 (188.8–262.2) | [claim:C121] [evidence:E009]
| J + 600/year, 75% to F | 61.4% (50.7%–66.6%) | 23.5% (22.5%–24.5%) | 0.0% | 218.6 (186.5–259.2) | [claim:C122] [evidence:E009]
| J + 1,000/year, 25% to F | 57.1% (46.3%–62.5%) | 25.5% (24.4%–26.5%) | 0.0% | 231.2 (197.2–273.7) | [claim:C123] [evidence:E009]
| J + 1,000/year, 50% to F | 62.4% (51.9%–67.6%) | 24.9% (23.9%–25.9%) | 0.0% | 227.0 (193.8–268.8) | [claim:C124] [evidence:E009]
| J + 1,000/year, 75% to F | 67.2% (57.0%–72.1%) | 24.3% (23.3%–25.4%) | 0.0% | 222.4 (189.9–263.6) | [claim:C125] [evidence:E009]
| J + 2,000/year, 25% to F | 62.4% (51.9%–67.6%) | 28.7% (27.4%–29.8%) | 0.0% | 249.4 (212.9–294.3) | [claim:C126] [evidence:E009]
| J + 2,000/year, 50% to F | 71.7% (62.0%–76.3%) | 27.5% (26.3%–28.6%) | 0.0% | 240.6 (205.3–284.0) | [claim:C127] [evidence:E009]
| J + 2,000/year, 75% to F | 79.3% (71.0%–83.4%) | 26.1% (25.0%–27.3%) | 0.0% | 230.2 (196.5–272.6) | [claim:C128] [evidence:E009]
| J + 5,000/year, 25% to F | 75.7% (66.6%–80.0%) | 37.1% (35.0%–38.7%) | 0.0% | 296.0 (252.4–348.2) | [claim:C129] [evidence:E009]
| J + 5,000/year, 50% to F | 90.0% (84.2%–92.9%) | 34.1% (32.3%–35.5%) | 0.0% | 272.9 (232.5–321.6) | [claim:C130] [evidence:E009]
| J + 5,000/year, 75% to F | 96.7% (93.7%–98.2%) | 30.2% (28.8%–31.5%) | 0.0% | 244.3 (208.5–288.3) | [claim:C131] [evidence:E009]
| J + 8,000/year, 25% to F | 85.4% (78.3%–88.9%) | 44.0% (41.4%–45.9%) | 0.0% | 332.7 (283.6–391.2) | [claim:C132] [evidence:E009]
| J + 8,000/year, 50% to F | 97.4% (94.8%–98.6%) | 39.4% (37.3%–41.0%) | 0.0% | 297.5 (253.6–350.3) | [claim:C133] [evidence:E009]
| J + 8,000/year, 75% to F | 99.7% (98.5%–99.9%) | 33.2% (31.7%–34.5%) | 0.0% | 255.8 (218.4–301.6) | [claim:C134] [evidence:E009]
| J + uncapped screening | 99.8% (98.7%–99.9%) | 71.4% (71.3%–71.6%) | 100.0% | 484.1 (408.9–583.6) | [claim:C135] [evidence:E009]

UI, uncertainty interval across posterior draws. Costs are 2026 Singapore dollars for J, F and X, discounted at 3% a year to 2026. [claim:C136] [evidence:E001]

**Table 2. Median falls in new infections and HCV deaths in 2030 by annual probability of starting treatment in prison from 2026.** [claim:C137] [evidence:E014,E009]

| Policy | 0% in prison | 50% in prison | 90% in prison | 99.995% in prison | [claim:C138] [evidence:E014,E009]
|---|---|---|---|---|
| J only | -0.0% / -0.1% | 23.9% / 14.0% | 41.8% / 19.7% | 51.5% / 22.1% | [claim:C139] [evidence:E014,E009]
| J + 600/year, 25% to F | 2.7% / 2.6% | 27.7% / 16.4% | 45.5% / 21.9% | 54.9% / 24.2% | [claim:C140] [evidence:E014,E009]
| J + 600/year, 50% to F | 5.6% / 2.9% | 31.6% / 16.3% | 49.1% / 21.7% | 58.2% / 23.8% | [claim:C141] [evidence:E014,E009]
| J + 600/year, 75% to F | 8.7% / 3.1% | 35.4% / 16.2% | 52.6% / 21.4% | 61.4% / 23.5% | [claim:C142] [evidence:E014,E009]
| J + 1,000/year, 25% to F | 4.6% / 4.4% | 30.3% / 18.0% | 47.9% / 23.3% | 57.1% / 25.5% | [claim:C143] [evidence:E014,E009]
| J + 1,000/year, 50% to F | 9.7% / 4.8% | 36.6% / 17.8% | 53.8% / 22.9% | 62.4% / 24.9% | [claim:C144] [evidence:E014,E009]
| J + 1,000/year, 75% to F | 15.1% / 5.1% | 42.8% / 17.6% | 59.3% / 22.4% | 67.2% / 24.3% | [claim:C145] [evidence:E014,E009]
| J + 2,000/year, 25% to F | 9.7% / 8.7% | 36.6% / 21.7% | 53.8% / 26.7% | 62.4% / 28.7% | [claim:C146] [evidence:E014,E009]
| J + 2,000/year, 50% to F | 20.7% / 9.3% | 48.7% / 21.3% | 64.4% / 25.8% | 71.7% / 27.5% | [claim:C147] [evidence:E014,E009]
| J + 2,000/year, 75% to F | 31.9% / 9.4% | 59.6% / 20.6% | 73.4% / 24.6% | 79.3% / 26.1% | [claim:C148] [evidence:E014,E009]
| J + 5,000/year, 25% to F | 26.4% / 20.2% | 54.3% / 31.5% | 69.2% / 35.6% | 75.7% / 37.1% | [claim:C149] [evidence:E014,E009]
| J + 5,000/year, 50% to F | 52.4% / 20.4% | 76.8% / 29.8% | 86.6% / 33.0% | 90.0% / 34.1% | [claim:C150] [evidence:E014,E009]
| J + 5,000/year, 75% to F | 71.3% / 18.8% | 89.6% / 26.8% | 95.2% / 29.4% | 96.7% / 30.2% | [claim:C151] [evidence:E014,E009]
| J + 8,000/year, 25% to F | 42.7% / 29.6% | 69.0% / 39.5% | 80.8% / 42.9% | 85.4% / 44.0% | [claim:C152] [evidence:E014,E009]
| J + 8,000/year, 50% to F | 74.2% / 28.8% | 91.3% / 36.4% | 96.2% / 38.7% | 97.4% / 39.4% | [claim:C153] [evidence:E014,E009]
| J + 8,000/year, 75% to F | 88.0% / 24.5% | 97.7% / 30.8% | 99.4% / 32.6% | 99.7% / 33.2% | [claim:C154] [evidence:E014,E009]
| J + uncapped screening | 89.5% / 65.4% | 98.2% / 70.3% | 99.6% / 71.3% | 99.8% / 71.4% | [claim:C155] [evidence:E014,E009]

Each cell gives the fall in new infections / the fall in HCV deaths (medians across posterior draws). The last column is the main analysis (Table 1). Targets: 80% and 65%. [claim:C156] [evidence:E014,E009]

**Table 3. Pareto-efficient strategies over 2026–2050: decision criteria, TOPSIS closeness with equal weights and probabilities across posterior draws.** [claim:C157] [evidence:E018,E020,E013]

| Policy | Infections averted | Liver deaths averted | Draws meeting both 2030 targets | Programme budget, S$ million | Treatment courses | TOPSIS closeness | Pareto-efficient, % of draws | Ranked first, % of draws | [claim:C158] [evidence:E018,E020,E013]
|---|---|---|---|---|---|---|---|---|
| J + 8,000/year, 75% to F | 12,950 (11,111–15,284) | 2,197 (1,853–2,613) | 0.0% | 311.1 (266.2–365.9) | 7,752 (6,596–9,155) | 0.565 | 100.0% | 99.2% | [claim:C159] [evidence:E018,E020,E013]
| J + 5,000/year, 75% to F | 12,593 (10,856–14,825) | 2,045 (1,730–2,433) | 0.0% | 287.9 (245.6–339.3) | 7,220 (6,132–8,552) | 0.559 | 100.0% | 0.8% | [claim:C160] [evidence:E018,E020,E013]
| J + uncapped screening | 12,981 (11,132–15,425) | 3,085 (2,574–3,766) | 100.0% | 543.0 (460.0–652.8) | 11,419 (9,630–13,783) | 0.551 | 100.0% | 0.0% | [claim:C161] [evidence:E018,E020,E013]
| J + 2,000/year, 75% to F | 11,678 (10,095–13,749) | 1,842 (1,561–2,196) | 0.0% | 268.1 (227.1–321.2) | 6,795 (5,733–8,172) | 0.524 | 99.9% | 0.0% | [claim:C162] [evidence:E018,E020,E013]
| J + 8,000/year, 50% to F | 12,651 (10,897–14,892) | 2,453 (2,066–2,922) | 0.0% | 365.8 (312.0–432.6) | 9,243 (7,842–11,000) | 0.509 | 100.0% | 0.0% | [claim:C163] [evidence:E018,E020,E013]
| J + 5,000/year, 50% to F | 12,198 (10,524–14,366) | 2,253 (1,900–2,676) | 0.0% | 337.3 (286.6–398.8) | 8,651 (7,318–10,279) | 0.501 | 100.0% | 0.0% | [claim:C164] [evidence:E018,E020,E013]
| J + 1,000/year, 75% to F | 11,019 (9,482–13,048) | 1,755 (1,488–2,101) | 0.0% | 265.6 (224.0–320.1) | 6,802 (5,707–8,244) | 0.493 | 86.5% | 0.0% | [claim:C165] [evidence:E018,E020,E013]
| J + 8,000/year, 25% to F | 11,969 (10,327–14,100) | 2,597 (2,183–3,103) | 0.0% | 409.6 (347.3–486.5) | 10,387 (8,766–12,434) | 0.428 | 100.0% | 0.0% | [claim:C166] [evidence:E018,E020,E013]

Medians (95% UI) across posterior draws. Outcomes averted are against no treatment, summed over 2026–2050 and not discounted. The programme budget covers testing, consultations and drugs in J, F and X, discounted at 3% a year to 2026; it excludes care for advanced liver disease. Strategies are ordered by TOPSIS closeness. [claim:C167] [evidence:E018,E020,E013,E001]

# Supporting information

**S1 Text. Technical appendix.** Model equations, parameters, prison rate estimation, calibration and its checks, the 2015 model population, scenario implementation, costs, additional results and the sensitivity analysis with prison treatment from 2019. [claim:C168] [evidence:E039]

**S1 Checklist. CHEERS 2022 checklist.** [claim:C169] [evidence:E042]
