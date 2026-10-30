---
title: "Prison treatment and post-release screening for hepatitis C among people who inject drugs in Singapore: A modelling study"
---

**Short title:** Hepatitis C treatment after prison release in Singapore

Yichen He\*

[Affiliation: to be provided by the authors]

\* Corresponding author. E-mail: [to be provided by the authors]

# Abstract

## Background

Hepatitis C virus (HCV) infection in Singapore is concentrated among people who inject drugs (PWID), many of whom cycle through prison. We asked whether treatment in prison, alone or with community screening of released people from 2026, can meet the World Health Organization 2030 targets (80% fewer new infections and 65% fewer HCV deaths than in 2015), and how limited screening should be shared after release.

## Methods and findings

We built an age-structured model of PWID never imprisoned (D), in prison (J), and released and still injecting (F) or stopped (X), with prison rates estimated from records of 20,720 people. We calibrated it by Markov chain Monte Carlo to 2015 prison seroprevalence, first arrests of new heroin and methamphetamine users and a survey estimate of active injectors. In the main analysis nothing changed from 2015 until 2026, when treatment in prison and 15 community screening policies (600 to 8,000 tests a year, 25%–75% to F) began; we varied prison treatment coverage and, in a sensitivity analysis, started prison treatment in 2019. Results are over 2,000 posterior draws.

Prison treatment alone from 2026 reduced new infections in 2030 by 49.8% (95% uncertainty interval [UI] 37.7%–56.9%) and HCV deaths by 32.6% (95% UI 30.6%–34.8%). Five of 15 screening policies met the incidence target, but none met the mortality target; the best, 8,000 tests a year with 25% to F, reduced deaths by 62.6% (95% UI 59.0%–65.6%), and only uncapped screening of all released people met both targets. Moving tests from F to X increased the fall in deaths and reduced the fall in infections. Lower prison treatment coverage reduced every effect. Had prison treatment started in 2019, it alone would have reduced new infections by 83.4% (95% UI 67.8%–89.5%), meeting the incidence target, and 5,000 tests a year with 25% to F would have met both targets in 90% of draws. Limitations include a modelled prison population (8,493) far above the 3,367 screened people recorded, weakly identified entry parameters, charge-based costs and parameter-only uncertainty.

## Conclusions

In this model, starting HCV treatment for PWID in Singapore in 2026 is too late for prison treatment and capped community screening to meet the 2030 target for HCV deaths; the incidence target needs thousands of community tests a year in addition to prison treatment.

# Author summary

## Why was this study done?

- Most hepatitis C infections in Singapore are among people who inject drugs, and many of them spend time in prison.

- It is unclear whether treatment in prison, with or without community screening after release, can meet the WHO 2030 elimination targets, which ask for large falls in both new infections and deaths.

- Community screening capacity is limited, and it is not known whether it should go to released people who still inject or to those who have stopped.

## What did the researchers do and find?

- We built a model of hepatitis C among people who inject drugs that follows them into and out of prison, and fitted it to prison and national data from 2015.

- If treatment starts in 2026, prison treatment alone cut new infections by about 50% and deaths by about 33% by 2030. With up to 8,000 community tests a year, several policies met the infection target but none met the death target.

- If prison treatment had started in 2019, it alone would have met the infection target, and 5,000 or more community tests a year, mostly for people who have stopped injecting, would have met the death target.

## What do these findings mean?

- Time matters: the later treatment starts, the less any policy can achieve by 2030.

- Prison treatment prevents many new infections, but most hepatitis C deaths are among released people who have stopped injecting, who need community testing.

- These results come from a model with important uncertainties, including the size of the community population and the costs, and should be checked against programme data.

# Abbreviations

CHEERS, Consolidated Health Economic Evaluation Reporting Standards; D, people who inject drugs and have never been imprisoned; F, released people who still inject; HCV, hepatitis C virus; J, people in prison; MCMC, Markov chain Monte Carlo; PWID, people who inject drugs; UI, uncertainty interval; WHO, World Health Organization; X, released people who have stopped injecting.

# Introduction

WHO set targets for 2030 to eliminate viral hepatitis as a public health threat, including large reductions in new HCV infections and in HCV deaths compared with 2015 [@who2022strategy; @who2023validation].
Direct-acting antivirals cure most people treated [@pawlotsky2020], and treating people who inject drugs can also prevent onward transmission [@martin2013; @dore2021].

In Singapore, people who inject drugs are a key population for HCV [@winslow2007; @cda_hcv], and their number has been estimated at about 11,000 [@teo2019].
Many of them pass through prison, where treatment with direct-acting antivirals has been offered [@prison2021; @tan2023].

Prison offers a setting where people can be found and treated, but linkage to care after release is often poor [@clink2022; @eckhardt2022].
Community testing programmes in Singapore have been described in policy statements [@moh2024ett; @moh2025ett; @cda_programme], and point-of-care testing can improve linkage to care [@hsiang2022].
A previous model examined what is needed to eliminate HCV in Singapore [@chaillon2021]. Our model adds the distinction between released people who return to injecting and those who stop.

This distinction matters for policy. People who still inject keep transmitting and are often re-arrested, so prison treatment reaches them again. People who have stopped do not transmit, but they may carry long-standing infections that progress to cirrhosis and liver cancer, and they are less likely to return to prison.

We used a calibrated transmission model to estimate whether prison treatment and community screening starting in 2026 can meet the WHO 2030 targets for PWID in Singapore, how a fixed annual number of community screening tests should be split between released people who still inject and those who have stopped, how the answer depends on prison treatment coverage, and what an earlier start in prison would have achieved.

# Methods

## Study design

This is a modelling study with a costing of the screening and treatment policies. We report it following the Consolidated Health Economic Evaluation Reporting Standards 2022 (CHEERS 2022) [@cheers2022] where they apply (S1 Checklist).

[Protocol registration: to be confirmed by the authors. No registered protocol is recorded in the project repository.]

## Model structure

We built a deterministic compartmental model of PWID in Singapore (Fig 1). People enter the model when they start injecting and are placed in group D (injecting, never imprisoned). First arrest moves them to J (in prison). On release, a share returns to injecting (F) and the rest stop (X). People in F can be re-arrested; people in X are not.

Each group was divided into six age groups (<20, 20–29, 30–39, 40–49, 50–59 and 60 years or older), four liver disease stages (non-cirrhotic, compensated cirrhosis, decompensated cirrhosis and hepatocellular carcinoma) and five HCV states: susceptible, acute, chronic, on treatment and recovered. Recovered people, after spontaneous clearance or cure, remain antibody positive and can be reinfected.

Transmission occurs only among people who inject (D and F). The force of infection in each age group depends on an age contact matrix and on the share of infectious (acute or chronic) people in each age group of D and F. There is no transmission in prison or after people stop injecting.

Natural history, treatment efficacy and mortality parameters came from the literature (S1 Text, Table A). Background mortality was the 2015 Singapore death rate by age, multiplied by a standardised mortality ratio for people who use drugs of 9.52 for viraemic people and 4.73 for others [@vansanten2014]. Decompensated cirrhosis and hepatocellular carcinoma add liver deaths, which are lower after cure.

The model was solved with a fourth-order Runge–Kutta method with a daily step. Each run started from a small seeded population and ran for 150 model years, so that 2015 was a steady state; runs whose last ten years changed by more than set tolerances were excluded (S1 Text, Section 1).

## Prison rates

We estimated release and re-arrest rates by age from Singapore Prison Service records of 20,720 people held for drug consumption with admissions from 2000 to 2017 (33,557 prison stays), using Bayesian survival models fitted in Stan (S1 Text, Section 3). Release rates were 0.33 to 0.64 per year by age. A cure-type model gave a 72.4% probability of eventual re-arrest (95% credible interval 71.4%–73.4%), and re-arrest rates of 0.27 to 0.55 per year among those re-arrested.

We assumed that people who are eventually re-arrested return to injecting (F) and that the others stop (X). The first-arrest rate from D was the re-arrest rate multiplied by a calibrated factor c.

## Entry and calibration targets

The Central Narcotics Bureau recorded 1,309 new drug abusers in 2015 (257, 644, 292, 94, 21, 1 by age group) [@singstat_drug]. Of the 1,309 new abusers that year, 113 mainly used heroin and 1,014 methamphetamine [@singstat_m891521], the drugs most often injected in Singapore; we took this share (86.1%) as the share of new abusers who inject, the same at every age, because the drug table has no age breakdown.
Entry into D by age was these counts multiplied by 86.1% (1,127 a year in total) and by a calibrated scale of at least 1, and was held constant in every year.

We calibrated the model to three 2015 targets (S1 Text, Table C): HCV antibody prevalence by age among people in prison at the end of 2015 with a result from the December 2014 to February 2016 screening exercise (3,367 people); first arrests by age, the same counts multiplied by 86.1% and rounded (221, 554, 251, 81, 18, 1); and 11,000 active injectors in the community (95% confidence interval 6,500–17,000), a network scale-up estimate [@teo2019], compared with the modelled number in D and F.

## Calibration

We estimated 13 parameters on the log scale: six multipliers of the contact matrix (one per age group), six entry scales and the first-arrest factor c. Priors were fixed before fitting: normal with mean 0 and standard deviation 1.5 for the contact multipliers, half-normal with standard deviation 1 for the entry scales, and normal with mean 0 and standard deviation 0.5 for log c, each truncated to wide bounds (S1 Text, Table D).

The likelihood combined a beta distribution for each seroprevalence (precision 100), a negative binomial distribution for each first-arrest count (size 50) and a log-normal distribution for the community estimate (log-scale standard deviation 0.2453, taken from its confidence interval).

A multistart Nelder–Mead search found the posterior mode. We then sampled the posterior with an adaptive random-walk Metropolis algorithm [@haario2001]: 8 chains started from independent prior draws, each with 60,000 model runs, of which 20,000 were warm-up. Convergence was checked with rank-normalised R-hat and effective sample sizes [@vehtari2021].

We checked the sampler with simulation-based calibration [@talts2018] on simulated data sets (S1 Text, Section 4).

## Policy scenarios

In the main analysis, nobody was treated before 2026, so the model stayed at its 2015 equilibrium until then. From 2026, people with chronic infection in J were treated at a rate that gives a 99.995% annual probability of starting treatment, and community screening of F and X began. We compared no treatment, prison treatment only (J only), 15 finite policies (600, 1,000, 2,000, 5,000, 8,000 tests per year, with 25%, 50% or 75% given to F and the rest to X), and uncapped screening of everyone in F and X each year as a reference.

To show how much depends on prison treatment, we repeated all 17 policies with an annual probability of starting treatment in J of 0%, 50%, 90% and 99.995% (68 combinations).

In a sensitivity analysis, prison treatment started in 2019 instead, at the main coverage, and community screening still started in 2026.

The smallest capacity, 600 tests per year, was set with reference to a public community testing programme [@moh2024ett]; the other capacities and shares are illustrative and are not announced policies.

At the start of each year, the tests given to F or X were applied to people not on treatment, in proportion to their HCV states. Each antibody-positive person had an RNA test, and the number of chronic infections found set that year's treatment capacity in the group. Within that capacity, people with chronic infection started treatment at a rate that gives a 99% annual probability. Sofosbuvir–velpatasvir without ribavirin was assumed, and people in D were never treated.

## Outcomes

The primary outcomes were the reductions in 2030, compared with 2015, of annual new HCV infections among PWID and of annual HCV deaths. New infections were the number of infections in D and F during the calendar year. HCV deaths were the liver deaths in decompensated cirrhosis and hepatocellular carcinoma during the year, plus the excess non-liver deaths of viraemic people over the rate of non-viraemic people. The 2015 values were the yearly numbers at the calibrated equilibrium.

A policy met a target if the median reduction was at least 80% for infections or 65% for deaths. We also report the share of posterior draws in which both targets were met.

## Costs

We costed antibody tests (S$25.31), RNA tests (S$307), one consultation per treatment start (S$25) and a 12-week sofosbuvir–velpatasvir course (S$40,405) in 2026 Singapore dollars, from published charges and a published price (S1 Text, Table H) [@wong2019; @nhg_dsc; @parkway_rna; @sgd_rate].

Costs for 2026–2030 in all three groups (J, F and X) were discounted at 3% a year to 2026. Unit prices were varied together by a factor of 0.5 and 1.5 in a sensitivity analysis. These are charges, not costs to the health system.

## Uncertainty

For each of 2,000 posterior draws we ran every scenario and computed every outcome. We report medians and 95% uncertainty intervals (UI; 2.5th to 97.5th percentiles) across draws. These intervals reflect parameter uncertainty only.

## Software

The model was written in C++ and called from R through Rcpp [@rcpp2011]; rate models used Stan. Code and generated results are in the project repository.

## Ethics statement

[Ethics statement: to be completed by the authors. Record-level Singapore Prison Service data were used to estimate prison rates and the 2015 seroprevalence targets; the approving ethics committee, the approval number and the consent arrangements are not recorded in the project repository.]

# Results

## Calibration and the 2015 population

The posterior chains converged (largest R-hat 1.004; smallest bulk effective sample size 2,137). The model reproduced the age pattern of seroprevalence in prison, and the 95% posterior intervals of all 13 modelled targets overlapped the target ranges (Fig 2).

At ages 30 and over, modelled first arrests were higher than the targets (for example 281 against 251 at ages 30–39 at the posterior mode), because entry could not fall below the scaled counts and the entry scales settled near their lower bound (posterior medians 1.09 to 1.17 below age 60).

The modelled number of active injectors in the community was 7,005 (6,125–8,319), below the survey estimate of 11,000 but within its range. The first-arrest factor c was 2.59 (95% credible interval 1.45–5.24), so first arrest was faster than re-arrest.

In 2015 the model had 50,494 people who had ever injected (95% credible interval 45,344–58,894): 1,040 in D, 8,493 in J, 5,934 in F and 34,964 in X. The modelled prison population was not a calibration target and was much larger than the 3,367 people with a screening result recorded in prison in 2015.

In 2015 the model had 552 (473–654) new infections and 266 (226–317) HCV deaths per year. Without treatment, both stayed at their 2015 level until 2030 (within 0.04% until 2025).

## Prison treatment alone

With treatment in prison from 2026, 3,047 people started treatment in J in 2026, as the infections accumulated since 2015 were treated, and fewer in later years. New infections fell by 49.8% (95% UI 37.7%–56.9%) by 2030 and HCV deaths by 32.6% (95% UI 30.6%–34.8%) (Fig 3). No posterior draw met either target.

In 2026, 71% of the 6,748 people with chronic infection in F and X were in X, and F and X together accounted for 91% of all HCV deaths. People in X are never re-arrested in the model, so prison treatment does not reach them.

People in D, who had never been imprisoned, accounted for 15% of new infections in 2030 under prison treatment only; no policy in this study reached them.

## Community screening policies

Five of the 15 finite screening policies met the incidence target at the median (Fig 4 and Table 1): 8,000 tests a year with any share to F, and 5,000 tests with 50% or 75% to F. None met the mortality target. The largest fall in deaths among finite policies was 62.6% (95% UI 59.0%–65.6%), with 8,000 tests and 25% to F; both targets were met in 5.7% of draws, and in no draw for any other finite policy.

At every capacity, giving a larger share of tests to F increased the fall in new infections and reduced the fall in deaths. All quotas were used in full, except the F quota with 8,000 tests and 75% to F (99.5% used over 2026–2030), where screening ran out of eligible people in F.

At the start of 2026, 5,951 people in F were eligible for screening, of whom 33.3% had chronic infection; in X, 35,052 were eligible and 13.6% had chronic infection. A test given to F therefore found more chronic infections, in people who could transmit, but people in F are re-arrested and can be treated in prison, whereas people in X can be reached only by community screening; this is a likely reason why tests given to X prevented more deaths.

Uncapped screening of everyone in F and X each year reduced new infections by 99.7% (95% UI 97.9%–99.9%) and deaths by 85.1% (95% UI 84.1%–86.1%), and was the only policy that met both targets. Even then, 12% of the remaining infections in 2030 were in D.

## Prison treatment coverage

Every reduction rose with the coverage of prison treatment (Fig 5 and Table 2). Without prison treatment, deaths fell by 43.2% (95% UI 38.7%–46.5%) with 8,000 tests and 25% to F, and only 8,000 tests with 75% to F met the incidence target, reducing new infections by 88.2% (95% UI 82.3%–90.8%). At 50% coverage, three finite policies met the incidence target, and at 90% the same five as at full coverage.

At no coverage did a finite policy meet the mortality target, while uncapped screening met it even without prison treatment, reducing deaths by 77.4% (95% UI 76.9%–77.9%).

## Sensitivity analysis: prison treatment from 2019

Had prison treatment started in 2019, prison treatment alone would have reduced new infections by 83.4% (95% UI 67.8%–89.5%) and deaths by 53.2% (95% UI 50.7%–55.9%) by 2030, meeting the incidence target at the median. All 15 finite policies would then have met the incidence target, and three would have met the mortality target: 5,000 tests with 25% to F, which reduced deaths by 66.9% (95% UI 64.1%–69.3%) and met both targets in 90% of draws, and 8,000 tests with 25% or 50% to F, which met both targets in 99.7% and 98.2% of draws (S1 Text, Section 9).

## Costs

Discounted costs for 2026–2030 were S$203.6 million (95% UI 172.8–243.6) with prison treatment only, mostly for drugs and mostly in 2026. The cheapest policy that met the incidence target, 5,000 tests a year with 75% to F, cost S$232.7 million (95% UI 198.2–276.9), or S$29.0 million (95% UI 24.5–34.2) more than prison treatment only.

Uncapped screening, the only policy that met both targets, cost S$419.6 million (95% UI 350.5–506.1). Drugs made up 98% of the cost of 5,000 tests a year with 25% to F, so totals scale almost directly with the drug price; S1 Text, Section 8 gives all policies.

# Discussion

In this model of PWID in Singapore, starting HCV treatment in 2026 left too little time to meet the WHO 2030 mortality target with prison treatment and capped community screening. Prison treatment from 2026 reduced new infections by about half and deaths by about a third. Up to 8,000 community tests a year could meet the incidence target, but only screening of every released person each year met the mortality target.

Timing was decisive. With the same prison treatment starting in 2019, prison treatment alone met the incidence target, and several thousand community tests a year met the mortality target. Each year without treatment adds to the pool of chronic infections that progress to cirrhosis and liver cancer, and deaths fall more slowly than infections after treatment begins.

The two WHO targets pull screening in different directions. Tests given to people who still inject prevent more infections, while tests given to people who have stopped prevent more deaths. Because the mortality target was harder to meet in every analysis, screening policy should be judged mainly by the deaths it prevents.

This suggests that the people who matter most for HCV mortality are not those most visible to drug services. People who have stopped injecting may no longer be in contact with prison or drug services, and reaching them may need testing in primary care or other community settings [@hsiang2022].

Part of the modelled HCV mortality is excess non-liver death among viraemic people, based on mortality ratios from a cohort of people who use drugs [@vansanten2014]. Evidence that viraemic infection raises mortality from causes other than liver disease supports this [@lee2012; @cepeda2022], but the excess is not found in every cohort [@piselli2021]; with a smaller excess, deaths would depend more on cirrhosis and liver cancer.

We extend earlier modelling for Singapore [@chaillon2021] by separating released people by injecting status and by modelling screening capacity explicitly. This separation gives a reason why the mortality target may be harder to meet than the incidence target.

## Strengths and limitations

The model was calibrated with a full Bayesian method that converged well, and every policy result carries posterior uncertainty from 2,000 draws. Prison flows were estimated from individual prison records rather than assumed, and screening capacity was modelled explicitly, so treatment depends on how many people are found.

The study has important limitations. First, the modelled prison population (8,493) was about 2.5 times the number of screened people recorded in prison in 2015. The prison size was not fitted, and the model may keep people in prison or in contact with it for too long, which would favour prison treatment.

Second, the main analysis assumes that nobody was treated before 2026. Treatment in prison has in fact been offered [@prison2021; @tan2023], but its start and coverage are not recorded in the project data. The sensitivity analysis with full coverage from 2019 shows how much an earlier start changes the results; the real history probably lies between the two analyses.

Third, the entry and first-arrest parameters were weakly identified. We took the share of new drug abusers who mainly use heroin or methamphetamine (86.1%) as the share who inject, the same at every age; not all of them inject, and some users of other drugs do. The entry scales settled near their lower bound, and the posterior of c was only modestly narrower than its prior (standard deviation 78% of the prior's).

Fourth, we assumed that released people who are eventually re-arrested return to injecting and that the others stop permanently, and that there is no transmission in prison. Relapse after a long gap would put more infections back into F.

Fifth, screening was an aggregate approximation: tests were spread across each group in proportion to its HCV states, with perfect linkage within capacity and no loss to follow-up. Real programmes will find and treat fewer people for the same number of tests.

Sixth, costs were published charges and a published drug price carried forward without inflation, not payer costs, and follow-up visits were excluded. Our uncertainty intervals cover model parameters only, not structural assumptions, prices or uptake.

Finally, the seroprevalence targets could not be reproduced exactly from the current processing of the prison records (largest difference 0.5 percentage points).

## Conclusions

In this model, treatment for PWID in Singapore that starts in 2026 can meet the 2030 target for new HCV infections only with prison treatment and thousands of community tests a year, and cannot meet the target for HCV deaths without screening every released person. Reducing deaths requires finding and treating people who have stopped injecting after release, and every year of delay makes the target harder to reach.

# Data availability

[Data availability statement: to be completed by the authors. Model code, aggregate inputs and all generated results are in the project repository (link and DOI to be added). The record-level Singapore Prison Service data cannot be shared by the authors; the access conditions and the contact for access requests must be stated.]

# Funding

[Funding statement: to be completed by the authors, including grant numbers and the role of any funder, or the PLOS statement for unfunded work.]

# Competing interests

[Competing interests: to be declared by each author.]

# Author contributions

[CRediT author contributions: to be completed by the authors.]

# Use of artificial intelligence

[AI use statement: to be completed and approved by the authors. Draft: Claude (Anthropic) was used to write and revise analysis code and to draft and revise the text, figures and supporting information, with procedures from Scientific Agent Skills [@kassis2026]. The authors are to state that they checked all content, data, code, references and conclusions, and take responsibility for them.]

# Acknowledgments

[Acknowledgments: to be completed by the authors.]

# References

::: {#refs}
:::

# Figure captions

**Fig 1. Model structure.** People who inject drugs enter D (never imprisoned) and move to J (prison) at first arrest. On release they move to F (still injecting) or X (stopped). People in F can be re-arrested. Transmission occurs only in D and F. In the main analysis, treatment is given in J and community screening is offered to F and X from 2026.

**Fig 2. Calibration to the 2015 targets.** (A) HCV antibody prevalence in prison by age. (B) First arrests by age (Central Narcotics Bureau new drug abusers, scaled to heroin and methamphetamine users). (C) Active injectors in the community (D and F). Black: target with the 95% range of its likelihood; blue: 95% posterior interval of the model value from the MCMC posterior.

**Fig 3. New infections and HCV deaths, 2015–2030.** Annual new HCV infections (A) and HCV deaths (B) relative to 2015, without treatment, with prison treatment only, and with 8,000 community tests a year given 25% or 75% to F, all from 2026. Yearly counts are drawn at mid-year and joined by a monotone curve; the vertical line marks 2026. Lines are medians and bands are 95% uncertainty intervals across 2,000 posterior draws. Dashed lines mark the WHO targets.

**Fig 4. Reductions in 2030 by screening capacity and share of tests given to F.** (A) New infections. (B) HCV deaths. Points are medians and bars are 95% uncertainty intervals. Dashed lines mark the targets; the dotted line is prison treatment only and the solid line is uncapped screening.

**Fig 5. Reductions in 2030 by prison treatment coverage.** Median reductions in new infections (A) and HCV deaths (B) compared with 2015, for each policy (rows) and annual probability of starting treatment in prison from 2026 (columns). Bold values meet the target.

# Tables

**Table 1. Reductions in 2030 compared with 2015, probability of meeting both WHO targets, and discounted costs for 2026–2030 (main analysis).**

| Policy | Fall in new infections, % (95% UI) | Fall in HCV deaths, % (95% UI) | Draws meeting both targets | Cost, S$ million (95% UI) |
|---|---|---|---|---|
| J only | 49.8% (37.7%–56.9%) | 32.6% (30.6%–34.8%) | 0.0% | 203.6 (172.8–243.6) |
| J + 600/year, 25% to F | 53.4% (41.1%–60.4%) | 35.7% (33.6%–37.9%) | 0.0% | 215.4 (182.8–256.9) |
| J + 600/year, 50% to F | 56.8% (44.6%–63.6%) | 35.3% (33.2%–37.6%) | 0.0% | 212.7 (180.8–254.4) |
| J + 600/year, 75% to F | 60.2% (48.0%–66.8%) | 34.9% (32.7%–37.3%) | 0.0% | 210.0 (178.5–251.1) |
| J + 1,000/year, 25% to F | 55.7% (43.4%–62.5%) | 37.6% (35.4%–39.8%) | 0.0% | 222.8 (189.6–265.6) |
| J + 1,000/year, 50% to F | 61.3% (49.2%–67.8%) | 37.0% (34.8%–39.3%) | 0.0% | 218.4 (185.5–260.6) |
| J + 1,000/year, 75% to F | 66.4% (54.6%–72.6%) | 36.3% (34.0%–38.7%) | 0.0% | 213.6 (181.6–255.6) |
| J + 2,000/year, 25% to F | 61.3% (49.2%–67.8%) | 42.2% (39.8%–44.6%) | 0.0% | 240.2 (205.3–285.3) |
| J + 2,000/year, 50% to F | 71.1% (59.6%–76.9%) | 41.0% (38.4%–43.4%) | 0.0% | 231.2 (197.2–275.5) |
| J + 2,000/year, 75% to F | 79.2% (68.6%–84.2%) | 39.2% (36.6%–41.9%) | 0.0% | 220.6 (187.5–263.8) |
| J + 5,000/year, 25% to F | 75.3% (64.2%–80.8%) | 53.8% (50.7%–56.7%) | 0.0% | 283.9 (242.0–335.3) |
| J + 5,000/year, 50% to F | 90.3% (82.8%–93.5%) | 50.3% (47.2%–53.1%) | 0.0% | 260.9 (222.8–309.6) |
| J + 5,000/year, 75% to F | 96.9% (93.1%–98.4%) | 45.1% (42.2%–48.0%) | 0.0% | 232.7 (198.2–276.9) |
| J + 8,000/year, 25% to F | 85.5% (76.4%–89.6%) | 62.6% (59.0%–65.6%) | 5.7% | 316.4 (268.7–373.8) |
| J + 8,000/year, 50% to F | 97.6% (94.3%–98.9%) | 57.2% (53.8%–60.1%) | 0.0% | 282.8 (240.8–334.7) |
| J + 8,000/year, 75% to F | 99.7% (97.8%–99.9%) | 49.0% (46.1%–51.9%) | 0.0% | 243.4 (208.1–288.8) |
| J + uncapped screening | 99.7% (97.9%–99.9%) | 85.1% (84.1%–86.1%) | 100.0% | 419.6 (350.5–506.1) |

UI, uncertainty interval across posterior draws. Costs are 2026 Singapore dollars for J, F and X, discounted at 3% a year to 2026.

**Table 2. Median falls in new infections and HCV deaths in 2030 by annual probability of starting treatment in prison from 2026.**

| Policy | 0% in prison | 50% in prison | 90% in prison | 99.995% in prison |
|---|---|---|---|---|
| J only | -0.0% / -0.1% | 22.6% / 21.7% | 40.1% / 29.5% | 49.8% / 32.6% |
| J + 600/year, 25% to F | 2.7% / 4.2% | 26.5% / 25.3% | 44.0% / 32.7% | 53.4% / 35.7% |
| J + 600/year, 50% to F | 5.6% / 4.9% | 30.5% / 25.4% | 47.8% / 32.5% | 56.8% / 35.3% |
| J + 600/year, 75% to F | 8.7% / 5.4% | 34.5% / 25.4% | 51.5% / 32.3% | 60.2% / 34.9% |
| J + 1,000/year, 25% to F | 4.6% / 7.0% | 29.2% / 27.6% | 46.6% / 34.8% | 55.7% / 37.6% |
| J + 1,000/year, 50% to F | 9.8% / 8.0% | 35.8% / 27.7% | 52.7% / 34.4% | 61.3% / 37.0% |
| J + 1,000/year, 75% to F | 15.3% / 8.7% | 42.2% / 27.6% | 58.5% / 33.9% | 66.4% / 36.3% |
| J + 2,000/year, 25% to F | 9.8% / 13.5% | 35.8% / 33.1% | 52.7% / 39.7% | 61.3% / 42.2% |
| J + 2,000/year, 50% to F | 21.0% / 15.1% | 48.4% / 33.0% | 63.9% / 38.8% | 71.1% / 41.0% |
| J + 2,000/year, 75% to F | 32.5% / 15.9% | 59.8% / 32.3% | 73.4% / 37.4% | 79.2% / 39.2% |
| J + 5,000/year, 25% to F | 26.8% / 30.3% | 54.2% / 46.9% | 68.8% / 52.0% | 75.3% / 53.8% |
| J + 5,000/year, 50% to F | 53.6% / 31.9% | 77.5% / 45.4% | 87.0% / 49.1% | 90.3% / 50.3% |
| J + 5,000/year, 75% to F | 72.7% / 30.0% | 90.3% / 41.3% | 95.6% / 44.1% | 96.9% / 45.1% |
| J + 8,000/year, 25% to F | 43.5% / 43.2% | 69.5% / 57.4% | 81.1% / 61.4% | 85.5% / 62.6% |
| J + 8,000/year, 50% to F | 75.5% / 43.5% | 91.9% / 54.0% | 96.5% / 56.4% | 97.6% / 57.2% |
| J + 8,000/year, 75% to F | 88.2% / 37.8% | 97.7% / 46.5% | 99.4% / 48.4% | 99.7% / 49.0% |
| J + uncapped screening | 88.9% / 77.4% | 97.9% / 84.1% | 99.5% / 85.0% | 99.7% / 85.1% |

Each cell gives the fall in new infections / the fall in HCV deaths (medians across posterior draws). The last column is the main analysis (Table 1). Targets: 80% and 65%.

# Supporting information

**S1 Text. Technical appendix.** Model equations, parameters, prison rate estimation, calibration and its checks, the 2015 model population, scenario implementation, costs, additional results and the sensitivity analysis with prison treatment from 2019.

**S1 Checklist. CHEERS 2022 checklist.**
