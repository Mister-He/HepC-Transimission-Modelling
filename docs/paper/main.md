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

Prison treatment alone from 2026 reduced new infections in 2030 by 51.5% (95% uncertainty interval [UI] 40.6%–57.3%) and HCV deaths by 22.1% (95% UI 21.2%–23.0%). Five of 15 screening policies met the incidence target, but none met the mortality target; the best, 8,000 tests a year with 25% to F, reduced deaths by 44.0% (95% UI 41.4%–45.9%). Only uncapped screening of all released people met both targets, reducing deaths by 71.4% (95% UI 71.3%–71.6%), but not if cure lowered the risk of liver death in decompensated cirrhosis less than assumed. Moving tests from F to X increased the fall in deaths and reduced the fall in infections. Lower prison treatment coverage reduced every effect. Had prison treatment started in 2019, it alone would have reduced new infections by 85.0% (95% UI 72.7%–89.7%), meeting the incidence target, but no capped policy would have met the mortality target; the best reduced deaths by 56.1% (95% UI 54.0%–57.7%). Limitations include treating new heroin and methamphetamine users as the people who inject, weakly identified entry parameters, liver disease inputs from other settings, to which the fall in deaths was sensitive, charge-based costs, and uncertainty intervals that cover the calibrated parameters only.

## Conclusions

In this model, prison treatment and capped community screening could not meet the 2030 target for HCV deaths among PWID in Singapore, even with prison treatment from 2019; with treatment from 2026, the incidence target needs thousands of community tests a year in addition to prison treatment.

# Author summary

## Why was this study done?

- Most hepatitis C infections in Singapore are among people who inject drugs, and many of them spend time in prison.

- It is unclear whether treatment in prison, with or without community screening after release, can meet the WHO 2030 elimination targets, which ask for large falls in both new infections and deaths.

- Community screening capacity is limited, and it is not known whether it should go to released people who still inject or to those who have stopped.

## What did the researchers do and find?

- We built a model of hepatitis C among people who inject drugs that follows them into and out of prison, and fitted it to prison and national data from 2015.

- If treatment starts in 2026, prison treatment alone cut new infections by about 52% and deaths by about 22% by 2030. With up to 8,000 community tests a year, several policies met the infection target but none met the death target.

- If prison treatment had started in 2019, it alone would have met the infection target, but only testing every released person every year would have met the death target.

## What do these findings mean?

- Time matters: the later treatment starts, the less any policy can achieve by 2030. Deaths from liver disease fall slowly after cure, because people with advanced liver disease remain at risk.

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

Natural history, treatment efficacy and mortality parameters came from the literature (S1 Text, Table A). Each cause of death had one route in the model. Deaths from causes other than liver disease occurred at the 2015 Singapore death rate by age [@singstat_m810141] (ages 60–74 for the oldest group) multiplied by 4.73, the standardised mortality ratio of people who use drugs without HCV or HIV infection [@vansanten2014], and for viraemic people also by 1.47, the hazard ratio of extrahepatic death with viraemic HCV infection [@lee2012]. Liver deaths occurred only in decompensated cirrhosis and hepatocellular carcinoma, at lower rates after cure.

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

The primary outcomes were the reductions in 2030, compared with 2015, of annual new HCV infections among PWID and of annual HCV deaths. New infections were the number of infections in D and F during the calendar year. HCV deaths were the liver deaths in decompensated cirrhosis and hepatocellular carcinoma during the year, the deaths from HCV-related cirrhosis and liver cancer that the WHO target counts. The 2015 values were the yearly numbers at the calibrated equilibrium. As a secondary outcome, we added the extrahepatic excess deaths of viraemic people (S1 Text, Section 10).

A policy met a target if the median reduction was at least 80% for infections or 65% for deaths. We also report the share of posterior draws in which both targets were met.

## Costs

We costed antibody tests (S$25.31), RNA tests (S$307), one consultation per treatment start (S$25) and a 12-week sofosbuvir–velpatasvir course (S$40,405) in 2026 Singapore dollars, from published charges and a published price (S1 Text, Table H) [@wong2019; @nhg_dsc; @parkway_rna; @sgd_rate].

Costs for 2026–2030 in all three groups (J, F and X) were discounted at 3% a year to 2026. Unit prices were varied together by a factor of 0.5 and 1.5 in a sensitivity analysis. These are charges, not costs to the health system.

## Uncertainty

For each of 2,000 posterior draws we ran every scenario and computed every outcome. We report medians and 95% uncertainty intervals (UI; 2.5th to 97.5th percentiles) across draws, and the share of draws meeting a target. These describe uncertainty in the 13 calibrated parameters given the fixed inputs, not the probability that a target will be met.

To show how much the fixed inputs matter, we held the calibrated parameters at the posterior mode and set each fixed epidemiological input in turn to the low and high end of its literature range, without refitting (S1 Text, Section 13). In further sensitivity analyses, entry of new injectors from 2016 to 2025 followed the observed yearly counts of new drug abusers instead of the 2015 level (S1 Text, Section 11), and the share of infections with genotype 3 was 10% or 50% instead of 78% (S1 Text, Section 12).

## Software

The model was written in C++ and called from R through Rcpp [@rcpp2011]; rate models used Stan. Code and generated results are in the project repository.

## Ethics statement

[Ethics statement: to be completed by the authors. Record-level Singapore Prison Service data were used to estimate prison rates and the 2015 seroprevalence targets; the approving ethics committee, the approval number and the consent arrangements are not recorded in the project repository.]

# Results

## Calibration and the 2015 population

The posterior chains converged (largest R-hat 1.005; smallest bulk effective sample size 1,641). The model reproduced the age pattern of seroprevalence in prison, and the 95% posterior intervals of all 13 modelled targets overlapped the target ranges (Fig 2).

At ages 30 and over, modelled first arrests were higher than the targets (for example 279 against 251 at ages 30–39 at the posterior mode), because entry could not fall below the scaled counts and the entry scales settled near their lower bound (posterior medians 1.08 to 1.15 below age 60).

The modelled number of active injectors in the community was 7,332 (6,513–8,610), below the survey estimate of 11,000 but within its range. The first-arrest factor c was 2.71 (95% credible interval 1.54–5.23), so first arrest was faster than re-arrest.

In 2015 the model had 59,599 people who had ever injected (95% credible interval 53,968–69,269): 1,000 in D, 8,768 in J, 6,303 in F and 43,510 in X. The modelled prison population was not a calibration target; it was close to the 8,725 people held in prison for drug consumption at the end of 2015, the population whose release and re-arrest rates the model uses.

In 2015 the model had 529 (454–628) new infections and 157 (130–191) HCV deaths per year. Without treatment, both stayed at their 2015 level until 2030 (within 0.09% until 2025).

## Prison treatment alone

With treatment in prison from 2026, 3,189 people started treatment in J in 2026, as the infections accumulated since 2015 were treated, and fewer in later years. New infections fell by 51.5% (95% UI 40.6%–57.3%) by 2030 and HCV deaths by 22.1% (95% UI 21.2%–23.0%) (Fig 3). No posterior draw met either target.

In 2026, 73% of the 8,027 people with chronic infection in F and X were in X, and F and X together accounted for 90% of all HCV deaths. People in X are never re-arrested in the model, so prison treatment does not reach them.

People in D, who had never been imprisoned, accounted for 15% of new infections in 2030 under prison treatment only; no policy in this study reached them.

## Community screening policies

Five of the 15 finite screening policies met the incidence target at the median (Fig 4 and Table 1): 8,000 tests a year with any share to F, and 5,000 tests with 50% or 75% to F. None met the mortality target. The largest fall in deaths among finite policies was 44.0% (95% UI 41.4%–45.9%), with 8,000 tests and 25% to F, and no finite policy met both targets in any posterior draw.

At every capacity, giving a larger share of tests to F increased the fall in new infections and reduced the fall in deaths. All quotas were used in full.

At the start of 2026, 6,318 people in F were eligible for screening, of whom 33.8% had chronic infection; in X, 43,683 were eligible and 13.5% had chronic infection. A test given to F therefore found more chronic infections, in people who could transmit, but people in F are re-arrested and can be treated in prison, whereas people in X can be reached only by community screening; this is a likely reason why tests given to X prevented more deaths.

Uncapped screening of everyone in F and X each year reduced new infections by 99.8% (95% UI 98.7%–99.9%) and deaths by 71.4% (95% UI 71.3%–71.6%), and was the only policy that met both targets. Even then, 12% of the remaining infections in 2030 were in D. In the one-way analysis of the fixed inputs, the only input that stopped it meeting the mortality target was the relative risk of liver death after cure in decompensated cirrhosis: at the less favourable end of its range, deaths fell by 59.4% (S1 Text, Section 13).

## Prison treatment coverage

Every reduction rose with the coverage of prison treatment (Fig 5 and Table 2). Without prison treatment, deaths fell by 29.6% (95% UI 26.3%–32.0%) with 8,000 tests and 25% to F, and only 8,000 tests with 75% to F met the incidence target, reducing new infections by 88.0% (95% UI 82.8%–90.7%). At 50% coverage, three finite policies met the incidence target, and at 90% the same five as at full coverage.

At no coverage did a finite policy meet the mortality target. Uncapped screening met it at every coverage, though only just when there was no prison treatment, reducing deaths by 65.4% (95% UI 65.1%–65.7%).

## Sensitivity analysis: prison treatment from 2019

Had prison treatment started in 2019, prison treatment alone would have reduced new infections by 85.0% (95% UI 72.7%–89.7%) and deaths by 40.8% (95% UI 39.5%–42.1%) by 2030, meeting the incidence target at the median. All 15 finite policies would then have met the incidence target, but none would have met the mortality target: the largest fall in deaths was 56.1% (95% UI 54.0%–57.7%), with 8,000 tests and 25% to F. Only uncapped screening met both targets, reducing deaths by 76.4% (95% UI 76.1%–76.7%) (S1 Text, Section 9).

## Costs

Discounted costs for 2026–2030 were S$211.6 million (95% UI 180.3–252.4) with prison treatment only, mostly for drugs and mostly in 2026. The cheapest policy that met the incidence target, 5,000 tests a year with 75% to F, cost S$244.3 million (95% UI 208.5–288.3), or S$32.5 million (95% UI 27.7–37.8) more than prison treatment only.

Uncapped screening, the only policy that met both targets, cost S$484.1 million (95% UI 408.9–583.6). Drugs made up 98% of the cost of 5,000 tests a year with 25% to F, so totals scale almost directly with the drug price; S1 Text, Section 8 gives all policies.

# Discussion

In this model of PWID in Singapore, starting HCV treatment in 2026 left too little time to meet the WHO 2030 mortality target with prison treatment and capped community screening. Prison treatment from 2026 reduced new infections by about half and deaths by about a fifth. Up to 8,000 community tests a year could meet the incidence target, but only screening of every released person each year met the mortality target, and that depended on how much cure lowers the risk of liver death in decompensated cirrhosis.

Timing mattered most for new infections. With the same prison treatment starting in 2019, prison treatment alone met the incidence target, but even then no capped policy met the mortality target. Each year without treatment adds to the pool of chronic infections that progress to cirrhosis and liver cancer. Cure lowers but does not remove the risk of liver death in people who already have decompensated cirrhosis or liver cancer, and people with compensated cirrhosis can still progress, so liver deaths fall much more slowly than infections after treatment begins.

The two WHO targets pull screening in different directions. Tests given to people who still inject prevent more infections, while tests given to people who have stopped prevent more deaths. Because the mortality target was harder to meet in every analysis, screening policy should be judged mainly by the deaths it prevents.

This suggests that the people who matter most for HCV mortality are not those most visible to drug services. People who have stopped injecting may no longer be in contact with prison or drug services, and reaching them may need testing in primary care or other community settings [@hsiang2022].

We counted only liver deaths, as the WHO target does. Viraemic infection may also raise mortality from other causes [@lee2012; @cepeda2022], although the excess is not found in every cohort [@piselli2021]. Cure removes this excess at once, so counting it as well gave larger falls in deaths (prison treatment only 25.8%, uncapped screening 80.6%; S1 Text, Section 10); still, no finite policy reached 65%.

We extend earlier modelling for Singapore [@chaillon2021] by separating released people by injecting status and by modelling screening capacity explicitly. This separation gives a reason why the mortality target may be harder to meet than the incidence target.

## Strengths and limitations

The model was calibrated with a full Bayesian method that converged well, and every policy result carries posterior uncertainty from 2,000 draws. Prison flows were estimated from individual prison records rather than assumed, and screening capacity was modelled explicitly, so treatment depends on how many people are found.

The study has important limitations. First, the seroprevalence targets came from the 3,367 people in prison at the end of 2015 who had a screening result, out of the 8,725 held for drug consumption. The modelled prison population (8,768) matched the number held, although it was not fitted, but the people screened may not represent everyone in prison.

Second, the main analysis assumes that nobody was treated before 2026. Treatment in prison has in fact been offered [@prison2021; @tan2023], but its start and coverage are not recorded in the project data. The sensitivity analysis with full coverage from 2019 shows how much an earlier start changes the results; the real history probably lies between the two analyses.

Third, we did not observe who injects. We treated new drug abusers whose main drug was heroin or methamphetamine (86.1% of 1,309 in 2015) as the people who start injecting, at the same share in every age group because the drug table has no age breakdown. The records give the main drug, not the route of use: users of these drugs who smoke rather than inject inflate entry, and people who inject other drugs are missed. Entry could rise above this share through the calibrated entry scales but could not fall below it, and it was held at its 2015 level in every year; a fall in new injectors after 2015 would lower infections without treatment, and a rise would raise them.

Fourth, the entry and first-arrest parameters were weakly identified. The entry scales settled near their lower bound, and the posterior of c was only modestly narrower than its prior (standard deviation 72% of the prior's), as was the contact multiplier for ages 50–59. The first-arrest rate was assumed to have the same age pattern as the re-arrest rate, because the records contain only people who had already been arrested.

Fifth, we assumed that released people who are eventually re-arrested return to injecting and that the others stop permanently and are never re-arrested, and that there is no transmission in prison. Relapse after a long gap would put more infections back into F.

Sixth, screening was an aggregate approximation: tests were spread across each group in proportion to its HCV states, with perfect linkage within capacity and no loss to follow-up, and people who had never been imprisoned were never screened or treated. Real programmes will find and treat fewer people for the same number of tests.

Seventh, the model assumed that the epidemic was in a steady state in 2015, with constant prison rates, mixing that depends only on age, and the same mortality ratios for people who have stopped injecting as for those who still inject. Costs were published charges and a published drug price carried forward without inflation, not payer costs, and follow-up visits were excluded. Our uncertainty intervals cover the calibrated parameters only, not the fixed inputs, structural assumptions, prices or uptake, so they are too narrow. The one-way analysis shows that the fall in deaths depends most on the liver disease inputs, such as the risk of liver death after cure in decompensated cirrhosis (with uncapped screening, 59.4% to 78.9% across its range), which come from cohorts outside Singapore. S1 Text, Section 14 (Table O) lists every assumption with its basis and its likely effect on the results.

Finally, the seroprevalence targets could not be reproduced exactly from the current processing of the prison records (largest difference 0.5 percentage points).

## Conclusions

In this model, treatment for PWID in Singapore that starts in 2026 can meet the 2030 target for new HCV infections only with prison treatment and thousands of community tests a year, and cannot meet the target for HCV deaths without screening every released person; even then the result depends on liver disease inputs from other settings. Reducing deaths requires finding and treating people who have stopped injecting after release, and every year of delay makes the target harder to reach.

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
| J only | 51.5% (40.6%–57.3%) | 22.1% (21.2%–23.0%) | 0.0% | 211.6 (180.3–252.4) |
| J + 600/year, 25% to F | 54.9% (44.0%–60.4%) | 24.2% (23.2%–25.1%) | 0.0% | 223.6 (190.8–264.8) |
| J + 600/year, 50% to F | 58.2% (47.4%–63.6%) | 23.8% (22.9%–24.8%) | 0.0% | 221.1 (188.8–262.2) |
| J + 600/year, 75% to F | 61.4% (50.7%–66.6%) | 23.5% (22.5%–24.5%) | 0.0% | 218.6 (186.5–259.2) |
| J + 1,000/year, 25% to F | 57.1% (46.3%–62.5%) | 25.5% (24.4%–26.5%) | 0.0% | 231.2 (197.2–273.7) |
| J + 1,000/year, 50% to F | 62.4% (51.9%–67.6%) | 24.9% (23.9%–25.9%) | 0.0% | 227.0 (193.8–268.8) |
| J + 1,000/year, 75% to F | 67.2% (57.0%–72.1%) | 24.3% (23.3%–25.4%) | 0.0% | 222.4 (189.9–263.6) |
| J + 2,000/year, 25% to F | 62.4% (51.9%–67.6%) | 28.7% (27.4%–29.8%) | 0.0% | 249.4 (212.9–294.3) |
| J + 2,000/year, 50% to F | 71.7% (62.0%–76.3%) | 27.5% (26.3%–28.6%) | 0.0% | 240.6 (205.3–284.0) |
| J + 2,000/year, 75% to F | 79.3% (71.0%–83.4%) | 26.1% (25.0%–27.3%) | 0.0% | 230.2 (196.5–272.6) |
| J + 5,000/year, 25% to F | 75.7% (66.6%–80.0%) | 37.1% (35.0%–38.7%) | 0.0% | 296.0 (252.4–348.2) |
| J + 5,000/year, 50% to F | 90.0% (84.2%–92.9%) | 34.1% (32.3%–35.5%) | 0.0% | 272.9 (232.5–321.6) |
| J + 5,000/year, 75% to F | 96.7% (93.7%–98.2%) | 30.2% (28.8%–31.5%) | 0.0% | 244.3 (208.5–288.3) |
| J + 8,000/year, 25% to F | 85.4% (78.3%–88.9%) | 44.0% (41.4%–45.9%) | 0.0% | 332.7 (283.6–391.2) |
| J + 8,000/year, 50% to F | 97.4% (94.8%–98.6%) | 39.4% (37.3%–41.0%) | 0.0% | 297.5 (253.6–350.3) |
| J + 8,000/year, 75% to F | 99.7% (98.5%–99.9%) | 33.2% (31.7%–34.5%) | 0.0% | 255.8 (218.4–301.6) |
| J + uncapped screening | 99.8% (98.7%–99.9%) | 71.4% (71.3%–71.6%) | 100.0% | 484.1 (408.9–583.6) |

UI, uncertainty interval across posterior draws. Costs are 2026 Singapore dollars for J, F and X, discounted at 3% a year to 2026.

**Table 2. Median falls in new infections and HCV deaths in 2030 by annual probability of starting treatment in prison from 2026.**

| Policy | 0% in prison | 50% in prison | 90% in prison | 99.995% in prison |
|---|---|---|---|---|
| J only | -0.0% / -0.1% | 23.9% / 14.0% | 41.8% / 19.7% | 51.5% / 22.1% |
| J + 600/year, 25% to F | 2.7% / 2.6% | 27.7% / 16.4% | 45.5% / 21.9% | 54.9% / 24.2% |
| J + 600/year, 50% to F | 5.6% / 2.9% | 31.6% / 16.3% | 49.1% / 21.7% | 58.2% / 23.8% |
| J + 600/year, 75% to F | 8.7% / 3.1% | 35.4% / 16.2% | 52.6% / 21.4% | 61.4% / 23.5% |
| J + 1,000/year, 25% to F | 4.6% / 4.4% | 30.3% / 18.0% | 47.9% / 23.3% | 57.1% / 25.5% |
| J + 1,000/year, 50% to F | 9.7% / 4.8% | 36.6% / 17.8% | 53.8% / 22.9% | 62.4% / 24.9% |
| J + 1,000/year, 75% to F | 15.1% / 5.1% | 42.8% / 17.6% | 59.3% / 22.4% | 67.2% / 24.3% |
| J + 2,000/year, 25% to F | 9.7% / 8.7% | 36.6% / 21.7% | 53.8% / 26.7% | 62.4% / 28.7% |
| J + 2,000/year, 50% to F | 20.7% / 9.3% | 48.7% / 21.3% | 64.4% / 25.8% | 71.7% / 27.5% |
| J + 2,000/year, 75% to F | 31.9% / 9.4% | 59.6% / 20.6% | 73.4% / 24.6% | 79.3% / 26.1% |
| J + 5,000/year, 25% to F | 26.4% / 20.2% | 54.3% / 31.5% | 69.2% / 35.6% | 75.7% / 37.1% |
| J + 5,000/year, 50% to F | 52.4% / 20.4% | 76.8% / 29.8% | 86.6% / 33.0% | 90.0% / 34.1% |
| J + 5,000/year, 75% to F | 71.3% / 18.8% | 89.6% / 26.8% | 95.2% / 29.4% | 96.7% / 30.2% |
| J + 8,000/year, 25% to F | 42.7% / 29.6% | 69.0% / 39.5% | 80.8% / 42.9% | 85.4% / 44.0% |
| J + 8,000/year, 50% to F | 74.2% / 28.8% | 91.3% / 36.4% | 96.2% / 38.7% | 97.4% / 39.4% |
| J + 8,000/year, 75% to F | 88.0% / 24.5% | 97.7% / 30.8% | 99.4% / 32.6% | 99.7% / 33.2% |
| J + uncapped screening | 89.5% / 65.4% | 98.2% / 70.3% | 99.6% / 71.3% | 99.8% / 71.4% |

Each cell gives the fall in new infections / the fall in HCV deaths (medians across posterior draws). The last column is the main analysis (Table 1). Targets: 80% and 65%.

# Supporting information

**S1 Text. Technical appendix.** Model equations, parameters, prison rate estimation, calibration and its checks, the 2015 model population, scenario implementation, costs, additional results and the sensitivity analysis with prison treatment from 2019.

**S1 Checklist. CHEERS 2022 checklist.**
