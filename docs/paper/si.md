---
title: "S1 Text. Technical appendix"
subtitle: "Prison treatment and post-release screening for hepatitis C among people who inject drugs in Singapore: A modelling study"
---

Tables and figures in this appendix are numbered Table A, Table B and so on, and Fig A, Fig B and so on.

# 1. Model

## 1.1 Compartments

The state $y_{s,k,h,i}$ counts people in group $s$ (D, J, F or X), liver stage $k$ (1 non-cirrhotic, 2 compensated cirrhosis, 3 decompensated cirrhosis, 4 hepatocellular carcinoma), HCV state $h$ (u susceptible, a acute, c chronic, t on treatment, r recovered) and age group $i$ (1–6). There are 480 compartments.

Recovered people (r), after spontaneous clearance or treatment, remain antibody positive and can be reinfected. Susceptible people are in stage 1 only. Recovered people keep the liver stage they had when they cleared.

## 1.2 Transmission

The force of infection in age group $i$ is $\gamma_i = q \sum_j \tilde C_{ij} \, \dfrac{a_j + c_j}{N_j}$, where $a_j$ and $c_j$ are the acute and chronic people and $N_j$ is everyone in D and F in age group $j$, and $\tilde C_{ij} = C_{ij} e^{\theta_i}$ is the base contact matrix with row $i$ multiplied by the calibrated factor $e^{\theta_i}$. Susceptible and recovered people in D and F become acute at rate $\gamma_i$. People in J and X are not infected.

## 1.3 Movement between groups

New injectors enter D (state u, stage 1) at rate $s_i \beta_i$ per year, where $\beta_i$ is the 2015 Central Narcotics Bureau count of new drug abusers in age group $i$ and $s_i \ge 1$ is a calibrated entry scale. Entry is the same in every year.

People leave D for J at the first-arrest rate $\lambda^{(1)}_i = c\,\lambda^{(3)}_i$ and leave F for J at the re-arrest rate $\lambda^{(3)}_i$. People leave J at the release rate $\lambda^{(2)}_i$; a share $\pi$ goes to F and $1-\pi$ to X. People in X do not return to J. Every group ages from age group $i$ to $i+1$ at rate 1/10 per year.

## 1.4 Natural history and treatment

Acute infection lasts on average $\iota_1$; a share $\kappa$ clears to r and the rest become chronic. Chronic people progress from stage 1 to 2 and from stage 2 to 3 or 4, and from 3 to 4, at rates $p^{*} = (\rho\, r_3 + 1 - \rho)\,p$ that weight the genotype 3 relative risk $r_3$ by the genotype 3 share $\rho$. Recovered people in stages 2 and 3 progress at $\phi\, p^{*}$.

Chronic people start treatment at rate $\tau_{s,k}$ (zero in D, and zero before prison treatment starts: 2026 in the main analysis, 2019 in the sensitivity analysis). Treatment lasts on average $\iota_2$; a share $\alpha_k$ (sustained virological response) moves to r and the rest return to c.

## 1.5 Mortality

Everyone dies at the background rate $\mu_i\,\sigma_h$, where $\mu_i$ is the 2015 Singapore death rate in age group $i$ and $\sigma_h$ is the standardised mortality ratio of people who use drugs: $\sigma_v$ for viraemic states (a, c, t) and $\sigma_n$ for others (u, r). People in stages 3 and 4 also die of liver disease at $\mu_{DC}$ or $\mu_{HCC}$, multiplied by $\psi_{DC}$ or $\psi_{HCC}$ after treatment (states t and r).

## 1.6 Numerical solution and equilibrium

The equations were integrated with a fourth-order Runge–Kutta method and a step of one day. Negative values were set to zero after each step, and ageing was applied after each step. Each run started from a seed population in D and ran for 150 model years; model year 150 is 2015.

A run was accepted as being at equilibrium in 2015 if, between model years 140 and 150, the prison population of every age group changed by at most 1% on the log scale, prison seroprevalence by at most 0.5 percentage points, total HCV, decompensated cirrhosis and carcinoma counts by at most 1%, and every compartment by at most 2%. Runs that failed had zero posterior density.

# 2. Fixed parameters

Table A lists the fixed parameters. They are the values used by the model code.

**Table A. Fixed parameters.**

| Parameter | Value | Source |
|---|---|---|
| Transmission probability per sharing event, $q$ | 0.0057 | [@boelen2014] |
| Spontaneous clearance, $\kappa$ | 0.26 | [@micallef2006] |
| Mean acute duration, $\iota_1$ | 180 days | [@fasano2024] |
| Mean treatment duration, $\iota_2$ | 84 days (12 weeks) | Regimen length |
| Genotype 3 share, $\rho$ | 0.78 | [@zhao2019] |
| Sustained response, stages 1–2 | 0.9787 | [@pawlotsky2020] |
| Sustained response, stage 3, without ribavirin (scenarios) | 0.549 | [@xiao2025] |
| Sustained response, stage 4 | 0.8818 | Assumed equal to stage 3 with ribavirin |
| Progression, stage 1 to 2 (per year) | 0.027 | [@thein2008] |
| Progression, stage 2 to 3; 2 to 4; 3 to 4 (per year) | 0.039; 0.014; 0.014 | [@lim2018] |
| Genotype 3 relative risk, to cirrhosis; to carcinoma | 1.36; 1.93 | [@kanwal2014] |
| Relative progression after cure, 2 to 3; 2 to 4; 3 to 4 | 0.07; 0.23; 1.0 | [@vandermeer2012; @morgan2013]; last value assumed |
| Liver mortality, stage 3; stage 4 (per year) | 0.13; 0.43 | [@lim2018] |
| Relative liver mortality after cure, stage 3; stage 4 | 0.5; 0.37 | [@premkumar2024; @dang2020] |
| Mortality ratio, viraemic; non-viraemic | 9.52; 4.73 | [@vansanten2014] |
| Background death rate by age (per 1,000 per year) | 0.20; 0.25; 0.45; 1.20; 3.35; 21.88 | Singapore 2015 death rates (source to be cited) |
| Release rate by age (per year) | 0.462; 0.636; 0.454; 0.337; 0.325; 0.394 | Prison records (Section 3) |
| Re-arrest rate by age (per year) | 0.548; 0.478; 0.469; 0.471; 0.384; 0.265 | Prison records (Section 3) |
| Probability of eventual re-arrest, $\pi$ | 0.724 | Prison records (Section 3) |
| Entry into D by age, $\beta_i$ (per year) | 221.3; 554.5; 251.4; 80.9; 18.1; 0.9 | [@singstat_drug; @singstat_m891521]; 2015 new drug abusers × 0.8610, the share mainly using heroin or methamphetamine |

The base contact matrix $C$ is given in the model code; its source is to be documented by the authors.

# 3. Prison rate estimation

Singapore Prison Service records were restricted to people ever held for drug consumption. Consecutive stays with no gap were merged. This gave 20,720 people and 33,557 stays with admissions from 2000 to 2017; follow-up ended on 31 December 2017.

Release: each stay was split at the boundaries of the 10-year age groups, and a constant release hazard per age group was estimated, with a random-walk prior on the log hazards across ages.

Re-arrest: each interval from release to the next admission, or to the end of follow-up, was split by age in the same way. A cure-type survival model estimated the probability $\pi$ of ever being re-arrested and, for people who would be re-arrested, an age-specific hazard.

Both models were fitted in Stan with four chains of 4,000 sampling iterations; the largest R-hat was 1.001. Table B gives the posterior means used in the model and their 95% credible intervals.

**Table B. Prison rates by age (per year; posterior mean and 95% credible interval).**

| Age group | Release, $\lambda^{(2)}$ | Re-arrest, $\lambda^{(3)}$ |
|---|---|---|
| <20 | 0.462 (0.422–0.503) | 0.548 (0.459–0.647) |
| 20–29 | 0.636 (0.622–0.651) | 0.478 (0.455–0.501) |
| 30–39 | 0.454 (0.444–0.465) | 0.469 (0.449–0.490) |
| 40–49 | 0.337 (0.329–0.345) | 0.471 (0.451–0.491) |
| 50–59 | 0.325 (0.316–0.334) | 0.384 (0.363–0.405) |
| 60+ | 0.394 (0.372–0.416) | 0.265 (0.235–0.298) |

The probability of eventual re-arrest was 0.724 (95% credible interval 0.714–0.734). In the model, people who will be re-arrested are assumed to return to injecting (F).

These rates come from everyone held for drug consumption, not only people who inject.

# 4. Calibration

## 4.1 Targets

**Table C. Calibration targets, 2015.**

| Target | <20 | 20–29 | 30–39 | 40–49 | 50–59 | 60+ |
|---|---|---|---|---|---|---|
| HCV antibody prevalence in prison | 12.7% | 19.2% | 30.0% | 48.3% | 53.1% | 42.9% |
| People screened in prison (end of 2015) | 122 | 995 | 803 | 747 | 557 | 143 |
| New drug abusers (Central Narcotics Bureau) | 257 | 644 | 292 | 94 | 21 | 1 |
| First arrests (new drug abusers × 0.8610, rounded) | 221 | 554 | 251 | 81 | 18 | 1 |

Community injectors (D + F, all ages): 11,000 (95% confidence interval 6,500–17,000) [@teo2019].

The share 0.8610 = (113 + 1,014) / 1,309 is the share of the 2015 new drug abusers whose main drug was heroin or methamphetamine [@singstat_m891521]. The drug table has no age breakdown, so the same share was used at every age. It scales both the entry into D (unrounded) and the first-arrest target (rounded).

The seroprevalence targets are stored in the targets file. Recomputing them from the current processing of the prison records (positive divided by screened, Table C) gives values that differ by up to 0.5 percentage points; the reason is not recorded and should be resolved by the authors.

The community estimate refers to 2017; because entry is constant and nobody is treated before 2019 in any analysis, the modelled 2015 value is used for that year.

## 4.2 Parameters and prior

**Table D. Calibrated parameters and prior.**

| Parameter | Meaning | Prior on $\theta$ | Bounds of the scale |
|---|---|---|---|
| $\theta_1$–$\theta_6$ | Log contact multiplier, one per age group | Normal(0, 1.5) | 0.01–1,000 |
| $\theta_7$–$\theta_{12}$ | Log entry scale $s_i$, one per age group | Half-normal(0, 1) | 1–100 |
| $\theta_{13}$ | Log first-arrest factor $c$ | Normal(0, 0.5) | 0.5–100 |

The prior was fixed before fitting and did not use the data. The same prior was used for the mode search and the MCMC.

## 4.3 Likelihood

Seroprevalence: each observed value $y_i$ follows Beta($p_i\phi$, $(1-p_i)\phi$) around the modelled value $p_i$, with $\phi$ = 100. First arrests: each count follows a negative binomial distribution with mean equal to the modelled flow $c\,\lambda^{(3)}_i D_i$ and size 50. Community injectors: the log estimate follows Normal(log($D + F$), 0.2453).

## 4.4 Posterior mode

A Nelder–Mead search started from $\theta$ = 0 and 7 random perturbations of it (standard deviation 0.9 on the log scale), each with up to 5,000 iterations and up to 5 restarts. Seven of the 8 starts reached the same mode (negative log posterior 18.77 to 18.84); the other stopped at 73.1.

At the mode, the largest absolute seroprevalence error was 2.4 percentage points, the first-arrest factor was 2.11 and the community population was 6,786.

## 4.5 MCMC

We used adaptive random-walk Metropolis [@haario2001] with 8 chains started from independent prior draws. The initial proposal covariance was 0.01 times the prior variance. Adaptation began after 5,000 proposals, used the recent half of the warm-up, was refreshed every 50 proposals and was frozen after 20,000 warm-up proposals; 5% of proposals always used the initial covariance.

Each chain ran 60,000 model evaluations, and every 10th retained iteration was stored. The run took 1.3 hours. Acceptance rates were 0.12 to 0.13, below the 0.15 lower limit used to flag chains, so the chains mixed slowly but converged.

All diagnostics used every retained iteration: rank-normalised split R-hat was at most 1.004, bulk effective sample size at least 2,137, tail effective sample size at least 2,219, and the Monte Carlo standard error of each mean at most 0.022 posterior standard deviations [@vehtari2021] (Table E; Fig A).

**Table E. Posterior summaries and convergence diagnostics.**

| Parameter | Median (95% credible interval), natural scale | Posterior SD / prior SD | R-hat | Bulk ESS | Tail ESS |
|---|---|---|---|---|---|
| Contact, <20 | 5.48 (2.83–8.75) | 0.19 | 1.003 | 2,426 | 3,542 |
| Contact, 20–29 | 0.0621 (0.0357–0.0963) | 0.17 | 1.004 | 2,219 | 3,208 |
| Contact, 30–39 | 0.0402 (0.0186–0.0676) | 0.22 | 1.004 | 2,137 | 2,219 |
| Contact, 40–49 | 0.41 (0.16–0.74) | 0.27 | 1.002 | 2,342 | 2,327 |
| Contact, 50–59 | 0.85 (0.0597–4.93) | 0.77 | 1.002 | 2,847 | 4,454 |
| Contact, 60+ | 0.55 (0.0466–3.44) | 0.75 | 1.003 | 2,919 | 4,847 |
| Entry scale, <20 | 1.17 (1.01–1.58) | 0.20 | 1.002 | 3,005 | 3,510 |
| Entry scale, 20–29 | 1.16 (1.01–1.53) | 0.19 | 1.004 | 2,901 | 4,533 |
| Entry scale, 30–39 | 1.09 (1.00–1.41) | 0.15 | 1.002 | 3,099 | 3,730 |
| Entry scale, 40–49 | 1.09 (1.00–1.42) | 0.15 | 1.001 | 3,474 | 4,225 |
| Entry scale, 50–59 | 1.13 (1.00–1.61) | 0.21 | 1.002 | 3,192 | 3,692 |
| Entry scale, 60+ | 1.50 (1.02–3.97) | 0.62 | 1.002 | 2,733 | 3,568 |
| First-arrest factor $c$ | 2.59 (1.45–5.24) | 0.78 | 1.003 | 2,605 | 3,964 |

SD, standard deviation; ESS, effective sample size. The prior SD is that of the truncated prior.

For $\theta_{5}$, $\theta_{6}$, $\theta_{13}$ the posterior standard deviation was more than 0.7 times the prior's, so these parameters were mainly informed by the prior: the contact multipliers for ages 50–59 and 60 or older, and the first-arrest factor. The posterior median of c (2.59) was nevertheless well above the prior median of 1.

The entry scales were close to their lower bound of 1 for ages under 60 (medians 1.09 to 1.17); for ages 60 and over the median was 1.50 (95% credible interval 1.02–3.97).

## 4.6 Fit to the targets

For each of 4,000 posterior draws we computed the modelled 2015 targets. The 95% posterior intervals of the model values (without observation noise) overlapped the target ranges for all 13 targets (main text Fig 2).

**Table F. Targets against the 95% posterior interval of each modelled value.**

| Target | Age group | Target (95% range) | Model, 95% posterior interval |
|---|---|---|---|
| Seroprevalence | <20 | 12.7% (7.0–19.9) | 5.5–18.4 |
| Seroprevalence | 20–29 | 19.2% (12.1–27.4) | 13.0–27.5 |
| Seroprevalence | 30–39 | 30.0% (21.5–39.3) | 22.9–39.9 |
| Seroprevalence | 40–49 | 48.3% (38.6–58.0) | 40.9–56.6 |
| Seroprevalence | 50–59 | 53.1% (43.4–62.8) | 44.2–57.1 |
| Seroprevalence | 60+ | 42.9% (33.4–52.6) | 38.1–51.8 |
| First arrests | <20 | 221 (158–293) | 206–326 |
| First arrests | 20–29 | 554 (405–725) | 529–804 |
| First arrests | 30–39 | 251 (180–332) | 274–373 |
| First arrests | 40–49 | 81 (55–111) | 93–130 |
| First arrests | 50–59 | 18 (9–28) | 22–35 |
| First arrests | 60+ | 1 (0–3) | 2–6 |
| Community injectors | All ages (D + F) | 11,000 (6,801–17,791) | 6,167–8,383 |

## 4.7 Checks of the posterior

Direct simulation-based calibration [@talts2018]: we reran the MCMC (one chain of 45,000 evaluations, 15,000 warm-up) on 8 data sets simulated from prior draws. Pooled over 104 ranks, 95% interval coverage was 0.92, 90% coverage 0.89 and the mean rank 0.49. With 8 data sets this guards only against gross miscalibration.

# 5. Model population in 2015

**Table G. Modelled population in 2015 (median and 95% credible interval across posterior draws).**

| Group | All people | Viraemic (a, c, t) |
|---|---|---|
| D, never imprisoned | 1,040 (490–1,903) | 73 (21–199) |
| J, in prison | 8,493 (7,687–9,783) | 2,312 (1,956–2,769) |
| F, released, injecting | 5,934 (5,353–6,886) | 2,166 (1,829–2,630) |
| X, released, stopped | 34,964 (31,198–41,117) | 4,758 (3,864–5,921) |
| D + F (active injectors outside prison) | 7,005 (6,125–8,319) | 2,250 (1,907–2,708) |
| All groups | 50,494 (45,344–58,894) | 9,320 (7,746–11,352) |

The modelled prison population was not a calibration target. In the prison records, 3,367 people held for drug consumption with a screening result were in prison at the end of 2015.

# 6. Scenario implementation

Treatment in J started on 1 January 2026 in the main analysis (2019 in the sensitivity analysis) at the rate $\tau_J = -\ln(1 - 0.99995) = 9.90$ per year for people with chronic infection in every stage and age. In the coverage analysis the annual probability $p$ was 0.0, 0.5, 0.9, 0.99995, with $\tau_J = -\ln(1 - p)$. From 2026, people in F and X started at $\tau_{FX} = -\ln(1 - 0.99) = 4.61$ per year, but only until that year's capacity was used.

At the start of each year, for each of F and X, $N = u + a + c + r$ people were eligible. The number screened was the smaller of the quota and $N$. Antibody tests equalled the number screened; RNA tests equalled screened times $(a + c + r)/N$; chronic infections found, and so the year's treatment capacity, equalled screened times $c/N$.

The solver located the moment in each year at which a group used its capacity and then stopped new starts in that group for the rest of the year. Unused capacity was not moved between F and X or carried over. The same person can be counted as screened in more than one year, and the model counts treatment courses, not people.

In J, every entrant not on treatment (u, a, c, r) was counted as having an antibody test and every antibody-positive entrant (a, c, r) an RNA test; each treatment start was costed as one course and one consultation.

New infections in a year were the infections (u and r to a, in D and F) during that calendar year, and HCV deaths the liver deaths in stages 3 and 4 plus $(\sigma_v - \sigma_n)\mu_i$ times the number of viraemic people, summed over the year (trapezoid rule over the daily solver output). The 2015 values were one year at the calibrated equilibrium; in the main analysis the 2016–2025 values, still without treatment, differed from them by at most 0.04%. Both were compared with their 2015 values.

For each of 2,000 posterior draws (sampling seed fixed), we ran the 2015 equilibrium, the shared 2016–2025 history and all 17 scenarios to 2030; 2,000 of 2,000 draws were valid.

# 7. Costs

**Table H. Unit costs (2026 Singapore dollars).**

| Item | Base | Range in sensitivity analysis | Basis |
|---|---|---|---|
| SOF/VEL 400/100 mg | S$40,405.30 | S$20,202.65–S$60,607.95 | Published price of US$2,577 per week for 12 weeks [@wong2019], converted at S$1.3066 per US$ [@sgd_rate]; not a procurement price |
| Anti-HCV antibody | S$25.31 | S$12.65–S$37.97 | Published clinic charge [@nhg_dsc] |
| HCV RNA PCR | S$307.00 | S$153.50–S$460.50 | Published laboratory charge [@parkway_rna] |
| Doctor consultation | S$25.00 | S$12.50–S$37.50 | Published clinic charge [@nhg_dsc] |

The source abstract for the drug price gives US$32,325 per 12 weeks, which differs from US$2,577 times 12 (US$30,924); we used the weekly figure. Prices were not inflated to 2026. Follow-up visits, tests of cure, care navigation and genotyping were not costed.

Costs were counted in the year they occurred and discounted at 3% a year to 2026; numbers of tests and courses were not discounted. Costs before 2026 were not included.

# 8. Additional results

Table I gives every scenario. Fig C shows which policies met each target at the median, Fig D the costs and treatment starts by group, and Fig E the yearly screening cascade in F and X.

**Table I. All scenarios, 2030 and 2026–2030 totals (median and 95% uncertainty interval).**

| Policy | Draws meeting incidence target | Draws meeting mortality target | Treatment starts in F and X, 2026–2030 | Extra cost over J only, S$ million |
|---|---|---|---|---|
| J only | 0.0% | 0.0% | 0 (0–0) | 0.0 (0.0–0.0) |
| J + 600/year, 25% to F | 0.0% | 0.0% | 422 (355–495) | 11.6 (9.7–13.8) |
| J + 600/year, 50% to F | 0.0% | 0.0% | 477 (407–554) | 9.1 (7.6–10.8) |
| J + 600/year, 75% to F | 0.0% | 0.0% | 522 (452–604) | 6.4 (5.4–7.6) |
| J + 1,000/year, 25% to F | 0.0% | 0.0% | 692 (582–811) | 19.0 (15.8–22.6) |
| J + 1,000/year, 50% to F | 0.0% | 0.0% | 773 (659–899) | 14.7 (12.4–17.4) |
| J + 1,000/year, 75% to F | 0.0% | 0.0% | 828 (716–960) | 10.0 (8.5–11.8) |
| J + 2,000/year, 25% to F | 0.0% | 0.0% | 1,326 (1,117–1,556) | 36.4 (30.2–43.1) |
| J + 2,000/year, 50% to F | 0.0% | 0.0% | 1,441 (1,227–1,675) | 27.5 (23.1–32.5) |
| J + 2,000/year, 75% to F | 40.1% | 0.0% | 1,467 (1,268–1,702) | 17.0 (14.5–20.0) |
| J + 5,000/year, 25% to F | 5.9% | 0.0% | 2,922 (2,467–3,421) | 79.8 (66.3–94.5) |
| J + 5,000/year, 50% to F | 99.2% | 0.0% | 2,957 (2,520–3,443) | 57.0 (47.6–67.2) |
| J + 5,000/year, 75% to F | 100.0% | 0.0% | 2,645 (2,271–3,087) | 29.0 (24.5–34.2) |
| J + 8,000/year, 25% to F | 92.5% | 6.1% | 4,127 (3,492–4,837) | 112.4 (93.5–133.0) |
| J + 8,000/year, 50% to F | 100.0% | 0.0% | 3,968 (3,371–4,642) | 78.9 (65.7–93.1) |
| J + 8,000/year, 75% to F | 100.0% | 0.0% | 3,250 (2,766–3,821) | 39.6 (32.6–46.9) |
| J + uncapped screening | 100.0% | 100.0% | 7,256 (6,016–8,828) | 216.0 (177.0–264.8) |

Cost sensitivity: with all unit prices halved or increased by half, the discounted cost of 5,000 tests a year with 25% to F was S$142.0 million or S$425.9 million; epidemiological results do not change.

# 9. Sensitivity analysis: prison treatment from 2019

This analysis repeats the main analysis with treatment in J from 1 January 2019 at the same coverage; community screening of F and X still starts in 2026. It uses the same 2,000 posterior draws (Table J; Fig F).

**Table J. Sensitivity analysis with prison treatment from 2019: reductions in 2030 compared with 2015 and probability of meeting both targets.**

| Policy | Fall in new infections, % (95% UI) | Fall in HCV deaths, % (95% UI) | Draws meeting both targets |
|---|---|---|---|
| J only | 83.4% (67.8%–89.5%) | 53.2% (50.7%–55.9%) | 0.0% |
| J + 600/year, 25% to F | 84.8% (70.0%–90.5%) | 55.1% (52.6%–57.7%) | 0.0% |
| J + 600/year, 50% to F | 86.1% (72.1%–91.3%) | 54.7% (52.1%–57.3%) | 0.0% |
| J + 600/year, 75% to F | 87.3% (74.2%–92.2%) | 54.2% (51.6%–56.9%) | 0.0% |
| J + 1,000/year, 25% to F | 85.7% (71.4%–91.1%) | 56.4% (53.9%–58.9%) | 0.0% |
| J + 1,000/year, 50% to F | 87.7% (74.9%–92.4%) | 55.6% (53.1%–58.3%) | 0.0% |
| J + 1,000/year, 75% to F | 89.6% (78.2%–93.6%) | 54.8% (52.2%–57.5%) | 0.0% |
| J + 2,000/year, 25% to F | 87.7% (74.9%–92.4%) | 59.3% (56.8%–61.9%) | 0.0% |
| J + 2,000/year, 50% to F | 91.2% (80.8%–94.7%) | 57.8% (55.2%–60.5%) | 0.0% |
| J + 2,000/year, 75% to F | 93.8% (85.4%–96.4%) | 56.1% (53.5%–58.9%) | 0.0% |
| J + 5,000/year, 25% to F | 92.6% (83.2%–95.6%) | 66.9% (64.1%–69.3%) | 90.0% |
| J + 5,000/year, 50% to F | 97.2% (92.4%–98.5%) | 63.6% (60.8%–66.3%) | 14.7% |
| J + 5,000/year, 75% to F | 99.1% (96.9%–99.6%) | 59.5% (56.7%–62.3%) | 0.0% |
| J + 8,000/year, 25% to F | 95.7% (89.3%–97.6%) | 72.7% (69.8%–75.2%) | 99.7% |
| J + 8,000/year, 50% to F | 99.3% (97.4%–99.7%) | 68.2% (65.3%–70.8%) | 98.2% |
| J + 8,000/year, 75% to F | 99.9% (99.2%–100.0%) | 62.2% (59.4%–65.0%) | 2.5% |
| J + uncapped screening | 99.9% (99.3%–100.0%) | 88.3% (87.5%–89.3%) | 100.0% |

With prison treatment from 2019, new infections under J only were 36% of the 2015 level by 2025, against 100% in the main analysis.

# 10. Software and reproducibility

The model is in C++ (scripts/model/sim.cpp), called from R through Rcpp [@rcpp2011]. Rate models used Stan and convergence diagnostics the R package posterior. Figures were drawn with matplotlib from the saved outputs.

Analysis commands, run times and expected outputs are listed in docs/training.md in the repository. This appendix and the main text were generated from the saved outputs by scripts/build_paper.py; drafting followed the procedures of Scientific Agent Skills [@kassis2026].

# Figures

**Fig A. MCMC trace plots.** Stored draws of the 13 parameters for each of the 8 chains after warm-up.

![](output/calibration/mcmc_reference/mcmc_traceplots.png){width=95%}

**Fig B. Marginal posterior densities of the 13 parameters from the MCMC.**

![](output/calibration/mcmc_reference/figures/marginals.png){width=95%}

**Fig C. Target attainment in 2030 by policy.** Median reductions in new infections and HCV deaths, and whether both targets were met.

![](output/scenarios_2026/scenario_target_heatmap.png){width=95%}

**Fig D. Costs and treatment starts, 2026–2030.** (A) Discounted costs by group. (B) Undiscounted treatment courses by group. Bars are medians of each group; whiskers are 95% uncertainty intervals of the total.

![](output/scenarios_2026/scenario_cumulative_costs.png){width=95%}

**Fig E. Screening cascade in F and X by year.** Number screened, RNA tests, RNA-positive people, chronic infections found and treatment starts, by capacity (columns) and share of tests given to F (rows).

![](output/scenarios_2026/scenario_care_cascade.png){width=95%}

**Fig F. Sensitivity analysis with prison treatment from 2019.** (A) New infections and HCV deaths relative to 2015 (medians), 2015–2030; the dotted line marks the start of prison treatment (2019) and the dashed line the start of community screening (2026). (B) Reductions in 2030 by policy, with 95% uncertainty intervals.

![](output/scenarios_2026_j2019/scenario_target_trajectories.png){width=95%}

# References

::: {#refs}
:::
