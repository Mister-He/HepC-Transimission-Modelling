---
title: "S1 Text. Technical appendix"
subtitle: "Prison treatment and post-release screening for hepatitis C among people who inject drugs in Singapore: A modelling study"
---

Tables and figures in this appendix are numbered Table A, Table B and so on, and Fig A, Fig B and so on.

# 1. Model

## 1.1 Compartments

The state $y_{s,k,h,i}$ counts people in group $s$ (D, J, F or X), liver stage $k$ (1 non-cirrhotic, 2 compensated cirrhosis, 3 decompensated cirrhosis, 4 hepatocellular carcinoma), HCV state $h$ (u susceptible, a acute, c chronic, t on treatment, r recovered) and age group $i$ (1–6). There are 480 compartments. [claim:C129] [evidence:E032]

Recovered people (r), after spontaneous clearance or treatment, remain antibody positive and can be reinfected. Susceptible people are in stage 1 only. Recovered people keep the liver stage they had when they cleared. [claim:C130] [evidence:E032]

## 1.2 Transmission

The force of infection in age group $i$ is $\gamma_i = q \sum_j \tilde C_{ij} \, \dfrac{a_j + c_j}{N_j}$, where $a_j$ and $c_j$ are the acute and chronic people and $N_j$ is everyone in D and F in age group $j$, and $\tilde C_{ij} = C_{ij} e^{\theta_i}$ is the base contact matrix with row $i$ multiplied by the calibrated factor $e^{\theta_i}$. Susceptible and recovered people in D and F become acute at rate $\gamma_i$. People in J and X are not infected. [claim:C131] [evidence:E032,E041]

## 1.3 Movement between groups

New injectors enter D (state u, stage 1) at rate $s_i \beta_i$ per year, where $\beta_i$ is the 2015 Central Narcotics Bureau count of new drug abusers in age group $i$ and $s_i \ge 1$ is a calibrated entry scale. Entry is the same in every year. [claim:C132] [evidence:E032,E004,E041]

People leave D for J at the first-arrest rate $\lambda^{(1)}_i = c\,\lambda^{(3)}_i$ and leave F for J at the re-arrest rate $\lambda^{(3)}_i$. People leave J at the release rate $\lambda^{(2)}_i$; a share $\pi$ goes to F and $1-\pi$ to X. People in X do not return to J. Every group ages from age group $i$ to $i+1$ at rate 1/10 per year. [claim:C133] [evidence:E032,E004]

## 1.4 Natural history and treatment

Acute infection lasts on average $\iota_1$; a share $\kappa$ clears to r and the rest become chronic. Chronic people progress from stage 1 to 2 and from stage 2 to 3 or 4, and from 3 to 4, at rates $p^{*} = (\rho\, r_3 + 1 - \rho)\,p$ that weight the genotype 3 relative risk $r_3$ by the genotype 3 share $\rho$. Recovered people in stages 2 and 3 progress at $\phi\, p^{*}$. [claim:C134] [evidence:E032,E004]

Chronic people start treatment at rate $\tau_{s,k}$ (zero in D, and zero before prison treatment starts: 2026 in the main analysis, 2019 in the sensitivity analysis). Treatment lasts on average $\iota_2$; a share $\alpha_k$ (sustained virological response) moves to r and the rest return to c. [claim:C135] [evidence:E032]

## 1.5 Mortality

Each cause of death has one route. Everyone dies of causes other than liver disease at $\mu_i\,\sigma_h$, where $\mu_i$ is the 2015 Singapore death rate in age group $i$ and $\sigma_n$ is the standardised mortality ratio of people who use drugs without HCV or HIV infection, which contains almost no liver deaths; viraemic people (a, c, t) have $\sigma_v = \sigma_n$ times the hazard ratio of extrahepatic death for viraemic infection. Liver deaths occur only in stages 3 and 4, at $\mu_{DC}$ or $\mu_{HCC}$, multiplied by $\psi_{DC}$ or $\psi_{HCC}$ after cure (state r only; people on treatment keep the untreated rate). $\psi$ are adjusted hazard ratios of liver-related death for sustained virological response against none. [claim:C136] [evidence:E032,E004]

## 1.6 Numerical solution and equilibrium

The equations were integrated with a fourth-order Runge–Kutta method and a step of one day. Negative values were set to zero after each step, and ageing was applied after each step. Each run started from a seed population in D and ran for 150 model years; model year 150 is 2015. [claim:C137] [evidence:E032,E004,E007]

A run was accepted as being at equilibrium in 2015 if, between model years 140 and 150, the prison population of every age group changed by at most 1% on the log scale, prison seroprevalence by at most 0.5 percentage points, total HCV, decompensated cirrhosis and carcinoma counts by at most 1%, and every compartment by at most 2%. Runs that failed had zero posterior density. [claim:C138] [evidence:E037,E041]

# 2. Fixed parameters

Table A lists the fixed parameters. They are the values used by the model code. [claim:C139] [evidence:E004]

**Table A. Fixed parameters.**

| Parameter | Value | Source |
|---|---|---|
| Transmission probability per sharing event, $q$ | 0.0057 | [@boelen2014] | [claim:C140] [evidence:E004,E069]
| Spontaneous clearance, $\kappa$ | 0.26 | [@micallef2006] | [claim:C141] [evidence:E004,E070]
| Mean acute duration, $\iota_1$ | 180 days | [@fasano2024] | [claim:C142] [evidence:E004,E071]
| Mean treatment duration, $\iota_2$ | 84 days (12 weeks) | Regimen length | [claim:C143] [evidence:E004]
| Genotype 3 share, $\rho$ | 0.78 | [@zhao2019] | [claim:C144] [evidence:E004,E072]
| Sustained response, stages 1–2 | 0.9787 | [@pawlotsky2020] | [claim:C145] [evidence:E004,E019]
| Sustained response, stage 3, without ribavirin (scenarios) | 0.549 | [@xiao2025] | [claim:C146] [evidence:E004,E073]
| Sustained response, stage 4 | 0.8818 | Assumed equal to stage 3 with ribavirin | [claim:C147] [evidence:E004]
| Progression, stage 1 to 2 (per year) | 0.027 | [@thein2008] | [claim:C148] [evidence:E004,E074]
| Progression, stage 2 to 3; 2 to 4; 3 to 4 (per year) | 0.039; 0.014; 0.014 | [@lim2018] | [claim:C149] [evidence:E004,E075]
| Genotype 3 relative risk, to cirrhosis; to carcinoma | 1.36; 1.93 | [@kanwal2014]; Derived from the genotype-specific hazard ratios: genotype 3 against a prevalence-weighted geometric mean of genotypes 1, 2 and 4 | [claim:C150] [evidence:E004,E076]
| Relative progression after cure, 2 to 3; 2 to 4; 3 to 4 | 0.07; 0.23; 1.0 | [@vandermeer2012; @morgan2013]; last value assumed | [claim:C151] [evidence:E004,E077,E078]
| Liver mortality, stage 3; stage 4 (per year) | 0.13; 0.43 | [@lim2018] | [claim:C152] [evidence:E004,E075]
| Relative liver mortality after cure, stage 3; stage 4 | 0.4; 0.34 | [@pageaux2022; @dang2020]; liver-related death | [claim:C153] [evidence:E004,E079,E080]
| Mortality ratio of people who use drugs, $\sigma_n$ (HCV and HIV uninfected) | 4.73 | [@vansanten2014]; All-cause; applied to everyone | [claim:C154] [evidence:E004,E035]
| Extrahepatic mortality hazard ratio, viraemic | 1.47 | [@lee2012]; HCV RNA positive against antibody negative, adjusted; $\sigma_v = \sigma_n \times$ this ratio | [claim:C155] [evidence:E004,E036]
| Background death rate by age (per 1,000 per year) | 0.20; 0.25; 0.45; 1.20; 3.35; 10.14 | [@singstat_m810141; @singstat_m810011]; Singapore residents, 2015; 60+ is the population-weighted rate for ages 60–74 | [claim:C156] [evidence:E004,E034,E081]
| Release rate by age (per year) | 0.462; 0.636; 0.454; 0.337; 0.325; 0.394 | Prison records (Section 3) | [claim:C157] [evidence:E004]
| Re-arrest rate by age (per year) | 0.548; 0.478; 0.469; 0.471; 0.384; 0.265 | Prison records (Section 3) | [claim:C158] [evidence:E004]
| Probability of eventual re-arrest, $\pi$ | 0.724 | Prison records (Section 3) | [claim:C159] [evidence:E004]
| Entry into D by age, $\beta_i$ (per year) | 221.3; 554.5; 251.4; 80.9; 18.1; 0.9 | [@singstat_drug; @singstat_m891521]; 2015 new drug abusers × 0.8610, the share mainly using heroin or methamphetamine | [claim:C160] [evidence:E004,E039,E040]

The base contact matrix $C$ is given in the model code; its source is to be documented by the authors.

# 3. Prison rate estimation

Singapore Prison Service records were restricted to people ever held for drug consumption. Consecutive stays with no gap were merged. This gave 20,720 people and 33,557 stays with admissions from 2000 to 2017; follow-up ended on 31 December 2017. [claim:C161] [evidence:E011,E082]

Release: each stay was split at the boundaries of the 10-year age groups, and a constant release hazard per age group was estimated, with a random-walk prior on the log hazards across ages. [claim:C162] [evidence:E083,E082]

Re-arrest: each interval from release to the next admission, or to the end of follow-up, was split by age in the same way. A cure-type survival model estimated the probability $\pi$ of ever being re-arrested and, for people who would be re-arrested, an age-specific hazard. [claim:C163] [evidence:E084,E082]

Both models were fitted in Stan with four chains of 4,000 sampling iterations; the largest R-hat was 1.001. Table B gives the posterior means used in the model and their 95% credible intervals. [claim:C164] [evidence:E010,E038,E082]

**Table B. Prison rates by age (per year; posterior mean and 95% credible interval).** [claim:C165] [evidence:E010]

| Age group | Release, $\lambda^{(2)}$ | Re-arrest, $\lambda^{(3)}$ | [claim:C166] [evidence:E010]
|---|---|---|
| <20 | 0.462 (0.422–0.503) | 0.548 (0.459–0.647) | [claim:C167] [evidence:E010,E038]
| 20–29 | 0.636 (0.622–0.651) | 0.478 (0.455–0.501) | [claim:C168] [evidence:E010,E038]
| 30–39 | 0.454 (0.444–0.465) | 0.469 (0.449–0.490) | [claim:C169] [evidence:E010,E038]
| 40–49 | 0.337 (0.329–0.345) | 0.471 (0.451–0.491) | [claim:C170] [evidence:E010,E038]
| 50–59 | 0.325 (0.316–0.334) | 0.384 (0.363–0.405) | [claim:C171] [evidence:E010,E038]
| 60+ | 0.394 (0.372–0.416) | 0.265 (0.235–0.298) | [claim:C172] [evidence:E010,E038]

The probability of eventual re-arrest was 0.724 (95% credible interval 0.714–0.734). In the model, people who will be re-arrested are assumed to return to injecting (F). [claim:C173] [evidence:E010,E004]

These rates come from everyone held for drug consumption, not only people who inject. [claim:C174] [evidence:E082]

# 4. Calibration

## 4.1 Targets

**Table C. Calibration targets, 2015.** [claim:C175] [evidence:E007]

| Target | <20 | 20–29 | 30–39 | 40–49 | 50–59 | 60+ | [claim:C176] [evidence:E007]
|---|---|---|---|---|---|---|
| HCV antibody prevalence in prison | 12.7% | 19.2% | 30.0% | 48.3% | 53.1% | 42.9% | [claim:C177] [evidence:E007]
| People screened in prison, end of 2015 (seroprevalence denominator) | 122 | 995 | 803 | 747 | 557 | 143 | [claim:C178] [evidence:E042]
| New drug abusers (Central Narcotics Bureau) | 257 | 644 | 292 | 94 | 21 | 1 | [claim:C179] [evidence:E005]
| First arrests (new drug abusers × 0.8610, rounded) | 221 | 554 | 251 | 81 | 18 | 1 | [claim:C180] [evidence:E005,E006,E007]

Community injectors (D + F, all ages): 11,000 (95% confidence interval 6,500–17,000) [@teo2019]. [claim:C181] [evidence:E007,E023]

The share 0.8610 = (113 + 1,014) / 1,309 is the share of the 2015 new drug abusers whose main drug was heroin or methamphetamine [@singstat_m891521]. The drug table has no age breakdown, so the same share was used at every age. It scales both the entry into D (unrounded) and the first-arrest target (rounded). [claim:C182] [evidence:E006,E004,E007,E040]

The seroprevalence targets are stored in the targets file. Recomputing them from the current processing of the prison records (positive divided by screened, Table C) gives values that differ by up to 0.5 percentage points; the reason is not recorded and should be resolved by the authors. [claim:C183] [evidence:E007,E042]

The community estimate refers to 2017; because entry is constant and nobody is treated before 2019 in any analysis, the modelled 2015 value is used for that year. [claim:C184] [evidence:E007,E032]

## 4.2 Parameters and prior

**Table D. Calibrated parameters and prior.**

| Parameter | Meaning | Prior on $\theta$ | Bounds of the scale | [claim:C185] [evidence:E041]
|---|---|---|---|
| $\theta_1$–$\theta_6$ | Log contact multiplier, one per age group | Normal(0, 1.5) | 0.01–1,000 | [claim:C186] [evidence:E041]
| $\theta_7$–$\theta_{12}$ | Log entry scale $s_i$, one per age group | Half-normal(0, 1) | 1–100 | [claim:C187] [evidence:E041]
| $\theta_{13}$ | Log first-arrest factor $c$ | Normal(0, 0.5) | 0.5–100 | [claim:C188] [evidence:E041]

The prior was fixed before fitting and did not use the data. The same prior was used for the mode search and the MCMC. [claim:C189] [evidence:E041]

## 4.3 Likelihood

Seroprevalence: each observed value $y_i$ follows Beta($p_i\phi$, $(1-p_i)\phi$) around the modelled value $p_i$, with $\phi$ = 100. First arrests: each count follows a negative binomial distribution with mean equal to the modelled flow $c\,\lambda^{(3)}_i D_i$ and size 50. Community injectors: the log estimate follows Normal(log($D + F$), 0.2453). [claim:C190] [evidence:E041,E007,E032]

## 4.4 Posterior mode

A Nelder–Mead search started from $\theta$ = 0 and 7 random perturbations of it (standard deviation 0.9 on the log scale), each with up to 5,000 iterations and up to 5 restarts. All 8 starts reached negative log posteriors within 1 of each other (19.25 to 20.13). [claim:C191] [evidence:E085,E043]

At the mode, the largest absolute seroprevalence error was 4.0 percentage points, the first-arrest factor was 2.20 and the community population was 7,096. [claim:C192] [evidence:E085]

## 4.5 MCMC

We used adaptive random-walk Metropolis [@haario2001] with 8 chains started from independent prior draws. The initial proposal covariance was 0.01 times the prior variance. Adaptation began after 5,000 proposals, used the recent half of the warm-up, was refreshed every 50 proposals and was frozen after 20,000 warm-up proposals; 5% of proposals always used the initial covariance. [claim:C193] [evidence:E008,E044]

Each chain ran 60,000 model evaluations, and every 10th retained iteration was stored. The run took 2.0 hours. Acceptance rates were 0.12 to 0.13, below the 0.15 lower limit used to flag chains, so the chains mixed slowly but converged. [claim:C194] [evidence:E008,E086]

All diagnostics used every retained iteration: rank-normalised split R-hat was at most 1.005, bulk effective sample size at least 1,641, tail effective sample size at least 1,490, and the Monte Carlo standard error of each mean at most 0.026 posterior standard deviations [@vehtari2021] (Table E; Fig A). [claim:C195] [evidence:E058,E045]

**Table E. Posterior summaries and convergence diagnostics.**

| Parameter | Median (95% credible interval), natural scale | Posterior SD / prior SD | R-hat | Bulk ESS | Tail ESS | [claim:C196] [evidence:E014]
|---|---|---|---|---|---|
| Contact, <20 | 5.62 (3.05–8.76) | 0.18 | 1.002 | 3,177 | 4,258 | [claim:C197] [evidence:E014,E058,E008]
| Contact, 20–29 | 0.062 (0.0353–0.0956) | 0.17 | 1.005 | 1,945 | 2,719 | [claim:C198] [evidence:E014,E058,E008]
| Contact, 30–39 | 0.0407 (0.0188–0.0672) | 0.21 | 1.005 | 1,641 | 1,490 | [claim:C199] [evidence:E014,E058,E008]
| Contact, 40–49 | 0.37 (0.16–0.65) | 0.24 | 1.005 | 2,106 | 2,232 | [claim:C200] [evidence:E014,E058,E008]
| Contact, 50–59 | 0.57 (0.0433–3.12) | 0.74 | 1.003 | 2,818 | 3,807 | [claim:C201] [evidence:E014,E058,E008]
| Contact, 60+ | 0.28 (0.0288–1.29) | 0.66 | 1.002 | 2,984 | 4,476 | [claim:C202] [evidence:E014,E058,E008]
| Entry scale, <20 | 1.15 (1.01–1.53) | 0.19 | 1.002 | 2,968 | 4,020 | [claim:C203] [evidence:E014,E058,E008]
| Entry scale, 20–29 | 1.15 (1.01–1.49) | 0.18 | 1.002 | 2,890 | 4,224 | [claim:C204] [evidence:E014,E058,E008]
| Entry scale, 30–39 | 1.09 (1.00–1.40) | 0.15 | 1.002 | 2,689 | 3,176 | [claim:C205] [evidence:E014,E058,E008]
| Entry scale, 40–49 | 1.08 (1.00–1.39) | 0.15 | 1.004 | 2,898 | 3,426 | [claim:C206] [evidence:E014,E058,E008]
| Entry scale, 50–59 | 1.13 (1.00–1.61) | 0.21 | 1.004 | 3,051 | 3,543 | [claim:C207] [evidence:E014,E058,E008]
| Entry scale, 60+ | 1.49 (1.02–4.05) | 0.63 | 1.004 | 2,314 | 3,074 | [claim:C208] [evidence:E014,E058,E008]
| First-arrest factor $c$ | 2.71 (1.54–5.23) | 0.72 | 1.003 | 2,822 | 4,344 | [claim:C209] [evidence:E014,E058,E008]

SD, standard deviation; ESS, effective sample size. The prior SD is that of the truncated prior. [claim:C210] [evidence:E008]

For $\theta_{5}$, $\theta_{13}$ the posterior standard deviation was more than 0.7 times the prior's, so these parameters were mainly informed by the prior (contact, 50–59; first-arrest factor $c$). The posterior median of c (2.71) was nevertheless well above the prior median of 1. [claim:C211] [evidence:E014,E008]

The entry scales were close to their lower bound of 1 for ages under 60 (medians 1.08 to 1.15); for ages 60 and over the median was 1.49 (95% credible interval 1.02–4.05). [claim:C212] [evidence:E014]

## 4.6 Fit to the targets

For each of 4,000 posterior draws we computed the modelled 2015 targets. The 95% posterior intervals of the model values (without observation noise) overlapped the target ranges for all 13 targets (main text Fig 2). [claim:C213] [evidence:E061,E059]

**Table F. Targets against the 95% posterior interval of each modelled value.** [claim:C214] [evidence:E059]

| Target | Age group | Target (95% range) | Model, 95% posterior interval | [claim:C215] [evidence:E059]
|---|---|---|---|
| Seroprevalence | <20 | 12.7% (7.0–19.9) | 6.0–18.7 | [claim:C216] [evidence:E059]
| Seroprevalence | 20–29 | 19.2% (12.1–27.4) | 13.0–27.6 | [claim:C217] [evidence:E059]
| Seroprevalence | 30–39 | 30.0% (21.5–39.3) | 22.8–39.9 | [claim:C218] [evidence:E059]
| Seroprevalence | 40–49 | 48.3% (38.6–58.0) | 42.0–54.8 | [claim:C219] [evidence:E059]
| Seroprevalence | 50–59 | 53.1% (43.4–62.8) | 43.4–55.4 | [claim:C220] [evidence:E059]
| Seroprevalence | 60+ | 42.9% (33.4–52.6) | 40.3–54.7 | [claim:C221] [evidence:E059]
| First arrests | <20 | 221 (158–293) | 206–316 | [claim:C222] [evidence:E059]
| First arrests | 20–29 | 554 (405–725) | 531–785 | [claim:C223] [evidence:E059]
| First arrests | 30–39 | 251 (180–332) | 272–373 | [claim:C224] [evidence:E059]
| First arrests | 40–49 | 81 (55–111) | 92–127 | [claim:C225] [evidence:E059]
| First arrests | 50–59 | 18 (9–28) | 22–35 | [claim:C226] [evidence:E059]
| First arrests | 60+ | 1 (0–3) | 2–7 | [claim:C227] [evidence:E059]
| Community injectors | All ages (D + F) | 11,000 (6,801–17,791) | 6,467–8,634 | [claim:C228] [evidence:E059]

## 4.7 Checks of the posterior

Direct simulation-based calibration [@talts2018]: we reran the MCMC (one chain of 45,000 evaluations, 15,000 warm-up) on 8 data sets simulated from prior draws. Pooled over 104 ranks, 95% interval coverage was 0.95, 90% coverage 0.89 and the mean rank 0.51. With 8 data sets this guards only against gross miscalibration. [claim:C229] [evidence:E087,E046]

# 5. Model population in 2015

**Table G. Modelled population in 2015 (median and 95% credible interval across posterior draws).** [claim:C230] [evidence:E061]

| Group | All people | Viraemic (a, c, t) |
|---|---|---|
| D, never imprisoned | 1,000 (511–1,767) | 67 (22–170) | [claim:C231] [evidence:E061]
| J, in prison | 8,768 (8,000–10,110) | 2,414 (2,067–2,860) | [claim:C232] [evidence:E061]
| F, released, injecting | 6,303 (5,736–7,296) | 2,300 (1,959–2,755) | [claim:C233] [evidence:E061]
| X, released, stopped | 43,510 (39,167–50,744) | 5,892 (4,856–7,216) | [claim:C234] [evidence:E061]
| D + F (active injectors outside prison) | 7,332 (6,513–8,610) | 2,373 (2,026–2,832) | [claim:C235] [evidence:E061]
| All groups | 59,599 (53,968–69,269) | 10,676 (8,953–12,867) | [claim:C236] [evidence:E061]

The modelled prison population was not a calibration target. In the prison records, 8,725 people held for drug consumption were in prison at the end of 2015, the population whose release and re-arrest rates the model uses (Section 3); 3,367 of them had a screening result and gave the seroprevalence targets. [claim:C237] [evidence:E061,E042]

# 6. Scenario implementation

Treatment in J started on 1 January 2026 in the main analysis (2019 in the sensitivity analysis) at the rate $\tau_J = -\ln(1 - 0.99995) = 9.90$ per year for people with chronic infection in every stage and age. In the coverage analysis the annual probability $p$ was 0.0, 0.5, 0.9, 0.99995, with $\tau_J = -\ln(1 - p)$. From 2026, people in F and X started at $\tau_{FX} = -\ln(1 - 0.99) = 4.61$ per year, but only until that year's capacity was used. [claim:C238] [evidence:E001,E002,E012,E013]

At the start of each year, for each of F and X, $N = u + a + c + r$ people were eligible. The number screened was the smaller of the quota and $N$. Antibody tests equalled the number screened; RNA tests equalled screened times $(a + c + r)/N$; chronic infections found, and so the year's treatment capacity, equalled screened times $c/N$. [claim:C239] [evidence:E002]

The solver located the moment in each year at which a group used its capacity and then stopped new starts in that group for the rest of the year. Unused capacity was not moved between F and X or carried over. The same person can be counted as screened in more than one year, and the model counts treatment courses, not people. [claim:C240] [evidence:E032,E002]

In J, every entrant not on treatment (u, a, c, r) was counted as having an antibody test and every antibody-positive entrant (a, c, r) an RNA test; each treatment start was costed as one course and one consultation. [claim:C241] [evidence:E002]

New infections in a year were the infections (u and r to a, in D and F) during that calendar year, and HCV deaths the liver deaths in stages 3 and 4 (deaths from HCV-related decompensated cirrhosis and liver cancer, the deaths counted by the WHO target), summed over the year (trapezoid rule over the daily solver output). The 2015 values were one year at the calibrated equilibrium; in the main analysis the 2016–2025 values, still without treatment, differed from them by at most 0.09%. Both were compared with their 2015 values. [claim:C242] [evidence:E002,E032,E062]

For each of 2,000 posterior draws (sampling seed fixed), we ran the 2015 equilibrium, the shared 2016–2025 history and all 17 scenarios to 2030; 2,000 of 2,000 draws were valid. [claim:C243] [evidence:E009,E052,E002]

# 7. Costs

**Table H. Unit costs (2026 Singapore dollars).** [claim:C244] [evidence:E047]

| Item | Base | Range in sensitivity analysis | Basis |
|---|---|---|---|
| SOF/VEL 400/100 mg | S$40,405.30 | S$20,202.65–S$60,607.95 | Published price of US$2,577 per week for 12 weeks [@wong2019], converted at S$1.3066 per US$ [@sgd_rate]; not a procurement price | [claim:C245] [evidence:E047,E048,E051]
| Anti-HCV antibody | S$25.31 | S$12.65–S$37.97 | Published clinic charge [@nhg_dsc] | [claim:C246] [evidence:E047,E049]
| HCV RNA PCR | S$307.00 | S$153.50–S$460.50 | Published laboratory charge [@parkway_rna] | [claim:C247] [evidence:E047,E050]
| Doctor consultation | S$25.00 | S$12.50–S$37.50 | Published clinic charge [@nhg_dsc] | [claim:C248] [evidence:E047,E049]

The source abstract for the drug price gives US$32,325 per 12 weeks, which differs from US$2,577 times 12 (US$30,924); we used the weekly figure. Prices were not inflated to 2026. Follow-up visits, tests of cure, care navigation and genotyping were not costed. [claim:C249] [evidence:E047]

Costs were counted in the year they occurred and discounted at 3% a year to 2026; numbers of tests and courses were not discounted. Costs before 2026 were not included. [claim:C250] [evidence:E001]

# 8. Additional results

Table I gives every scenario. Fig C shows which policies met each target at the median, Fig D the costs and treatment starts by group, and Fig E the yearly screening cascade in F and X. [claim:C251] [evidence:E009]

**Table I. All scenarios, 2030 and 2026–2030 totals (median and 95% uncertainty interval).** [claim:C252] [evidence:E009]

| Policy | Draws meeting incidence target | Draws meeting mortality target | Treatment starts in F and X, 2026–2030 | Extra cost over J only, S$ million | [claim:C253] [evidence:E009]
|---|---|---|---|---|
| J only | 0.0% | 0.0% | 0 (0–0) | 0.0 (0.0–0.0) | [claim:C254] [evidence:E009]
| J + 600/year, 25% to F | 0.0% | 0.0% | 429 (363–498) | 11.9 (10.0–14.0) | [claim:C255] [evidence:E009]
| J + 600/year, 50% to F | 0.0% | 0.0% | 485 (416–556) | 9.5 (8.0–11.1) | [claim:C256] [evidence:E009]
| J + 600/year, 75% to F | 0.0% | 0.0% | 533 (460–606) | 6.9 (5.9–8.0) | [claim:C257] [evidence:E009]
| J + 1,000/year, 25% to F | 0.0% | 0.0% | 705 (596–818) | 19.6 (16.4–23.0) | [claim:C258] [evidence:E009]
| J + 1,000/year, 50% to F | 0.0% | 0.0% | 788 (674–903) | 15.4 (13.0–17.9) | [claim:C259] [evidence:E009]
| J + 1,000/year, 75% to F | 0.0% | 0.0% | 847 (732–962) | 10.8 (9.2–12.4) | [claim:C260] [evidence:E009]
| J + 2,000/year, 25% to F | 0.0% | 0.0% | 1,359 (1,149–1,575) | 37.7 (31.6–44.3) | [claim:C261] [evidence:E009]
| J + 2,000/year, 50% to F | 0.0% | 0.0% | 1,474 (1,264–1,692) | 28.9 (24.4–33.7) | [claim:C262] [evidence:E009]
| J + 2,000/year, 75% to F | 38.2% | 0.0% | 1,506 (1,302–1,713) | 18.6 (15.8–21.5) | [claim:C263] [evidence:E009]
| J + 5,000/year, 25% to F | 2.6% | 0.0% | 3,043 (2,580–3,529) | 84.5 (70.8–99.1) | [claim:C264] [evidence:E009]
| J + 5,000/year, 50% to F | 99.5% | 0.0% | 3,064 (2,620–3,545) | 61.1 (51.5–71.4) | [claim:C265] [evidence:E009]
| J + 5,000/year, 75% to F | 100.0% | 0.0% | 2,743 (2,361–3,164) | 32.5 (27.7–37.8) | [claim:C266] [evidence:E009]
| J + 8,000/year, 25% to F | 94.8% | 0.0% | 4,367 (3,708–5,085) | 121.4 (102.3–142.8) | [claim:C267] [evidence:E009]
| J + 8,000/year, 50% to F | 100.0% | 0.0% | 4,157 (3,557–4,826) | 85.5 (72.3–100.4) | [claim:C268] [evidence:E009]
| J + 8,000/year, 75% to F | 100.0% | 0.0% | 3,414 (2,938–3,967) | 43.9 (37.2–51.2) | [claim:C269] [evidence:E009]
| J + uncapped screening | 100.0% | 100.0% | 8,637 (7,214–10,498) | 272.5 (227.2–334.1) | [claim:C270] [evidence:E009]

Cost sensitivity: with all unit prices halved or increased by half, the discounted cost of 5,000 tests a year with 25% to F was S$148.0 million or S$444.0 million; epidemiological results do not change. [claim:C271] [evidence:E009,E001]

# 9. Sensitivity analysis: prison treatment from 2019

This analysis repeats the main analysis with treatment in J from 1 January 2019 at the same coverage; community screening of F and X still starts in 2026. It uses the same 2,000 posterior draws (Table J; Fig F). [claim:C272] [evidence:E013,E002]

**Table J. Sensitivity analysis with prison treatment from 2019: reductions in 2030 compared with 2015 and probability of meeting both targets.** [claim:C273] [evidence:E013]

| Policy | Fall in new infections, % (95% UI) | Fall in HCV deaths, % (95% UI) | Draws meeting both targets | [claim:C274] [evidence:E013]
|---|---|---|---|
| J only | 85.0% (72.7%–89.7%) | 40.8% (39.5%–42.1%) | 0.0% | [claim:C275] [evidence:E013]
| J + 600/year, 25% to F | 86.3% (74.6%–90.5%) | 42.2% (40.9%–43.6%) | 0.0% | [claim:C276] [evidence:E013]
| J + 600/year, 50% to F | 87.5% (76.5%–91.3%) | 41.8% (40.5%–43.2%) | 0.0% | [claim:C277] [evidence:E013]
| J + 600/year, 75% to F | 88.6% (78.2%–92.1%) | 41.4% (40.1%–42.8%) | 0.0% | [claim:C278] [evidence:E013]
| J + 1,000/year, 25% to F | 87.1% (75.9%–91.0%) | 43.1% (41.8%–44.5%) | 0.0% | [claim:C279] [evidence:E013]
| J + 1,000/year, 50% to F | 88.9% (78.7%–92.4%) | 42.5% (41.2%–43.8%) | 0.0% | [claim:C280] [evidence:E013]
| J + 1,000/year, 75% to F | 90.5% (81.4%–93.5%) | 41.8% (40.5%–43.2%) | 0.0% | [claim:C281] [evidence:E013]
| J + 2,000/year, 25% to F | 88.9% (78.7%–92.4%) | 45.3% (43.9%–46.6%) | 0.0% | [claim:C282] [evidence:E013]
| J + 2,000/year, 50% to F | 91.9% (83.9%–94.6%) | 44.0% (42.7%–45.4%) | 0.0% | [claim:C283] [evidence:E013]
| J + 2,000/year, 75% to F | 94.3% (88.2%–96.2%) | 42.7% (41.4%–44.1%) | 0.0% | [claim:C284] [evidence:E013]
| J + 5,000/year, 25% to F | 93.2% (86.2%–95.5%) | 51.2% (49.4%–52.6%) | 0.0% | [claim:C285] [evidence:E013]
| J + 5,000/year, 50% to F | 97.3% (94.0%–98.4%) | 48.3% (46.7%–49.7%) | 0.0% | [claim:C286] [evidence:E013]
| J + 5,000/year, 75% to F | 99.1% (97.7%–99.5%) | 45.1% (43.7%–46.5%) | 0.0% | [claim:C287] [evidence:E013]
| J + 8,000/year, 25% to F | 96.0% (91.5%–97.4%) | 56.1% (54.0%–57.7%) | 0.0% | [claim:C288] [evidence:E013]
| J + 8,000/year, 50% to F | 99.3% (98.1%–99.7%) | 52.0% (50.2%–53.5%) | 0.0% | [claim:C289] [evidence:E013]
| J + 8,000/year, 75% to F | 99.9% (99.5%–100.0%) | 47.2% (45.7%–48.6%) | 0.0% | [claim:C290] [evidence:E013]
| J + uncapped screening | 99.9% (99.6%–100.0%) | 76.4% (76.1%–76.7%) | 100.0% | [claim:C291] [evidence:E013]

With prison treatment from 2019, new infections under J only were 34% of the 2015 level by 2025, against 100% in the main analysis. [claim:C292] [evidence:E088,E062]

# 10. Secondary outcome: liver and extrahepatic HCV deaths

The WHO mortality target counts deaths from HCV-related cirrhosis and liver cancer, so the primary outcome counts liver deaths only. As a secondary outcome we added the extrahepatic excess deaths of viraemic people, $(\sigma_v - \sigma_n)\mu_i$ times the number of viraemic people, summed over each year (Table K). [claim:C293] [evidence:E002,E032,E004]

With extrahepatic deaths included, zero of the 15 finite policies reduced HCV deaths by 65% or more at the median; prison treatment only reduced them by 25.8% (95% UI 24.4%–27.3%). [claim:C294] [evidence:E009]

**Table K. Fall in HCV deaths in 2030 compared with 2015 by death definition (main analysis).** [claim:C295] [evidence:E009]

| Policy | Liver deaths (primary), % (95% UI) | Liver and extrahepatic deaths, % (95% UI) | [claim:C296] [evidence:E009]
|---|---|---|
| J only | 22.1% (21.2%–23.0%) | 25.8% (24.4%–27.3%) | [claim:C297] [evidence:E009]
| J + 600/year, 25% to F | 24.2% (23.2%–25.1%) | 28.3% (26.8%–29.8%) | [claim:C298] [evidence:E009]
| J + 600/year, 50% to F | 23.8% (22.9%–24.8%) | 27.9% (26.5%–29.6%) | [claim:C299] [evidence:E009]
| J + 600/year, 75% to F | 23.5% (22.5%–24.5%) | 27.6% (26.1%–29.2%) | [claim:C300] [evidence:E009]
| J + 1,000/year, 25% to F | 25.5% (24.4%–26.5%) | 29.9% (28.3%–31.5%) | [claim:C301] [evidence:E009]
| J + 1,000/year, 50% to F | 24.9% (23.9%–25.9%) | 29.3% (27.7%–31.0%) | [claim:C302] [evidence:E009]
| J + 1,000/year, 75% to F | 24.3% (23.3%–25.4%) | 28.6% (27.0%–30.4%) | [claim:C303] [evidence:E009]
| J + 2,000/year, 25% to F | 28.7% (27.4%–29.8%) | 33.7% (31.8%–35.4%) | [claim:C304] [evidence:E009]
| J + 2,000/year, 50% to F | 27.5% (26.3%–28.6%) | 32.4% (30.6%–34.2%) | [claim:C305] [evidence:E009]
| J + 2,000/year, 75% to F | 26.1% (25.0%–27.3%) | 30.9% (29.2%–32.8%) | [claim:C306] [evidence:E009]
| J + 5,000/year, 25% to F | 37.1% (35.0%–38.7%) | 43.6% (40.9%–45.8%) | [claim:C307] [evidence:E009]
| J + 5,000/year, 50% to F | 34.1% (32.3%–35.5%) | 40.2% (37.8%–42.5%) | [claim:C308] [evidence:E009]
| J + 5,000/year, 75% to F | 30.2% (28.8%–31.5%) | 35.8% (33.7%–37.8%) | [claim:C309] [evidence:E009]
| J + 8,000/year, 25% to F | 44.0% (41.4%–45.9%) | 51.6% (48.4%–54.0%) | [claim:C310] [evidence:E009]
| J + 8,000/year, 50% to F | 39.4% (37.3%–41.0%) | 46.3% (43.6%–48.7%) | [claim:C311] [evidence:E009]
| J + 8,000/year, 75% to F | 33.2% (31.7%–34.5%) | 39.2% (37.0%–41.2%) | [claim:C312] [evidence:E009]
| J + uncapped screening | 71.4% (71.3%–71.6%) | 80.6% (80.0%–81.3%) | [claim:C313] [evidence:E009]

# 11. Sensitivity analysis: observed entry of new drug abusers, 2016–2025

The main analysis holds entry into D at its 2015 level in every year. Here entry in each age group from 2016 to 2025 followed the Central Narcotics Bureau count of new drug abusers in that year relative to 2015 [@singstat_drug], and the 2025 level was held to 2030; the 2015 equilibrium, the calibration and the posterior draws are unchanged (Table L; Fig H). [claim:C314] [evidence:E005,E055,E002,E039]

Prison treatment only reduced new infections by 56.0% (95% UI 45.1%–61.9%) and liver deaths by 22.2% (95% UI 21.4%–23.1%) (main analysis 51.5% and 22.1%). Six finite policies met the incidence target (main analysis 5) and zero met the mortality target (main analysis 0). [claim:C315] [evidence:E055,E009]

**Table L. Sensitivity analysis with observed entry, 2016–2025: reductions in 2030 compared with 2015 and probability of meeting both targets.** [claim:C316] [evidence:E055]

| Policy | Fall in new infections, % (95% UI) | Fall in liver deaths, % (95% UI) | Draws meeting both targets | [claim:C317] [evidence:E055]
|---|---|---|---|
| J only | 56.0% (45.1%–61.9%) | 22.2% (21.4%–23.1%) | 0.0% | [claim:C318] [evidence:E055]
| J + 600/year, 25% to F | 59.4% (48.5%–65.0%) | 24.3% (23.4%–25.2%) | 0.0% | [claim:C319] [evidence:E055]
| J + 600/year, 50% to F | 62.7% (51.8%–67.9%) | 24.0% (23.1%–24.9%) | 0.0% | [claim:C320] [evidence:E055]
| J + 600/year, 75% to F | 65.7% (55.0%–70.8%) | 23.6% (22.8%–24.6%) | 0.0% | [claim:C321] [evidence:E055]
| J + 1,000/year, 25% to F | 61.6% (50.7%–66.9%) | 25.7% (24.7%–26.6%) | 0.0% | [claim:C322] [evidence:E055]
| J + 1,000/year, 50% to F | 66.8% (56.1%–71.7%) | 25.1% (24.2%–26.1%) | 0.0% | [claim:C323] [evidence:E055]
| J + 1,000/year, 75% to F | 71.5% (61.3%–76.2%) | 24.5% (23.6%–25.5%) | 0.0% | [claim:C324] [evidence:E055]
| J + 2,000/year, 25% to F | 66.8% (56.1%–71.7%) | 28.9% (27.7%–30.0%) | 0.0% | [claim:C325] [evidence:E055]
| J + 2,000/year, 50% to F | 75.7% (66.2%–80.0%) | 27.8% (26.6%–28.8%) | 0.0% | [claim:C326] [evidence:E055]
| J + 2,000/year, 75% to F | 82.9% (74.9%–86.6%) | 26.4% (25.3%–27.5%) | 0.0% | [claim:C327] [evidence:E055]
| J + 5,000/year, 25% to F | 79.6% (70.7%–83.6%) | 37.4% (35.5%–39.0%) | 0.0% | [claim:C328] [evidence:E055]
| J + 5,000/year, 50% to F | 92.6% (87.8%–94.9%) | 34.4% (32.7%–35.8%) | 0.0% | [claim:C329] [evidence:E055]
| J + 5,000/year, 75% to F | 97.9% (95.6%–98.9%) | 30.5% (29.1%–31.6%) | 0.0% | [claim:C330] [evidence:E055]
| J + 8,000/year, 25% to F | 88.5% (82.1%–91.5%) | 44.4% (41.9%–46.3%) | 0.0% | [claim:C331] [evidence:E055]
| J + 8,000/year, 50% to F | 98.4% (96.4%–99.2%) | 39.7% (37.7%–41.3%) | 0.0% | [claim:C332] [evidence:E055]
| J + 8,000/year, 75% to F | 99.8% (98.8%–99.9%) | 33.4% (31.9%–34.6%) | 0.0% | [claim:C333] [evidence:E055]
| J + uncapped screening | 99.8% (98.8%–99.9%) | 71.5% (71.3%–71.6%) | 100.0% | [claim:C334] [evidence:E055]

# 12. Sensitivity analysis: share of genotype 3

The main analysis assumes that 78% of infections are genotype 3 [@zhao2019]. Here the share was 10% or 50%. The share moves the genotype-weighted progression rates and the cure rates (Section 1.4): with RBV-free treatment, cure in decompensated cirrhosis was 85.5% at 10% and 67.5% at 50%, against 54.9% in the main analysis. The calibrated parameters and the posterior draws were not refitted (Table M). [claim:C335] [evidence:E056,E004,E072]

Prison treatment only reduced liver deaths by 19.7% (95% UI 18.8%–20.7%) at 10% and 21.2% (95% UI 20.3%–22.1%) at 50% (main analysis 22.1%); uncapped screening reduced them by 71.1% and 71.3%. The same policies met each target as in the main analysis. [claim:C336] [evidence:E056,E009]

**Table M. Sensitivity analysis with a genotype 3 share of 10% and 50%: falls in 2030 compared with 2015 (median, 95% UI).** [claim:C337] [evidence:E056]

| Policy | New infections, 10% | Liver deaths, 10% | New infections, 50% | Liver deaths, 50% | [claim:C338] [evidence:E056]
|---|---|---|---|---|
| J only | 51.2% (40.2%–56.9%) | 19.7% (18.8%–20.7%) | 51.4% (40.4%–57.1%) | 21.2% (20.3%–22.1%) | [claim:C339] [evidence:E056]
| J + 600/year, 25% to F | 54.6% (43.6%–60.1%) | 22.0% (21.1%–23.0%) | 54.8% (43.9%–60.3%) | 23.3% (22.4%–24.3%) | [claim:C340] [evidence:E056]
| J + 600/year, 50% to F | 57.9% (47.1%–63.3%) | 21.6% (20.6%–22.6%) | 58.1% (47.3%–63.5%) | 22.9% (22.0%–24.0%) | [claim:C341] [evidence:E056]
| J + 600/year, 75% to F | 61.1% (50.4%–66.4%) | 21.1% (20.2%–22.1%) | 61.3% (50.6%–66.5%) | 22.5% (21.6%–23.6%) | [claim:C342] [evidence:E056]
| J + 1,000/year, 25% to F | 56.8% (46.0%–62.2%) | 23.5% (22.5%–24.6%) | 57.0% (46.2%–62.4%) | 24.7% (23.7%–25.7%) | [claim:C343] [evidence:E056]
| J + 1,000/year, 50% to F | 62.2% (51.6%–67.4%) | 22.8% (21.8%–23.8%) | 62.3% (51.8%–67.5%) | 24.1% (23.1%–25.1%) | [claim:C344] [evidence:E056]
| J + 1,000/year, 75% to F | 67.0% (56.8%–71.9%) | 22.0% (21.0%–23.0%) | 67.2% (56.9%–72.1%) | 23.4% (22.4%–24.4%) | [claim:C345] [evidence:E056]
| J + 2,000/year, 25% to F | 62.2% (51.6%–67.4%) | 27.1% (25.8%–28.3%) | 62.3% (51.8%–67.5%) | 28.1% (26.8%–29.2%) | [claim:C346] [evidence:E056]
| J + 2,000/year, 50% to F | 71.5% (61.8%–76.1%) | 25.6% (24.3%–26.8%) | 71.6% (61.9%–76.2%) | 26.8% (25.5%–27.9%) | [claim:C347] [evidence:E056]
| J + 2,000/year, 75% to F | 79.3% (70.9%–83.4%) | 23.9% (22.7%–25.0%) | 79.3% (71.0%–83.4%) | 25.3% (24.1%–26.4%) | [claim:C348] [evidence:E056]
| J + 5,000/year, 25% to F | 75.6% (66.5%–80.0%) | 36.6% (34.4%–38.3%) | 75.7% (66.6%–80.0%) | 36.9% (34.8%–38.5%) | [claim:C349] [evidence:E056]
| J + 5,000/year, 50% to F | 90.1% (84.3%–93.0%) | 32.9% (31.0%–34.3%) | 90.1% (84.3%–93.0%) | 33.6% (31.8%–35.0%) | [claim:C350] [evidence:E056]
| J + 5,000/year, 75% to F | 96.8% (93.8%–98.2%) | 28.2% (26.8%–29.5%) | 96.8% (93.7%–98.2%) | 29.4% (28.0%–30.7%) | [claim:C351] [evidence:E056]
| J + 8,000/year, 25% to F | 85.4% (78.3%–88.9%) | 44.2% (41.4%–46.3%) | 85.4% (78.3%–88.9%) | 44.1% (41.5%–46.1%) | [claim:C352] [evidence:E056]
| J + 8,000/year, 50% to F | 97.5% (94.9%–98.7%) | 38.7% (36.5%–40.5%) | 97.5% (94.9%–98.7%) | 39.2% (37.0%–40.8%) | [claim:C353] [evidence:E056]
| J + 8,000/year, 75% to F | 99.7% (98.6%–99.9%) | 31.5% (29.9%–32.9%) | 99.7% (98.6%–99.9%) | 32.5% (31.0%–33.8%) | [claim:C354] [evidence:E056]
| J + uncapped screening | 99.8% (98.8%–99.9%) | 71.1% (70.9%–71.2%) | 99.8% (98.7%–99.9%) | 71.3% (71.2%–71.5%) | [claim:C355] [evidence:E056]

# 13. One-way sensitivity analysis of the fixed inputs

The uncertainty intervals cover the 13 calibrated parameters only. To show how much the fixed inputs matter, we held the calibrated parameters at the posterior mode (Section 4.4) and set each fixed input in turn to the low and the high end of its range, without refitting, running all 17 policies each time (Table N; Fig G). Ranges are the intervals of the technical parameter table where given, and otherwise from the literature; prison rates use the 95% credible intervals of the rate models (Table B). Inputs from one estimate were moved together (the two genotype 3 relative risks for cirrhosis, the two for carcinoma, and progression to carcinoma from stages 2 and 3). Treatment duration and background death rates were not varied. [claim:C356] [evidence:E053,E054,E089,E079]

**Table N. Fixed inputs varied in the one-way sensitivity analysis.** [claim:C357] [evidence:E054]

| Input | Base | Low | High | Source of range | [claim:C358] [evidence:E054]
|---|---|---|---|---|
| Transmission probability q | 0.0057 | 0.0032 | 0.0105 | PDF Table 1 (Boelen 2014) | [claim:C359] [evidence:E054]
| Spontaneous clearance | 0.26 | 0.22 | 0.29 | PDF Table 1 (Micallef 2006) | [claim:C360] [evidence:E054]
| Acute duration (years) | 0.4932 | 0.25 | 1 | Grebely 2014 (3-12 months) | [claim:C361] [evidence:E054]
| GT3 share | 0.78; 0.9787; 0.549; 0.8818; 0.8818 | 0.74; 0.9792; 0.567; 0.8838; 0.8838 | 0.82; 0.9782; 0.531; 0.8798; 0.8798 | PDF Table 1 (Zhao 2019) | [claim:C362] [evidence:E054]
| NC to CC | 0.027 | 0.02543 | 0.02857 | PDF Table 2 (Thein 2008) | [claim:C363] [evidence:E054]
| CC to DC | 0.039 | 0.02183 | 0.06072 | PDF Table 2 (Lim 2018) | [claim:C364] [evidence:E054]
| CC/DC to HCC | 0.014; 0.014 | 0.00155; 0.00155 | 0.03875; 0.03875 | PDF Table 2 (Lim 2018, one parameter) | [claim:C365] [evidence:E054]
| GT3 RR, cirrhosis | 1.36; 1.36 | 1.27; 1.27 | 1.46; 1.46 | PDF Table 2 (derived from Kanwal 2014) | [claim:C366] [evidence:E054]
| GT3 RR, HCC | 1.93; 1.93 | 1.71; 1.71 | 2.17; 2.17 | PDF Table 2 (derived from Kanwal 2014) | [claim:C367] [evidence:E054]
| Post-SVR RR, CC to DC | 0.07 | 0.03 | 0.2 | PDF Table 2 (van der Meer 2012) | [claim:C368] [evidence:E054]
| Post-SVR RR, CC to HCC | 0.23 | 0.16 | 0.35 | PDF Table 2 (Morgan 2013) | [claim:C369] [evidence:E054]
| Post-SVR RR, DC to HCC | 1 | 0.52 | 1 | Pageaux 2022 (HR 0.52); base 1 assumed | [claim:C370] [evidence:E054]
| Liver mortality, DC | 0.13 | 0.111 | 0.1502 | PDF Table 2 (Lim 2018) | [claim:C371] [evidence:E054]
| Liver mortality, HCC | 0.43 | 0.3719 | 0.4892 | PDF Table 2 (Lim 2018) | [claim:C372] [evidence:E054]
| Post-SVR mortality RR, DC | 0.4 | 0.17 | 0.96 | Pageaux 2022 (liver death HR 0.40, 0.17-0.96) | [claim:C373] [evidence:E054]
| Post-SVR mortality RR, HCC | 0.34 | 0.13 | 0.88 | Dang 2020 (liver death HR 0.34, 0.13-0.88) | [claim:C374] [evidence:E054]
| PWID SMR | 4.73; 6.953 | 3.47; 5.101 | 6.45; 9.482 | van Santen 2014 Table 3 | [claim:C375] [evidence:E054]
| Extrahepatic HR, viraemic | 6.953 | 5.818 | 8.372 | Lee 2012 Table 3 | [claim:C376] [evidence:E054]
| SVR GT3, NC/CC | 0.9787 | 0.9725 | 0.9834 | PDF Table 1 (EASL 2020) | [claim:C377] [evidence:E054]
| SVR other GT, NC/CC | 0.9787 | 0.9778 | 0.9794 | PDF Table 1 (EASL 2020) | [claim:C378] [evidence:E054]
| SVR GT3, DC (RBV-free) | 0.549 | 0.2994 | 0.8142 | PDF Table 1 (EASL 2020; Xiao 2025) | [claim:C379] [evidence:E054]
| SVR other GT, DC (RBV-free) | 0.549 | 0.5314 | 0.5622 | PDF Table 1 (EASL 2020; Xiao 2025) | [claim:C380] [evidence:E054]
| SVR GT3, HCC (RBV+) | 0.8818; 0.8818 | 0.757; 0.757 | 0.9676; 0.9676 | PDF Table 1 (EASL 2020; Xiao 2025) | [claim:C381] [evidence:E054]
| SVR other GT, HCC (RBV+) | 0.8818; 0.8818 | 0.8708; 0.8708 | 0.8884; 0.8884 | PDF Table 1 (EASL 2020; Xiao 2025) | [claim:C382] [evidence:E054]
| Release rates | 0.462; 0.636; 0.454; 0.337; 0.325; 0.394 | 0.4224; 0.6215; 0.4441; 0.329; 0.3161; 0.3722 | 0.5033; 0.6513; 0.4648; 0.3449; 0.3343; 0.4158 | Stan fit 95% interval | [claim:C383] [evidence:E054]
| Re-arrest rates | 0.548; 0.478; 0.469; 0.471; 0.384; 0.265 | 0.4594; 0.4555; 0.4489; 0.451; 0.3632; 0.2349 | 0.6473; 0.5006; 0.4897; 0.4908; 0.405; 0.298 | Stan fit 95% interval | [claim:C384] [evidence:E054]
| Return to injecting after release | 0.724 | 0.7136 | 0.7344 | Stan fit 95% interval | [claim:C385] [evidence:E054]

For j only, the fall in new infections at the mode was 53.9% and ranged from 8.5% (transmission probability q, high) to 57.5% (re-arrest rates, high); the fall in liver deaths was 22.1% and ranged from 18.8% (CC/DC to HCC, low) to 26.1% (CC/DC to HCC, high). [claim:C386] [evidence:E015]

For J + 8,000/year, 25% to F, the fall in new infections at the mode was 87.8% and ranged from 57.6% (transmission probability q, high) to 89.9% (acute duration, low); the fall in liver deaths was 45.0% and ranged from 37.7% (post-SVR mortality RR, DC, high) to 49.3% (post-SVR mortality RR, HCC, low). [claim:C387] [evidence:E015]

Pass or fail at the mode changed for at least one policy when one of these inputs was moved: acute duration (both ends); re-arrest rates (both ends); post-SVR mortality RR, DC (high end); transmission probability q (high end) (output/tornado/pass_flips.csv). With transmission probability q at its low end or PWID SMR at its low end, the 2015 snapshot was not at equilibrium after 150 years and the run was excluded. [claim:C388] [evidence:E065,E015]

# 14. Assumptions and limitations

Table O lists the main assumptions of the model and analysis, the basis for each and its likely effect on the results. Where the direction of an effect is not clear, the table says so.

**Table O. Main assumptions, their basis and likely effect on the results.**

| Area | Assumption | Basis | Likely effect |
|---|---|---|---|
| Who injects | New drug abusers whose main drug was heroin or methamphetamine (86.1%) are the people who start injecting | 113 heroin and 1,014 methamphetamine of 1,309 new abusers in 2015 [@singstat_m891521]; the table records the main drug, not the route of use | Non-injecting users of these drugs raise entry; injectors of other drugs are missed. Net direction unknown | [claim:C389] [evidence:E006,E040]
| Who injects | The same share at every age | The drug table has no age breakdown | Direction unknown; the entry scales can raise entry in any age group | [claim:C390] [evidence:E006,E041]
| Entry | Entry into D is at least the scaled 2015 count and is the same in every year, before and after 2015 | Only 2015 counts were used; entry scales are bounded below by 1 | A fall in new injectors after 2015 would make the reductions larger than modelled; a rise would make them smaller | [claim:C391] [evidence:E004,E032,E041]
| Steady state | The epidemic was at equilibrium in 2015 after 150 model years with constant parameters | Needed to start the model without historical data (Section 1.6) | Past trends are ignored; the 2016–2025 history without treatment stays at the 2015 level | [claim:C392] [evidence:E032,E037]
| First arrest | The first-arrest rate is $c$ times the re-arrest rate, with the same age pattern | The records contain only people who were arrested; $c$ is calibrated | $c$ is mainly informed by its prior (posterior SD 72% of the prior's) | [claim:C393] [evidence:E004,E014]
| Release | Released people who will be re-arrested (72.4%) return to injecting (F); the others stop permanently (X) and are never re-arrested | Cure-type survival model of re-arrest (Section 3) | Relapse after a long gap would put more infections into F; people in X can be reached only by community screening | [claim:C394] [evidence:E010,E004]
| Prison rates | Release and re-arrest rates are constant over time and the same in every HCV state | Estimated from records of people held for drug consumption, admissions 2000–2017 | Direction unknown | [claim:C395] [evidence:E010,E011,E004]
| Prison population | J represents everyone held for drug consumption, and its seroprevalence is that of the screened subset | Modelled 8,768 against 8,725 people held for drug consumption at the end of 2015, of whom 3,367 were screened | If people screened differ from the other people held, the prison seroprevalence targets, and so transmission, are biased; direction unknown | [claim:C396] [evidence:E061,E042]
| Transmission | Only people in D and F transmit and are infected; there is no transmission in prison or after people stop injecting | Model structure (Section 1.2) | Any transmission in prison would add infections that prison treatment could prevent | [claim:C397] [evidence:E032]
| Mixing | Contact depends only on age, through a base matrix with a calibrated multiplier per age group, and is the same for D and F; reinfection occurs at the same rate as first infection | Model structure (Section 1.2) | High-risk networks are not represented; direction unknown | [claim:C398] [evidence:E032,E004]
| Natural history | Progression, clearance and treatment efficacy come from the literature, with a fixed genotype mix | Table A | Not included in the uncertainty intervals; Section 13 varies each value one at a time, and Section 12 the genotype mix | [claim:C399] [evidence:E004,E015]
| Mortality | All PWID, including those who stopped injecting, die of causes other than liver disease at the 2015 Singapore rate (ages 60–74 for the oldest group) times 4.73, and viraemic people times 1.47 more; liver deaths occur only in decompensated cirrhosis and liver cancer | Mortality ratio of people who use drugs without HCV [@vansanten2014]; extrahepatic hazard ratio of viraemic infection from a general-population cohort [@lee2012] | The extrahepatic excess is not found in every cohort [@piselli2021]; it does not enter the WHO outcome (Section 10) | [claim:C400] [evidence:E004,E032,E035,E036,E067]
| Seroprevalence targets | Antibody prevalence among people in prison with a screening result represents people in J by age | Prison records, December 2014 to February 2016 (Section 4.1) | Targets differ from the current processing of the records by up to 0.5 percentage points | [claim:C401] [evidence:E007,E042]
| Community target | The 2017 network scale-up estimate applies to 2015 | Entry is constant and there is no treatment before 2019 in any analysis | Small if the injecting population changed little between 2015 and 2017 | [claim:C402] [evidence:E007]
| Treatment history | Nobody was treated before 2026 in the main analysis | Start and coverage of treatment in prison are not recorded in the project data [@prison2021; @tan2023] | Understates past treatment; the sensitivity analysis (Section 9) gives the other bound | [claim:C403] [evidence:E001,E013,E017,E024]
| Prison treatment | People with chronic infection in J start treatment with annual probability 0.99995 in the main analysis | Policy assumption; varied from 0 in the coverage analysis | Close to the most that prison treatment can achieve in the model | [claim:C404] [evidence:E001,E012]
| Screening | Tests are spread over each group in proportion to its HCV states, with perfect linkage within capacity, no loss to follow-up and repeat screening allowed; D is never screened or treated | Aggregate model (Section 6) | Real programmes will find and treat fewer people for the same number of tests | [claim:C405] [evidence:E002]
| Targets | Reductions are measured against the modelled 2015 yearly numbers, and a target is met if the median reduction reaches it | WHO 2030 targets for incidence and mortality | Results depend on the modelled 2015 level, not on observed counts | [claim:C406] [evidence:E002,E001]
| Costs | Published charges and a published drug price, not inflated and not payer costs; follow-up visits, tests of cure and care navigation are excluded | Section 7 | Totals scale almost directly with the drug price | [claim:C407] [evidence:E047,E001]
| Uncertainty | Intervals reflect posterior uncertainty in the 13 calibrated parameters only, given the fixed inputs; the model is deterministic | Section 4 | Fixed inputs, structural, price and uptake uncertainty are not included, so intervals are likely too narrow; Section 13 shows the effect of each fixed input | [claim:C408] [evidence:E009,E008,E015]

# 15. Software and reproducibility

The model is in C++ (scripts/model/sim.cpp), called from R through Rcpp [@rcpp2011]. Rate models used Stan and convergence diagnostics the R package posterior. Figures were drawn with matplotlib from the saved outputs. [claim:C409] [evidence:E032,E090,E057]

Analysis commands, run times and expected outputs are listed in docs/training.md in the repository. This appendix and the main text were generated from the saved outputs by scripts/build_paper.py; drafting followed the procedures of Scientific Agent Skills [@kassis2026]. [claim:C410] [evidence:E090,E068]

# Figures

**Fig A. MCMC trace plots.** Stored draws of the 13 parameters for each of the 8 chains after warm-up. [claim:C411] [evidence:E008]

![](output/calibration/mcmc_reference/mcmc_traceplots.png){width=95%}

**Fig B. Marginal posterior densities of the 13 parameters from the MCMC.** [claim:C412] [evidence:E008]

![](output/calibration/mcmc_reference/figures/marginals.png){width=95%}

**Fig C. Target attainment in 2030 by policy.** Median reductions in new infections and HCV deaths, and whether both targets were met. [claim:C413] [evidence:E009]

![](output/scenarios_2026/scenario_target_heatmap.png){width=95%}

**Fig D. Costs and treatment starts, 2026–2030.** (A) Discounted costs by group. (B) Undiscounted treatment courses by group. Bars are medians of each group; whiskers are 95% uncertainty intervals of the total. [claim:C414] [evidence:E091,E009]

![](output/scenarios_2026/scenario_cumulative_costs.png){width=95%}

**Fig E. Screening cascade in F and X by year.** Number screened, RNA tests, RNA-positive people, chronic infections found and treatment starts, by capacity (columns) and share of tests given to F (rows).

![](output/scenarios_2026/scenario_care_cascade.png){width=95%}

**Fig F. Sensitivity analysis with prison treatment from 2019.** (A) New infections and HCV deaths relative to 2015 (medians), 2015–2030; the dotted line marks the start of prison treatment (2019) and the dashed line the start of community screening (2026). (B) Reductions in 2030 by policy, with 95% uncertainty intervals. [claim:C415] [evidence:E013,E088]

![](output/scenarios_2026_j2019/scenario_target_trajectories.png){width=95%}

**Fig G. One-way sensitivity of the 2030 reductions to the fixed inputs.** Prison treatment only (A, B) and the finite policy with the largest fall in liver deaths (C, D). Bars show the result with each input at the low (blue) or high (orange) end of its range; the solid line is the result at the posterior mode and the dashed line the WHO target. [claim:C416] [evidence:E015]

![](output/tornado/tornado.png){width=95%}

**Fig H. Sensitivity analysis with observed entry of new drug abusers, 2016–2025.** (A) New infections and HCV deaths relative to 2015 (medians), 2015–2030; the dashed line marks the start of prison treatment and community screening (2026). (B) Reductions in 2030 by policy, with 95% uncertainty intervals. [claim:C417] [evidence:E055,E092]

![](output/scenarios_2026_inflow/scenario_target_trajectories.png){width=95%}

# References

::: {#refs}
:::
