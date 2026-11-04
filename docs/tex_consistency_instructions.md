## 1. Changes that apply to every equation in 3.2–3.5

### 1.1 Background mortality: replace `\omega\mu_i` by state-specific SMRs (most important)

The model no longer uses one PWID SMR `\omega` (14.68). Each cause of death has one route
(`sim.cpp:74-85`):

| HCV state | Non-liver death rate | Symbol to use |
|---|---|---|
| u, r (non-viraemic) | `mu_i * smr_nonviraemic` | `\sigma_n\,\mu_i` |
| a, c, t (viraemic) | `mu_i * smr_viraemic` | `\sigma_v\,\mu_i` |

with `\sigma_n = 4.73` and `\sigma_v = \sigma_n \times 1.47 = 6.95`.

Edit rule, in every `d/dt` line of D, J, F and X:
- rows for `u` and `r`: `\omega\mu_i` → `\sigma_n\mu_i`
- rows for `a`, `c` and `t`: `\omega\mu_i` → `\sigma_v\mu_i`

### 1.2 Liver deaths in stages 3 and 4

Liver deaths occur only in stage 3 (DC) and stage 4 (HCC). They are added to the background rate. The
post-SVR modifier `\psi` multiplies **only the liver death rate**, and **only for cured people (r)**:
people on treatment (t) keep the untreated liver death rate (author decision 2026-11-04).

| Stage | u, a, c, t | r |
|---|---|---|
| 3 (DC) | `+ \mu^{\mathrm{DC}}` | `+ \psi^{\mathrm{DC}}\mu^{\mathrm{DC}}` |
| 4 (HCC) | `+ \mu^{\mathrm{HCC}}` | `+ \psi^{\mathrm{HCC}}\mu^{\mathrm{HCC}}` |

Edit rule:
- `\omega\mu^{\mathrm{DC}}` → `\mu^{\mathrm{DC}}` and `\omega\mu^{\mathrm{HCC}}` → `\mu^{\mathrm{HCC}}` (u, a, c rows).
  `\mu^{\mathrm{DC}}`, `\mu^{\mathrm{HCC}}` are not multiplied by any SMR.
- **t rows of stages 3 and 4** (`dD_{t,3,i}`, `dD_{t,4,i}`, and the same for J, F, X):
  `\omega\mu^{\mathrm{DC}}_{\mathrm{SVR}}` → `\mu^{\mathrm{DC}}`;
  `\omega\mu^{\mathrm{HCC}}_{\mathrm{SVR}}` → `\mu^{\mathrm{HCC}}`.
- **r rows of stages 3 and 4**: `\omega\mu^{\mathrm{DC}}_{\mathrm{SVR}}` → `\psi^{\mathrm{DC}}\mu^{\mathrm{DC}}`;
  `\omega\mu^{\mathrm{HCC}}_{\mathrm{SVR}}` → `\psi^{\mathrm{HCC}}\mu^{\mathrm{HCC}}`.

Example, D stage 3:

```latex
\frac{dD_{t,3,i}}{dt}
  &= \tau_{D,3}\,D_{c,3,i}
     - \left(\frac{1}{\iota_2} + \lambda_i^{(1)}
             + \sigma_v\mu_i + \mu^{\mathrm{DC}}\right)D_{t,3,i}
\\[0.8em]
\frac{dD_{r,3,i}}{dt}
  &= \frac{\kappa}{\iota_1}\,D_{a,3,i}
     + \frac{\bar{\alpha}_3^{-}}{\iota_2}\,D_{t,3,i}
     + \tilde{p}_{r,2\to3}\,D_{r,2,i}
     - \bigl(\gamma_i(t) + \lambda_i^{(1)} + \sigma_n\mu_i
             + \psi^{\mathrm{DC}}\mu^{\mathrm{DC}}
             + \tilde{p}_{r,3\to4}\bigr)\,D_{r,3,i}
```

### 1.3 Force of infection: `\gamma_{i,j}(t)` → `\gamma_i(t)`

The force of infection depends only on the susceptible person's age group `i` (it already sums over
`j`). Replace every `\gamma_{i,j}(t)` by `\gamma_i(t)` and define it once in 3.1 (Section 3).

### 1.4 First arrest: `\lambda^{(1)}_i = c\,\lambda^{(3)}_i`

`\lambda^{(1)}_i` is not fitted to prison data. It equals the calibrated multiplier `c` times the
re-arrest rate (`sim.cpp:114`, `l1 = c_first * lambda3`). Keep the symbol `\lambda^{(1)}_i` in the
equations and define it in 3.1.

### 1.5 Treatment rates: `\tau_k` → `\tau_{S,k}`

Treatment initiation is stratum-specific (`tau_by_stratum`, 16 values). Replace `\tau_k` by `\tau_{D,k}`,
`\tau_{J,k}`, `\tau_{F,k}`, `\tau_{X,k}` in the D, J, F, X subsections respectively (both the outflow
from `c` and the inflow to `t`). Values are defined in 3.1 (Section 3).

### 1.6 SVR in stages 3 and 4: resolve `\bar{\alpha}_3^{\pm}`

The costed regimen is RBV-free SOF/VEL, so every reported analysis uses (`scenarios_2026.R`,
`policy_params$alpha_DC <- alpha_DC_neg`):
- stage 3 rows (`t → r` and `t → c`): `\bar{\alpha}_3^{\pm}` → `\bar{\alpha}_3^{-}`;
- stage 4 rows: `\bar{\alpha}_3^{\pm}` → `\bar{\alpha}_4 \equiv \bar{\alpha}_3^{+}` (`alpha_HCC = 0.8818`,
  RBV+ composite; no HCC-specific trial data).

### 1.7 Ageing (missing from every equation)

Every compartment ages at rate 1/10 per year into the next 10-year group; nobody leaves the oldest
group (`sim.cpp`, `AGE_RATE = 0.1`, applied after each daily RK4 step). Do **not** add the term to all
480 equations; add one statement in 3.1 (Section 3) and keep the equations as written.

## 2. Already consistent — leave as is

- Entry paragraph at the top of `equations.tex` (constant entry `\beta_i = b_i s_i` into
  `D_{u,1,i}` only; observed-inflow sensitivity). Updated and authorised earlier.
- Genotype-weighted progression `\tilde{p}_{c,\cdot}` and post-SVR progression
  `\tilde{p}_{r,\cdot} = \phi\,\tilde{p}_{c,\cdot}`; recovered people do not progress from stage 1.
- The HCV state transitions (`a → r` at `\kappa/\iota_1`, `a → c` at `(1-\kappa)/\iota_1`,
  `c → t` at `\tau`, `t → r` at `\bar\alpha/\iota_2`, `t → c` at `(1-\bar\alpha)/\iota_2`,
  `r → a` at `\gamma`), and that only D and F are infected; u and r at every stage can be infected.
- Movements between strata: D → J at `\lambda^{(1)}_i`, F → J at `\lambda^{(3)}_i`, J → F at
  `\pi\lambda^{(2)}_i`, J → X at `(1-\pi)\lambda^{(2)}_i`; X is never re-arrested. All HCV states and
  stages move.
- `\bar{\alpha}_1` for stages 1–2.

Optional notation clean-up: Section 3.1 writes the baseline progression as `p_3^{\cdot}` for both
genotypes (`\rho r_3 p_3 + (1-\rho) p_3`). Because the baseline is the **other**-genotype rate, rename it
`p_o^{\cdot}` (or `p^{\cdot}`) in 3.1 and Table 2.

## 3. Section 3.1: replacement text

1. **Delete** the block

   ```latex
   \noindent The effective disease-specific mortality rates incorporating
   SVR benefit are:
   \begin{align*}
   \omega\mu^{\mathrm{DC}}_{\mathrm{SVR}} &= \psi^{\mathrm{DC}}\,(\omega\mu^{\mathrm{DC}}+\omega\mu_i) \\
   \omega\mu^{\mathrm{HCC}}_{\mathrm{SVR}} &= \psi^{\mathrm{HCC}}\,(\omega\mu^{\mathrm{HCC}}+\omega\mu_i)
   \end{align*}
   ```

   It is wrong for the current model in three ways: `\psi` also scaled the background rate, the
   background rate was then added a second time in each equation, and `\omega` multiplied the liver
   death rates.

2. **Insert** in its place (adapt wording, keep the content):

   ```latex
   \noindent\textbf{Mortality.} Each cause of death has one route. Deaths from causes other than
   liver disease occur at $\sigma_n\mu_i$ for non-viraemic people ($u$, $r$) and at
   $\sigma_v\mu_i$ for viraemic people ($a$, $c$, $t$), where $\mu_i$ is the 2015 Singapore resident
   death rate of age group $i$, $\sigma_n$ is the standardised mortality ratio of people who use
   drugs without HCV or HIV infection, and $\sigma_v=\sigma_n\,\mathrm{HR}_v$ with $\mathrm{HR}_v$
   the hazard ratio of extrahepatic death for viraemic infection. Liver deaths occur only in
   stages 3 and 4, at $\mu^{\mathrm{DC}}$ and $\mu^{\mathrm{HCC}}$, added to the background rate.
   After cure ($r$), the liver death rate is multiplied by $\psi^{\mathrm{DC}}$ or
   $\psi^{\mathrm{HCC}}$; people on treatment ($t$) keep the untreated rate. The death rate of a
   person in stratum $S$, stage $k$, state $h$ and age $i$ is
   \begin{align*}
   m_{k,h,i} &= \sigma_h\,\mu_i
     + \mathbb{1}[k=3]\,\psi_h^{\mathrm{DC}}\mu^{\mathrm{DC}}
     + \mathbb{1}[k=4]\,\psi_h^{\mathrm{HCC}}\mu^{\mathrm{HCC}},\\
   \sigma_h &= \begin{cases}\sigma_v & h\in\{a,c,t\}\\ \sigma_n & h\in\{u,r\}\end{cases},\qquad
   \psi_h^{\cdot} = \begin{cases}\psi^{\cdot} & h=r\\ 1 & \text{otherwise.}\end{cases}
   \end{align*}
   The WHO mortality outcome counts liver deaths only
   ($\mu^{\mathrm{DC}}$ and $\mu^{\mathrm{HCC}}$ terms); the secondary outcome adds the
   extrahepatic excess $(\sigma_v-\sigma_n)\mu_i$ of viraemic people.
   ```

3. **Insert** the definitions used by 3.2–3.5:

   ```latex
   \noindent\textbf{Force of infection.} For susceptible or recovered people in D and F of age
   group $i$,
   \[
   \gamma_i(t) = q\sum_{j=1}^{6} C_{ij}\,e^{\theta_i}\,
     \frac{\sum_{S\in\{D,F\}}\sum_k\bigl(S_{a,k,j}+S_{c,k,j}\bigr)}
          {\sum_{S\in\{D,F\}}\sum_k\sum_h S_{h,k,j}},
   \]
   where $C$ is the base contact matrix and $e^{\theta_i}$ the calibrated contact scale of row $i$.
   People in J and X are not infected and do not transmit.

   \noindent\textbf{Arrest.} The first-arrest rate is $\lambda^{(1)}_i = c\,\lambda^{(3)}_i$, with
   $c$ calibrated; $\lambda^{(2)}_i$ (release), $\lambda^{(3)}_i$ (re-arrest) and $\pi$ (probability
   of eventual re-arrest) are estimated from prison records.

   \noindent\textbf{Treatment.} $\tau_{S,k}$ is the treatment initiation rate in stratum $S$ and
   stage $k$; it is the same at every stage and age. $\tau_{D,k}=0$ always. Before treatment starts
   (calendar 2026 in the main analysis, 2019 in a sensitivity analysis) every $\tau_{S,k}=0$; after,
   $\tau_{J,k}=-\ln(1-0.99995)$ per year (varied: annual probability 0, 0.5, 0.9), and
   $\tau_{F,k}=\tau_{X,k}=-\ln(1-0.99)$ per year from 2026 until the year's treatment capacity
   (chronic infections found by that year's screening quota) is used, then 0.

   \noindent\textbf{Ageing.} Every compartment moves to the next age group at rate $1/10$ per year
   (10-year groups; the oldest group, 60+, has no outflow). The term is omitted from the equations
   below.
   ```

## 4. Section 4 (Notation): row changes

Values below are from `setup.R`, `likelihood.R`, `data/calibration_targets.json` and the MCMC posterior
(`output/calibration/mcmc_reference/`). Keep every other row as it is (q, κ, ι1, ι2, ρ, the six SVR12
rates, the four progression rates, the four GT3 relative risks and the three φ were checked and match
`setup.R`).

### Table 1

| Row | Current | Change to |
|---|---|---|
| HCV states note | "Recovered individuals can be reinfected" | Add: "r is antibody-reactive and non-infectious; only people in D and F are reinfected." |
| `\gamma_{i,j}(t)` | force of infection "from infectious age group j" | `\gamma_i(t)`, force of infection on age group i (defined in 3.1); model-derived |
| ρ | 0.78, U(0.74–0.82) [20] | Keep; add "sensitivity analysis 0.10 and 0.50" |
| `\tau_k` | Treatment initiation rate under scenario k = 1..4 | `\tau_{S,k}`: rate in stratum S and stage k. J: −ln(1−0.99995)/yr; F, X: −ln(1−0.99)/yr until annual capacity used; D: 0. Policy assumption |
| `\lambda^{(1)}_i` | First-arrest rate, i = 1..9, fitted to SPS data | `c\,\lambda^{(3)}_i`; derived. Age groups i = 1..6 (<20, 20–29, 30–39, 40–49, 50–59, 60+) |
| `c` (new row) | — | First-arrest multiplier; calibrated (θ13). Prior log c ~ N(0, 0.5) truncated to c ∈ [0.5, 100]; posterior median 2.71 (95% CrI 1.54–5.23) |
| `\lambda^{(2)}_i` | Release rate, age-varying, fitted to SPS data | 0.462, 0.636, 0.454, 0.337, 0.325, 0.394 per year (posterior means; 95% CrIs in S1 Text Table B). Bayesian survival model of 33,557 stays of 20,720 people held for drug consumption, admissions 2000–2017 |
| `\lambda^{(3)}_i` | Re-arrest rate, age-varying, fitted to SPS data | 0.548, 0.478, 0.469, 0.471, 0.384, 0.265 per year; cure-type survival model, same records |
| π | Recidivism probability 0.65, U(0.54–0.70) | Probability of eventual re-arrest 0.724 (95% CrI 0.714–0.734), same model; re-arrested people return to injecting (F), the rest stop (X) |
| `\beta_i(t)` | Entry rate, time- and age-varying | `\beta_i = b_i\,s_i`, constant in time (main analysis). `b_i` = CNB 2015 new drug abusers by age (257, 644, 292, 94, 21, 1) × 0.861 (heroin + methamphetamine share, SingStat M891521). `s_i = e^{\theta_{6+i}}` calibrated, prior N⁺(0, 1) on log scale, s_i ∈ [1, 100]. Sensitivity: observed CNB counts 2016–2025 |
| `C_{ij}`, `e^{\theta_i}` (new rows) | — | Base contact matrix (source to be cited by the authors) and calibrated contact scales θ1–θ6, prior N(0, 1.5), scale ∈ [0.01, 1000] |
| `N_0` | Initial population of ever-PWID in 2017, 11,000 (6,500–17,000) [15] | Community injectors (D + F), 11,000 (95% CI 6,500–17,000) [15]. **Calibration target** for the modelled 2015 D + F (log-normal, SD 0.2453). Not an initial condition: the model starts from an arbitrary `y_0` and runs 150 years to equilibrium (model year 150 = 2015) |
| "HCV seroprevalence among ever-PWIDs in 2017, 46.5%, Bernoulli" | — | **Remove.** Replace with "HCV antibody prevalence in prison by age, 2015: 12.7%, 19.2%, 30.0%, 48.3%, 53.1%, 42.9% (3,367 people screened); calibration target, Beta likelihood with precision φ = 100" |
| First-arrest target (new row) | — | First arrests by age, 2015: 221, 554, 251, 81, 18, 1 (CNB new drug abusers × 0.861, rounded); calibration target, negative binomial, size 50 |

### Table 2

| Row | Current | Change to |
|---|---|---|
| `p_3^{\cdot}` symbols | "(other GT)" with subscript 3 | Optional: rename to `p_o^{\cdot}` (Section 2) |
| `\mu_i` | Background mortality, fitted to SingStat data | 2015 Singapore resident age-specific death rates (SingStat), fixed: 0.20, 0.25, 0.45, 1.20, 3.35, 10.14 per 1,000 per year. 60+ = population-weighted mean of ages 60–74 (SingStat M810141, M810011); ages 75+ excluded (author decision) |
| ω | SMR of ever-PWID 14.68 (13.01–16.35) [8] | **Remove.** Replace with two rows below |
| `\sigma_n` (new) | — | All-cause SMR of people who use drugs without HCV or HIV infection, 4.73 (95% CI 3.47–6.45); van Santen et al. 2014, AIDS 28:2589, Table 3. Applies to u, r |
| `\mathrm{HR}_v` (new) | — | Hazard ratio of extrahepatic death, viraemic (anti-HCV+, HCV RNA detectable) vs anti-HCV−, 1.47 (1.23–1.77); Lee et al. 2012, J Infect Dis 206:469, Table 3. `\sigma_v = \sigma_n\mathrm{HR}_v = 6.95` applies to a, c, t |
| `\mu^{\mathrm{DC}}` | Additional disease mortality rate in DC, 0.13 | Keep value; interpretation "liver death rate in DC, added to the background rate"; note Lim 2018 reports annual probabilities, used here as rates |
| `\mu^{\mathrm{HCC}}` | 0.43 | Same as above |
| `\psi^{\mathrm{DC}}` | 0.45, U(0.31, 0.63) [13] (Premkumar 2024) | 0.40 (95% CI 0.17–0.96): adjusted hazard ratio of liver-related death, SVR vs no SVR, in patients with previous decompensated cirrhosis; Pageaux et al. 2022, BMC Infect Dis 22:94 (ANRS CO22 HEPATHER). Applies to cured people (r) only. Add the reference; drop [13] here (Premkumar 2024 reports no such ratio in its abstract) |
| `\psi^{\mathrm{HCC}}` | 0.37 (0.16–0.83) [2] | 0.34 (95% CI 0.13–0.88): adjusted hazard ratio of **liver-related** death, SVR vs untreated, HCV-related HCC; Dang et al. 2020, Hepatology 71:1910 (0.37 is their all-cause HR). Applies to r only |

Reference to add:
`Pageaux G-P, Nzinga CL, Ganne N, et al. Clinical outcomes after treatment with direct antiviral agents: beyond the virological response in patients with previous HCV-related decompensated cirrhosis. BMC Infect Dis. 2022;22:94. doi:10.1186/s12879-022-07076-0`
(also in `docs/references.bib` as `pageaux2022`; van Santen 2014 and Lee 2012 are there too).

## 5. Related passages outside 3–4 that now contradict the model (report, do not edit unless asked)

- §1.6 Mortality: "applied to the treatment compartment (t) only" (now r only), the use of ω, and
  "fitted to Singapore life-table data" (now fixed 2015 rates, 60+ = 60–74).
- §5 schematic: `\mu_i + \psi^{\mathrm{DC}}\mu^{\mathrm{DC}}` drawn on the treatment compartment.
- §5 schematic labels `D_{u,k,i}` in stages 2–4 "post-SVR"; `u` is susceptible (cured people are `r`).

## 6. Checks after editing

1. `grep -n '\\omega' docs/equations.tex` returns nothing; `grep -n 'SVR}' docs/equations.tex`
   returns nothing; `grep -n 'gamma_{i,j}' docs/equations.tex` returns nothing.
2. In every subsection, each of the 20 equations per stratum uses `\sigma_n` in u/r rows and `\sigma_v`
   in a/c/t rows; stage 3–4 t rows contain `\mu^{\mathrm{DC}}`/`\mu^{\mathrm{HCC}}` without `\psi`;
   stage 3–4 r rows contain `\psi^{\cdot}\mu^{\cdot}`.
3. Stage 3 rows use `\bar\alpha_3^{-}`, stage 4 rows `\bar\alpha_4`; no `\pm` remains.
4. D equations use `\tau_{D,k}` (stated as 0), J `\tau_{J,k}`, F `\tau_{F,k}`, X `\tau_{X,k}`.
5. Cross-check one equation per stratum against `rhs()` in `sim.cpp` (D/F loop lines 123–172, J
   174–212, X 214–252) and against the independent R implementation
   `tests/integration/test_independent_rhs.R`.
6. The document compiles; every new symbol (`\sigma_n`, `\sigma_v`, `\mathrm{HR}_v`, `c`, `C_{ij}`,
   `\theta_i`, `\tau_{S,k}`, `\bar\alpha_4`) appears in Section 4.
7. Report to the user: files changed, anything in this list you could not apply, and any further
   inconsistency you found.
