# =============================================================================
# si_text.py — S1 Text (technical appendix) for the PLOS Medicine article. Numbers come from
# Facts; parameter values are read from setup.R.
# =============================================================================
from fmt import F, pct, ui, col, red, money, n0, interval, cap_word, scen_label, WORDS

AGES = ["<20", "20–29", "30–39", "40–49", "50–59", "60+"]


def build(reg, f):
    from registry import Doc
    d = Doc(reg, "si")
    s = f.summ; caps = f.cfg["annual_screening_totals"]; shares = [round(100 * x) for x in f.cfg["f_allocation_fractions"]]
    st = f.setup
    d.raw("---\ntitle: \"S1 Text. Technical appendix\"\nsubtitle: \"Prison treatment and post-release screening for hepatitis C among people who inject drugs in Singapore: A modelling study\"\n---")
    d.raw("Tables and figures in this appendix are numbered Table A, Table B and so on, and Fig A, Fig B and so on.")

    # ---------------- 1 model ----------------
    d.h(1, "1. Model")
    d.h(2, "1.1 Compartments")
    d.p("The state $y_{s,k,h,i}$ counts people in group $s$ (D, J, F or X), liver stage $k$ (1 non-cirrhotic, 2 compensated cirrhosis, 3 decompensated cirrhosis, 4 hepatocellular carcinoma), HCV state $h$ (u susceptible, a acute, c chronic, t on treatment, r recovered) and age group $i$ (1–6). There are 480 compartments.", ("sim",), kind="method")
    d.p("Recovered people (r), after spontaneous clearance or treatment, remain antibody positive and can be reinfected. Susceptible people are in stage 1 only. Recovered people keep the liver stage they had when they cleared.", ("sim",), kind="method")
    d.h(2, "1.2 Transmission")
    d.p("The force of infection in age group $i$ is $\\gamma_i = q \\sum_j \\tilde C_{ij} \\, \\dfrac{a_j + c_j}{N_j}$, where $a_j$ and $c_j$ are the acute and chronic people and $N_j$ is everyone in D and F in age group $j$, and $\\tilde C_{ij} = C_{ij} e^{\\theta_i}$ is the base contact matrix with row $i$ multiplied by the calibrated factor $e^{\\theta_i}$. Susceptible and recovered people in D and F become acute at rate $\\gamma_i$. People in J and X are not infected.",
        ("sim", "likelihood"), kind="method")
    d.h(2, "1.3 Movement between groups")
    d.p("New injectors enter D (state u, stage 1) at rate $s_i \\beta_i$ per year, where $\\beta_i$ is the 2015 Central Narcotics Bureau count of new drug abusers in age group $i$ and $s_i \\ge 1$ is a calibrated entry scale. Entry is the same in every year.", ("sim", "setup", "likelihood"), kind="method")
    d.p("People leave D for J at the first-arrest rate $\\lambda^{(1)}_i = c\\,\\lambda^{(3)}_i$ and leave F for J at the re-arrest rate $\\lambda^{(3)}_i$. People leave J at the release rate $\\lambda^{(2)}_i$; a share $\\pi$ goes to F and $1-\\pi$ to X. People in X do not return to J. Every group ages from age group $i$ to $i+1$ at rate 1/10 per year.",
        ("sim", "setup"), kind="method")
    d.h(2, "1.4 Natural history and treatment")
    d.p("Acute infection lasts on average $\\iota_1$; a share $\\kappa$ clears to r and the rest become chronic. Chronic people progress from stage 1 to 2 and from stage 2 to 3 or 4, and from 3 to 4, at rates $p^{*} = (\\rho\\, r_3 + 1 - \\rho)\\,p$ that weight the genotype 3 relative risk $r_3$ by the genotype 3 share $\\rho$. Recovered people in stages 2 and 3 progress at $\\phi\\, p^{*}$.",
        ("sim", "setup"), kind="method")
    d.p("Chronic people start treatment at rate $\\tau_{s,k}$ (zero in D and before 2019). Treatment lasts on average $\\iota_2$; a share $\\alpha_k$ (sustained virological response) moves to r and the rest return to c.", ("sim",), kind="method")
    d.h(2, "1.5 Mortality")
    d.p("Everyone dies at the background rate $\\mu_i\\,\\sigma_h$, where $\\mu_i$ is the 2015 Singapore death rate in age group $i$ and $\\sigma_h$ is the standardised mortality ratio of people who use drugs: $\\sigma_v$ for viraemic states (a, c, t) and $\\sigma_n$ for others (u, r). People in stages 3 and 4 also die of liver disease at $\\mu_{DC}$ or $\\mu_{HCC}$, multiplied by $\\psi_{DC}$ or $\\psi_{HCC}$ after treatment (states t and r).",
        ("sim", "setup"), kind="method")
    d.h(2, "1.6 Numerical solution and equilibrium")
    d.p("The equations were integrated with a fourth-order Runge–Kutta method and a step of one day. Negative values were set to zero after each step, and ageing was applied after each step. Each run started from a seed population in D and ran for 150 model years; model year 150 is 2015.", ("sim", "setup", "targets"), kind="method")
    d.p("A run was accepted as being at equilibrium in 2015 if, between model years 140 and 150, the prison population of every age group changed by at most 1% on the log scale, prison seroprevalence by at most 0.5 percentage points, total HCV, decompensated cirrhosis and carcinoma counts by at most 1%, and every compartment by at most 2%. Runs that failed had zero posterior density.",
        ("equilibrium", "likelihood"), kind="method")

    # ---------------- 2 parameters ----------------
    d.h(1, "2. Fixed parameters")
    d.p("Table A lists the fixed parameters. They are the values used by the model code.", ("setup",), kind="method")
    d.raw("**Table A. Fixed parameters.**")
    mu = st("mu"); lam2 = st("lambda2"); lam3 = st("lambda3")
    rows = [  # (parameter, value, citation keys, note)
        ["Transmission probability per sharing event, $q$", f"{st('q')}", ["boelen2014"], ""],
        ["Spontaneous clearance, $\\kappa$", f"{st('kappa')}", ["micallef2006"], ""],
        ["Mean acute duration, $\\iota_1$", "180 days", ["fasano2024"], ""],
        ["Mean treatment duration, $\\iota_2$", "84 days (12 weeks)", [], "Regimen length"],
        ["Genotype 3 share, $\\rho$", f"{st('rho')}", ["zhao2019"], ""],
        ["Sustained response, stages 1–2", f"{st('alpha_NC')}", ["pawlotsky2020"], ""],
        ["Sustained response, stage 3, without ribavirin (scenarios)", f"{st('alpha_DC_neg')}", ["xiao2025"], ""],
        ["Sustained response, stage 4", f"{st('alpha_HCC')}", [], "Assumed equal to stage 3 with ribavirin"],
        ["Progression, stage 1 to 2 (per year)", f"{st('p_NC_CC')}", ["thein2008"], ""],
        ["Progression, stage 2 to 3; 2 to 4; 3 to 4 (per year)", f"{st('p_CC_DC')}; {st('p_CC_HCC')}; {st('p_DC_HCC')}", ["lim2018"], ""],
        ["Genotype 3 relative risk, to cirrhosis; to carcinoma", f"{st('r3_NC_CC')}; {st('r3_CC_HCC')}", ["kanwal2014"], ""],
        ["Relative progression after cure, 2 to 3; 2 to 4; 3 to 4", f"{st('phi_CC_DC')}; {st('phi_CC_HCC')}; {st('phi_DC_HCC')}", ["vandermeer2012", "morgan2013"], "last value assumed"],
        ["Liver mortality, stage 3; stage 4 (per year)", f"{st('mu_DC')}; {st('mu_HCC')}", ["lim2018"], ""],
        ["Relative liver mortality after cure, stage 3; stage 4", f"{st('psi_DC')}; {st('psi_HCC')}", ["premkumar2024", "dang2020"], ""],
        ["Mortality ratio, viraemic; non-viraemic", f"{st('smr_viraemic')}; {st('smr_nonviraemic')}", ["vansanten2014"], ""],
        ["Background death rate by age (per 1,000 per year)", "; ".join(f"{1000 * m:.2f}" for m in mu), [], "Singapore 2015 death rates (source to be cited)"],
        ["Release rate by age (per year)", "; ".join(f"{x:.3f}" for x in lam2), [], "Prison records (Section 3)"],
        ["Re-arrest rate by age (per year)", "; ".join(f"{x:.3f}" for x in lam3), [], "Prison records (Section 3)"],
        ["Probability of eventual re-arrest, $\\pi$", f"{st('pi_recid')}", [], "Prison records (Section 3)"],
        ["Entry into D by age, $\\beta_i$ (per year)", "; ".join(str(x) for x in f.cnb_counts), ["singstat_drug"], ""],
    ]
    d.raw("| Parameter | Value | Source |")
    d.raw("|---|---|---|")
    for name, val, keys, note in rows:
        src = "; ".join(x for x in ([f"[{'; '.join('@' + k for k in keys)}]"] if keys else []) + ([note] if note else []))
        d.p(f"| {name} | {val} | {src} |", ("setup",) + tuple(keys), kind="method")
    d.raw("The base contact matrix $C$ is given in the model code; its source is to be documented by the authors.")

    # ---------------- 3 prison rates ----------------
    d.h(1, "3. Prison rate estimation")
    d.p(f"Singapore Prison Service records were restricted to people ever held for drug consumption. Consecutive stays with no gap were merged. This gave {n0(f.sps_people)} people and {n0(f.sps_stays)} stays with admissions from {f.sps_first} to {f.sps_last}; follow-up ended on 31 December 2017.",
        ("sps", "preprocess"), kind="method")
    d.p("Release: each stay was split at the boundaries of the 10-year age groups, and a constant release hazard per age group was estimated, with a random-walk prior on the log hazards across ages.", ("stan_release", "preprocess"), kind="method")
    d.p("Re-arrest: each interval from release to the next admission, or to the end of follow-up, was split by age in the same way. A cure-type survival model estimated the probability $\\pi$ of ever being re-arrested and, for people who would be re-arrested, an age-specific hazard.", ("stan_rearrest", "preprocess"), kind="method")
    d.p(f"Both models were fitted in Stan with four chains of 4,000 sampling iterations; the largest R-hat was {max(f.rate_rhat):.3f}. Table B gives the posterior means used in the model and their 95% credible intervals.", ("rates", "release", "preprocess"), kind="method")
    d.p("**Table B. Prison rates by age (per year; posterior mean and 95% credible interval).**", ("rates",), kind="result")
    d.p("| Age group | Release, $\\lambda^{(2)}$ | Re-arrest, $\\lambda^{(3)}$ |", ("rates",), kind="result")
    d.raw("|---|---|---|")
    for i in range(1, 7):
        a, b = f.rate_ci[("lambda2", i)], f.rate_ci[("lambda3", i)]
        d.p(f"| {AGES[i - 1]} | {a[0]:.3f} ({a[1]:.3f}–{a[2]:.3f}) | {b[0]:.3f} ({b[1]:.3f}–{b[2]:.3f}) |", ("rates", "release"))
    pv = f.rate_ci[("pi", 1)]
    d.p(f"The probability of eventual re-arrest was {pv[0]:.3f} (95% credible interval {pv[1]:.3f}–{pv[2]:.3f}). In the model, people who will be re-arrested are assumed to return to injecting (F).", ("rates", "setup"))
    d.p("These rates come from everyone held for drug consumption, not only people who inject.", ("preprocess",), kind="method")

    # ---------------- 4 calibration ----------------
    d.h(1, "4. Calibration")
    d.h(2, "4.1 Targets")
    d.p("**Table C. Calibration targets, 2015.**", ("targets",), kind="result")
    d.p("| Target | <20 | 20–29 | 30–39 | 40–49 | 50–59 | 60+ |", ("targets",), kind="result")
    d.raw("|---|---|---|---|---|---|---|")
    d.p("| HCV antibody prevalence in prison | " + " | ".join(f"{100 * x:.1f}%" for x in f.tg["seroprevalence"]) + " |", ("targets",))
    d.p("| People screened in prison (end of 2015) | " + " | ".join(x["screened"] for x in f.prison) + " |", ("prison",))
    d.p("| First arrests (new drug abusers) | " + " | ".join(str(x) for x in f.cnb_counts) + " |", ("cnb", "targets"))
    d.p(f"Community injectors (D + F, all ages): {f.tg['community_pwid']:,} (95% confidence interval 6,500–17,000) [@teo2019].", ("targets",))
    d.p(f"The seroprevalence targets are stored in the targets file. Recomputing them from the current processing of the prison records (positive divided by screened, Table C) gives values that differ by up to {100 * f.sero_max_diff:.1f} percentage points; the reason is not recorded and should be resolved by the authors.",
        ("targets", "prison"), kind="result", unc="moderate")
    d.p(f"The community estimate refers to {f.tg['community_pwid_year']}; because entry is constant and treatment starts in 2019, the modelled 2015 value is used for that year.", ("targets", "sim"), kind="method")
    d.h(2, "4.2 Parameters and prior")
    d.raw("**Table D. Calibrated parameters and prior.**")
    d.p("| Parameter | Meaning | Prior on $\\theta$ | Bounds of the scale |", ("likelihood",), kind="result")
    d.raw("|---|---|---|---|")
    d.p("| $\\theta_1$–$\\theta_6$ | Log contact multiplier, one per age group | Normal(0, 1.5) | 0.01–1,000 |", ("likelihood",), kind="method")
    d.p("| $\\theta_7$–$\\theta_{12}$ | Log entry scale $s_i$, one per age group | Half-normal(0, 1) | 1–100 |", ("likelihood",), kind="method")
    d.p("| $\\theta_{13}$ | Log first-arrest factor $c$ | Normal(0, 0.5) | 0.5–100 |", ("likelihood",), kind="method")
    d.p("The prior was fixed before fitting and did not use the data. The same prior was used for the mode search, the MCMC and the neural posterior.", ("likelihood",), kind="method")
    d.h(2, "4.3 Likelihood")
    d.p(f"Seroprevalence: each observed value $y_i$ follows Beta($p_i\\phi$, $(1-p_i)\\phi$) around the modelled value $p_i$, with $\\phi$ = {f.tg['beta_precision']}. "
        f"First arrests: each count follows a negative binomial distribution with mean equal to the modelled flow $c\\,\\lambda^{{(3)}}_i D_i$ and size {f.tg['first_arrests_nb_size']}. "
        f"Community injectors: the log estimate follows Normal(log($D + F$), {f.tg['community_pwid_log_sd']}).", ("likelihood", "targets", "sim"), kind="method")
    d.h(2, "4.4 Posterior mode")
    best = f.nm_best; same = f.nm_same
    d.p(f"A Nelder–Mead search started from $\\theta$ = 0 and {len(f.nm) - 1} random perturbations of it (standard deviation {f.nm_cfg['sd_perturb']} on the log scale), each with up to {f.nm_cfg['maxit']:,} iterations and up to {f.nm_cfg['nm_restarts']} restarts. "
        f"{cap_word(len(same))} of the {len(f.nm)} starts reached the same mode (negative log posterior {min(F(x['objective']) for x in same):.2f} to {max(F(x['objective']) for x in same):.2f}); the other stopped at {F(f.nm_other[0]['objective']):.1f}.",
        ("nm", "nm_cfg"), unc="low")
    d.p(f"At the mode, the largest absolute seroprevalence error was {100 * F(best['prevalence_max_abs_err']):.1f} percentage points, the first-arrest factor was {F(best['c_first']):.2f} and the community population was {n0(best['community_pwid'])}.", ("nm",), unc="low")
    d.h(2, "4.5 MCMC")
    acc = f.acc
    d.p(f"We used adaptive random-walk Metropolis [@haario2001] with {f.mcfg['n_chains']} chains started from independent prior draws. The initial proposal covariance was {f.mcfg['init_scale'] ** 2:g} times the prior variance. Adaptation began after {f.mcfg['adapt_start']:,} proposals, used the recent half of the warm-up, was refreshed every {f.mcfg['adapt_every']} proposals and was frozen after {f.mcfg['n_warmup']:,} warm-up proposals; {pct(f.mcfg['safety_p'], 0)} of proposals always used the initial covariance.",
        ("mcmc_cfg",), kind="method")
    d.p(f"Each chain ran {f.mcfg['max_evals_per_chain']:,} model evaluations, and every {f.mcfg['thin_storage']}th retained iteration was stored. The run took {F(f.mcfg['method_sec']) / 3600:.1f} hours. Acceptance rates were {min(acc):.2f} to {max(acc):.2f}, below the 0.15 lower limit used to flag chains, so the chains mixed slowly but converged.",
        ("mcmc_cfg", "mcmc_chain"), unc="low")
    d.p(f"All diagnostics used every retained iteration: rank-normalised split R-hat was at most {f.rhat_max:.3f}, bulk effective sample size at least {n0(f.ess_bulk_min)}, tail effective sample size at least {n0(f.ess_tail_min)}, and the Monte Carlo standard error of each mean at most {f.mcse_max:.3f} posterior standard deviations [@vehtari2021] (Table E; Fig A).",
        ("mcmc_diag",), unc="low")
    d.raw("**Table E. Posterior summaries and convergence diagnostics.**")
    d.p("| Parameter | Median (95% credible interval), natural scale | Posterior SD / prior SD | R-hat | Bulk ESS | Tail ESS |", ("mcmc_par",), kind="result")
    d.raw("|---|---|---|---|---|---|")
    diag = {x["parameter"]: x for x in f.mdiag}
    names = [f"Contact, {a}" for a in AGES] + [f"Entry scale, {a}" for a in AGES] + ["First-arrest factor $c$"]
    for i in range(1, 14):
        p = f.mpar[f"theta{i}"]; g = diag[f"theta{i}"]
        import math
        m, lo, hi = (math.exp(F(p[k])) for k in ("median", "q025", "q975"))
        g3 = lambda v: f"{v:.3g}" if v < 0.1 else f"{v:.2f}"
        d.p(f"| {names[i - 1]} | {g3(m)} ({g3(lo)}–{g3(hi)}) | {f.sd_ratio[i - 1]:.2f} | {F(g['rhat']):.3f} | {n0(g['ess_bulk'])} | {n0(g['ess_tail'])} |", ("mcmc_par", "mcmc_diag", "mcmc_cfg"), unc="low")
    d.p("SD, standard deviation; ESS, effective sample size. The prior SD is that of the truncated prior.", ("mcmc_cfg",), kind="method")
    pdm = ", ".join(f"$\\theta_{{{i}}}$" for i in f.prior_dominated)
    d.p(f"For {pdm} the posterior standard deviation was more than 0.7 times the prior's, so these parameters were mainly informed by the prior: the contact multipliers for ages 50–59 and 60 or older, and the first-arrest factor. The posterior median of c ({f.c_post[0]:.2f}) was nevertheless well above the prior median of 1.",
        ("mcmc_par", "mcmc_cfg"), unc="moderate")
    d.p(f"The entry scales were close to their lower bound of 1 for ages under 60 (medians {min(v[0] for k, v in f.entry_post.items() if k < 12):.2f} to {max(v[0] for k, v in f.entry_post.items() if k < 12):.2f}); for ages 60 and over the median was {f.entry_post[12][0]:.2f} (95% credible interval {f.entry_post[12][1]:.2f}–{f.entry_post[12][2]:.2f}).", ("mcmc_par",), unc="moderate")
    d.h(2, "4.6 Fit to the targets")
    d.p(f"For each of {F(f.strata['population_D']['n_draws_used']):,.0f} posterior draws we computed the modelled 2015 targets. The 95% posterior intervals of the model values (without observation noise) overlapped the target ranges for all {f.n_ovl} targets (main text Fig 2).", ("strata", "overlap"), unc="moderate")
    d.p("**Table F. Targets against the 95% posterior interval of each modelled value.**", ("overlap",), kind="result")
    d.p("| Target | Age group | Target (95% range) | Model, 95% posterior interval |", ("overlap",), kind="result")
    d.raw("|---|---|---|---|")
    for o in f.ovl:
        sc = 100 if o["measure"] == "seroprevalence" else 1
        dd = 1 if sc == 100 else 0
        u = "%" if sc == 100 else ""
        meas = {"seroprevalence": "Seroprevalence", "first_arrests": "First arrests", "community_pwid": "Community injectors"}[o["measure"]]
        d.p(f"| {meas} | {o['age_group'].replace('-', '–')} | {F(o['target']) * sc:,.{dd}f}{u} ({F(o['target_lower']) * sc:,.{dd}f}–{F(o['target_upper']) * sc:,.{dd}f}) | {F(o['MCMC_lower']) * sc:,.{dd}f}–{F(o['MCMC_upper']) * sc:,.{dd}f} |", ("overlap",), unc="moderate")
    d.h(2, "4.7 Checks of the posterior")
    m = f.msbc["all"]
    d.p(f"Direct simulation-based calibration [@talts2018]: we reran the MCMC (one chain of 45,000 evaluations, 15,000 warm-up) on {m['n_datasets']} data sets simulated from prior draws. Pooled over {13 * int(m['n_datasets'])} ranks, 95% interval coverage was {F(m['coverage_95']):.2f}, 90% coverage {F(m['coverage_90']):.2f} and the mean rank {F(m['mean_rank']):.2f}. With {m['n_datasets']} data sets this guards only against gross miscalibration.",
        ("mcmc_sbc",), unc="moderate")
    ns = f.npe_sim; acc_ = f.npe_acc["npe_lhs"]
    sbc = f.npe_sbc
    d.p(f"Neural posterior estimation [@greenberg2019; @sbi2020]: we simulated {ns['n_attempted']:,} parameter sets from the prior by Latin hypercube sampling ({ns['n_valid']:,} valid), added observation noise from the likelihood, and trained an ensemble of {f.npe_cfg['n_ensemble']} masked autoregressive flows [@papamakarios2017].",
        ("npe_sim", "npe_cfg"), kind="method")
    bad = ", ".join(f"$\\theta_{{{x[5:]}}}$" for x in f.npe_bad_shape)
    d.p(f"Its simulation-based calibration on {int(F(sbc[0]['n_sbc'])):,} data sets gave 95% coverage of {min(F(x['coverage_95']) for x in sbc):.3f} to {max(F(x['coverage_95']) for x in sbc):.3f}, but the rank histograms of {bad} were not uniform, so the neural posterior failed its predeclared quality check (Fig B).",
        ("npe_sbc",), unc="moderate")
    d.p(f"Of {f.npe_checked:,} neural posterior draws, {pct(1 - f.npe_pass_eq / f.npe_checked, 0)} failed the equilibrium check. Compared with the MCMC, the neural posterior intervals were on average {F(acc_['mean_width95_ratio']):.2f} times as wide, and {pct(acc_['share_pass_all'], 0)} of parameters met all agreement criteria (Fig C). We therefore used the MCMC posterior for every result and report the neural posterior only as a failed cross-check.",
        ("npe_support", "npe_acc"), unc="moderate")

    # ---------------- 5 population ----------------
    d.h(1, "5. Model population in 2015")
    d.p("**Table G. Modelled population in 2015 (median and 95% credible interval across posterior draws).**", ("strata",), kind="result")
    d.raw("| Group | All people | Viraemic (a, c, t) |")
    d.raw("|---|---|---|")
    lab = {"D": "D, never imprisoned", "J": "J, in prison", "F": "F, released, injecting", "X": "X, released, stopped",
           "current_PWID": "D + F (active injectors outside prison)", "ever_PWID": "All groups"}
    for g in ["D", "J", "F", "X", "current_PWID", "ever_PWID"]:
        a, v = f.strata[f"population_{g}"], f.strata[f"active_hcv_{g}"]
        d.p(f"| {lab[g]} | {interval(a['median'], a['lower_95'], a['upper_95'])} | {interval(v['median'], v['lower_95'], v['upper_95'])} |", ("strata",), unc="moderate")
    d.p(f"The modelled prison population was not a calibration target. In the prison records, {f.j_obs_2015:,} people held for drug consumption with a screening result were in prison at the end of 2015.", ("strata", "prison"), unc="moderate")

    # ---------------- 6 scenarios ----------------
    d.h(1, "6. Scenario implementation")
    d.p(f"Treatment in J from 2019 used the rate $\\tau_J = -\\ln(1 - {f.cfg['j_annual_start_probability']}) = {f.tau_j:.2f}$ per year for people with chronic infection in every stage and age. From 2026, people in F and X started at $\\tau_{{FX}} = -\\ln(1 - {f.cfg['fx_annual_start_probability']}) = {f.tau_fx:.2f}$ per year, but only until that year's capacity was used.",
        ("config", "scen_code"), kind="method")
    d.p("At the start of each year, for each of F and X, $N = u + a + c + r$ people were eligible. The number screened was the smaller of the quota and $N$. Antibody tests equalled the number screened; RNA tests equalled screened times $(a + c + r)/N$; chronic infections found, and so the year's treatment capacity, equalled screened times $c/N$.",
        ("scen_code",), kind="method")
    d.p("The solver located the moment in each year at which a group used its capacity and then stopped new starts in that group for the rest of the year. Unused capacity was not moved between F and X or carried over. The same person can be counted as screened in more than one year, and the model counts treatment courses, not people.", ("sim", "scen_code"), kind="method")
    d.p("In J, every entrant not on treatment (u, a, c, r) was counted as having an antibody test and every antibody-positive entrant (a, c, r) an RNA test; each treatment start was costed as one course and one consultation.", ("scen_code",), kind="method")
    d.p("New infections in a year were the infections (u and r to a, in D and F) during that calendar year, and HCV deaths the liver deaths in stages 3 and 4 plus $(\\sigma_v - \\sigma_n)\\mu_i$ times the number of viraemic people, summed over the year (trapezoid rule over the daily solver output). The 2015 values were one year at the calibrated equilibrium; the 2016 values, still without treatment, differed from them by less than 0.1%. Both were compared with their 2015 values.", ("scen_code", "sim"), kind="method")
    d.p(f"For each of {f.ndraw:,} posterior draws (sampling seed fixed), we ran the 2015 equilibrium, the shared 2016–2025 history and all {len(s)} scenarios to 2030; {int(f.valid['draws_valid']):,} of {int(f.valid['draws_requested']):,} draws were valid.", ("summary", "draws_valid", "scen_code"))

    # ---------------- 7 costs ----------------
    d.h(1, "7. Costs")
    d.p("**Table H. Unit costs (2026 Singapore dollars).**", ("costev",), kind="result")
    d.raw("| Item | Base | Range in sensitivity analysis | Basis |")
    d.raw("|---|---|---|---|")
    basis = {"SOF/VEL 400/100 mg": "Published price of US$2,577 per week for 12 weeks [@wong2019], converted at S$1.3066 per US$ [@sgd_rate]; not a procurement price",
             "Anti-HCV antibody": "Published clinic charge [@nhg_dsc]", "HCV RNA PCR": "Published laboratory charge [@parkway_rna]",
             "Doctor consultation": "Published clinic charge [@nhg_dsc]"}
    for x in f.costev:
        d.p(f"| {x['item']} | S${F(x['base_cost_sgd']):,.2f} | S${F(x['low_sgd']):,.2f}–S${F(x['high_sgd']):,.2f} | {basis[x['item']]} |", ("costev",), kind="method")
    d.p("The source abstract for the drug price gives US$32,325 per 12 weeks, which differs from US$2,577 times 12 (US$30,924); we used the weekly figure. Prices were not inflated to 2026. Follow-up visits, tests of cure, care navigation and genotyping were not costed.", ("costev",), kind="method")
    d.p(f"Costs were counted in the year they occurred and discounted at {pct(f.cfg['discount_rate'], 0)} a year to 2026; numbers of tests and courses were not discounted. Costs before 2026 were not included.", ("config",), kind="method")

    # ---------------- 8 additional results ----------------
    d.h(1, "8. Additional results")
    d.p(f"Table I gives every scenario. Fig D shows which policies met each target at the median, Fig E the costs and treatment starts by group, and Fig F the yearly screening cascade in F and X.", ("summary",))
    d.p("**Table I. All scenarios, 2030 and 2026–2030 totals (median and 95% uncertainty interval).**", ("summary",), kind="result")
    d.p("| Policy | Draws meeting incidence target | Draws meeting mortality target | Treatment starts in F and X, 2026–2030 | Extra cost over J only, S$ million |", ("summary",), kind="result")
    d.raw("|---|---|---|---|---|")
    for k in ["J_only"] + [f"S{c}_F{x}" for c in caps for x in shares] + ["SInf"]:
        r = s[k]
        d.p(f"| {scen_label(r)} | {pct(r['who_incidence_pass_prob'], 1)} | {pct(r['who_mortality_pass_prob'], 1)} | {interval(r['fx_treatment_starts'], r['fx_treatment_starts_lo'], r['fx_treatment_starts_hi'])} | "
            f"{interval(r['incremental_discounted_cost_vs_j_only_sgd'], r['incremental_discounted_cost_vs_j_only_sgd_lo'], r['incremental_discounted_cost_vs_j_only_sgd_hi'], d=1, scale=1e-6)} |", ("summary",), unc="moderate")
    p5 = s["S5000_F25"]
    d.p(f"Cost sensitivity: with all unit prices halved or increased by half, the discounted cost of 5,000 tests a year with 25% to F was S${F(p5['cost_x0_5_sgd']) / 1e6:.1f} million or S${F(p5['cost_x1_5_sgd']) / 1e6:.1f} million; epidemiological results do not change.", ("summary", "config"), unc="moderate")

    # ---------------- 9 software ----------------
    d.h(1, "9. Software and reproducibility")
    d.p("The model is in C++ (scripts/model/sim.cpp), called from R through Rcpp [@rcpp2011]. Rate models used Stan, convergence diagnostics the R package posterior, and the neural posterior the sbi package [@sbi2020]. Figures were drawn with matplotlib from the saved outputs.", ("sim", "figures"), kind="method")
    d.p("Analysis commands, run times and expected outputs are listed in docs/training.md in the repository. This appendix and the main text were generated from the saved outputs by scripts/build_paper.py; drafting followed the procedures of Scientific Agent Skills [@kassis2026].", ("figures",), kind="declaration")

    # ---------------- figures ----------------
    d.h(1, "Figures")
    d.p(f"**Fig A. MCMC trace plots.** Stored draws of the 13 parameters for each of the {f.mcfg['n_chains']} chains after warm-up.", ("mcmc_cfg",), kind="result")
    d.raw(f"![]({dict(SI_FIGURES)['Fig A']}){{width=95%}}")
    d.p(f"**Fig B. Rank histograms of the neural posterior's simulation-based calibration.** {int(F(sbc[0]['n_sbc'])):,} data sets, 20 bins; the band is the 99% range under uniformity.", ("npe_sbc",), kind="result")
    d.raw(f"![]({dict(SI_FIGURES)['Fig B']}){{width=95%}}")
    d.raw("**Fig C. Marginal posteriors from the MCMC (grey) and the neural posterior (orange).**")
    d.raw(f"![]({dict(SI_FIGURES)['Fig C']}){{width=95%}}")
    d.p("**Fig D. Target attainment in 2030 by policy.** Median reductions in new infections and HCV deaths, and whether both targets were met.", ("summary",), kind="result")
    d.raw(f"![]({dict(SI_FIGURES)['Fig D']}){{width=95%}}")
    d.p("**Fig E. Costs and treatment starts, 2026–2030.** (A) Discounted costs by group. (B) Undiscounted treatment courses by group. Bars are medians of each group; whiskers are 95% uncertainty intervals of the total.", ("stratum_tot","summary"), kind="result")
    d.raw(f"![]({dict(SI_FIGURES)['Fig E']}){{width=95%}}")
    d.raw("**Fig F. Screening cascade in F and X by year.** Number screened, RNA tests, RNA-positive people, chronic infections found and treatment starts, by capacity (columns) and share of tests given to F (rows).")
    d.raw(f"![]({dict(SI_FIGURES)['Fig F']}){{width=95%}}")
    d.h(1, "References")
    d.raw("::: {#refs}\n:::")
    return d


SI_FIGURES = [  # (label, repository path)
    ("Fig A", "output/calibration/mcmc_reference/mcmc_traceplots.png"),
    ("Fig B", "output/calibration/npe_lhs/sbc_rank_histograms.png"),
    ("Fig C", "output/calibration/npe_vs_mcmc/figures/marginals.png"),
    ("Fig D", "output/scenarios_2026/scenario_target_heatmap.png"),
    ("Fig E", "output/scenarios_2026/scenario_cumulative_costs.png"),
    ("Fig F", "output/scenarios_2026/scenario_care_cascade.png"),
]
