# =============================================================================
# main_text.py — the PLOS Medicine research article (main text). Every number comes from
# Facts; every factual line is bound to its evidence through Doc.p().
# =============================================================================
from fmt import F, pct, ui, col, red, money, n0, interval, cap_word, scen_label, WORDS

OUT = ("summary",)          # scenario summary file
TRJ = ("traj",)
CAL = ("mcmc_par", "mcmc_diag", "mcmc_cfg")


def checks(f):
    """The narrative below is written for these results; stop the build if they change."""
    s = f.summ
    assert len(f.inc_pass) == 15, "not every finite policy meets the incidence target"
    assert set(f.mort_pass) == {"S5000_F25", "S8000_F25", "S8000_F50"}, "mortality-passing set changed"
    assert s["J_only"]["who_incidence_pass"] == "TRUE" and s["J_only"]["who_mortality_pass"] != "TRUE"
    assert s["SInf"]["who_joint_pass"] == "TRUE"
    assert f.inc_order_ok and f.mort_order_ok, "F/X ordering of the reductions changed"
    assert f.x_chronic_share > 0.85 and f.fx_death_share_2026 > 0.95, "X / F+X shares changed"
    assert f.x_chronic / f.x_eligible > f.f_chronic / f.f_eligible, "X no longer has the higher chronic share"
    assert all(abs(v[st] - 1) < 1e-6 for v in f.quota_use.values() for st in "FX"), "a quota is not fully used"
    assert min(f.mort_pass, key=lambda k: F(s[k]["total_discounted_cost_2026_2030_sgd"])) == "S5000_F25"
    assert F(s["S5000_F50"]["joint_pass_probability"]) < 0.2
    assert f.rhat_max < 1.01 and f.n_ovl == 13
    assert F(f.strata["population_J"]["median"]) > 2 * f.j_obs_2015, "modelled prison size no longer far above observed"
    assert f.prior_dominated == [5, 6, 13]
    assert all(abs(F(f.tr("no_treatment", 2030, f.notx)[c]) - 1) < 0.01 for c in ("incidence_ratio_2015", "mortality_ratio_2015"))


def build(reg, f):
    from registry import Doc
    checks(f)
    d = Doc(reg, "main")
    s = f.summ; J, P5, P8, P85, P550, INF = (s[k] for k in ("J_only", "S5000_F25", "S8000_F25", "S8000_F50", "S5000_F50", "SInf"))
    caps = f.cfg["annual_screening_totals"]; shares = [round(100 * x) for x in f.cfg["f_allocation_fractions"]]
    cnb_total = sum(f.cnb_counts)
    jpred = f.strata["population_J"]; df = f.strata["population_current_PWID"]; ever = f.strata["population_ever_PWID"]
    xs, fx26 = pct(f.x_chronic_share, 0), pct(f.fx_death_share_2026, 0)
    fchr, xchr = pct(f.f_chronic / f.f_eligible, 1), pct(f.x_chronic / f.x_eligible, 1)
    b15 = f.tr("J_only", 2015)

    # ---------------- front matter ----------------
    d.raw("---\ntitle: \"Prison treatment and post-release screening for hepatitis C among people who inject drugs in Singapore: A modelling study\"\n---")
    d.raw("**Short title:** Hepatitis C treatment after prison release in Singapore")
    d.raw("Yichen He\\*")
    d.raw("[Affiliation: to be provided by the authors]")
    d.raw("\\* Corresponding author. E-mail: [to be provided by the authors]")

    # ---------------- abstract ----------------
    d.h(1, "Abstract")
    d.h(2, "Background")
    d.p("Hepatitis C virus (HCV) infection in Singapore is concentrated among people who inject drugs (PWID), many of whom cycle through prison, where treatment expanded from 2019. "
        "We asked whether prison treatment alone can meet the World Health Organization 2030 targets (80% fewer new infections and 65% fewer HCV deaths than in 2015), and how limited community screening should be shared after release.",
        ("config", "scen_code", "who2022strategy"), kind="factual")
    d.h(2, "Methods and findings")
    d.p(f"We built an age-structured model of PWID never imprisoned (D), in prison (J), and released and still injecting (F) or stopped (X), with prison rates estimated from records of {n0(f.sps_people)} people. "
        f"We calibrated it by Markov chain Monte Carlo to 2015 prison seroprevalence, first arrests and a survey estimate of active injectors, and projected {len(f.finite)} policies screening {caps[0]:,} to {caps[-1]:,} released people a year from 2026 ({shares[0]}%–{shares[2]}% of tests to F) over {f.ndraw:,} posterior draws.",
        ("setup", "cnb", "targets", "mcmc_cfg", "config", "summary", "rates", "sps"), kind="method")
    d.p(f"Prison treatment alone reduced new infections in 2030 by {red(J, 'incidence_reduction', first=True)} but HCV deaths by only {red(J, 'mortality_reduction')}. "
        f"In 2026, {xs} of chronic infections among released people were in X, and F and X accounted for {fx26} of HCV deaths. "
        f"All screening policies met the incidence target; the mortality target was met only with {int(F(P5['screening_total'])):,} tests a year and 25% to F ({red(P5, 'mortality_reduction').replace(' (95% UI ', ', 95% UI ').rstrip(')')}; both targets in {pct(P5['joint_pass_probability'], 0)} of draws) or {caps[-1]:,} tests with 25% or 50% to F. "
        f"Moving tests from F to X increased the fall in deaths and reduced the fall in infections. "
        f"The cheapest passing policy cost {money(P5, 'total_discounted_cost_2026_2030_sgd')} over 2026–2030, mostly for drugs. "
        f"Limitations include a modelled prison population ({n0(jpred['median'])}) far above the {f.j_obs_2015:,} screened people recorded, weakly identified entry parameters, charge-based costs and parameter-only uncertainty.",
        ("summary", "annual", "traj", "strata", "prison", "mcmc_par", "costev"), unc="moderate", intent="exploratory")
    d.h(2, "Conclusions")
    d.p("In this model, prison treatment is likely to meet the 2030 target for new HCV infections among PWID in Singapore but not the target for deaths, which required screening several thousand released people a year, mostly people who have stopped injecting.",
        ("summary",), kind="interpretive", unc="moderate")

    # ---------------- author summary ----------------
    d.h(1, "Author summary")
    d.h(2, "Why was this study done?")
    d.p("- Most hepatitis C infections in Singapore are among people who inject drugs, and many of them spend time in prison, where treatment expanded from 2019.", ("config", "winslow2007", "prison2021"), kind="factual")
    d.p("- It is unclear whether prison treatment alone can meet the WHO 2030 elimination targets, which ask for large falls in both new infections and deaths.", ("who2022strategy",), kind="factual")
    d.raw("- Community screening capacity is limited, and it is not known whether it should go to released people who still inject or to those who have stopped.")
    d.h(2, "What did the researchers do and find?")
    d.p(f"- We built a model of hepatitis C among people who inject drugs that follows them into and out of prison, and fitted it to prison and national data from 2015.", ("targets",), kind="method")
    d.p(f"- In the model, prison treatment alone cut new infections by about {round(100 * F(J['incidence_reduction']))}% by 2030 but deaths by only about {round(100 * F(J['mortality_reduction']))}%, because most deaths occur among released people who have stopped injecting, whom prison treatment does not reach again.",
        ("summary", "traj"), unc="moderate")
    d.p(f"- The death target was met only when {int(F(P5['screening_total'])):,} or more released people were screened each year, with most tests going to people who have stopped injecting.", ("summary",), unc="moderate")
    d.h(2, "What do these findings mean?")
    d.raw("- Prison treatment can prevent most new infections, but it cannot by itself prevent most hepatitis C deaths.")
    d.raw("- Finding and treating people who have stopped injecting, through community testing outside prison, is needed to reduce deaths.")
    d.raw("- These results come from a model with important uncertainties, including the size of the community population and the costs, and should be checked against programme data.")

    # ---------------- abbreviations (PLOS: non-standard ones, alphabetical) ----------------
    d.h(1, "Abbreviations")
    d.raw("CHEERS, Consolidated Health Economic Evaluation Reporting Standards; D, people who inject drugs and have never been imprisoned; "
          "F, released people who still inject; HCV, hepatitis C virus; J, people in prison; MCMC, Markov chain Monte Carlo; "
          "PWID, people who inject drugs; UI, uncertainty interval; WHO, World Health Organization; X, released people who have stopped injecting.")

    # ---------------- introduction ----------------
    d.h(1, "Introduction")
    d.p("WHO set targets for 2030 to eliminate viral hepatitis as a public health threat, including large reductions in new HCV infections and in HCV deaths compared with 2015 [@who2022strategy; @who2023validation].", kind="factual")
    d.p("Direct-acting antivirals cure most people treated [@pawlotsky2020], and treating people who inject drugs can also prevent onward transmission [@martin2013; @dore2021].", kind="factual", join=True)
    d.p("In Singapore, people who inject drugs are a key population for HCV [@winslow2007; @cda_hcv], and their number has been estimated at about 11,000 [@teo2019].", ("targets",), kind="factual")
    d.p("Many of them pass through prison, and treatment in prison with direct-acting antivirals expanded from 2019 [@prison2021; @tan2023].", ("config",), kind="factual", join=True)
    d.p("Prison offers a setting where people can be found and treated, but linkage to care after release is often poor [@clink2022; @eckhardt2022].", kind="factual")
    d.p("Community testing programmes in Singapore have been described in policy statements [@moh2024ett; @moh2025ett; @cda_programme], and point-of-care testing can improve linkage to care [@hsiang2022].", kind="factual", join=True)
    d.p("A previous model examined what is needed to eliminate HCV in Singapore [@chaillon2021]. Our model adds the distinction between released people who return to injecting and those who stop.", kind="interpretive", join=True)
    d.p("This distinction matters for policy. People who still inject keep transmitting and are often re-arrested, so prison treatment reaches them again. People who have stopped do not transmit, but they may carry long-standing infections that progress to cirrhosis and liver cancer, and they are less likely to return to prison.", ("rates", "sim"), kind="interpretive")
    d.p("We used a calibrated transmission model to estimate whether prison treatment alone can meet the WHO 2030 targets for PWID in Singapore, and how a fixed annual number of community screening tests should be split between released people who still inject and those who have stopped.", ("config", "scen_code"), kind="method")

    # ---------------- methods ----------------
    d.h(1, "Methods")
    d.h(2, "Study design")
    d.p("This is a modelling study with a costing of the screening and treatment policies. We report it following the Consolidated Health Economic Evaluation Reporting Standards 2022 (CHEERS 2022) [@cheers2022] where they apply (S1 Checklist).", kind="method")
    d.raw("[Protocol registration: to be confirmed by the authors. No registered protocol is recorded in the project repository.]")
    d.h(2, "Model structure")
    d.p(f"We built a deterministic compartmental model of PWID in Singapore (Fig 1). People enter the model when they start injecting and are placed in group D (injecting, never imprisoned). "
        f"First arrest moves them to J (in prison). On release, a share returns to injecting (F) and the rest stop (X). People in F can be re-arrested; people in X are not.", ("sim", "setup"), kind="method")
    d.p(f"Each group was divided into six age groups (<20, 20–29, 30–39, 40–49, 50–59 and 60 years or older), four liver disease stages (non-cirrhotic, compensated cirrhosis, decompensated cirrhosis and hepatocellular carcinoma) and five HCV states: susceptible, acute, chronic, on treatment and recovered. "
        f"Recovered people, after spontaneous clearance or cure, remain antibody positive and can be reinfected.", ("sim", "targets"), kind="method")
    d.p("Transmission occurs only among people who inject (D and F). The force of infection in each age group depends on an age contact matrix and on the share of infectious (acute or chronic) people in each age group of D and F. There is no transmission in prison or after people stop injecting.", ("sim",), kind="method")
    d.p(f"Natural history, treatment efficacy and mortality parameters came from the literature (S1 Text, Table A). Background mortality was the 2015 Singapore death rate by age, multiplied by a standardised mortality ratio for people who use drugs of {f.setup('smr_viraemic')} for viraemic people and {f.setup('smr_nonviraemic')} for others [@vansanten2014]. "
        "Decompensated cirrhosis and hepatocellular carcinoma add liver deaths, which are lower after cure.", ("setup", "sim"), kind="method")
    d.p("The model was solved with a fourth-order Runge–Kutta method with a daily step. Each run started from a small seeded population and ran for 150 model years, so that 2015 was a steady state; runs whose last ten years changed by more than set tolerances were excluded (S1 Text, Section 1).", ("sim", "equilibrium", "targets"), kind="method")

    d.h(2, "Prison rates")
    lam2 = [f.rate_ci[("lambda2", i)][0] for i in range(1, 7)]; lam3 = [f.rate_ci[("lambda3", i)][0] for i in range(1, 7)]
    piv = f.rate_ci[("pi", 1)]
    d.p(f"We estimated release and re-arrest rates by age from Singapore Prison Service records of {n0(f.sps_people)} people held for drug consumption with admissions from {f.sps_first} to {f.sps_last} ({n0(f.sps_stays)} prison stays), using Bayesian survival models fitted in Stan (S1 Text, Section 3). "
        f"Release rates were {min(lam2):.2f} to {max(lam2):.2f} per year by age. A cure-type model gave a {pct(piv[0])} probability of eventual re-arrest (95% credible interval {pct(piv[1])}–{pct(piv[2])}), and re-arrest rates of {min(lam3):.2f} to {max(lam3):.2f} per year among those re-arrested.",
        ("rates", "release", "setup", "sps"), kind="method")
    d.p("We assumed that people who are eventually re-arrested return to injecting (F) and that the others stop (X). The first-arrest rate from D was the re-arrest rate multiplied by a calibrated factor c.", ("setup", "sim"), kind="method")

    d.h(2, "Entry and calibration targets")
    d.p(f"Entry into D by age was the number of new drug abusers recorded by the Central Narcotics Bureau in 2015 ({', '.join(str(x) for x in f.cnb_counts)} by age group; {cnb_total:,} in total) [@singstat_drug], multiplied by a calibrated scale of at least 1, and was held constant in every year.", ("cnb", "likelihood"), kind="method")
    d.p(f"We calibrated the model to three 2015 targets (S1 Text, Table C): HCV antibody prevalence by age among people in prison at the end of 2015 with a result from the December 2014 to February 2016 screening exercise ({f.screened_2015:,} people); the same Central Narcotics Bureau counts as the annual number of first arrests by age; "
        f"and {f.tg['community_pwid']:,} active injectors in the community (95% confidence interval 6,500–17,000), a network scale-up estimate [@teo2019], compared with the modelled number in D and F.", ("targets", "prison", "cnb"), kind="method")
    d.h(2, "Calibration")
    d.p(f"We estimated 13 parameters on the log scale: six multipliers of the contact matrix (one per age group), six entry scales and the first-arrest factor c. Priors were fixed before fitting: normal with mean 0 and standard deviation 1.5 for the contact multipliers, half-normal with standard deviation 1 for the entry scales, and normal with mean 0 and standard deviation 0.5 for log c, each truncated to wide bounds (S1 Text, Table D).",
        ("likelihood",), kind="method")
    d.p(f"The likelihood combined a beta distribution for each seroprevalence (precision {f.tg['beta_precision']}), a negative binomial distribution for each first-arrest count (size {f.tg['first_arrests_nb_size']}) and a log-normal distribution for the community estimate (log-scale standard deviation {f.tg['community_pwid_log_sd']}, taken from its confidence interval).", ("likelihood", "targets"), kind="method")
    d.p(f"A multistart Nelder–Mead search found the posterior mode. We then sampled the posterior with an adaptive random-walk Metropolis algorithm [@haario2001]: {f.mcfg['n_chains']} chains started from independent prior draws, each with {f.mcfg['max_evals_per_chain']:,} model runs, of which {f.mcfg['n_warmup']:,} were warm-up. "
        f"Convergence was checked with rank-normalised R-hat and effective sample sizes [@vehtari2021].", ("nm_cfg", "mcmc_cfg"), kind="method")
    d.p("We checked the sampler with simulation-based calibration [@talts2018] on simulated data sets, and compared the posterior with neural posterior estimation [@greenberg2019; @sbi2020] as a cross-check (S1 Text, Section 4).", kind="method")

    d.h(2, "Policy scenarios")
    d.p(f"All scenarios share the same history: no treatment until 2018 and treatment in prison from 2019, at a rate that gives a {pct(f.cfg['j_annual_start_probability'], 3)} annual probability of starting treatment for people with chronic infection in J. "
        f"From 2026, we added community screening for F and X. We compared prison treatment only (J only), {len(f.finite)} finite policies ({', '.join(f'{c:,}' for c in caps)} tests per year, with {shares[0]}%, {shares[1]}% or {shares[2]}% given to F and the rest to X), and uncapped screening of everyone in F and X each year as a reference.",
        ("config", "scen_code"), kind="method")
    d.p(f"The smallest capacity, {caps[0]} tests per year, was set with reference to a public community testing programme [@moh2024ett]; the other capacities and shares are illustrative and are not announced policies.", ("config",), kind="method")
    d.p("At the start of each year, the tests given to F or X were applied to people not on treatment, in proportion to their HCV states. Each antibody-positive person had an RNA test, and the number of chronic infections found set that year's treatment capacity in the group. "
        f"Within that capacity, people with chronic infection started treatment at a rate that gives a {pct(f.cfg['fx_annual_start_probability'], 0)} annual probability. Sofosbuvir–velpatasvir without ribavirin was assumed, and people in D were never treated.", ("config", "scen_code", "setup"), kind="method")
    d.h(2, "Outcomes")
    d.p("The primary outcomes were the reductions in 2030, compared with 2015, of annual new HCV infections among PWID and of annual HCV deaths. "
        "New infections were the annual number of infections in D and F. HCV deaths were liver deaths in decompensated cirrhosis and hepatocellular carcinoma, plus the excess non-liver deaths of viraemic people over the rate of non-viraemic people.", ("scen_code",), kind="method")
    d.p("A policy met a target if the median reduction was at least 80% for infections or 65% for deaths. We also report the share of posterior draws in which both targets were met.", ("scen_code", "config"), kind="method")
    d.h(2, "Costs")
    d.p(f"We costed antibody tests (S${f.cfg['unit_cost_sgd']['antibody']}), RNA tests (S${f.cfg['unit_cost_sgd']['rna']}), one consultation per treatment start (S${f.cfg['unit_cost_sgd']['consultation']}) and a 12-week sofosbuvir–velpatasvir course (S${F(f.costev[0]['base_cost_sgd']):,.0f}) in 2026 Singapore dollars, from published charges and a published price (S1 Text, Table H) [@wong2019; @nhg_dsc; @parkway_rna; @sgd_rate].",
        ("costev", "config"), kind="method")
    d.p(f"Costs for 2026–2030 in all three groups (J, F and X) were discounted at {pct(f.cfg['discount_rate'], 0)} a year to 2026. Unit prices were varied together by a factor of 0.5 and 1.5 in a sensitivity analysis. These are charges, not costs to the health system.", ("config", "costev"), kind="method")
    d.h(2, "Uncertainty")
    d.p(f"For each of {f.ndraw:,} posterior draws we ran every scenario and computed every outcome. We report medians and 95% uncertainty intervals (UI; 2.5th to 97.5th percentiles) across draws. These intervals reflect parameter uncertainty only.", ("summary", "draws_valid"), kind="method")
    d.h(2, "Software")
    d.p("The model was written in C++ and called from R through Rcpp [@rcpp2011]; rate models used Stan. Code and generated results are in the project repository.", kind="method")
    d.h(2, "Ethics statement")
    d.p("[Ethics statement: to be completed by the authors. Record-level Singapore Prison Service data were used to estimate prison rates and the 2015 seroprevalence targets; the approving ethics committee, the approval number and the consent arrangements are not recorded in the project repository.]", ("sps",), kind="declaration")

    # ---------------- results ----------------
    d.h(1, "Results")
    d.h(2, "Calibration and the 2015 population")
    d.p(f"The posterior chains converged (largest R-hat {f.rhat_max:.3f}; smallest bulk effective sample size {n0(f.ess_bulk_min)}). The model reproduced the age pattern of seroprevalence in prison, and the 95% posterior intervals of all {f.n_ovl} modelled targets overlapped the target ranges (Fig 2).",
        ("mcmc_diag", "overlap"), unc="moderate")
    pred = f.nm_pred
    d.p(f"At ages 30 and over, modelled first arrests were higher than the recorded counts (for example {F(pred[2]['first_arrests_hat']):.0f} against {pred[2]['first_arrests_target']} at ages 30–39 at the posterior mode), because entry could not fall below the recorded counts and the entry scales settled at their lower bound (posterior medians {min(v[0] for v in f.entry_post.values()):.2f} to {max(v[0] for k, v in f.entry_post.items() if k < 12):.2f} below age 60).",
        ("nm_pred", "mcmc_par"), unc="moderate")
    d.p(f"The modelled number of active injectors in the community was {interval(df['median'], df['lower_95'], df['upper_95'])}, below the survey estimate of {f.tg['community_pwid']:,} but within its range. The first-arrest factor c was {f.c_post[0]:.2f} (95% credible interval {f.c_post[1]:.2f}–{f.c_post[2]:.2f}), so first arrest was faster than re-arrest.",
        ("strata", "mcmc_par", "targets"), unc="moderate")
    d.p(f"In 2015 the model had {n0(ever['median'])} people who had ever injected (95% credible interval {n0(ever['lower_95'])}–{n0(ever['upper_95'])}): {n0(f.strata['population_D']['median'])} in D, {n0(jpred['median'])} in J, {n0(f.strata['population_F']['median'])} in F and {n0(f.strata['population_X']['median'])} in X. "
        f"The modelled prison population was not a calibration target and was much larger than the {f.j_obs_2015:,} people with a screening result recorded in prison in 2015.", ("strata", "prison"), unc="moderate")
    d.p(f"In 2015 the model had {interval(b15['incidence'], b15['incidence_lo'], b15['incidence_hi'])} new infections and {interval(b15['mortality'], b15['mortality_lo'], b15['mortality_hi'])} HCV deaths per year. Without treatment, both stayed at their 2015 level until 2030.",
        ("traj", "notx"), unc="moderate")

    d.h(2, "Prison treatment alone")
    j25 = f.tr("J_only", 2025)
    d.p(f"With treatment in prison from 2019, new infections fell steadily (Fig 3). By 2025 they were {pct(j25['incidence_ratio_2015'], 0)} of the 2015 level, and by 2030 the reduction was {red(J, 'incidence_reduction')}, meeting the incidence target at the median and in {pct(J['who_incidence_pass_prob'], 0)} of draws.",
        TRJ + OUT, unc="moderate")
    d.p(f"HCV deaths fell less: by {red(J, 'mortality_reduction')} in 2030. No posterior draw met the mortality target with prison treatment alone.", OUT, unc="moderate")
    d.p(f"The reason is where infections remain. Released people who still inject are often re-arrested and treated again in prison, but people who have stopped injecting are never re-arrested in the model. In 2026, under prison treatment only, {xs} of the {n0(f.x_chronic + f.f_chronic)} people with chronic infection in F and X were in X, and F and X together accounted for {fx26} of all HCV deaths.",
        ("annual", "traj"), unc="moderate")
    d.p(f"People in D, who had never been imprisoned, accounted for {pct(f.d_inc_share_j, 0)} of new infections in 2030 under prison treatment only; no policy in this study reached them.", OUT, unc="moderate")

    d.h(2, "Community screening policies")
    d.p(f"All {len(f.finite)} finite screening policies met the incidence target at the median, but only {cap_word(len(f.mort_pass)).lower()} met the mortality target (Fig 4 and Table 1). "
        f"With {caps[-1]:,} tests a year, the fall in deaths was {red(P8, 'mortality_reduction')} with 25% to F and {red(P85, 'mortality_reduction')} with 50% to F; both targets were met in {pct(P8['joint_pass_probability'], 1)} and {pct(P85['joint_pass_probability'], 0)} of draws.",
        OUT, unc="moderate")
    d.p(f"With 5,000 tests a year and 25% to F, deaths fell by {red(P5, 'mortality_reduction')}, and both targets were met in {pct(P5['joint_pass_probability'], 0)} of draws. With 50% to F at the same capacity, deaths fell by {red(P550, 'mortality_reduction')} and both targets were met in only {pct(P550['joint_pass_probability'], 0)} of draws.",
        OUT, unc="moderate")
    d.p("At every capacity, giving a larger share of tests to F increased the fall in new infections and reduced the fall in deaths. All quotas were used in full in every year and group.", ("summary", "annual"), unc="moderate")
    d.p(f"This trade-off follows from the two groups. At the start of 2026, {n0(f.f_eligible)} people in F were eligible for screening, of whom {fchr} had chronic infection; in X, {n0(f.x_eligible)} were eligible and {xchr} had chronic infection. "
        "A test given to X therefore found more chronic infections, whereas a test given to F found fewer infections but in people who could transmit.", ("annual",), unc="moderate")
    d.p(f"Uncapped screening of everyone in F and X each year reduced new infections by {red(INF, 'incidence_reduction')} and deaths by {red(INF, 'mortality_reduction')}. Even then, {pct(f.d_inc_share_sinf, 0)} of the remaining infections in 2030 were in D.", OUT, unc="moderate")

    d.h(2, "Costs")
    d.p(f"Discounted costs for 2026–2030 were {money(J, 'total_discounted_cost_2026_2030_sgd')} with prison treatment only. The cheapest policy that met both targets at the median, 5,000 tests a year with 25% to F, cost {money(P5, 'total_discounted_cost_2026_2030_sgd')}, or {money(P5, 'incremental_discounted_cost_vs_j_only_sgd')} more than prison treatment only.",
        OUT, unc="moderate")
    d.p(f"Drugs made up {pct(f.drug_share_5000, 0)} of the cost of this policy, so the total scales almost directly with the drug price. With {caps[-1]:,} tests a year and 25% to F, the cost was {money(P8, 'total_discounted_cost_2026_2030_sgd')}; S1 Text, Section 8 gives all policies.", ("annual", "summary"), unc="moderate")

    # ---------------- discussion ----------------
    d.h(1, "Discussion")
    d.p(f"In this model of PWID in Singapore, prison treatment from 2019 was enough to reduce new HCV infections by more than four fifths by 2030, but it reduced HCV deaths by only about half. "
        f"The gap arises because released people account for nearly all HCV deaths, and most of their chronic infections are among people who have stopped injecting, whom prison treatment does not reach again. Meeting the mortality target required community screening of several thousand released people a year, with most tests given to people who have stopped injecting.",
        OUT + ("annual",), kind="interpretive", unc="moderate")
    d.p("The two WHO targets pull screening in different directions. Tests given to people who still inject prevent more infections, while tests given to people who have stopped prevent more deaths. Because prison treatment already achieves the incidence target at the median, the binding target in Singapore is mortality, and screening policy should be judged mainly by the deaths it prevents.",
        OUT, kind="interpretive", unc="moderate")
    d.p("This suggests that the people who matter most for HCV mortality are not those most visible to drug services. People who have stopped injecting may no longer be in contact with prison or drug services, and reaching them may need testing in primary care or other community settings [@hsiang2022].", kind="interpretive")
    d.p("Part of the modelled HCV mortality is excess non-liver death among viraemic people, based on mortality ratios from a cohort of people who use drugs [@vansanten2014]. Evidence that viraemic infection raises mortality from causes other than liver disease supports this [@lee2012; @cepeda2022], but the excess is not found in every cohort [@piselli2021]; with a smaller excess, deaths would depend more on cirrhosis and liver cancer.", kind="interpretive")
    d.p("We extend earlier modelling for Singapore [@chaillon2021] by separating released people by injecting status and by modelling screening capacity explicitly. This separation gives a reason why the mortality target may be harder to meet than the incidence target.", kind="interpretive")
    d.h(2, "Strengths and limitations")
    d.p(f"The model was calibrated with a full Bayesian method that converged well, and every policy result carries posterior uncertainty from {f.ndraw:,} draws. Prison flows were estimated from individual prison records rather than assumed, and screening capacity was modelled explicitly, so treatment depends on how many people are found.",
        ("mcmc_diag", "summary", "rates"), kind="interpretive", unc="moderate")
    d.p(f"The study has important limitations. First, the modelled prison population ({n0(jpred['median'])}) was about {F(jpred['median']) / f.j_obs_2015:.0f} times the number of screened people recorded in prison in 2015. The prison size was not fitted, and the model may keep people in prison or in contact with it for too long, which would favour prison treatment.",
        ("strata", "prison"), kind="interpretive", unc="high")
    d.p(f"Second, the entry and first-arrest parameters were weakly identified. Central Narcotics Bureau counts include people who do not inject, the entry scales settled at their lower bound, and the posterior of c was only modestly narrower than its prior (standard deviation {pct(f.sd_ratio[12], 0)} of the prior's). The size of the community population that never enters prison is therefore uncertain.",
        ("mcmc_par", "mcmc_cfg", "cnb"), kind="interpretive", unc="high")
    d.p("Third, we assumed that released people who are eventually re-arrested return to injecting and that the others stop permanently, and that there is no transmission in prison. Relapse after a long gap would put more infections back into F.", ("setup", "sim"), kind="interpretive", unc="high")
    d.p("Fourth, screening was an aggregate approximation: tests were spread across each group in proportion to its HCV states, with perfect linkage within capacity and no loss to follow-up. Real programmes will find and treat fewer people for the same number of tests.", ("scen_code",), kind="interpretive", unc="high")
    d.p("Fifth, costs were published charges and a published drug price carried forward without inflation, not payer costs, and follow-up visits were excluded. Our uncertainty intervals cover model parameters only, not structural assumptions, prices or uptake.", ("costev", "config"), kind="interpretive", unc="high")
    d.p(f"Finally, the seroprevalence targets could not be reproduced exactly from the current processing of the prison records (largest difference {100 * f.sero_max_diff:.1f} percentage points), and the neural posterior used as a cross-check failed its calibration check (S1 Text, Section 4).",
        ("targets", "prison", "npe_sbc"), kind="result", unc="moderate")
    d.h(2, "Conclusions")
    d.p("In this model, treatment in prison is likely to meet the 2030 target for new HCV infections among people who inject drugs in Singapore, but not the target for HCV deaths. "
        "Reducing deaths requires finding and treating people who have stopped injecting after release, which needs community testing on a scale of thousands of people a year.",
        OUT, kind="interpretive", unc="moderate")

    # ---------------- statements ----------------
    d.h(1, "Data availability")
    d.raw("[Data availability statement: to be completed by the authors. Model code, aggregate inputs and all generated results are in the project repository (link and DOI to be added). The record-level Singapore Prison Service data cannot be shared by the authors; the access conditions and the contact for access requests must be stated.]")
    d.h(1, "Funding")
    d.raw("[Funding statement: to be completed by the authors, including grant numbers and the role of any funder, or the PLOS statement for unfunded work.]")
    d.h(1, "Competing interests")
    d.raw("[Competing interests: to be declared by each author.]")
    d.h(1, "Author contributions")
    d.raw("[CRediT author contributions: to be completed by the authors.]")
    d.h(1, "Use of artificial intelligence")
    d.p("[AI use statement: to be completed and approved by the authors. Draft: Claude (Anthropic) was used to write and revise analysis code and to draft and revise the text, figures and supporting information, with procedures from Scientific Agent Skills [@kassis2026]. The authors are to state that they checked all content, data, code, references and conclusions, and take responsibility for them.]", ("kassis2026",), kind="declaration")
    d.h(1, "Acknowledgments")
    d.raw("[Acknowledgments: to be completed by the authors.]")
    d.h(1, "References")
    d.raw("::: {#refs}\n:::")

    # ---------------- figure and table legends ----------------
    d.h(1, "Figure captions")
    d.p("**Fig 1. Model structure.** People who inject drugs enter D (never imprisoned) and move to J (prison) at first arrest. On release they move to F (still injecting) or X (stopped). People in F can be re-arrested. Transmission occurs only in D and F. Treatment is given in J from 2019 and community screening is offered to F and X from 2026.", ("sim", "config"), kind="method")
    d.p("**Fig 2. Calibration to the 2015 targets.** (A) HCV antibody prevalence in prison by age. (B) First arrests by age (Central Narcotics Bureau new drug abusers). (C) Active injectors in the community (D and F). Black: target with the 95% range of its likelihood; blue: 95% posterior interval of the model value from the MCMC posterior; orange: neural posterior estimate (cross-check only).", ("overlap",), kind="result")
    d.p(f"**Fig 3. New infections and HCV deaths, 2015–2030.** Annual new HCV infections (A) and HCV deaths (B) relative to 2015, without treatment, with prison treatment only, and with {caps[-1]:,} community tests a year from 2026 given 25% or 75% to F. Yearly counts are drawn at mid-year and joined by a monotone curve; vertical lines mark the start of prison treatment (2019) and community screening (2026). Lines are medians and bands are 95% uncertainty intervals across {f.ndraw:,} posterior draws. Dashed lines mark the WHO targets.", OUT + TRJ, kind="result")
    d.p("**Fig 4. Reductions in 2030 by screening capacity and share of tests given to F.** (A) New infections. (B) HCV deaths. Points are medians and bars are 95% uncertainty intervals. Dashed lines mark the targets; the dotted line is prison treatment only and the solid line is uncapped screening.", OUT, kind="result")

    # ---------------- table 1 ----------------
    d.h(1, "Tables")
    d.p("**Table 1. Reductions in 2030 compared with 2015, probability of meeting both WHO targets, and discounted costs for 2026–2030.**", OUT, kind="result")
    rows = []
    for k in ["J_only"] + [f"S{c}_F{x}" for c in caps for x in shares] + ["SInf"]:
        r = s[k]
        rows.append([scen_label(r), col(r, "incidence_reduction", scale=100, d=1, suf="%").replace(" (95% UI ", " ("),
                     col(r, "mortality_reduction", scale=100, d=1, suf="%").replace(" (95% UI ", " ("),
                     pct(r["joint_pass_probability"], 1), col(r, "total_discounted_cost_2026_2030_sgd", scale=1e-6, d=1).replace(" (95% UI ", " (")])
    d.table(["Policy", "Fall in new infections, % (95% UI)", "Fall in HCV deaths, % (95% UI)", "Draws meeting both targets", "Cost, S$ million (95% UI)"], rows, OUT)
    d.p("UI, uncertainty interval across posterior draws. Costs are 2026 Singapore dollars for J, F and X, discounted at 3% a year to 2026.", ("config",), kind="method")

    # ---------------- supporting information ----------------
    d.h(1, "Supporting information")
    d.p("**S1 Text. Technical appendix.** Model equations, parameters, prison rate estimation, calibration and its checks, the 2015 model population, scenario implementation, costs and additional results.", ("sim",), kind="declaration")
    d.p("**S1 Checklist. CHEERS 2022 checklist.**", ("cheers2022",), kind="declaration")
    return d
