# =============================================================================
# main_text.py — the PLOS Medicine research article (main text). Every number comes from
# Facts; every factual line is bound to its evidence through Doc.p().
# =============================================================================
from fmt import F, pct, ui, col, red, money, n0, interval, cap_word, scen_label, WORDS, lc

OUT = ("summary",)          # scenario summary file
TRJ = ("traj",)
CAL = ("mcmc_par", "mcmc_diag", "mcmc_cfg")


def checks(f):
    """The narrative below is written for these results; stop the build if they change."""
    s, sens = f.summ, f.j2019
    assert f.j_start == "2026" and f.j2019_valid() == "2019"
    assert set(f.inc_pass) == {"S8000_F25", "S5000_F50", "S8000_F50", "S5000_F75", "S8000_F75"}, "incidence-passing set changed"
    assert f.mort_pass == [], "a finite policy now meets the mortality target"
    assert F(s["J_only"]["who_incidence_pass_prob"]) == 0 and F(s["J_only"]["who_mortality_pass_prob"]) == 0
    assert 0.45 < F(s["J_only"]["incidence_reduction"]) < 0.55 and 0.15 < F(s["J_only"]["mortality_reduction"]) < 0.25, "'about half' / 'about a fifth' changed"
    assert s["SInf"]["who_joint_pass"] == "TRUE" and 0.70 < F(s["SInf"]["mortality_reduction"]) < 0.75, "uncapped screening margin changed"
    assert {x["input"] for x in f.tor_flips if x["scenario_id"] == "SInf"} == {"psi_DC"}, "inputs that flip uncapped screening changed"
    best = max(f.finite, key=lambda k: F(s[k]["mortality_reduction"]))
    assert best == "S8000_F25", "best finite policy changed"
    assert all(F(s[k]["joint_pass_probability"]) == 0 for k in f.finite)
    assert f.inc_order_ok and f.mort_order_ok, "F/X ordering of the reductions changed"
    assert 0.65 < f.x_chronic_share < 0.75 and f.fx_death_share_2026 > 0.85, "X / F+X shares changed"
    assert f.f_chronic / f.f_eligible > f.x_chronic / f.x_eligible, "F no longer has the higher chronic share"
    assert all(abs(v[st] - 1) < 1e-6 for v in f.quota_use.values() for st in "FX"), "a screening quota is no longer used in full"
    assert f.history_dev < 1e-3, "the main run is not flat before 2026"
    assert f.rhat_max < 1.01 and f.n_ovl == 13
    assert abs(F(f.strata["population_J"]["median"]) / f.j_held_2015 - 1) < 0.1, "modelled prison size no longer close to the number held"
    assert f.prior_dominated == [5, 13]
    assert all(abs(F(f.tr("no_treatment", 2030, f.notx)[c]) - 1) < 0.01 for c in ("incidence_ratio_2015", "mortality_ratio_2015"))
    # J coverage grid: reductions rise with coverage; only uncapped screening meets the mortality target
    lv = sorted(f.jcov)
    assert all([F(f.jcov[p][k][c]) for p in lv] == sorted(F(f.jcov[p][k][c]) for p in lv)
               for k in s for c in ("incidence_reduction", "mortality_reduction"))
    assert all([k for k in f.jcov[p] if f.jcov[p][k]["who_mortality_pass"] == "TRUE"] == ["SInf"] for p in lv)
    assert 0.65 <= F(f.jcov[0.0]["SInf"]["mortality_reduction"]) < 0.67, "uncapped screening without prison treatment no longer just meets the target"
    assert [k for k in f.jcov[0.0] if f.jcov[0.0][k]["who_incidence_pass"] == "TRUE"] == ["S8000_F75", "SInf"]
    assert abs(F(f.jcov[0.0]["J_only"]["incidence_reduction"])) < 0.005
    assert [p for p in lv if sum(f.jcov[p][k]["who_incidence_pass"] == "TRUE" for k in f.finite) == len(f.inc_pass)] == lv[2:]
    # sensitivity: prison treatment from 2019
    assert sens["J_only"]["who_incidence_pass"] == "TRUE" and sens["J_only"]["who_mortality_pass"] != "TRUE"
    assert not any(sens[k]["who_mortality_pass"] == "TRUE" for k in f.finite) and sens["SInf"]["who_joint_pass"] == "TRUE"
    assert max(f.finite, key=lambda k: F(sens[k]["mortality_reduction"])) == "S8000_F25"
    assert all(sens[k]["who_incidence_pass"] == "TRUE" for k in f.finite)
    # genotype 3 share: same pass sets as the main analysis
    for g in f.gt3.values():
        assert [k for k in g if g[k]["who_incidence_pass"] == "TRUE"] == [k for k in s if s[k]["who_incidence_pass"] == "TRUE"]
        assert [k for k in g if g[k]["who_mortality_pass"] == "TRUE"] == ["SInf"]
    # long horizon, cost per outcome averted and strategy selection (decision horizon 2026-end year)
    hz, E, c = f.hz, f.end_year, f.cum
    assert E == 2050, "the text describes a projection to 2050"
    yr = first_years(f)
    assert yr[("J_only", "mortality_ratio_2015")] > 2045 and yr[("SInf", "mortality_ratio_2015")] <= 2030
    assert yr[("S8000_F25", "mortality_ratio_2015")] < yr[("S8000_F75", "mortality_ratio_2015")] < yr[("J_only", "mortality_ratio_2015")]
    caps = f.cfg["annual_screening_totals"]
    assert all(F(c[f"S{n}_F25"][f"liver_deaths_averted_vs_j_only_2026_{E}"]) > F(c[f"S{n}_F75"][f"liver_deaths_averted_vs_j_only_2026_{E}"]) for n in caps), \
        "F 25% no longer averts more liver deaths than F 75% at every capacity"
    lowF75 = [f"S{n}_F75" for n in caps[:3]]
    assert all(F(c[k][f"incremental_programme_cost_vs_j_only_sgd_2026_{E}_lo"]) < 0 < F(c[k][f"incremental_programme_cost_vs_j_only_sgd_2026_{E}_hi"]) for k in lowF75), \
        "low-capacity F 75% policies are no longer budget-neutral against J only"
    top = lambda st, h: next(x for x in f.da if x["set"] == st and x["horizon"] == h and x["topsis_rank"] == "1")
    assert top("s17", hz)["strategy"] == "S8000_F75" and F(top("s17", hz)["p_topsis_top"]) > 0.95
    assert top("s17", "2026-2030")["strategy"] == "SInf" and top("s68", hz)["strategy"] == "SInf @ J 99.995%"
    assert frontier(f, "liver_deaths_averted") == ["S1000_F75", "S2000_F75", "S5000_F75", "S8000_F75", "S8000_F50", "SInf"]
    assert frontier(f, "infections_averted") == ["S1000_F75", "S2000_F75", "S5000_F75", "S8000_F75", "SInf"]
    d17 = [x for x in f.da if x["set"] == "s17" and x["horizon"] == hz]
    assert {x["strategy"] for x in d17 if x["rank_cpa_infection"] == "1" or x["rank_cpa_liver_death"] == "1"} == {"S5000_F75"}
    assert max(d17, key=lambda x: F(x["p_top_random_weights"]))["strategy"] == "SInf"
    fr = {(x["benefit"], x["strategy"]): x for x in f.da_frontier if x["set"] == "s17" and x["horizon"] == hz}
    assert all(fr[(b, k)]["status"] != "frontier" for b in ("infections_averted", "liver_deaths_averted") for k in ["J_only"] + [f"S{n}_F25" for n in caps])
    assert F(fr[("liver_deaths_averted", "S1000_F75")]["cost"]) <= F(fr[("liver_deaths_averted", "J_only")]["cost"])
    d68 = [x for x in f.da if x["set"] == "s68" and x["horizon"] == hz]
    assert {x["strategy"].split(" @ ")[0] for x in d68 if x["topsis_rank"] in ("1", "2", "3", "4")} == {"SInf"}
    by_pol = {}
    for x in d68:
        if x["pareto_efficient"] == "TRUE":
            pol, cov = x["strategy"].split(" @ J "); by_pol.setdefault(pol, []).append((float(cov.strip("%")), F(x["topsis_closeness"])))
    assert all([v for _, v in sorted(cl)] == sorted(v for _, v in cl) for cl in by_pol.values()), "closeness no longer rises with prison coverage"
    # linkage: lower linkage lowers every finite effect; still only uncapped screening meets the mortality target
    for l, L in f.link.items():
        assert [k for k in L if L[k]["who_mortality_pass"] == "TRUE"] == ["SInf"] and L["SInf"]["who_joint_pass"] == "TRUE"
        assert all(F(L[k]["mortality_reduction"]) < F(s[k]["mortality_reduction"]) for k in f.finite)
    assert all(F(f.link[0.6][k]["mortality_reduction"]) < F(f.link[0.8][k]["mortality_reduction"]) for k in f.finite)


def first_years(f):
    """First year in which the median yearly count is at or below the WHO target share of 2015."""
    out = {}
    for x in f.long:
        for col, thr in (("incidence_ratio_2015", 0.2), ("mortality_ratio_2015", 0.35)):
            if F(x[col]) <= thr: out.setdefault((x["scenario_id"], col), int(F(x["year"])))
    return out


def frontier(f, benefit, st="s17"):
    """Strategies on the cost-effectiveness frontier (in order of cost) at the decision horizon."""
    fr = [x for x in f.da_frontier if x["set"] == st and x["horizon"] == f.hz and x["benefit"] == benefit and x["status"] == "frontier"]
    return [x["strategy"] for x in sorted(fr, key=lambda x: F(x["cost"]))]


def kS(row, c):  # a cost per outcome in S$, rounded to the nearest thousand, with its UI
    r = lambda v: f"S${round(F(v), -3):,.0f}"
    return f"{r(row[c])} (95% UI {r(row[c + '_lo'])}–{r(row[c + '_hi'])})"


def build(reg, f):
    from registry import Doc
    checks(f)
    d = Doc(reg, "main")
    s = f.summ; J, P8, P85, P875, INF = (s[k] for k in ("J_only", "S8000_F25", "S8000_F50", "S8000_F75", "SInf"))
    sens = f.j2019; SJ, S8 = sens["J_only"], sens["S8000_F25"]
    caps = f.cfg["annual_screening_totals"]; shares = [round(100 * x) for x in f.cfg["f_allocation_fractions"]]
    levels = sorted(f.jcov); J0 = f.jcov[0.0]
    beta_total = sum(f.beta)
    jpred = f.strata["population_J"]; df = f.strata["population_current_PWID"]; ever = f.strata["population_ever_PWID"]
    xs, fx26 = pct(f.x_chronic_share, 0), pct(f.fx_death_share_2026, 0)
    fchr, xchr = pct(f.f_chronic / f.f_eligible, 1), pct(f.x_chronic / f.x_eligible, 1)
    b15 = f.tr("J_only", 2015)
    inc_cheapest = min(f.inc_pass, key=lambda k: F(s[k]["total_discounted_cost_2026_2030_sgd"]))
    jstarts26 = F(f.an("J_only", 2026, "J")["treatment_starts"])
    JC = ("jcov",) + OUT
    lab = {x["input"]: x["label"] for x in f.tor_ranges}
    plain = {"p_HCC": "progression to liver cancer", "phi_CC_DC": "the relative risk of decompensation after cure",
             "psi_DC": "the relative risk of liver death after cure in decompensated cirrhosis",
             "psi_HCC": "the relative risk of liver death after cure in liver cancer"}
    sinf_flip = "; ".join(plain.get(i, lc(lab[i])) for i in sorted({x["input"] for x in f.tor_flips if x["scenario_id"] == "SInf"}))
    sinf_psi = sorted(F(x["mortality_reduction"]) for x in f.tor if x["scenario_id"] == "SInf" and x["input"] == "psi_DC")
    SENS = ("j2019",)
    # long horizon, decision analysis and linkage
    E, hz, C = f.end_year, f.hz, f.cum
    H = f"2026–{E}"
    cu = lambda k, n, h=E: col(C[k], f"{n}_2026_{h}", d=0)
    yr = first_years(f)
    LONG = ("long", "cum"); DA = ("da", "da_frontier", "da_code"); LK = ("linkage",)
    d17 = {x["strategy"]: x for x in f.da if x["set"] == "s17" and x["horizon"] == hz}
    d17_30 = {x["strategy"]: x for x in f.da if x["set"] == "s17" and x["horizon"] == "2026-2030"}
    d68 = {x["strategy"]: x for x in f.da if x["set"] == "s68" and x["horizon"] == hz}
    top17 = next(k for k, x in d17.items() if x["topsis_rank"] == "1"); top17_30 = next(k for k, x in d17_30.items() if x["topsis_rank"] == "1")
    top68 = next(k for k, x in d68.items() if x["topsis_rank"] == "1")
    par17 = sorted((k for k, x in d17.items() if x["pareto_efficient"] == "TRUE"), key=lambda k: int(d17[k]["topsis_rank"]))
    rw = sorted(d17, key=lambda k: -F(d17[k]["p_top_random_weights"]))
    fr_ld, fr_inf = frontier(f, "liver_deaths_averted"), frontier(f, "infections_averted")
    lab17 = lambda k: scen_label(s[k])
    budget = [x for x in f.da_budget if x["set"] == "s17" and x["horizon"] == hz and x["benefit"] == "liver_deaths_averted"]
    wts = {x["weights"]: x["top_strategy"] for x in f.da_weights if x["set"] == "s17" and x["horizon"] == hz}
    lowF75 = [f"S{n}_F75" for n in caps[:3]]
    L6, L8 = f.link[0.6], f.link[0.8]; LC6 = f.link_cum[0.6]

    # ---------------- front matter ----------------
    d.raw("---\ntitle: \"Whom to test, where and how intensively: prison treatment and post-release screening to eliminate hepatitis C among people who inject drugs in Singapore. A modelling study with multi-criteria decision analysis\"\n---")
    d.raw("**Short title:** Hepatitis C elimination after prison release in Singapore")
    d.raw("Yichen He\\*")
    d.raw("[Affiliation: to be provided by the authors]")
    d.raw("\\* Corresponding author. E-mail: [to be provided by the authors]")

    # ---------------- abstract ----------------
    d.h(1, "Abstract")
    d.h(2, "Background")
    d.p("Hepatitis C virus (HCV) infection in Singapore is concentrated among people who inject drugs (PWID), many of whom cycle through prison. "
        "Eliminating HCV among PWID needs decisions on whom to test and treat, through which delivery point and how intensively. We asked whether treatment in prison, alone or with community screening of released people from 2026, can meet the World Health Organization 2030 elimination targets (80% fewer new infections and 65% fewer HCV deaths than in 2015), and which strategy is preferred when health benefits, the WHO targets and the programme budget are weighed together.",
        ("config", "scen_code", "who2022strategy"), kind="factual")
    d.h(2, "Methods and findings")
    d.p(f"We built an age-structured model of PWID never imprisoned (D), in prison (J), and released and still injecting (F) or stopped (X), with prison rates estimated from records of {n0(f.sps_people)} people. "
        f"We calibrated it by Markov chain Monte Carlo to 2015 prison seroprevalence, first arrests of new heroin and methamphetamine users and a survey estimate of active injectors. "
        f"In the main analysis nothing changed from 2015 until 2026, when treatment in prison and {len(f.finite)} community screening policies ({caps[0]:,} to {caps[-1]:,} tests a year, {shares[0]}%–{shares[2]}% to F) began; we varied prison treatment coverage and, in a sensitivity analysis, started prison treatment in 2019. We projected to {E} because liver deaths lag infection by decades, counted infections and liver deaths averted against prison treatment only (the status quo) and against no treatment, and chose among strategies with Pareto efficiency and TOPSIS, as in a multi-criteria decision analysis. Results are over {f.ndraw:,} posterior draws.",
        ("setup", "cnb", "drug", "targets", "mcmc_cfg", "config", "summary", "rates", "sps", "cum", "da_code"), kind="method")
    d.p(f"Prison treatment alone from 2026 reduced new infections in 2030 by {red(J, 'incidence_reduction', first=True)} and HCV deaths by {red(J, 'mortality_reduction')}. "
        f"{cap_word(len(f.inc_pass))} of {len(f.finite)} screening policies met the incidence target, but none met the mortality target; the best, {caps[-1]:,} tests a year with 25% to F, reduced deaths by {red(P8, 'mortality_reduction')}. Only uncapped screening of all released people met both targets, reducing deaths by {red(INF, 'mortality_reduction')}, but not if cure lowered the risk of liver death in decompensated cirrhosis less than assumed. "
        f"Moving tests from F to X increased the fall in deaths and reduced the fall in infections. Lower prison treatment coverage reduced every effect. "
        f"Had prison treatment started in 2019, it alone would have reduced new infections by {red(SJ, 'incidence_reduction')}, meeting the incidence target, but no capped policy would have met the mortality target; the best reduced deaths by {red(S8, 'mortality_reduction')}. "
        f"Over {H}, {caps[-1]:,} tests a year with 75% to F averted {cu('S8000_F75', 'liver_deaths_averted_vs_j_only')} liver deaths more than prison treatment alone for an extra programme budget of {money(C['S8000_F75'], f'incremental_programme_cost_vs_j_only_sgd_2026_{E}')}, and was ranked first by TOPSIS with equal weights in {pct(d17[top17]['p_topsis_top'], 0)} of draws; uncapped screening averted the most liver deaths and was preferred when the WHO targets or liver deaths were weighted more. "
        f"Limitations include treating new heroin and methamphetamine users as the people who inject, weakly identified entry parameters, liver disease inputs from other settings, to which the fall in deaths was sensitive, charge-based costs, and uncertainty intervals that cover the calibrated parameters only.",
        ("summary", "jcov", "j2019", "mcmc_par", "drug", "tornado", "cum", "da", "da_weights"), unc="moderate", intent="exploratory")
    d.h(2, "Conclusions")
    d.p("In this model, prison treatment and capped community screening could not meet the 2030 target for HCV deaths among PWID in Singapore, even with prison treatment from 2019; with treatment from 2026, the incidence target needs thousands of community tests a year in addition to prison treatment. "
        "Over a longer horizon, several thousand community tests a year after release, most of them for people who still inject, gave the best balance of infections and liver deaths averted against the programme budget, while screening every released person each year is needed if meeting the 2030 mortality target is the priority.",
        OUT + SENS + DA, kind="interpretive", unc="moderate")

    # ---------------- author summary ----------------
    d.h(1, "Author summary")
    d.h(2, "Why was this study done?")
    d.p("- Most hepatitis C infections in Singapore are among people who inject drugs, and many of them spend time in prison.", ("winslow2007", "prison2021"), kind="factual")
    d.p("- It is unclear whether treatment in prison, with or without community screening after release, can meet the WHO 2030 elimination targets, which ask for large falls in both new infections and deaths.", ("who2022strategy",), kind="factual")
    d.raw("- Community screening capacity is limited, and it is not known whether it should go to released people who still inject or to those who have stopped, or how many tests a year are worth their cost.")
    d.h(2, "What did the researchers do and find?")
    d.p("- We built a model of hepatitis C among people who inject drugs that follows them into and out of prison, and fitted it to prison and national data from 2015.", ("targets",), kind="method")
    d.p(f"- If treatment starts in 2026, prison treatment alone cut new infections by about {round(100 * F(J['incidence_reduction']))}% and deaths by about {round(100 * F(J['mortality_reduction']))}% by 2030. With up to {caps[-1]:,} community tests a year, several policies met the infection target but none met the death target.",
        OUT, unc="moderate")
    d.p("- If prison treatment had started in 2019, it alone would have met the infection target, but only testing every released person every year would have met the death target.", SENS, unc="moderate")
    d.p(f"- Looking to {E}, and weighing infections and liver deaths prevented, the WHO targets and the programme budget together, {caps[-1]:,} community tests a year with most given to people who still inject was the preferred strategy; testing every released person each year prevented the most liver deaths but cost much more.", DA + LONG, unc="moderate")
    d.h(2, "What do these findings mean?")
    d.p("- Time matters: the later treatment starts, the less any policy can achieve by 2030. Deaths from liver disease fall slowly after cure, because people with advanced liver disease remain at risk.", OUT + SENS, kind="interpretive", unc="moderate")
    d.raw("- Prison treatment prevents many new infections, but most hepatitis C deaths are among released people who have stopped injecting, who need community testing.")
    d.p("- Because liver disease takes decades to develop, most of the deaths prevented by treatment now occur after 2030, so strategies should be judged over a longer horizon as well as against the 2030 targets.", LONG + ("thein2008",), kind="interpretive", unc="moderate")
    d.raw("- These results come from a model with important uncertainties, including the size of the community population and the costs, and should be checked against programme data.")

    # ---------------- abbreviations (PLOS: non-standard ones, alphabetical) ----------------
    d.h(1, "Abbreviations")
    d.raw("CHEERS, Consolidated Health Economic Evaluation Reporting Standards; D, people who inject drugs and have never been imprisoned; DAA, direct-acting antiviral; "
          "F, released people who still inject; HCV, hepatitis C virus; J, people in prison; MCMC, Markov chain Monte Carlo; "
          "PWID, people who inject drugs; TOPSIS, Technique for Order of Preference by Similarity to Ideal Solution; UI, uncertainty interval; WHO, World Health Organization; X, released people who have stopped injecting.")

    # ---------------- introduction ----------------
    d.h(1, "Introduction")
    d.p("WHO set targets for 2030 to eliminate viral hepatitis as a public health threat, including large reductions in new HCV infections and in HCV deaths compared with 2015 [@who2022strategy; @who2023validation].", kind="factual")
    d.p("Direct-acting antivirals cure most people treated [@pawlotsky2020], and treating people who inject drugs can also prevent onward transmission [@martin2013; @dore2021].", kind="factual", join=True)
    d.p("In Singapore, people who inject drugs are a key population for HCV [@winslow2007; @cda_hcv], and their number has been estimated at about 11,000 [@teo2019].", ("targets",), kind="factual")
    d.p("Many of them pass through prison, where treatment with direct-acting antivirals has been offered [@prison2021; @tan2023].", kind="factual", join=True)
    d.p("Prison offers a setting where people can be found and treated, but linkage to care after release is often poor [@clink2022; @eckhardt2022].", kind="factual")
    d.p("Community testing programmes in Singapore have been described in policy statements [@moh2024ett; @moh2025ett; @cda_programme], and point-of-care testing can improve linkage to care [@hsiang2022].", kind="factual", join=True)
    d.p("A previous model examined what is needed to eliminate HCV in Singapore [@chaillon2021]. Our model adds the distinction between released people who return to injecting and those who stop.", kind="interpretive", join=True)
    d.p("This distinction matters for policy. People who still inject keep transmitting and are often re-arrested, so prison treatment reaches them again. People who have stopped do not transmit, but they may carry long-standing infections that progress to cirrhosis and liver cancer, and they are less likely to return to prison.", ("rates", "sim"), kind="interpretive")
    d.p("Elimination policy for PWID therefore turns on three questions: whom to test and treat (people who have never been imprisoned, people in prison, or released people who still inject or have stopped), through which delivery point (prison or the community after release), and how intensively (how many community tests a year). "
        "These choices trade off several criteria at once: the two WHO targets, the infections and deaths averted, and the programme budget. Multi-criteria decision analysis makes such trade-offs explicit [@thokala2016]: strategies that are worse on every criterion than another can be set aside (Pareto efficiency), and the remainder ranked by their closeness to an ideal strategy (TOPSIS) [@hwang1981].", kind="interpretive")
    d.p("The time horizon also matters. Chronic HCV infection progresses to cirrhosis over decades [@thein2008], so most of the liver deaths that treatment in 2026 prevents will occur after 2030, and a horizon that ends in 2030 understates the value of treating people with long-standing infection.", kind="interpretive")
    d.p(f"We used a calibrated transmission model to estimate whether prison treatment and community screening starting in 2026 can meet the WHO 2030 targets for PWID in Singapore, how a fixed annual number of community screening tests should be split between released people who still inject and those who have stopped, how the answer depends on prison treatment coverage, and what an earlier start in prison would have achieved. "
        f"We then projected the strategies to {E} and selected among them with Pareto efficiency and TOPSIS, weighing infections and liver deaths averted, the WHO targets and the programme budget.", ("config", "scen_code", "da_code"), kind="method")

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
    d.p(f"Natural history, treatment efficacy and mortality parameters came from the literature (S1 Text, Table A). Each cause of death had one route in the model. Deaths from causes other than liver disease occurred at the 2015 Singapore death rate by age [@singstat_m810141] (ages 60–74 for the oldest group) multiplied by {f.setup('smr_nonviraemic')}, the standardised mortality ratio of people who use drugs without HCV or HIV infection [@vansanten2014], and for viraemic people also by {f.setup('smr_viraemic') / f.setup('smr_nonviraemic'):.2f}, the hazard ratio of extrahepatic death with viraemic HCV infection [@lee2012]. "
        "Liver deaths occurred only in decompensated cirrhosis and hepatocellular carcinoma, at lower rates after cure.", ("setup", "sim"), kind="method")
    d.p("The model was solved with a fourth-order Runge–Kutta method with a daily step. Each run started from a small seeded population and ran for 150 model years, so that 2015 was a steady state; runs whose last ten years changed by more than set tolerances were excluded (S1 Text, Section 1).", ("sim", "equilibrium", "targets"), kind="method")

    d.h(2, "Prison rates")
    lam2 = [f.rate_ci[("lambda2", i)][0] for i in range(1, 7)]; lam3 = [f.rate_ci[("lambda3", i)][0] for i in range(1, 7)]
    piv = f.rate_ci[("pi", 1)]
    d.p(f"We estimated release and re-arrest rates by age from Singapore Prison Service records of {n0(f.sps_people)} people held for drug consumption with admissions from {f.sps_first} to {f.sps_last} ({n0(f.sps_stays)} prison stays), using Bayesian survival models fitted in Stan (S1 Text, Section 3). "
        f"Release rates were {min(lam2):.2f} to {max(lam2):.2f} per year by age. A cure-type model gave a {pct(piv[0])} probability of eventual re-arrest (95% credible interval {pct(piv[1])}–{pct(piv[2])}), and re-arrest rates of {min(lam3):.2f} to {max(lam3):.2f} per year among those re-arrested.",
        ("rates", "release", "setup", "sps"), kind="method")
    d.p("We assumed that people who are eventually re-arrested return to injecting (F) and that the others stop (X). The first-arrest rate from D was the re-arrest rate multiplied by a calibrated factor c.", ("setup", "sim"), kind="method")

    d.h(2, "Entry and calibration targets")
    d.p(f"The Central Narcotics Bureau recorded {sum(f.cnb_counts):,} new drug abusers in 2015 ({', '.join(str(x) for x in f.cnb_counts)} by age group) [@singstat_drug]. "
        f"Of the {f.new_total:,} new abusers that year, {f.new_heroin} mainly used heroin and {f.new_meth:,} methamphetamine [@singstat_m891521], the drugs most often injected in Singapore; we took this share ({pct(f.pwid_share)}) as the share of new abusers who inject, the same at every age, because the drug table has no age breakdown.",
        ("cnb", "drug", "setup"), kind="method")
    d.p(f"Entry into D by age was these counts multiplied by {pct(f.pwid_share)} ({beta_total:,.0f} a year in total) and by a calibrated scale of at least 1, and was held constant in every year.", ("setup", "likelihood", "drug"), kind="method", join=True)
    d.p(f"We calibrated the model to three 2015 targets (S1 Text, Table C): HCV antibody prevalence by age among people in prison at the end of 2015 with a result from the December 2014 to February 2016 screening exercise ({f.screened_2015:,} people); first arrests by age, the same counts multiplied by {pct(f.pwid_share)} and rounded ({', '.join(str(x) for x in f.tg['first_arrests'])}); "
        f"and {f.tg['community_pwid']:,} active injectors in the community (95% confidence interval 6,500–17,000), a network scale-up estimate [@teo2019], compared with the modelled number in D and F.", ("targets", "prison", "cnb", "drug"), kind="method")
    d.h(2, "Calibration")
    d.p(f"We estimated 13 parameters on the log scale: six multipliers of the contact matrix (one per age group), six entry scales and the first-arrest factor c. Priors were fixed before fitting: normal with mean 0 and standard deviation 1.5 for the contact multipliers, half-normal with standard deviation 1 for the entry scales, and normal with mean 0 and standard deviation 0.5 for log c, each truncated to wide bounds (S1 Text, Table D).",
        ("likelihood",), kind="method")
    d.p(f"The likelihood combined a beta distribution for each seroprevalence (precision {f.tg['beta_precision']}), a negative binomial distribution for each first-arrest count (size {f.tg['first_arrests_nb_size']}) and a log-normal distribution for the community estimate (log-scale standard deviation {f.tg['community_pwid_log_sd']}, taken from its confidence interval).", ("likelihood", "targets"), kind="method")
    d.p(f"A multistart Nelder–Mead search found the posterior mode. We then sampled the posterior with an adaptive random-walk Metropolis algorithm [@haario2001]: {f.mcfg['n_chains']} chains started from independent prior draws, each with {f.mcfg['max_evals_per_chain']:,} model runs, of which {f.mcfg['n_warmup']:,} were warm-up. "
        f"Convergence was checked with rank-normalised R-hat and effective sample sizes [@vehtari2021].", ("nm_cfg", "mcmc_cfg"), kind="method")
    d.p("We checked the sampler with simulation-based calibration [@talts2018] on simulated data sets (S1 Text, Section 4).", kind="method")
    d.p("To validate the model, we checked that the 95% posterior predictive interval of every calibration target overlapped the target's range, and compared the modelled 2015 prison population, which was not a calibration target, with the number of people held in prison for drug consumption at the end of 2015. "
        "The model code was checked by unit and integration tests, including a comparison with an independent R implementation of the model equations. No observed series of new infections or HCV deaths among PWID in Singapore was available for external validation.", ("overlap", "prison", "strata"), kind="method")

    d.h(2, "Policy scenarios")
    d.p(f"In the main analysis, nobody was treated before 2026, so the model stayed at its 2015 equilibrium until then. From 2026, people with chronic infection in J were treated at a rate that gives a {pct(f.cfg['j_annual_start_probability'], 3)} annual probability of starting treatment, and community screening of F and X began. "
        f"We compared no treatment, prison treatment only (J only), {len(f.finite)} finite policies ({', '.join(f'{c:,}' for c in caps)} tests per year, with {shares[0]}%, {shares[1]}% or {shares[2]}% given to F and the rest to X), and uncapped screening of everyone in F and X each year as a reference.",
        ("config", "scen_code"), kind="method")
    d.p("Two comparators were used. No treatment, in which nobody is treated in any group, is the counterfactual for the infections and deaths averted by each strategy. Prison treatment only (J only) is the status quo comparator: it represents continuing the treatment offered in prison without community screening after release, and all incremental costs and outcomes of community screening are measured against it.", ("config", "scen_code", "cum"), kind="method")
    d.p(f"To show how much depends on prison treatment, we repeated all {len(s)} policies with an annual probability of starting treatment in J of {', '.join(pct(p, 0) for p in levels[:-1])} and {pct(levels[-1], 3)} ({len(levels) * len(s)} combinations).", ("config", "jcov", "scen_code"), kind="method")
    d.p("In a sensitivity analysis, prison treatment started in 2019 instead, at the main coverage, and community screening still started in 2026.", ("j2019", "scen_code"), kind="method")
    d.p(f"The smallest capacity, {caps[0]} tests per year, was set with reference to a public community testing programme [@moh2024ett]; the other capacities and shares are illustrative and are not announced policies.", ("config",), kind="method")
    d.p("At the start of each year, the tests given to F or X were applied to people not on treatment, in proportion to their HCV states. Each antibody-positive person had an RNA test, and the number of chronic infections found set that year's treatment capacity in the group. "
        f"Within that capacity, people with chronic infection started treatment at a rate that gives a {pct(f.cfg['fx_annual_start_probability'], 0)} annual probability. Sofosbuvir–velpatasvir without ribavirin was assumed, and people in D were never treated.", ("config", "scen_code", "setup"), kind="method")
    d.p("The care cascade was ideal in the main analysis: every chronic infection found by community screening was linked to treatment within the year's capacity, every person who started treatment completed it, and cure followed at the efficacy reported in trials. In a sensitivity analysis, only 80% or 60% of the chronic infections found in F and X were linked to treatment; treatment completion was not varied separately.", ("scen_code", "linkage", "setup"), kind="method")
    d.h(2, "Time horizon")
    d.p(f"The WHO targets were assessed in 2030. Because HCV-related cirrhosis and liver cancer develop over decades [@thein2008], we also projected every strategy to {E} and report cumulative outcomes over 2026–2030 and {H}; {H} was the main horizon for strategy selection. Entry of new injectors was held at its 2015 level to {E}, as before 2026, and the prison rates and model parameters were held constant.",
        ("config", "scen_code", "long", "setup"), kind="method")
    d.h(2, "Outcomes")
    d.p("The primary outcomes were the reductions in 2030, compared with 2015, of annual new HCV infections among PWID and of annual HCV deaths. "
        "New infections were the number of infections in D and F during the calendar year. HCV deaths were the liver deaths in decompensated cirrhosis and hepatocellular carcinoma during the year, the deaths from HCV-related cirrhosis and liver cancer that the WHO target counts. The 2015 values were the yearly numbers at the calibrated equilibrium. "
        "As a secondary outcome, we added the extrahepatic excess deaths of viraemic people (S1 Text, Section 10).", ("scen_code", "who2022strategy"), kind="method")
    d.p("A policy met a target if the median reduction was at least 80% for infections or 65% for deaths. We also report the share of posterior draws in which both targets were met.", ("scen_code", "config"), kind="method")
    d.p(f"For each posterior draw we also summed new infections and liver deaths over 2026–2030 and {H}, and subtracted them from the totals without treatment and with J only to give the infections and liver deaths averted. Averted outcomes were not discounted. Incident cases of decompensated cirrhosis and liver cancer were not counted separately.",
        ("scen_code", "cum"), kind="method")
    d.h(2, "Costs")
    d.p("Costs are the programme budget only: testing, consultations and direct-acting antiviral (DAA) drugs. They exclude the costs of caring for people with decompensated cirrhosis or liver cancer and any savings from the disease that treatment averts, so they are not the net costs of a strategy to the health system.", ("scen_code", "config"), kind="method")
    d.p(f"We costed antibody tests (S${f.cfg['unit_cost_sgd']['antibody']}), RNA tests (S${f.cfg['unit_cost_sgd']['rna']}), one consultation per treatment start (S${f.cfg['unit_cost_sgd']['consultation']}) and a 12-week sofosbuvir–velpatasvir course (S${F(f.costev[0]['base_cost_sgd']):,.0f}) in 2026 Singapore dollars, from published charges and a published price (S1 Text, Table H) [@wong2019; @nhg_dsc; @parkway_rna; @sgd_rate].",
        ("costev", "config"), kind="method")
    d.p(f"Programme budgets for 2026–2030 and {H} in all three groups (J, F and X) were discounted at {pct(f.cfg['discount_rate'], 0)} a year to 2026. Unit prices were varied together by a factor of 0.5 and 1.5 in a sensitivity analysis. These are charges, not costs to the health system. "
        f"The cost per infection or liver death averted was the incremental programme budget against J only divided by the outcomes averted against J only, both discounted at the same rate; it was not computed when fewer than one outcome was averted.", ("config", "costev", "scen_code", "cum"), kind="method")
    d.h(2, "Strategy selection")
    d.p(f"We compared strategies on five criteria at each horizon: three benefits (infections averted and liver deaths averted against no treatment, and the probability of meeting both WHO targets in 2030) and two costs (the discounted programme budget and the number of treatment courses). "
        f"We used two strategy sets: the {len(s)} policies of the main analysis, and the {len(s) * len(levels)} combinations of policy and prison treatment coverage. Each criterion was rescaled to 0–1 across the set, with the costs reverse-coded so that 1 was best.", DA + ("cum",), kind="method")
    d.p("We applied three decision rules: (a) the strategy with the largest health benefit whose programme budget was within a cap; (b) the strategy with the lowest cost per outcome averted against no treatment; and (c) TOPSIS among the Pareto-efficient strategies, ranking each by its relative closeness to the ideal (best on every criterion) and away from the anti-ideal (worst on every criterion) [@hwang1981]. "
        "We also drew the cost-effectiveness frontier for each benefit, flagging strategies that were dominated (another cost no more and gained no less) or extendedly dominated (a mix of two other strategies would do better) [@cantor1994].", DA, kind="method")
    d.p(f"The main analysis used posterior medians and equal weights. To carry the uncertainty through, we repeated Pareto efficiency, TOPSIS and the budget rule within each of the {f.ndraw:,} posterior draws and report the probability that each strategy was Pareto-efficient and that it was ranked first. "
        "We varied the weights by doubling or removing each criterion in turn, by using only the benefits or only the costs, and by drawing 1,000 weight vectors uniformly from the simplex.", DA + ("da_weights",), kind="method")
    d.h(2, "Uncertainty")
    d.p(f"For each of {f.ndraw:,} posterior draws we ran every scenario and computed every outcome. We report medians and 95% uncertainty intervals (UI; 2.5th to 97.5th percentiles) across draws, and the share of draws meeting a target. These describe uncertainty in the 13 calibrated parameters given the fixed inputs, not the probability that a target will be met. They do not include uncertainty in screening uptake, linkage to care, treatment completion, drug prices or model structure.", ("summary", "draws_valid"), kind="method")
    st_rho = f.setup("rho")
    d.p(f"To show how much the fixed inputs matter, we held the calibrated parameters at the posterior mode and set each fixed epidemiological input in turn to the low and high end of its literature range, without refitting (S1 Text, Section 13). In further sensitivity analyses, entry of new injectors from 2016 to 2025 followed the observed yearly counts of new drug abusers instead of the 2015 level (S1 Text, Section 11), and the share of infections with genotype 3 was 10% or 50% instead of {pct(st_rho, 0)} (S1 Text, Section 12). These analyses were run to 2030; the coverage and linkage analyses were run to {E} (S1 Text, Section 17).",
        ("tornado_code", "tornado_ranges", "inflow", "gt3", "scen_code", "linkage"), kind="method")
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
    d.p(f"At ages 30 and over, modelled first arrests were higher than the targets (for example {F(pred[2]['first_arrests_hat']):.0f} against {pred[2]['first_arrests_target']} at ages 30–39 at the posterior mode), because entry could not fall below the scaled counts and the entry scales settled near their lower bound (posterior medians {min(v[0] for v in f.entry_post.values()):.2f} to {max(v[0] for k, v in f.entry_post.items() if k < 12):.2f} below age 60).",
        ("nm_pred", "mcmc_par"), unc="moderate")
    d.p(f"The modelled number of active injectors in the community was {interval(df['median'], df['lower_95'], df['upper_95'])}, below the survey estimate of {f.tg['community_pwid']:,} but within its range. The first-arrest factor c was {f.c_post[0]:.2f} (95% credible interval {f.c_post[1]:.2f}–{f.c_post[2]:.2f}), so first arrest was faster than re-arrest.",
        ("strata", "mcmc_par", "targets"), unc="moderate")
    d.p(f"In 2015 the model had {n0(ever['median'])} people who had ever injected (95% credible interval {n0(ever['lower_95'])}–{n0(ever['upper_95'])}): {n0(f.strata['population_D']['median'])} in D, {n0(jpred['median'])} in J, {n0(f.strata['population_F']['median'])} in F and {n0(f.strata['population_X']['median'])} in X. "
        f"The modelled prison population was not a calibration target; it was close to the {f.j_held_2015:,} people held in prison for drug consumption at the end of 2015, the population whose release and re-arrest rates the model uses.", ("strata", "prison"), unc="moderate")
    d.p(f"In 2015 the model had {interval(b15['incidence'], b15['incidence_lo'], b15['incidence_hi'])} new infections and {interval(b15['mortality'], b15['mortality_lo'], b15['mortality_hi'])} HCV deaths per year. Without treatment, both stayed at their 2015 level until 2030 (within {100 * f.history_dev:.2f}% until 2025).",
        ("traj", "notx"), unc="moderate")

    d.h(2, "Prison treatment alone")
    d.p(f"With treatment in prison from 2026, {n0(jstarts26)} people started treatment in J in 2026, as the infections accumulated since 2015 were treated, and fewer in later years. "
        f"New infections fell by {red(J, 'incidence_reduction')} by 2030 and HCV deaths by {red(J, 'mortality_reduction')} (Fig 3). No posterior draw met either target.",
        ("annual",) + TRJ + OUT, unc="moderate")
    d.p(f"In 2026, {xs} of the {n0(f.x_chronic + f.f_chronic)} people with chronic infection in F and X were in X, and F and X together accounted for {fx26} of all HCV deaths. People in X are never re-arrested in the model, so prison treatment does not reach them.",
        ("annual", "traj", "sim"), unc="moderate")
    d.p(f"People in D, who had never been imprisoned, accounted for {pct(f.d_inc_share_j, 0)} of new infections in 2030 under prison treatment only; no policy in this study reached them.", OUT, unc="moderate")

    d.h(2, "Community screening policies")
    d.p(f"{cap_word(len(f.inc_pass))} of the {len(f.finite)} finite screening policies met the incidence target at the median (Fig 4 and Table 1): {caps[-1]:,} tests a year with any share to F, and {caps[-2]:,} tests with 50% or 75% to F. None met the mortality target. "
        f"The largest fall in deaths among finite policies was {red(P8, 'mortality_reduction')}, with {caps[-1]:,} tests and 25% to F, and no finite policy met both targets in any posterior draw.",
        OUT, unc="moderate")
    d.p(f"At every capacity, giving a larger share of tests to F increased the fall in new infections and reduced the fall in deaths. All quotas were used in full.",
        ("summary", "annual"), unc="moderate")
    d.p(f"At the start of 2026, {n0(f.f_eligible)} people in F were eligible for screening, of whom {fchr} had chronic infection; in X, {n0(f.x_eligible)} were eligible and {xchr} had chronic infection. "
        "A test given to F therefore found more chronic infections, in people who could transmit, but people in F are re-arrested and can be treated in prison, whereas people in X can be reached only by community screening; this is a likely reason why tests given to X prevented more deaths.", ("annual", "rates", "sim"), kind="interpretive", unc="moderate")
    d.p(f"Uncapped screening of everyone in F and X each year reduced new infections by {red(INF, 'incidence_reduction')} and deaths by {red(INF, 'mortality_reduction')}, and was the only policy that met both targets. Even then, {pct(f.d_inc_share_sinf, 0)} of the remaining infections in 2030 were in D. "
        f"In the one-way analysis of the fixed inputs, the only input that stopped it meeting the mortality target was {sinf_flip}: at the less favourable end of its range, deaths fell by {pct(sinf_psi[0])} (S1 Text, Section 13).", OUT + ("tornado_flips", "tornado"), unc="moderate")

    d.h(2, "Prison treatment coverage")
    d.p(f"Every reduction rose with the coverage of prison treatment (Fig 5 and Table 2). Without prison treatment, deaths fell by {red(J0['S8000_F25'], 'mortality_reduction')} with {caps[-1]:,} tests and 25% to F, and only {caps[-1]:,} tests with 75% to F met the incidence target, reducing new infections by {red(J0['S8000_F75'], 'incidence_reduction')}. "
        f"At {pct(levels[1], 0)} coverage, {cap_word(sum(f.jcov[levels[1]][k]['who_incidence_pass'] == 'TRUE' for k in f.finite)).lower()} finite policies met the incidence target, and at {pct(levels[2], 0)} the same {cap_word(len(f.inc_pass)).lower()} as at full coverage.",
        JC, unc="moderate")
    d.p(f"At no coverage did a finite policy meet the mortality target. Uncapped screening met it at every coverage, though only just when there was no prison treatment, reducing deaths by {red(J0['SInf'], 'mortality_reduction')}.", JC, unc="moderate")

    d.h(2, "Sensitivity analysis: prison treatment from 2019")
    d.p(f"Had prison treatment started in 2019, prison treatment alone would have reduced new infections by {red(SJ, 'incidence_reduction')} and deaths by {red(SJ, 'mortality_reduction')} by 2030, meeting the incidence target at the median. "
        f"All {len(f.finite)} finite policies would then have met the incidence target, but none would have met the mortality target: the largest fall in deaths was {red(S8, 'mortality_reduction')}, with {caps[-1]:,} tests and 25% to F. Only uncapped screening met both targets, reducing deaths by {red(sens['SInf'], 'mortality_reduction')} (S1 Text, Section 9).",
        SENS, unc="moderate")

    d.h(2, "Programme budget, 2026–2030")
    d.p(f"Discounted programme budgets for 2026–2030 were {money(J, 'total_discounted_cost_2026_2030_sgd')} with prison treatment only, mostly for drugs and mostly in 2026. The cheapest policy that met the incidence target, {int(F(s[inc_cheapest]['screening_total'])):,} tests a year with {round(100 * F(s[inc_cheapest]['f_share']))}% to F, cost {money(s[inc_cheapest], 'total_discounted_cost_2026_2030_sgd')}, or {money(s[inc_cheapest], 'incremental_discounted_cost_vs_j_only_sgd')} more than prison treatment only.",
        OUT + ("annual",), unc="moderate")
    d.p(f"Uncapped screening, the only policy that met both targets, cost {money(INF, 'total_discounted_cost_2026_2030_sgd')}. Drugs made up {pct(f.drug_share_5000, 0)} of the cost of 5,000 tests a year with 25% to F, so totals scale almost directly with the drug price; S1 Text, Section 8 gives all policies.", ("annual", "summary"), unc="moderate")

    sh = lambda k: 1 - F(C[k]["liver_deaths_averted_vs_no_treatment_2026_2030"]) / F(C[k][f"liver_deaths_averted_vs_no_treatment_2026_{E}"])
    yv = lambda k: yr[(k, "mortality_ratio_2015")]
    d.h(2, f"Long-term outcomes, {H}")
    d.p(f"Liver deaths kept falling after 2030 under every strategy (Fig 6). The median yearly number of liver deaths fell to 35% of the 2015 level in {yv('SInf')} with uncapped screening, {yv('S8000_F25')} with {caps[-1]:,} tests and 25% to F, {yv('S8000_F75')} with {caps[-1]:,} tests and 75% to F, and {yv('J_only')} with prison treatment only; new infections under prison treatment only fell to 20% of the 2015 level in {yr[('J_only', 'incidence_ratio_2015')]}.",
        LONG, unc="moderate")
    d.p(f"Against no treatment, prison treatment only averted {cu('J_only', 'infections_averted_vs_no_treatment')} new infections and {cu('J_only', 'liver_deaths_averted_vs_no_treatment')} liver deaths over {H}; {pct(sh('J_only'), 0)} of these liver deaths averted occurred after 2030. "
        f"Against prison treatment only, uncapped screening averted a further {cu('SInf', 'liver_deaths_averted_vs_j_only')} liver deaths, {caps[-1]:,} tests with 25% to F {cu('S8000_F25', 'liver_deaths_averted_vs_j_only')} and with 75% to F {cu('S8000_F75', 'liver_deaths_averted_vs_j_only')} (S1 Text, Section 15).",
        LONG, unc="moderate")
    d.p(f"At every capacity, giving 25% of tests to F averted more liver deaths over {H} than giving 75%, while giving 75% averted more infections and needed fewer treatment courses, because infections prevented in people who still inject do not need treatment later.",
        LONG, unc="moderate")

    d.h(2, "Cost per outcome averted")
    d.p(f"Over {H}, the extra programme budget against prison treatment only was {money(C['S8000_F75'], f'incremental_programme_cost_vs_j_only_sgd_2026_{E}')} with {caps[-1]:,} tests and 75% to F, {money(C['S8000_F25'], f'incremental_programme_cost_vs_j_only_sgd_2026_{E}')} with 25% to F and {money(C['SInf'], f'incremental_programme_cost_vs_j_only_sgd_2026_{E}')} with uncapped screening. "
        f"The cost per liver death averted was {kS(C['S8000_F75'], f'cost_per_liver_death_averted_vs_j_only_sgd_2026_{E}')}, {kS(C['S8000_F25'], f'cost_per_liver_death_averted_vs_j_only_sgd_2026_{E}')} and {kS(C['SInf'], f'cost_per_liver_death_averted_vs_j_only_sgd_2026_{E}')}, respectively.",
        LONG, unc="moderate")
    low_caps = [f"{n:,}" for n in caps[:3]]
    d.p(f"With {', '.join(low_caps[:-1])} or {low_caps[-1]} tests a year and 75% to F, the programme budget over {H} was about the same as with prison treatment only (95% UI of the difference including zero), because the infections they prevented did not need treatment later; these policies averted more infections and liver deaths at no clear extra cost.",
        LONG, unc="moderate")

    d.h(2, "Strategy selection")
    d.p(f"Over {H}, the cost-effectiveness frontier for liver deaths averted ran, in order of cost, through {'; '.join(lab17(k).replace('J + ', '') for k in fr_ld)} (Fig 7). Prison treatment only was dominated: {lab17(fr_ld[0]).replace('J + ', '')} cost no more and averted more infections and liver deaths. "
        f"Every policy giving 25% of tests to F was dominated or extendedly dominated; for infections averted the frontier ran through the 75%-to-F policies from {caps[1]:,} tests and uncapped screening. By rule (b), the lowest cost per outcome averted against no treatment was with {lab17('S5000_F75').replace('J + ', '')} for both infections and liver deaths.",
        DA, unc="moderate")
    d.p(f"{cap_word(len(par17))} of the {len(s)} strategies were Pareto-efficient at the median. With equal weights, TOPSIS ranked {lab17(top17)} first (closeness {F(d17[top17]['topsis_closeness']):.3f}; first in {pct(d17[top17]['p_topsis_top'], 1)} of posterior draws), followed by {lab17(par17[1])} ({F(d17[par17[1]]['topsis_closeness']):.3f}) and {lab17(par17[2])} ({F(d17[par17[2]]['topsis_closeness']):.3f}) (Fig 8 and Table 3). "
        f"Over 2026–2030 alone, uncapped screening, the only strategy that met both targets, ranked first.",
        DA, unc="moderate")
    d.p(f"The ranking depended on the weights. Doubling the weight on liver deaths averted or on meeting both targets made uncapped screening first; doubling the weight on the budget or on treatment courses made {lab17(wts['double_programme_cost']).replace('J + ', '')} first. "
        f"Across 1,000 random weight vectors, the strategy ranked first most often was {lab17(rw[0])} ({pct(d17[rw[0]]['p_top_random_weights'], 0)}), then {lab17(rw[1])} ({pct(d17[rw[1]]['p_top_random_weights'], 0)}), {lab17(rw[2])} ({pct(d17[rw[2]]['p_top_random_weights'], 0)}) and {lab17(rw[3])} ({pct(d17[rw[3]]['p_top_random_weights'], 0)}).",
        DA + ("da_weights",), unc="moderate")
    runs = []  # consecutive caps with the same best strategy are merged into one range
    for x in budget:
        m = F(x["budget_cap_sgd"]) / 1e6
        if runs and runs[-1][0] == x["best_strategy"]: runs[-1][2] = m
        else: runs.append([x["best_strategy"], m, m])
    caps_txt = "; ".join(f"{lab17(k).replace('J + ', '')} with a cap of S${a:,.0f}" + (f"–{b:,.0f}" if b > a else "") + " million" for k, a, b in runs)
    d.p(f"With a cap on the programme budget over {H}, the strategy averting the most liver deaths was {caps_txt} (rule a; S1 Text, Section 16).",
        DA + ("da_budget",), unc="moderate")
    d.p(f"When prison treatment coverage was also a choice, TOPSIS ranked uncapped screening with full prison coverage first in {pct(d68[top68]['p_topsis_top'], 0)} of draws, and the top four were uncapped screening at each prison coverage. For every screening policy with more than one Pareto-efficient coverage level, TOPSIS closeness rose with prison coverage.",
        DA, unc="moderate")

    d.h(2, "Linkage to care")
    d.p(f"If only 80% or 60% of the chronic infections found in F and X were linked to treatment, every finite policy had smaller effects. With 60% linkage, {caps[-1]:,} tests and 25% to F reduced liver deaths in 2030 by {red(L6['S8000_F25'], 'mortality_reduction')} (main analysis {pct(P8['mortality_reduction'])}) and averted {col(LC6['S8000_F25'], f'liver_deaths_averted_vs_j_only_2026_{E}', d=0)} liver deaths over {H} against prison treatment only. "
        f"Uncapped screening still met both targets at 80% and 60% linkage (fall in liver deaths {red(L8['SInf'], 'mortality_reduction')} and {red(L6['SInf'], 'mortality_reduction')}), because people missed in one year were screened again the next (S1 Text, Section 17).",
        LK + OUT, unc="moderate")

    # ---------------- discussion ----------------
    d.h(1, "Discussion")
    d.p(f"In this model of PWID in Singapore, starting HCV treatment in 2026 left too little time to meet the WHO 2030 mortality target with prison treatment and capped community screening. "
        f"Prison treatment from 2026 reduced new infections by about half and deaths by about a fifth. Up to {caps[-1]:,} community tests a year could meet the incidence target, but only screening of every released person each year met the mortality target, and that depended on how much cure lowers the risk of liver death in decompensated cirrhosis.",
        OUT, kind="interpretive", unc="moderate")
    d.p("Timing mattered most for new infections. With the same prison treatment starting in 2019, prison treatment alone met the incidence target, but even then no capped policy met the mortality target. Each year without treatment adds to the pool of chronic infections that progress to cirrhosis and liver cancer. Cure lowers but does not remove the risk of liver death in people who already have decompensated cirrhosis or liver cancer, and people with compensated cirrhosis can still progress, so liver deaths fall much more slowly than infections after treatment begins.",
        SENS + OUT, kind="interpretive", unc="moderate")
    d.p("The two WHO targets pull screening in different directions. Tests given to people who still inject prevent more infections, while tests given to people who have stopped prevent more deaths. Because the mortality target was harder to meet in every analysis, a policy aimed at the 2030 targets should be judged mainly by the deaths it prevents.",
        OUT + SENS, kind="interpretive", unc="moderate")
    d.p(f"Over a longer horizon the choice looks different. Most of the liver deaths that treatment averts occur after 2030, and preventing infections in people who still inject also saves the treatment they would later need. Over {H}, policies giving most community tests to F lay on the cost-effectiveness frontier, the smallest of them at no clear extra budget over prison treatment only, and {caps[-1]:,} tests a year with 75% to F was ranked first with equal weights in almost every posterior draw. "
        "Giving more tests to X averted more liver deaths but cost more, because it found and treated more long-standing infections; uncapped screening averted the most liver deaths and was the only strategy to meet the 2030 targets, but at the highest budget.",
        LONG + DA, kind="interpretive", unc="moderate")
    d.p(f"In terms of whom, where and how intensively, the analysis gives a consistent picture. Whom: people released from prison, both those who still inject (to prevent infections) and those who have stopped (to prevent liver deaths); people who have never been imprisoned were not reached by any strategy. Where: prison treatment is the base on which community screening builds, and higher prison coverage was preferred in every comparison. How intensively: several thousand community tests a year, well above the {caps[0]} a year of the smallest policy, which was set with reference to a public programme, and screening every released person if the 2030 mortality target is to be met. "
        "Which strategy is preferred depends on the weight placed on meeting the 2030 targets and on averting liver deaths relative to the budget, and the decision analysis makes that dependence explicit rather than resolving it.",
        OUT + DA + JC, kind="interpretive", unc="moderate")
    d.p("This suggests that the people who matter most for HCV mortality are not those most visible to drug services. People who have stopped injecting may no longer be in contact with prison or drug services, and reaching them may need testing in primary care or other community settings [@hsiang2022].", kind="interpretive")
    d.p(f"We counted only liver deaths, as the WHO target does. Viraemic infection may also raise mortality from other causes [@lee2012; @cepeda2022], although the excess is not found in every cohort [@piselli2021]. Cure removes this excess at once, so counting it as well gave larger falls in deaths (prison treatment only {pct(J['hcv_mortality_reduction'])}, uncapped screening {pct(INF['hcv_mortality_reduction'])}; S1 Text, Section 10); still, no finite policy reached 65%.",
        OUT, kind="interpretive", unc="moderate")
    d.p("We extend earlier modelling for Singapore [@chaillon2021] by separating released people by injecting status and by modelling screening capacity explicitly. This separation gives a reason why the mortality target may be harder to meet than the incidence target.", kind="interpretive")
    d.h(2, "Strengths and limitations")
    d.p(f"The model was calibrated with a full Bayesian method that converged well, and every policy result carries posterior uncertainty from {f.ndraw:,} draws. Prison flows were estimated from individual prison records rather than assumed, and screening capacity was modelled explicitly, so treatment depends on how many people are found.",
        ("mcmc_diag", "summary", "rates"), kind="interpretive", unc="moderate")
    d.p(f"The study has important limitations. First, the seroprevalence targets came from the {f.j_obs_2015:,} people in prison at the end of 2015 who had a screening result, out of the {f.j_held_2015:,} held for drug consumption. The modelled prison population ({n0(jpred['median'])}) matched the number held, although it was not fitted, but the people screened may not represent everyone in prison.",
        ("strata", "prison"), kind="interpretive", unc="high")
    d.p(f"Second, the main analysis assumes that nobody was treated before 2026. Treatment in prison has in fact been offered [@prison2021; @tan2023], but its start and coverage are not recorded in the project data. The sensitivity analysis with full coverage from 2019 shows how much an earlier start changes the results; the real history probably lies between the two analyses.",
        ("j2019", "config"), kind="interpretive", unc="high")
    d.p(f"Third, we did not observe who injects. We treated new drug abusers whose main drug was heroin or methamphetamine ({pct(f.pwid_share)} of {f.new_total:,} in 2015) as the people who start injecting, at the same share in every age group because the drug table has no age breakdown. The records give the main drug, not the route of use: users of these drugs who smoke rather than inject inflate entry, and people who inject other drugs are missed. "
        "Entry could rise above this share through the calibrated entry scales but could not fall below it, and it was held at its 2015 level in every year; a fall in new injectors after 2015 would lower infections without treatment, and a rise would raise them.",
        ("cnb", "drug", "setup", "likelihood", "sim"), kind="interpretive", unc="high")
    d.p(f"Fourth, the entry and first-arrest parameters were weakly identified. The entry scales settled near their lower bound, and the posterior of c was only modestly narrower than its prior (standard deviation {pct(f.sd_ratio[12], 0)} of the prior's), as was the contact multiplier for ages 50–59. "
        "The first-arrest rate was assumed to have the same age pattern as the re-arrest rate, because the records contain only people who had already been arrested.",
        ("mcmc_par", "mcmc_cfg", "setup"), kind="interpretive", unc="high")
    d.p("Fifth, we assumed that released people who are eventually re-arrested return to injecting and that the others stop permanently and are never re-arrested, and that there is no transmission in prison. Relapse after a long gap would put more infections back into F.", ("setup", "sim"), kind="interpretive", unc="high")
    d.p("Sixth, screening was an aggregate approximation: tests were spread across each group in proportion to its HCV states, with perfect linkage within capacity and no loss to follow-up, and people who had never been imprisoned were never screened or treated. Real programmes will find and treat fewer people for the same number of tests. "
        "The linkage sensitivity analysis shows that lower linkage reduces the effect of every finite policy; treatment completion and reinfection after cure in the community were not varied separately.", ("scen_code", "linkage"), kind="interpretive", unc="high")
    d.p("Seventh, the model assumed that the epidemic was in a steady state in 2015, with constant prison rates, mixing that depends only on age, and the same mortality ratios for people who have stopped injecting as for those who still inject. Costs were published charges and a published drug price carried forward without inflation, not payer costs, and follow-up visits were excluded. Our uncertainty intervals cover the calibrated parameters only, not the fixed inputs, structural assumptions, prices or uptake, so they are too narrow. "
        f"The one-way analysis shows that the fall in deaths depends most on the liver disease inputs, such as the risk of liver death after cure in decompensated cirrhosis (with uncapped screening, {pct(sinf_psi[0])} to {pct(sinf_psi[1])} across its range), which come from cohorts outside Singapore. S1 Text, Section 14 (Table O) lists every assumption with its basis and its likely effect on the results.",
        ("sim", "equilibrium", "setup", "costev", "config", "tornado"), kind="interpretive", unc="high")
    d.p(f"Eighth, the projection to {E} holds entry of new injectors, prison rates, prices and treatment efficacy at their current values for 25 years, so the long-horizon results are scenarios rather than forecasts. "
        "The programme budget excludes the cost of caring for advanced liver disease and so understates the value of strategies that avert liver deaths. "
        "The decision analysis depends on the criteria and weights chosen: infections averted and treatment courses are correlated, the benefit criteria were not discounted, and TOPSIS ranks can change when strategies are added to or removed from the set, as the different first choices in the two strategy sets show. We did not value outcomes in quality-adjusted life years, by design.",
        LONG + DA + ("da_weights",), kind="interpretive", unc="high")
    d.p(f"Finally, the seroprevalence targets could not be reproduced exactly from the current processing of the prison records (largest difference {100 * f.sero_max_diff:.1f} percentage points).",
        ("targets", "prison"), kind="result", unc="moderate")
    d.h(2, "Conclusions")
    d.p("In this model, treatment for PWID in Singapore that starts in 2026 can meet the 2030 target for new HCV infections only with prison treatment and thousands of community tests a year, and cannot meet the target for HCV deaths without screening every released person; even then the result depends on liver disease inputs from other settings. "
        "Reducing deaths requires finding and treating people who have stopped injecting after release, and every year of delay makes the target harder to reach. "
        f"Judged over {H} on infections and liver deaths averted, the WHO targets and the programme budget together, prison treatment combined with several thousand community tests a year, most of them for released people who still inject, was the preferred strategy unless meeting the 2030 mortality target is given priority.",
        OUT + SENS + DA, kind="interpretive", unc="moderate")

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
    d.p("**Fig 1. Model structure.** People who inject drugs enter D (never imprisoned) and move to J (prison) at first arrest. On release they move to F (still injecting) or X (stopped). People in F can be re-arrested. Transmission occurs only in D and F. In the main analysis, treatment is given in J and community screening is offered to F and X from 2026.", ("sim", "config"), kind="method")
    d.p("**Fig 2. Calibration to the 2015 targets.** (A) HCV antibody prevalence in prison by age. (B) First arrests by age (Central Narcotics Bureau new drug abusers, scaled to heroin and methamphetamine users). (C) Active injectors in the community (D and F). Black: target with the 95% range of its likelihood; blue: 95% posterior interval of the model value from the MCMC posterior.", ("overlap",), kind="result")
    d.p(f"**Fig 3. New infections and HCV deaths, 2015–2030.** Annual new HCV infections (A) and HCV deaths (B) relative to 2015, without treatment, with prison treatment only, and with {caps[-1]:,} community tests a year given 25% or 75% to F, all from 2026. Yearly counts are drawn at mid-year and joined by a monotone curve; the vertical line marks 2026. Lines are medians and bands are 95% uncertainty intervals across {f.ndraw:,} posterior draws. Dashed lines mark the WHO targets.", OUT + TRJ, kind="result")
    d.p("**Fig 4. Reductions in 2030 by screening capacity and share of tests given to F.** (A) New infections. (B) HCV deaths. Points are medians and bars are 95% uncertainty intervals. Dashed lines mark the targets; the dotted line is prison treatment only and the solid line is uncapped screening.", OUT, kind="result")
    d.p("**Fig 5. Reductions in 2030 by prison treatment coverage.** Median reductions in new infections (A) and HCV deaths (B) compared with 2015, for each policy (rows) and annual probability of starting treatment in prison from 2026 (columns). Bold values meet the target.", JC, kind="result")
    d.p(f"**Fig 6. New infections and liver deaths, 2015–{E}.** Annual new HCV infections (A) and liver deaths from decompensated cirrhosis and hepatocellular carcinoma (B) relative to 2015, for no treatment, prison treatment only, {caps[-1]:,} community tests a year with 25% or 75% to F, and uncapped screening. Lines are medians and bands are 95% uncertainty intervals across {f.ndraw:,} posterior draws. Dashed lines mark the WHO 2030 targets; dotted lines mark 2026 and 2030. Entry of new injectors is held at its 2015 level throughout.", LONG, kind="result")
    d.p(f"**Fig 7. Cost-effectiveness frontiers, {H}.** Discounted programme budget (testing, consultations and drugs; 2026 Singapore dollars) against new infections averted (A) and liver deaths averted (B) compared with no treatment, for the {len(s)} policies (posterior medians). Numbered points are on the frontier; squares are extendedly dominated and crosses dominated. Costs exclude care for advanced liver disease.", DA, kind="result")
    d.p(f"**Fig 8. TOPSIS ranking of the Pareto-efficient strategies, {H}.** Relative closeness to the ideal strategy with equal weights on infections averted, liver deaths averted, the probability of meeting both WHO targets, the programme budget and treatment courses. (A) The {len(s)} policies at full prison treatment coverage. (B) The {len(s) * len(levels)} combinations of policy and prison treatment coverage (at most 15 shown). Labels give the share of posterior draws in which each strategy ranked first and was Pareto-efficient.", DA, kind="result")

    # ---------------- tables ----------------
    d.h(1, "Tables")
    d.p("**Table 1. Reductions in 2030 compared with 2015, probability of meeting both WHO targets, and discounted costs for 2026–2030 (main analysis).**", OUT, kind="result")
    rows = []
    for k in ["J_only"] + [f"S{c}_F{x}" for c in caps for x in shares] + ["SInf"]:
        r = s[k]
        rows.append([scen_label(r), col(r, "incidence_reduction", scale=100, d=1, suf="%").replace(" (95% UI ", " ("),
                     col(r, "mortality_reduction", scale=100, d=1, suf="%").replace(" (95% UI ", " ("),
                     pct(r["joint_pass_probability"], 1), col(r, "total_discounted_cost_2026_2030_sgd", scale=1e-6, d=1).replace(" (95% UI ", " (")])
    d.table(["Policy", "Fall in new infections, % (95% UI)", "Fall in HCV deaths, % (95% UI)", "Draws meeting both targets", "Cost, S$ million (95% UI)"], rows, OUT)
    d.p("UI, uncertainty interval across posterior draws. Costs are 2026 Singapore dollars for J, F and X, discounted at 3% a year to 2026.", ("config",), kind="method")
    d.p("**Table 2. Median falls in new infections and HCV deaths in 2030 by annual probability of starting treatment in prison from 2026.**", JC, kind="result")
    rows = []
    for k in ["J_only"] + [f"S{c}_F{x}" for c in caps for x in shares] + ["SInf"]:
        rows.append([scen_label(s[k])] + [f"{pct(f.jcov[p][k]['incidence_reduction'])} / {pct(f.jcov[p][k]['mortality_reduction'])}" for p in levels])
    d.table(["Policy"] + [f"{pct(p, 0) if p < 0.999 else pct(p, 3)} in prison" for p in levels], rows, JC)
    d.p("Each cell gives the fall in new infections / the fall in HCV deaths (medians across posterior draws). The last column is the main analysis (Table 1). Targets: 80% and 65%.", JC, kind="method")
    d.p(f"**Table 3. Pareto-efficient strategies over {H}: decision criteria, TOPSIS closeness with equal weights and probabilities across posterior draws.**", DA, kind="result")
    rows = []
    for k in par17:
        x = d17[k]
        rows.append([lab17(k), interval(x["infections_averted"], x["infections_averted_lo"], x["infections_averted_hi"]),
                     interval(x["liver_deaths_averted"], x["liver_deaths_averted_lo"], x["liver_deaths_averted_hi"]), pct(x["who_joint_pass"], 1),
                     interval(x["programme_cost"], x["programme_cost_lo"], x["programme_cost_hi"], d=1, scale=1e-6),
                     interval(x["treatment_courses"], x["treatment_courses_lo"], x["treatment_courses_hi"]),
                     f"{F(x['topsis_closeness']):.3f}", pct(x["p_pareto_efficient"], 1), pct(x["p_topsis_top"], 1)])
    d.table(["Policy", "Infections averted", "Liver deaths averted", "Draws meeting both 2030 targets", "Programme budget, S$ million", "Treatment courses", "TOPSIS closeness", "Pareto-efficient, % of draws", "Ranked first, % of draws"], rows, DA)
    d.p(f"Medians (95% UI) across posterior draws. Outcomes averted are against no treatment, summed over {H} and not discounted. The programme budget covers testing, consultations and drugs in J, F and X, discounted at 3% a year to 2026; it excludes care for advanced liver disease. Strategies are ordered by TOPSIS closeness.", DA + ("config",), kind="method")

    # ---------------- supporting information ----------------
    d.h(1, "Supporting information")
    d.p("**S1 Text. Technical appendix.** Model equations, parameters, prison rate estimation, calibration and its checks, the 2015 model population, scenario implementation, costs, additional results and the sensitivity analysis with prison treatment from 2019.", ("sim",), kind="declaration")
    d.p("**S1 Checklist. CHEERS 2022 checklist.**", ("cheers2022",), kind="declaration")
    return d
