# Key assumptions (state of 2026-11-12)

Summary of the assumptions behind the HCV model of people who inject drugs (PWID) in Singapore, the
policy scenarios and the decision analysis. The full table with basis and likely effect is S1 Text,
Section 14 (Table O); limitations are in the main text Discussion. Numbers are from the current outputs.

## Population and entry
- **Who injects.** New drug abusers whose main drug was heroin or methamphetamine are treated as the
  people who start injecting: 86.1% of the 1,309 new abusers in 2015 (113 heroin + 1,014
  methamphetamine; SingStat M891521). The table records the main drug, not the route, so smokers of
  these drugs inflate entry and injectors of other drugs are missed. Net direction unknown.
- **Same share at every age.** The drug table has no age breakdown.
- **Entry is constant.** Entry into D = calibrated scale (≥ 1) × 0.861 × CNB 2015 new abusers by age,
  held at the 2015 level in every year, before 2015 and to 2050. The observed-inflow sensitivity
  (2016–2025 CNB series) holds the 2025 level afterwards (run to 2030 only).
- **Steady state in 2015.** 150 model years with constant parameters; 2016–2025 stays at the 2015
  level without treatment (within 0.09%).

## Prison and movement between groups
- Groups: D (never imprisoned, injecting), J (prison), F (released, injecting), X (released, stopped).
- **First arrest** = c × re-arrest rate, same age pattern (records only contain arrested people);
  c is calibrated and weakly identified (posterior SD 72% of the prior's).
- **Release.** People eventually re-arrested (π ≈ 72%) go to F; the others stop permanently (X) and
  are never re-arrested. No relapse from X.
- Prison rates constant over time and the same in every HCV state.
- J represents everyone held for drug consumption (8,725 at end 2015); its seroprevalence is that of
  the 3,367 screened. The modelled J (8,768) was not fitted.

## Transmission and natural history
- Transmission only in D and F; none in prison or after stopping. Age-only mixing with a calibrated
  multiplier per age group; reinfection at the same rate as first infection.
- Progression, clearance and treatment efficacy from the literature with a fixed genotype 3 share
  (varied 10% / 50% in a sensitivity analysis); not in the uncertainty intervals.
- **Mortality.** Non-liver deaths = 2015 Singapore rate × SMR 4.73, × 1.47 more for viraemic people
  (extrahepatic HR, general-population cohort). Liver deaths only in DC and HCC; after cure (state r
  only) multiplied by ψ_DC 0.40 and ψ_HCC 0.34. Liver deaths are the WHO outcome.

## Policies and care cascade
- **Main analysis: nobody treated before 2026.** Prison treatment and community screening both start
  in 2026 (sensitivity: prison treatment from 2019).
- **Status quo comparator = J only** (prison treatment, annual start probability 99.995%, no community
  screening). **No treatment** is the counterfactual for outcomes averted.
- **Ideal cascade.** Every chronic infection found by community screening is linked to treatment within
  the year's capacity; everyone who starts completes; cure at trial efficacy. Sensitivity: 80% / 60%
  linkage (completion not varied separately).
- Screening is an aggregate approximation: tests spread over each group in proportion to its HCV
  states, repeat screening allowed, unused quota not moved between F and X. D is never screened or
  treated.
- Capacities (600–8,000 tests a year) and F/X shares (25/50/75%) are illustrative; 600 was set with
  reference to a public programme.

## Outcomes, horizon and costs
- WHO targets assessed in **2030** against the modelled 2015 yearly numbers: new infections −80%,
  liver (DC/HCC) deaths −65%. A target is met if the median reduction reaches it.
- **Horizon to 2050** because liver deaths lag infection by decades. Entry, prison rates, prices and
  treatment efficacy are held at current values for 25 years, so these are scenarios, not forecasts.
- **Costs = programme budget only:** antibody and RNA tests, one consultation per treatment start and a
  12-week sofosbuvir–velpatasvir course, at published charges and a published price (2026 SGD, not
  inflated), discounted 3% a year to 2026. They exclude care for decompensated cirrhosis and liver cancer
  and savings from disease averted, follow-up visits, tests of cure and navigation.
- Outcomes averted are not discounted; cost per outcome averted uses outcomes discounted at the same rate.
- No QALYs or DALYs (by design).

## Decision analysis
- Five criteria with equal weights in the main analysis: infections averted, liver deaths averted,
  P(meeting both 2030 targets) (benefits); programme budget, treatment courses (costs).
- Min–max normalisation within each strategy set; TOPSIS among Pareto-efficient strategies.
- Main decision horizon 2026–2050. Ranks depend on the weights, on correlated criteria (infections
  averted and courses) and on the set of strategies compared (rank reversal).

## Uncertainty
- Intervals reflect posterior uncertainty in the 13 calibrated parameters only (2,000 reference-MCMC
  draws; the model is deterministic). Not included: fixed epidemiological inputs (one-way analysis
  instead), model structure, prices, screening uptake, linkage, completion. Intervals are too narrow.
