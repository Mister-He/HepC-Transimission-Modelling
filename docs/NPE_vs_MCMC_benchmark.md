# NPE / SNPE versus MCMC: Agent Benchmark Protocol

## 1. Objective and scope

Determine which inference method produces sufficiently accurate parameter posteriors and policy estimates at lower end-to-end cost for this project. Do not assume that NPE outperforms MCMC. The user reports inexpensive ODE evaluations and is interested in inference for one fixed observed dataset.

This document specifies an executable benchmark plan, not verified results. Inspect the current repository before implementing it. Read AGENTS.md, preserve user changes, and write isolated benchmark scripts and outputs without overwriting production posterior files. Do not change the scientific model to improve benchmark results.

Use two complementary comparisons:

1. **Approximately matched forward-model budget:** around 26,000 forward-evaluation opportunities per method and replicate. This measures simulation efficiency, not equal total computation.
2. **Matched end-to-end wall-clock budget:** the same hardware and the measured total SNPE runtime per replicate. This is the primary practical comparison because neural-network training is not free.

Never equate MCMC iterations with ODE calls without checking the sampler. HMC/NUTS, slice sampling, delayed rejection, numerical gradients, and initialization searches can each use multiple forward evaluations per iteration.

## 2. Fixed shared settings

| Setting | Value |
| --- | --- |
| Full independent replicates | 3 |
| Base seeds | 20260929, 20260930, 20260931 |
| Maximum CPU resources | 4 physical cores, or 4 available cores if physical-core control is unavailable |
| Primary hardware track | CPU only; no GPU/MPS |
| Concurrent methods | None; run methods sequentially to avoid contention |
| Simulation workers | 4, one numerical-library thread per worker |
| MCMC workers | 4 independent chains, one thread per chain |
| Neural training | One process, up to 4 CPU threads; no concurrent simulation workers |
| Numerical precision | Preserve existing ODE precision and tolerance; neural training float32 unless existing compatibility requires otherwise |
| NPE exported posterior draws | 20,000 per replicate |
| MCMC thinning | 1: retain all post-warmup states, including repeated states after rejection |
| Predictive evaluation subset | 1,000 posterior states per method and replicate, chain-balanced for MCMC |
| Initial audit timing sample | 128 prior parameter draws; optionally 128 high-density-region draws if already available |

The seed numbers are integers, not calendar dates. Derive separate reproducible RNG streams for simulation, network initialization, validation split, posterior sampling, chain initialization, and each chain. Same integer seeds do not imply identical random experiments across methods. Repeat NPE simulation and training, not merely posterior sampling.

Pin package versions, hardware details, thread settings, Git commit, input hashes, and configuration. Record any unavoidable hardware deviation. A separate GPU track is optional and must not be pooled with CPU-only results.

## 3. Statistical-target audit: a mandatory gate

Write model_audit.md documenting actual source files and functions for:

- Parameter dimension, units, prior, support, fixed parameters, and any joint prior dependence.
- Parameter transformations and Jacobians. Evaluate all reported posteriors in the same original parameter coordinates.
- ODE initial conditions, solver, tolerances, integration horizon, equilibrium rules, and numerical-failure handling.
- Observed data vector: order, denominators, years, missingness, transformations, and summary statistics.
- Exact observation model p(x | theta), including noise and dependencies.
- Whether the NPE simulator reproduces this observation model rather than returning only deterministic ODE means.
- Whether MCMC and NPE condition on the same information. A lossy summary and the full dataset do not generally define the same posterior.
- Penalties, filtering, truncation, equilibrium rejection, or parameter-dependent simulator failures that may change the target.
- Actual installed SNPE implementation and proposal correction. Never enable an uncorrected first-round loss for posterior-proposal data without a separately justified algorithm.
- What existing posterior/checkpoint files actually contain and how they were generated.

Correct implementation discrepancies without silently changing the scientific model. If the two methods cannot be made to target the same posterior, report that limitation before drawing an algorithm-performance conclusion. If the exact likelihood cannot be evaluated, do not present an arbitrary synthetic-likelihood MCMC as an exact same-target baseline.

## 4. Fixed NPE/SNPE configuration

### 4.1 Primary neural method: two-round SNPE-C / APT

| Parameter | Fixed value |
| --- | --- |
| Algorithm | SNPE-C / APT; explicit class for the installed sbi version |
| Density estimator | Masked autoregressive flow (MAF) |
| Hidden features | 50 |
| Flow transforms | 5 |
| Round 1 prior simulations | 20,000 |
| Round 2 posterior-proposal simulations | 6,000 |
| Rounds | 2 |
| Proposal for round 2 | Round-1 posterior conditioned on x_obs |
| Reuse round-1 training simulations | Yes, using the implementation's supported corrected loss |
| Retrain network from scratch in round 2 | No; continue network weights using the supported round-training API |
| Training batch size | 256 |
| Optimizer | Adam |
| Learning rate | 0.0005 |
| Validation fraction | 0.10 |
| Early-stopping patience | 20 epochs without validation improvement |
| Maximum epochs | 250 per round |
| Gradient clipping norm | 5.0 |
| Atomic candidate count | 10, when atomic loss is used and configurable |
| Input standardization | Per-feature standardization for theta and x, using the installed estimator's supported training-data mechanism |
| Posterior draws | 20,000 |
| Sampling | Direct sampling from the trained posterior with supported prior-support rejection |

Map these conceptual settings to the installed API; do not assume keyword names. Freeze unspecified architecture options to the installed version's defaults and serialize the resolved settings. Respect constant features and structured input requirements. Do not manually re-standardize round-2 inputs under incompatible statistics while reusing a round-1 network; use the package-supported mechanism and record its behavior.

Keep proposal state frozen as required by the training API. Use its supported multi-round validation and data-reuse behavior, and document actual splits and loss variant. Do not compare raw validation losses across different objectives as if they measured the same quantity.

These are benchmark settings, not a claim that the current project already uses this architecture. Record differences from existing project defaults; preserve existing production code. An existing-configuration result may be reported separately if it materially differs, without silently substituting it for the standardized benchmark.

Count all actual ODE work, including failed simulations and proposal-generation work. The 26,000 figure assumes one forward trajectory per simulation. Do not silently replace failed simulations until 26,000 successes are reached: stop and audit the failure mechanism, then apply the same scientifically justified rule to all methods. Track failures and retries explicitly.

### 4.2 Neural control: single-round NPE

Use the same MAF architecture, batch size, learning rate, validation fraction, stopping rule, and 20,000 exported draws, but:

- Round count: 1.
- Prior simulations: 26,000.
- Standard prior NLL loss.
- Maximum training epochs: 250 total.

This isolates the practical benefit of spending 6,000 simulations adaptively rather than from the prior. It matches simulation budget, not neural optimization steps: SNPE can train up to 250 epochs twice. Include training cost in both results. Within each seed, reuse the first 20,000 prior simulations where valid, but charge the full attributed simulation cost to each cold-start method. Report actual shared-cache savings separately.

Do not add SNPE-B as a mandatory fourth method; the user's main question is NPE versus MCMC. Add it only to investigate a specific SNPE-C failure or a separately requested B-versus-C comparison.

## 5. Fixed MCMC configuration and sampler-selection rule

### 5.1 Concrete forward-budget baseline: adaptive random-walk Metropolis

Use a reliable existing implementation where possible. This baseline needs no ODE gradients and makes the forward budget transparent; it is not assumed to be the strongest MCMC method for every posterior.

| Parameter | Fixed value |
| --- | --- |
| Independent chains | 4 |
| Initial target-density evaluation slots | 1 per chain |
| Warmup/adaptation proposals | 1,499 per chain |
| Retained sampling proposals | 5,000 per chain |
| Maximum target-density evaluations | 6,500 per chain; 26,000 total |
| Saved post-warmup states | 5,000 per chain; 20,000 total |
| Proposal | Symmetric Gaussian random walk in suitable unconstrained coordinates |
| Initial covariance | 0.1^2 times diagonal transformed-prior variances |
| Begin covariance adaptation | After 200 warmup proposals |
| Update covariance | Every 50 warmup proposals |
| Adaptive covariance component | (2.38^2 / d) times regularized empirical covariance |
| Regularization | Add 1e-6 times the diagonal transformed-prior covariance |
| Fixed safety component | With probability 0.05, use the initial proposal covariance |
| Remaining proposal probability | 0.95, adaptive covariance component |
| Freeze adaptation | Before the first retained proposal |
| Thinning | None |

Estimate transformed-prior variances from 2,048 prior draws requiring no ODE solves; record this initialization overhead. For priors without stable finite transformed variances, use a documented robust IQR-based scale and record the change. Use log/logit transformations only when valid; add the correct Jacobian to the transformed target. Preserve multivariate prior dependencies in density evaluation.

Initialize chains independently from valid prior draws, not from an NPE-trained posterior. Do not select only favorable initializations. Any additional ODE evaluations needed to find a valid initial state come out of the chain's 6,500-slot budget, reducing retained length accordingly. Record this explicitly. A chain with inadequate retained output or failed diagnostics is a failed budgeted result, not grounds for silent extra computation.

Cache the current log density; each random-walk proposal should require at most one new target evaluation. A rejection retains the current state and does not require re-evaluating it. Prior-support rejections may avoid an ODE solve, so 26,000 target evaluations is an upper bound on ODE calls rather than an exact equality. Record target evaluations, forward trajectories, solver retries, and CPU time separately. Multiple ODE trajectories per likelihood require separate accounting.

The default scale 2.38^2/d is a starting configuration, not an optimality guarantee. Inspect acceptance and mixing; do not claim that SNPE beats MCMC generally if this particular random walk mixes poorly.

### 5.2 Strong existing MCMC comparator

If the repository already has a reliable, reasonably tuned direct-likelihood HMC/NUTS, slice, or other sampler, also benchmark it instead of relying exclusively on the random-walk baseline. Preserve its supported parameterization and settings, documenting them. Do not build a new differentiable ODE implementation solely for this benchmark.

For available NUTS, use these starting settings:

- Independent chains: 4.
- Target acceptance: 0.85.
- Maximum tree depth: 10.
- Mass matrix: diagonal initially; existing justified dense-metric configuration may be retained and explicitly recorded.
- Initial warmup request: 1,000 transitions per chain.
- Initial retained request: 1,000 transitions per chain.

These transition counts are not a 26,000-ODE budget. Run this comparator primarily under the matched-time protocol below. Instrument all value/gradient evaluations and forward/sensitivity solves. Do not equate one automatic-differentiation sensitivity solve to one plain trajectory in a compute claim. If an implementation cannot safely stop at a completed transition, use measured transition costs to choose a bounded run, report actual runtime, and call it approximately time-matched.

For a secondary forward-call-capped HMC/NUTS run, cap at 6,500 actual forward solves per chain where meaningful and supported; complete the in-flight transition and report any overrun. If the cap is exhausted during warmup, record that outcome; do not present nonstationary warmup draws as posterior samples. Different solver sensitivity costs make this track descriptive rather than definitive.

## 6. Primary fairness track: matched end-to-end wall-clock time

For replicate r, let T_SNPE,r be the measured cold-start runtime of the full two-round SNPE pipeline, including simulation, preprocessing, training, posterior construction, and 20,000-draw sampling. Exclude only common data loading already shared by all methods; list compilation and startup separately and include method-specific costs in totals.

Use T_SNPE,r as the MCMC deadline on the same four-core CPU allocation. Reuse the same MCMC trajectory as the forward-budget run when technically possible: checkpoint at the forward-budget endpoint and at the time deadline, and attribute elapsed times correctly. Do not restart warmup merely to obtain the second endpoint, and do not provide the time-matched run with free warmup.

For random-walk Metropolis:

1. Run the same 1,499 warmup proposals per chain.
2. Freeze adaptation, then continue sampling until the shared wall-clock deadline.
3. Preserve all completed states. Record draws and actual forward calls instead of forcing exactly 20,000 samples.
4. If the deadline arrives before usable post-warmup output, report failure at this budget.

For existing NUTS or another competent sampler, include its warmup and initialization in the same deadline. Warmup completion is mandatory; never shorten it after inspecting results merely to create a favorable comparison. Record incomplete or overrun runs transparently.

A budget deviation up to 10% may be labeled approximately matched if unavoidable due to transition granularity. Larger deviations require reporting results at the actual time rather than asserting equal compute. Do not pad faster methods with idle time. Use the full MCMC output for diagnostics/statistics; drawing-count differences do not invalidate a time-matched comparison.

Report paired per-seed results and median/range across the three repetitions. Run each method without concurrent unrelated benchmark tasks. Do not generalize a CPU-only result to GPU-equipped environments.

## 7. Reference posterior: validation infrastructure, not a budgeted competitor

Build a separate, more thoroughly checked numerical reference using the competent direct-likelihood sampler. It is not analytic ground truth.

Default reference for the random-walk baseline:

- 4 independent chains with different seeds from all benchmark runs.
- 1 initialization evaluation + 4,999 adaptive warmup proposals per chain.
- 45,000 retained proposals per chain, with adaptation frozen.
- At most 200,000 target evaluations total, including warmup.
- No automatic extension beyond this cap.

Prefer an already available, well-diagnosed strong sampler or archived same-model reference when provenance is verified. Reference generation may stop at a saved checkpoint once key statistics are demonstrably precise and stable. If the 200,000-evaluation cap is inadequate or random-walk mixing fails, report that accuracy conclusions are unresolved rather than silently escalating computation. A stronger existing sampler can be used with its actual costs documented.

Reference checks:

- Rank-normalized split R-hat < 1.01 as a screening criterion.
- Bulk/tail ESS and MCSE for means and relevant quantiles.
- Reference mean MCSE <= 0.02 posterior SD for key parameters, so it is below the default 0.1-SD comparison tolerance.
- Trace plots, dispersed-start agreement, and investigation of possible missed modes.
- For HMC/NUTS: divergences, tree-depth saturation, and energy diagnostics.

R-hat alone does not establish global mode exploration. Ensemble walkers are not independent chains; use independent ensembles if relevant. Reference construction and validation costs must appear separately in the report, not be charged only to the MCMC competitor or hidden entirely.

## 8. Accuracy and uncertainty metrics

### Parameter marginals

For every parameter and replicate, on original units, report mean, median, SD, and 50%, 90%, and 95% equal-tailed intervals. Compute:

- Mean difference divided by reference SD.
- 95% interval endpoint differences divided by reference SD.
- SD ratio and interval-width ratio relative to reference.
- Reference MCSE, finite-draw variability within a fixed NPE network, and variation across independently trained networks, separately.

Use a documented alternative scale when reference SD is effectively zero. Do not claim NPE is accurate because independent draws have high ESS; independence does not measure neural approximation bias.

### Joint structure

Compare marginal plots, representative pair plots, correlation-matrix differences, boundary behavior, and mode masses. If these leave a material joint-distribution question unresolved, use one additional metric such as C2ST or MMD with fixed scaling and comparable sample sizes. Account for MCMC dependence in train/test splits and reference-versus-reference controls. A nonsignificant test is not proof of equality.

### Predictive and policy outcomes

Use 1,000 states per method and replicate for initial evaluation, sampled across all chains/runs without treating correlated states as independent evidence. Apply identical scenarios and horizons, and distinguish latent ODE outputs from noisy posterior predictive observations. Report Monte Carlo uncertainty; increase this subset only if it prevents a decision.

Select a few existing scientifically important outcomes before inspecting method differences: for example cumulative infections, deaths, costs, or elimination-target probabilities. Do not invent policy definitions. Good predictive fit alone cannot establish a correct parameter posterior.

## 9. Predeclared decision criteria

These are provisional engineering tolerances, not universal statistical standards:

- Absolute standardized mean difference <= 0.1 for each key parameter.
- Absolute standardized 95% interval endpoint difference <= 0.2.
- 95% interval-width ratio between 0.9 and 1.1.
- No demonstrated omission of a materially weighted mode.
- Policy differences below the project's existing decision tolerance. If no tolerance exists, show sensitivity at 1%, 5%, and 10% relative error for nonzero continuous outcomes, with absolute errors for near-zero values; use 1 and 5 percentage points for probabilities. Do not retrospectively choose whichever threshold favors one method.

Show worst-case and per-parameter failures, not only averages. If MCMC sampling uncertainty is comparable to a tolerance, mark that comparison inconclusive. Freeze the key-parameter list and tolerances in the manifest before formal comparisons.

Conclude one of: MCMC preferred; NPE/SNPE faster at acceptable error; comparable performance; or insufficient evidence. State which MCMC implementation was actually tested. Neither a budget-exhausted weak MCMC run nor a very fast biased neural posterior establishes general superiority.

## 10. Optional extensions and stopping rules

Do not automatically run a large hyperparameter search, SNPE-B experiment, or full SBC campaign.

If fixed-observation comparisons remain unresolved, add only the experiment that addresses the identified risk. For research-level coverage claims, use complete synthetic-data calibration experiments, generating theta from the prior and x from the full observation model. Re-run sequential inference for each synthetic observation; a single network specialized to the real x_obs is not an appropriate replacement. Determine the number of datasets from desired binomial precision, and report uncertainty. A 5-10-dataset smoke test does not establish calibrated coverage.

Keep tuning budgets explicit and balanced; never report only the best seed. Stop optional testing once the practical method-selection question is answered.

## 11. Required deliverables

Use an isolated benchmark/npe_vs_mcmc/ directory or the repository's equivalent:

- model_audit.md and an English report.md.
- A machine-readable resolved configuration containing all values above and any justified deviations.
- Reproducible execution and analysis scripts without machine-specific hard-coded paths.
- run_manifest.json with seeds, source/data hashes, versions, resources, budgets, actual counts, runtimes, and status.
- timings.csv, parameter_comparison.csv, policy_comparison.csv, and diagnostics.csv.
- Raw chains, posterior samples, trained-network checkpoints, and logs.
- Marginal/joint posterior plots and a table comparing accuracy at matched forward budget and matched time.

Report timing columns separately: initialization/compilation, ODE simulation, preprocessing, neural training, MCMC warmup, posterior sampling, and validation. Record allocated resources and peak memory. Count failed and incomplete runs. Clearly label everything not executed.

Execution order: audit -> timing/implementation smoke test -> primary SNPE and fixed-budget MCMC with time checkpoints -> single-round NPE control -> independent repetitions -> checked reference and accuracy analysis -> only necessary follow-up -> final report.

## 12. Sources

- APT / SNPE-C: https://proceedings.mlr.press/v97/greenberg19a.html
- SNPE-B original method: https://arxiv.org/abs/1711.01861
- sbi multi-round tutorial (verify compatibility with installed version): https://sbi-dev.github.io/sbi/v0.24.0/tutorials/02_multiround_inference/
- SBI benchmarking: https://proceedings.mlr.press/v130/lueckmann21a.html
- Stan posterior diagnostics: https://mc-stan.org/docs/reference-manual/analysis.html

These sources support methodological choices; they do not demonstrate superiority for this project.
