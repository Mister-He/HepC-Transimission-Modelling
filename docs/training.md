# Training

## Deterministic calibration

```bash
Rscript scripts/run_calibration.R \
  --run-id r_YYYYMMDD_agepop --maxit 1000 --n-starts 5 --seed 20260929
```

Targets are read from `data/calibration_targets.json`. Every objective
evaluation runs model time 0--150, checks the final five-year equilibrium, and
uses the final state as the 2015 snapshot. The likelihood contains only six
binomial seroprevalence terms and six log-normal population terms.

The six fitted contact scales and six age-specific inflow scales multiply the
2015 SingStat/CNB baseline in `data/singstat_new_drug_abusers.csv`. Values from
2016 onward use the same fitted scales; the latest observed year is held
constant through 2030.

## Neural posterior estimation

```bash
Rscript scripts/run_npe.R --step all \
  --fit output/calibration/r_YYYYMMDD_agepop/fit.rds \
  --out-dir output/calibration/npe_agepop \
  --n-sims 600 --n-proposal 100 --n-draws 1000 \
  --n-cores 4 --n-rounds 1 --max-epochs 250 --seed 20260929 \
  --python /path/to/python
```

NPE uses binomial prevalence noise and log-normal population noise
(`sigma = 0.10`). SBC coverage and posterior spread are quality gates. No
Laplace approximation is run or stored. The current accepted run uses
`/tmp/bayes-venv/bin/python`; use an equivalent environment containing
PyTorch and `sbi`.

## 2026 treatment scenarios

```bash
Rscript scripts/scenarios_2026.R \
  --fit output/calibration/r_YYYYMMDD_agepop/fit.rds
```

The script projects 2015--2025 without treatment and starts each policy in
2026. Treatment is prohibited in D in every scenario. The run evaluates 75
combinations: J and F/X treatment initiation probabilities of 0%, 25%, 50%,
75%, and 100%, crossed with age settings 20--50, 30--60, and all age groups.
These map to the model's 20--49 and 30--59 five-year age bands.
The WHO impact screen requires 80% incidence reduction and 65% HCV-related
mortality reduction versus the 2015 baseline. The script converts initiation
probability to the ODE hazard using `tau = -log(1 - probability)`; 100% is
represented numerically by a 10/y hazard (99.995%). It writes age-specific
Nature-style figures, one combined heatmap, treatment-demand figures only for
age settings with at least one passing policy, one combined trade-off figure,
endpoint tables and an optional Excel workbook.

Treatment-demand counts are integrated treatment starts over 2026--2030. They
are treatment courses rather than unique people; because the model permits
reinfection, one person may contribute more than one course. The primary
trade-off is reported using treatment starts versus cumulative viremic
person-years averted relative to status quo, with a Pareto frontier and an
efficiency ratio. The exact individual-level number of people treated at least
once is not identifiable from the aggregate compartment model.
