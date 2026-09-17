# Fastest treatment designs by outcome stratum

The revised analysis evaluates all 12 combinations of community uptake
(`C = 1%, 5%, 10%, 15%/y`) and prison uptake (`J = 40%, 70%, 100%/y`). All
columns use the same red sequential palette and shared 12-entry legend.

The selected summary figure is
`figures/analysis_ver2/fig_ver2_best_treatment_by_stratum.png`. It has five
columns ordered D, J, F, X, and Total. Each column contains the status-quo
trajectory plus the single treatment scenario selected as fastest for that
column, shown across all five outcome rows.

## Selection rule

For each outcome column (D, J, F, X, and Total), the pipeline compares all
non-status-quo treatment designs except the `target_coverage` family. The
`target_coverage` scenarios remain available in their own figure but are not
used to define the best-treatment figure. Designs are ranked by the lowest
posterior median infection count at year 20, with the year-10 count as the
tie-breaker.
The selected scenario is then displayed across all five outcome rows in its
5-column figure. If the same scenario is selected for multiple columns, it
uses the same colour and is listed once in the right-hand legend.

## Selected designs

| Outcome column | Fastest design | Best level | Median infections at 10 years | Median infections at 20 years |
|---|---|---|---:|---:|
| D | J + D/F | C=15%/y; J=100%/y | 295.1 | 111.5 |
| J | J + D/F | C=15%/y; J=100%/y | 234.5 | 91.1 |
| F | J + D/F | C=15%/y; J=100%/y | 224.2 | 92.3 |
| X | J + D/F/X | C=15%/y; J=100%/y | 1,013.2 | 197.1 |
| Total | J + D/F/X | C=15%/y; J=100%/y | 1,774.4 | 492.7 |

After excluding `target_coverage`, the selected designs are unchanged: treating D/F together with J is sufficient to reduce
infection most rapidly in D, J, and F, whereas X and the total population
benefit from extending community treatment to X. The common winning uptake
level is the highest explored community coverage and the prison coverage cap;
it should therefore be interpreted as an elimination-oriented upper-bound
scenario rather than an immediately feasible policy.

The selected trajectories remain non-zero at year 20. They should be described
as substantial reduction or near-elimination, not literal elimination. Costs,
capacity, retention, reinfection, and operational feasibility are not included
in this ranking.
