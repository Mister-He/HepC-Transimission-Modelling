Please conduct a comprehensive scientific and implementation audit of this HCV modelling project.

The primary files are:

- `Hepc_modelling_technical_details.pdf`
- `sim.cpp`
- `setup.R`
- `preprocess.R`

The objective is to determine whether the scientific specification is defensible and whether the preprocessing, configuration, and simulation implementation are mutually consistent.

## 1. Working principles

- Read the repository’s `AGENTS.md` and relevant project documentation before starting.
- Locate the actual files and inspect their dependencies, input data, and generated outputs where necessary.
- Record the repository commit and relevant file versions or hashes.
- Treat this as an audit. Do not modify the scientific document, model code, configuration, or original data. You may create separate audit reports and focused diagnostic scripts.
- Do not assume that the PDF is scientifically correct merely because it is the specification. Assess scientific validity separately from implementation consistency.
- Do not assume that code is incorrect merely because it differs from the PDF. Identify whether the discrepancy reflects an implementation error, outdated documentation, or an unresolved modelling decision.
- Distinguish confirmed errors, justified approximations, ambiguous descriptions, and unverifiable claims.
- Complete all accessible checks. If evidence or data are unavailable, explain precisely what remains unresolved and why.

## 2. Audit the scientific specification

Review `Hepc_modelling_technical_details.pdf` in full, including equations, diagrams, parameter tables, footnotes, and references. Inspect rendered pages where text extraction could misread mathematical notation or table structure.

### 2.1 Assumptions

For every substantive assumption:

- Explain its role in the model.
- Assess its epidemiological, clinical, demographic, and mathematical plausibility.
- Determine whether it is supported by evidence, a standard modelling convention, or an explicit pragmatic approximation.
- Assess applicability to the modelled population, setting, and time period.
- Identify contradictions between sections, equations, diagrams, and tables.
- Explain the likely consequences of questionable assumptions and whether sensitivity analysis is warranted.

Where applicable, examine population entry and exit, ageing, injecting behaviour, incarceration and release, transmission, disease progression, treatment, cure, reinfection, mortality, and equilibrium assumptions.

Do not label a simplifying assumption as an error solely because it omits real-world complexity.

### 2.2 Equations and mathematical definitions

Check whether:

- Every state variable, parameter, index, function, and initial condition is defined.
- Units and time scales are consistent.
- Rates, hazards, probabilities, proportions, counts, and durations are distinguished correctly.
- Conversions between these quantities are mathematically appropriate.
- Compartment flows have the correct signs, origins, and destinations.
- Population accounting is consistent with specified recruitment, ageing, migration, and mortality.
- Competing events and mortality are not inadvertently double-counted.
- Denominators and boundary cases are well defined, including empty populations.
- The formulation preserves the intended non-negativity and probability constraints.
- Initial conditions, equilibrium criteria, and observation equations are sufficiently specified to reproduce the model.

Distinguish properties of the mathematical model from properties of its numerical implementation. Identify under-specified equations rather than silently supplying missing assumptions.

### 2.3 Parameter provenance and use

Audit **every parameter-table entry**, including fixed, derived, calibrated, and assumed parameters.

For each entry:

1. Identify the stated value, uncertainty, units, population, and source.
2. Retrieve and inspect the original source where accessible.
3. Locate the exact supporting table, page, section, or supplementary material.
4. Verify that the source supports the claimed quantity and numerical value.
5. Check population, setting, period, endpoint, and subgroup applicability.
6. Reconstruct any conversion, averaging, extrapolation, or other derivation.
7. Assess whether the resulting quantity is appropriate for its actual role in the model.

Pay particular attention to:

- Annual versus daily quantities.
- Hazard-to-probability and duration-to-rate conversions.
- Relative versus absolute risks.
- Seroprevalence versus viraemic prevalence.
- Incident versus prevalent cases.
- Conditional versus unconditional probabilities.
- Age-specific versus aggregate estimates.
- Treatment eligibility, uptake, completion, and response.
- Background versus excess or disease-specific mortality.

If an entry is calibrated or assumed, verify that it is labelled accordingly; do not demand a literature source for a quantity explicitly estimated within the model.

### 2.4 References

Check:

- Whether each cited source exists.
- Bibliographic accuracy: authors, title, year, journal, DOI, and identifiers where provided.
- Agreement between in-text citations, parameter-table citations, and the reference list.
- Whether the cited source actually supports the associated claim.
- Misattributions, inappropriate secondary citations, duplicates, and missing references.

Use primary publications and authoritative official sources where possible. Search snippets alone are insufficient to verify a numerical claim. If the full source is inaccessible, mark the claim as unverified rather than incorrect.

## 3. Audit `sim.cpp` against the specification

Build an explicit mapping between the PDF and the implementation.

Check:

- State-vector ordering, compartment definitions, and indexing.
- Age-group definitions and boundary handling.
- Every implemented transition and its mathematical expression.
- Force of infection, mixing assumptions, and denominator populations.
- Recruitment, ageing, cessation, incarceration, release, and other applicable demographic or behavioural flows.
- Disease progression, treatment, cure, reinfection, and mortality.
- Parameter meanings, units, transformations, and time scales.
- Initial conditions and equilibrium-related logic.
- Outcome accumulation and aggregation.
- Numerical safeguards, clipping, and handling of invalid states.
- Omitted, duplicated, or additional flows.

Trace parameters through the actual call path, including any R-to-C++ interface, rather than comparing names alone.

For each discrepancy, show:

1. What the PDF specifies.
2. What the code implements.
3. The exact source locations.
4. Whether the difference changes model behaviour.
5. The likely scientific or numerical consequence.
6. A recommended resolution, without applying it.

Where useful, run focused diagnostics such as population-balance checks, zero-transmission checks, zero-treatment checks, and boundary-state checks. Test only properties implied by the specification; do not assume a closed population if recruitment or exits are intended.

## 4. Audit `setup.R` against the specification and simulation

Check whether `setup.R` correctly implements the PDF’s:

- Parameter values and uncertainty distributions.
- Units and conversion rules.
- Population and age-group structure.
- Initial conditions.
- Simulation horizon and time step.
- Fixed versus estimated parameter assignments.
- Relevant intervention settings and calibration targets.

Trace the **effective runtime configuration**, including imported files, defaults, overrides, and transformations before values reach `sim.cpp`.

Identify:

- Incorrect values or conversions.
- Swapped indices or mismatched age groups.
- Missing parameters.
- Unused configuration entries.
- Hard-coded values in `sim.cpp` that bypass `setup.R`.
- Duplicate definitions or overrides that make the effective value differ from the documented value.

## 5. Audit `preprocess.R` outputs against `setup.R`

Trace each relevant preprocessed quantity from raw input to its final use.

For every output consumed by `setup.R`:

- Identify the source data and preprocessing operations.
- Check filters, exclusions, missing-value handling, aggregation, weighting, denominators, and time periods.
- Check age bins, labels, factor levels, ordering, and joins.
- Reconstruct derived quantities and transformations.
- Verify output filenames, schemas, units, and precision.
- Confirm that `setup.R` reads or copies the intended output and applies any additional transformations correctly.
- Check for stale generated files, outdated hard-coded copies, and path-dependent behaviour.

When dependencies and data are available, rerun preprocessing in an isolated output directory and compare the regenerated results with the values actually used by `setup.R`.

Use exact equality for quantities that should match exactly. For floating-point differences, report the absolute and relative differences and justify the tolerance.

If raw data are unavailable, distinguish between:

- Code-level consistency that can be checked.
- Existing-output consistency that can be checked.
- Reproducibility that cannot be established.

Do not claim successful end-to-end reproduction without executing it.

## 6. Required audit outputs

Create a structured Markdown report, supported by tables where appropriate.

### A. Executive assessment

Summarize:

- Overall scientific defensibility.
- Consistency between the PDF and implementation.
- Consistency between preprocessing and runtime configuration.
- The most consequential findings.
- The limits of the audit.

Avoid an unconditional “validated” verdict when material checks remain unresolved.

### B. Findings register

For each finding, include:

| Field | Required content |
| --- | --- |
| ID | Unique identifier |
| Category | Scientific assumption, equation, parameter, citation, implementation, configuration, or preprocessing |
| Status | Confirmed error, potential issue, ambiguity, justified deviation, or unverified |
| Severity | Critical, major, minor, or informational |
| Location | PDF page/equation/table or source file and line numbers |
| Expected | What the specification or supporting evidence indicates |
| Observed | What the document, source, or implementation actually contains |
| Evidence | Reproducible calculation, code trace, diagnostic result, or source citation |
| Impact | Likely effect on model behaviour, estimates, or interpretation |
| Recommendation | Specific proposed correction or clarification |
| Confidence | High, moderate, or low, with a brief reason |

Prioritize severity by consequences, not by presentation or coding style.

### C. Traceability tables

Provide:

1. An assumption/equation-to-code mapping.
2. A complete parameter provenance table.
3. A PDF-to-`setup.R`-to-`sim.cpp` parameter mapping.
4. A preprocessing lineage table linking raw inputs, transformations, outputs, and runtime use.
5. A reference-verification table.

Include verified items as well as problematic ones so the coverage of the audit is visible.

### D. Reproducibility and unresolved items

List:

- Files and versions inspected.
- Commands and diagnostic scripts executed.
- Tests passed, failed, or not performed.
- Missing data, inaccessible sources, and environmental limitations.
- Questions requiring a scientific decision from the project owner.

### E. Prioritized action plan

Separate:

- Corrections needed before relying on model results.
- Documentation changes needed for reproducibility.
- Assumptions requiring justification or sensitivity analysis.
- Optional improvements.

## 7. Execution expectations

Start by inventorying the relevant files and tracing the data/configuration/simulation pipeline. Then audit the scientific specification and compare it systematically with the implementation.

Do not stop after discovering the first few discrepancies. Account for every assumption, equation, and parameter-table entry within scope, marking its verification status.

Keep conclusions proportional to the evidence. Do not fabricate citations, infer inaccessible source contents, silently repair discrepancies, or present proposed diagnostics as completed tests.
