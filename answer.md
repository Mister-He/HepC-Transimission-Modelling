# What this study can say about HCV elimination in Singapore

## Short answer

This study is designed to evaluate HCV elimination strategies among people who
inject drugs (PWID) and people moving through the prison--community pathway. It
can therefore provide a quantitative answer about treatment coverage, delivery
setting, age targeting, and the residual burden in the modelled PWID strata
`D`, `J`, `F`, and `X`.

It cannot, by itself, establish elimination across the whole Singapore
population. The model does not include people with no PWID history, and it does
not directly calibrate national HCV incidence, mortality, diagnosis, treatment
uptake, or linkage-to-care data.

## Why the PWID result is relevant nationally

The Singapore modelling literature describes HCV as concentrated among people
with a history of injecting drug use and calibrates the 2017 ever-PWID
population and seroprevalence. It also reports that the population moves between
the community, halfway houses, and correctional facilities. This supports using
the present model to study the high-burden transmission network and the
prison-community treatment cascade, but not to claim complete national
coverage. See [Chaillon et al., Singapore modelling study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8174139/).

The 47% estimate from Hsiang et al. is a halfway-house observation among
ever-PWID, not a national prevalence target or the prevalence of every annual
release event. It is therefore not used in the current calibration. The model
uses only the age-specific prison seroprevalence and population targets in
`data/calibration_targets.json`.

The time-varying inflow is anchored to the official [SingStat/CNB annual
age-specific new-drug-abuser dataset](https://data.gov.sg/datasets/d_f174679381e9085e32869c32bdb01b87/view).
This improves representation of the observed drug-abuser inflow, but those
counts remain an enforcement-based proxy rather than a complete census of all
PWID.

## What can be claimed from the new scenarios

The scenario analysis now compares two WHO-style impact indicators among the
modelled PWID network: an 80% reduction in HCV incidence and a 65% reduction
in HCV-related mortality by 2030 versus 2015. In the model, incidence is the
new-infection flow into acute infection, and HCV-related mortality is the
disease-specific DC/HCC mortality component. The former 1% viremic-prevalence
criterion is retained only as a secondary operational diagnostic.

These indicators are still modelled-PWID proxies, not national validation
measures. The model does not include non-PWID populations and does not have
national incidence or mortality surveillance calibration.

The most defensible interpretation is therefore:

> A prison-centred strategy can estimate the reduction achievable in the
> high-burden PWID network; a mixed prison-plus-community strategy is needed to
> test whether residual infection in formerly incarcerated community strata is
> also controlled.

## Evidence needed for a national-elimination conclusion

To extend the conclusion from “among modelled PWID” to “Singapore overall”, the
study would need at least:

1. national or population-representative anti-HCV and HCV-RNA surveillance;
2. incidence and mortality estimates for current and former PWID;
3. treatment initiation, completion, SVR, reinfection, and linkage-to-care data;
4. estimates for non-PWID groups and other recognised high-risk groups; and
5. a model extension that links those groups to the PWID transmission network.

Until those data are added, the study can cover the PWID-driven component of
Singapore's epidemic and compare intervention strategies within that component,
but it should not be presented as a national elimination verification study.
