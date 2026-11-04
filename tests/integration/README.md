# Integration tests guide

Integration tests compile and exercise the simulator, including recovery,
reinfection, treatment/SVR flows, constant entry, output thinning, and likelihood wiring.
`test_independent_rhs.R` checks a full `sim.cpp` step against an independently derived R right-hand
side, the liver/excess split of HCV deaths and the whole-step count.
