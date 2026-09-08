# Interpretation of the survival curves

The figure shows posterior survival curves for two age-specific events over
the first 10 years after an index event:

- **Recidivism-event survival:** the probability that a person has not had a
  recidivism/re-arrest event by time `t`.
- **Release-event survival:** the probability that a person has not yet been
  released by time `t`.

All curves start at 1 at time 0 because no event has occurred at the beginning
of follow-up. A steeper decline indicates that the event is occurring more
quickly in that age group. Thus, at any given time, a lower curve indicates a
higher accumulated probability of the corresponding event.

## Model quantities

For age group `i`, the release-event survival curve is based on the fitted
release rate `lambda_release[i]`:

```text
S_release(t | i) = exp(-lambda_release[i] * t)
```

The rearrest model has an age-specific rearrest hazard `lambda_rearrest[i]`
and a common recidivism propensity `pi`. Its survival curve is:

```text
S_recidivism(t | i) = (1 - pi) + pi * exp(-lambda_rearrest[i] * t)
```

This mixture form means that the recidivism curves approach `1 - pi` rather
than zero: the model assigns probability `1 - pi` to people who do not enter
the recidivism process. By contrast, release survival approaches zero when
the release rate is positive.

The curves treat the displayed age group as fixed throughout the 10-year
interval. They are therefore a direct comparison of the fitted age-specific
rates, rather than a projection that moves people between age groups as they
age.

## Uncertainty and colours

The solid line is the posterior median. The translucent band is the 95%
credible interval, representing posterior uncertainty in the fitted rates and
recidivism propensity. The `Age Group` legend uses a sequential gradient from
younger groups (`<20`) to older groups (`60+`).

The recidivism panel should be read separately from the release panel: the
vertical scales describe survival probabilities, but the two panels use
different fitted hazards and represent different events.
