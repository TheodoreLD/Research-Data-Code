# Final Model Details

This note gives a compact technical description of the two final 2024 wolf
relative encounter-frequency models.

## Common Modelling Target

The response is an independent wolf event count. The exposure is camera effort
in camera-days. INLA receives effort through `E`, so the linear predictor
describes the expected daily encounter rate and map outputs are converted to
expected events per 100 camera-days.

The models estimate relative encounter frequency, not abundance, density,
occupancy, or population size.

## Data Units

Both final models use camera-month rows:

1. A camera deployment is split when it crosses a calendar-month boundary.
2. The effort assigned to a row is the number of active camera-days inside that
   camera-month interval.
3. Wolf events are assigned by `eventStart` month.
4. The model includes month fixed effects and a spatial SPDE field.

This avoids mixing September and October effort/events in deployments that
cross month boundaries.

## Small/Local 2024 Model

Final script:

```text
scripts/wolf_small_2024_month_refit.R
```

Final output folder:

```text
results/small_2024/
```

Model:

```text
y_i ~ NegativeBinomial(mu_i, size)
log(mu_i) = log(effort_i) + intercept + month_i + spatial(s_i)
```

Key settings:

- 356 camera-month rows;
- 53 cameras;
- 46 independent wolf events;
- 4423.0 camera-days;
- reference month: 2024-06;
- prediction month: 2024-06;
- INLA-SPDE spatial random field;
- negative-binomial likelihood.

Priors:

- intercept: Gaussian centered on crude observed daily rate, SD 2.5 on log
  scale;
- month log-rate ratios: Gaussian(0, SD 1);
- negative-binomial log(size): Gaussian(log(2), SD 2);
- spatial range: PC prior, `P(range < 1000 m) = 0.5`;
- spatial marginal SD: PC prior, `P(SD > 1.5) = 0.05`.

Main diagnostics:

- posterior predictive total events: pass;
- posterior predictive zero fraction: pass;
- posterior predictive maximum camera count: pass;
- row Pearson dispersion: 0.550;
- camera Pearson dispersion: 0.567;
- residual Moran's I: -0.041, p = 0.656;
- PIT KS p-value: 0.952;
- required diagnostics pass: TRUE;
- spatial block CV mean log predictive density: -0.415;
- spatial block CV 90 percent coverage: 0.978;
- month-level residual lag-1 ACF: 0.050.

Main limitation:

The small/local dataset contains only 46 independent wolf events. The final
model is valid for relative encounter-frequency mapping, but month and spatial
effects should be interpreted with wide uncertainty.

## Large 2024 Model

Final script:

```text
scripts/wolf_2024_zinb_month_split_workflow.R
```

Final output folder:

```text
results/large_2024/
```

Model:

```text
y_i ~ ZeroInflatedNegativeBinomial1(mu_i, size, pi)
log(mu_i) = log(effort_i) + intercept + month_i + spatial(s_i)
```

Key settings:

- 344 camera-month rows;
- 60 cameras;
- 479 independent wolf events;
- 3574.0 camera-days;
- reference month: 2024-09;
- prediction month: 2024-09;
- INLA-SPDE spatial random field;
- zero-inflated negative-binomial type 1 likelihood.

Priors:

- intercept: Gaussian centered on crude observed daily rate, SD 2.5 on log
  scale;
- month log-rate ratios: Gaussian(0, SD 1);
- zero-inflation logit probability: Gaussian centered at 5 percent structural
  zero probability, SD 1.5 on logit scale;
- negative-binomial log(size): Gaussian(log(2), SD 2);
- spatial range: PC prior, `P(range < 5000 m) = 0.5`;
- spatial marginal SD: PC prior, `P(SD > 2.5) = 0.05`.

Model comparison:

| Model | WAIC | Delta WAIC |
| --- | ---: | ---: |
| ZINB spatial-month | 933.67 | 0.00 |
| NB spatial-month | 937.27 | 3.60 |
| Poisson spatial-month | 997.23 | 63.56 |
| NB non-spatial month | 1042.52 | 108.85 |
| Poisson non-spatial month | 1327.08 | 393.41 |

Main diagnostics:

- posterior predictive camera total events: pass;
- posterior predictive camera zero fraction: pass;
- posterior predictive camera maximum count: pass;
- row Pearson dispersion: 0.583;
- camera Pearson dispersion: 0.236;
- residual Moran's I: -0.034, p = 0.336;
- row PIT KS p-value: 0.0905;
- camera PIT KS p-value: 0.000342;
- required diagnostics pass: TRUE;
- spatial block CV row mean log predictive density: -1.513;
- spatial block CV camera 90 percent coverage: 0.933.

Temporal checks:

- deployment-order lag-1 residual correlation remains detectable:
  r = -0.180, p = 0.00239, median gap 12.7 days;
- formal equal-time 7-day lag-1 check: r = -0.153, p = 0.118;
- formal equal-time 14-day lag-1 check: r = -0.136, p = 0.094.

The deployment-order check is retained as a warning-style diagnostic because
the gaps are not fixed. The equal-time weekly and biweekly diagnostics are the
formal checks, and they are not significant at lag 1 after the camera-month
split.

Temporal-bin model sensitivity:

| Sensitivity model | Delta WAIC vs month model | 14-day lag-1 p-value |
| --- | ---: | ---: |
| Month fixed effects | 0.00 | 0.094 |
| Biweekly time-bin fixed effects | 2.58 | 0.105 |
| Weekly time-bin fixed effects | 7.18 | 0.102 |

The month model is retained because it has the best WAIC and gives comparable
temporal residual behaviour with fewer fixed effects.

## Final Interpretation

Both models are final for 2024 relative encounter-frequency mapping. The small
model is data-limited but passes diagnostics and prior sensitivity checks. The
large model is better supported statistically and remains robust after the
camera-month temporal correction, prior sensitivity, mesh sensitivity, spatial
block cross-validation, and formal equal-time temporal diagnostics.

