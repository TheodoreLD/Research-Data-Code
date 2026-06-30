# Final Model Details

This note gives a compact technical description of the three final 2023-2024
wolf relative encounter-frequency models.

## Common Modelling Target

The response is an independent wolf event count. The exposure is camera effort
in camera-days. INLA receives effort through `E`, so the linear predictor
describes the expected daily encounter rate and map outputs are converted to
expected events per 100 camera-days.

The models estimate relative encounter frequency, not abundance, density,
occupancy, or population size.

For camera-month row `i`, the common model structure is:

```text
log(mu_i) = log(E_i) + beta_0 + gamma[m_i] + u(s_i)
```

where `mu_i` is the expected event count, `E_i` is active camera-days,
`beta_0` is the baseline log encounter rate, `gamma[m_i]` is the fixed effect
for calendar month, and `u(s_i)` is the INLA-SPDE spatial random field at the
camera location.

## Data Units

All final models use camera-month rows:

1. A camera deployment is split when it crosses a calendar-month boundary.
2. The effort assigned to a row is the number of active camera-days inside that
   camera-month interval.
3. Wolf events are assigned by `eventStart` month.
4. The model includes calendar-month fixed effects as a temporal control and a
   spatial SPDE field.

This avoids mixing September and October effort/events in deployments that
cross month boundaries.

Month enters the model as a fixed effect, but the final mapped quantity is not
a single-month prediction. The maps report an effort-weighted annualized
survey-year surface:

```text
lambda_year(s) = sum_m w_m * 100 * exp(beta_0 + gamma[m] + u(s))
```

where `w_m` is the proportion of sampled camera-days in month `m`. This keeps
month in the model as a temporal control while reporting the spatial pattern
for the sampled survey-year period as a whole. The reference month remains only
the coding baseline for month-rate ratios.

Month is treated as a fixed effect. INLA can fit a random month effect, for
example with `f(month, model = "iid")`, but that was not used because each
dataset contains only a few sampled months and the aim was direct temporal
control rather than estimating a month-level variance.

## Road-Camera 2023 Model

Final script:

```text
scripts/wolf_2023_nb_month_split_workflow.R
```

Final output folder:

```text
results/road_2023/
```

Model:

```text
y_i ~ NegativeBinomial(mu_i, size)
log(mu_i) = log(effort_i) + intercept + month_i + spatial(s_i)
```

Key settings:

- 490 camera-month rows;
- 60 cameras;
- 586 independent wolf events;
- 5222.2 camera-days;
- reference month for coefficients: 2023-08;
- map target: effort-weighted annualized 2023 surface;
- annualization factor relative to the baseline month: 1.117;
- INLA-SPDE spatial random field;
- negative-binomial likelihood.

Priors:

- intercept: Gaussian centered on crude observed daily rate, SD 2.5 on log
  scale;
- month log-rate ratios: Gaussian(0, SD 1);
- negative-binomial log(size): Gaussian(log(2), SD 2);
- spatial range: PC prior, `P(range < 5000 m) = 0.5`;
- spatial marginal SD: PC prior, `P(SD > 2.0) = 0.05`.

Model comparison:

| Model | WAIC | Delta WAIC |
| --- | ---: | ---: |
| ZINB spatial-month | 1162.72 | 0.00 |
| NB spatial-month | 1162.97 | 0.25 |
| Poisson spatial-month | 1303.16 | 140.44 |

ZINB is only marginally lower by WAIC and has low estimated zero inflation
(`p = 0.032`), so the negative-binomial spatial-month model is retained for
parsimony.

Main diagnostics:

- posterior predictive camera total events: pass;
- posterior predictive camera zero fraction: pass;
- posterior predictive camera maximum count: pass;
- row Pearson dispersion: 0.665;
- camera Pearson dispersion: 0.290;
- residual Moran's I: -0.008, p = 0.638;
- row PIT KS p-value: 0.268;
- camera PIT KS p-value: 0.000520;
- deployment-order lag-1 residual correlation: r = 0.018, p = 0.711;
- formal equal-time 7-day lag-1 check: r = -0.072, p = 0.573;
- formal equal-time 14-day lag-1 check: r = -0.046, p = 0.412;
- required diagnostics pass: TRUE;
- spatial block CV row mean log predictive density: -1.407;
- spatial block CV row 90 percent coverage: 0.933;
- spatial block CV camera 90 percent coverage: 0.900;
- prior sensitivity: retained variants pass required diagnostics;
- mesh sensitivity: final, finer, and coarser mesh variants pass required diagnostics.

## Forest-Camera 2024 Model

Final script:

```text
scripts/wolf_forest_month_refit.R
```

Final output folder:

```text
results/forest_2024/
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
- reference month for coefficients: 2024-08;
- map target: effort-weighted annualized 2024 surface;
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
- spatial block CV mean log predictive density: -0.412;
- spatial block CV 90 percent coverage: 0.978;
- month-level residual lag-1 ACF: 0.050;
- prior sensitivity: all retained variants pass required diagnostics;
- spatial prior/range sensitivity: fixed and estimated spatial-range variants
  are reported in `results/forest_2024/wolf_forest_month_prior_sensitivity.csv`.

Main limitation:

The forest-camera dataset contains only 46 independent wolf events. The final
model is valid for relative encounter-frequency mapping, but month and spatial
effects should be interpreted with wide uncertainty.

## Road-Camera 2024 Model

Final script:

```text
scripts/wolf_2024_zinb_month_split_workflow.R
```

Final output folder:

```text
results/road_2024/
```

Model:

```text
y_i ~ ZeroInflatedNegativeBinomial1(mu_i, size, pi)
log(mu_i) = log(effort_i) + intercept + month_i + spatial(s_i)
```

`ZeroInflatedNegativeBinomial1` is INLA's type-1 zero-inflated negative-binomial
parameterization. The parameter `pi` represents additional structural-zero
probability, while the negative-binomial component models overdispersed event
counts through `mu_i` and `size`.

Key settings:

- 344 camera-month rows;
- 60 cameras;
- 479 independent wolf events;
- 3574.0 camera-days;
- reference month for coefficients: 2024-08;
- map target: effort-weighted annualized 2024 surface;
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
| NB spatial-month | 937.32 | 3.65 |
| Poisson spatial-month | 997.19 | 63.52 |

Main diagnostics:

- posterior predictive camera total events: pass;
- posterior predictive camera zero fraction: pass;
- posterior predictive camera maximum count: pass;
- row Pearson dispersion: 0.583;
- camera Pearson dispersion: 0.236;
- residual Moran's I: -0.033, p = 0.370;
- row PIT KS p-value: 0.118;
- camera PIT KS p-value: 0.000492;
- required diagnostics pass: TRUE;
- spatial block CV row mean log predictive density: -1.524;
- spatial block CV row 90 percent coverage: 0.962;
- spatial block CV camera 90 percent coverage: 0.933;
- prior sensitivity: retained variants are stable and pass required diagnostics;
- mesh sensitivity: final, finer, and coarser mesh variants pass required diagnostics.

Temporal checks:

- deployment-order lag-1 residual correlation remains detectable:
  r = -0.178, p = 0.00267;
- formal equal-time 7-day lag-1 check: r = -0.153, p = 0.118;
- formal equal-time 14-day lag-1 check: r = -0.136, p = 0.094.

The deployment-order check is retained as a warning-style diagnostic because
the gaps are not fixed. The equal-time 7-day and 14-day diagnostics are the
formal checks, and they are not significant at lag 1 after the camera-month
split.

## Final Interpretation

All three models are final for relative encounter-frequency mapping. The
road-camera 2023 model passes diagnostics after the camera-month temporal
correction and is retained as a parsimonious NB model. The forest-camera 2024
model passes diagnostics and prior sensitivity checks. The road-camera 2024
model passes diagnostics after the camera-month temporal correction, prior
sensitivity, mesh sensitivity, spatial block cross-validation, and formal
equal-time temporal diagnostics.
