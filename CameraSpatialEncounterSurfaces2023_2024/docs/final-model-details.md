# Final Model Details

This note is the technical reference for the three final 2023-2024 wolf
relative encounter-frequency models: full methodology, priors, and complete
diagnostic numbers. `README.md` gives the short summary and points here for
detail, so numbers are stated once, in this file. All numbers below are from
a `WOLF_RUN_PROFILE=final` rerun of the three scripts against the private
camera-trap data, and match the files committed under `results/`.

## Common Modelling Target

The response is an independent wolf event count per camera-month row. The
exposure is camera effort in camera-days, passed to INLA through the `E`
argument, so the linear predictor describes the expected daily encounter
rate and map outputs are converted to expected events per 100 camera-days.

The models estimate relative encounter frequency, not abundance, density,
occupancy, or population size. Camera-trap encounter rates are a standard
relative-abundance index in ecology when detection probability cannot be
estimated directly (Rowcliffe et al. 2008; O'Brien 2011).

For camera-month row `i`, the shared model structure is:

```text
log(mu_i) = log(E_i) + beta_0 + gamma[m_i] + u(s_i)
```

where `mu_i` is the expected event count, `E_i` is active camera-days,
`beta_0` is the model intercept on the log encounter-rate scale, `gamma[m_i]`
is the fixed effect for calendar month, and `u(s_i)` is a spatial random
field estimated with the INLA-SPDE method: the field is represented on a
triangulated mesh over the study area and given a Matérn covariance via a
stochastic partial differential equation, fitted by integrated nested
Laplace approximation rather than MCMC (Rue, Martino & Chopin 2009; Lindgren,
Rue & Lindström 2011). Two of the three models use a negative-binomial
likelihood for `y_i` to absorb overdispersion beyond what a Poisson count
model allows (Hilbe 2011); the road-camera 2024 model additionally uses a
zero-inflated negative-binomial likelihood (see that survey's section).

## Data Units: Camera-Month Rows

All final models are fit on camera-month rows, not raw per-day or
per-deployment records:

1. A camera deployment is split at every calendar-month boundary it crosses.
2. Each resulting row's exposure is the number of active camera-days inside
   that specific camera-month interval.
3. Wolf events are assigned to the row whose month contains the event's
   `eventStart` timestamp.

For example, a camera deployed 2023-07-20 to 2023-08-15 produces two rows:
one for July (12 active days, 2023-07-20 to 2023-07-31) and one for August
(15 active days, 2023-08-01 to 2023-08-15). Any wolf event recorded in that
window is assigned to whichever row's month contains it. This keeps effort
and events aligned and avoids mixing two different months' exposure into one
row.

The calendar month itself enters the model as a fixed effect (`gamma[m_i]`
above), estimated jointly with the spatial field, not fit separately.

## Annualized Map Surface

Month is a fixed effect in the model, but the published maps are not a
single-month prediction. Each survey's sampled months are combined into one
effort-weighted, annualized surface:

```text
lambda_year(s) = sum_m w_m * 100 * exp(beta_0 + gamma[m] + u(s))
```

`w_m` is the share of that survey's total camera-days that fell in month
`m`. Equivalently, the reference-month prediction is scaled by an
annualization factor: the effort-weighted average of `exp(gamma[m] -
gamma[m_ref])` over all sampled months. Concretely, for road-camera 2023
(reference month 2023-08):

| Month | Share of camera-days | Rate ratio vs. 2023-08 |
| --- | ---: | ---: |
| 2023-07 | 0.296 | 0.969 |
| 2023-08 (reference) | 0.341 | 1.000 |
| 2023-09 | 0.334 | 1.371 |
| 2023-10 | 0.030 | 1.078 |

The effort-weighted average of the rate-ratio column is 1.117 — the
annualization factor reported for that survey below. It means the
sampled-year average rate runs about 11.7% above the August-only rate,
mostly because September (33% of camera-days) has a higher fitted rate. Each
survey's full weight table is in
`results/<survey>/wolf_<survey>_annualization_weights.csv`.

The mapped central estimate is the posterior mean of this annualized
quantity; the accompanying posterior-SD map is its matching uncertainty
surface, on the same units.

## Diagnostic Gate

A model is called final only if it clears a specific gate: the camera-level
posterior predictive checks (total events, zero fraction, max count; Gelman,
Meng & Stern 1996) and the residual Moran's I spatial-autocorrelation test
(a two-sided permutation test, scaled by `WOLF_RUN_PROFILE`, in all three
scripts; Moran 1950, applied to model residuals following Dormann et al.
2007).

PIT (probability integral transform) KS p-values (Czado, Gneiting & Held
2009, for discrete/count outcomes) are also computed and reported for every
model as a supporting calibration diagnostic, but are not part of the gate.
Camera-level PIT is sensitive to small-count discreteness and camera-level
clustering in ways that do not track whether the mapped spatial surface
itself is distorted, so a low PIT KS p-value is reported but does not on its
own fail a model. For example, the road-camera 2023 model reports "required
diagnostics pass: TRUE" alongside a camera PIT KS p-value near 0.0015: the
low value is shown, not dropped, but does not gate the pass/fail call.

Model comparison across candidate likelihoods (Poisson / NB / ZINB) uses
WAIC (Watanabe 2010), cross-checked against DIC (Spiegelhalter et al. 2002).
Out-of-sample predictive performance is checked with spatial block
cross-validation: cameras are grouped into spatial folds by k-means
clustering, each fold is held out in turn, and the SPDE mesh for that fold
is rebuilt from the training cameras only, so no held-out location leaks
into the training mesh (following the general spatially-blocked
cross-validation approach of Roberts et al. 2017). Held-out counts are
simulated from full joint posterior draws of the fitted model.

**Sensitivity checks.** All three scripts refit the final model under
perturbed priors and perturbed SPDE mesh resolution (finer/coarser) and
report whether WAIC, DIC, and posterior hyperparameters stay stable. The
forest-camera script additionally recomputes the full PPC/Moran's I gate for
every sensitivity variant, independently re-verifying "passes required
diagnostics" at each one; the two road-camera scripts check WAIC/DIC/
hyperparameter stability only and do not recompute the gate per variant.
This is a difference in how much each script's sensitivity loop checks, not
a difference in the final fitted models themselves.

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
- observed mean encounter frequency: 11.221 events per 100 camera-days;
- map target: effort-weighted annualized 2023 surface (annualization factor
  1.117, see above).

Weakly informative priors:

- intercept: Gaussian(mean = -2.187, SD 2.5 on log scale), centered on the
  crude observed daily rate;
- month log-rate ratios: Gaussian(0, SD 1);
- negative-binomial log(size): Gaussian(log(2), SD 2);
- spatial range: PC prior, `P(range < 5000 m) = 0.5` (Fuglstad et al. 2019);
- spatial marginal SD: PC prior, `P(SD > 2.0) = 0.05` (Simpson et al. 2017).

Model comparison:

| Model | WAIC | Delta WAIC |
| --- | ---: | ---: |
| ZINB spatial-month | 1162.65 | 0.00 |
| NB spatial-month | 1162.90 | 0.25 |
| Poisson spatial-month | 1303.06 | 140.41 |

ZINB is only marginally lower by WAIC and has a low estimated zero-inflation
probability (mean 0.032), so the negative-binomial spatial-month model is
retained for parsimony.

Main diagnostics:

- posterior predictive camera total events / zero fraction / maximum count:
  all pass;
- row Pearson dispersion: 0.665; camera Pearson dispersion: 0.292;
- residual Moran's I: -0.004 (expected -0.017), two-sided p = 0.492;
- row PIT KS p-value: 0.1019; camera PIT KS p-value: 0.001507
  (supporting diagnostic, not part of the gate — see above);
- negative-binomial size posterior mean: 1.708;
- required diagnostics pass: TRUE;
- temporal residual autocorrelation: within-camera lag-1 r = 0.017,
  p = 0.7281 (n = 430 pairs); no evidence of residual temporal
  autocorrelation; date-ordered mean-residual lag-1 ACF: -0.176;
- spatial block cross-validation: row 90 percent coverage = 0.96, camera 90
  percent coverage = 0.95;
- prior sensitivity: WAIC, DIC, and posterior hyperparameters (NB size,
  spatial range, spatial SD) are stable across 6 prior variants (WAIC 1162.80
  to 1163.56; delta WAIC 0.00 to 0.76; stability checked, gate not
  recomputed per variant — see "Sensitivity checks" above);
- mesh sensitivity: WAIC and hyperparameters are stable across the final,
  finer, and coarser mesh variants (WAIC 1162.81 to 1163.05; delta WAIC 0.00
  to 0.24; same basis as above).

Exploratory check: a Spearman correlation of deployment start day-of-year
against UTM northing (`results/road_2023/wolf_2023_exploratory_timing_vs_northing.csv`)
gives rho = 0.042, p = 0.357 (n = 490) — deployment timing is not
meaningfully correlated with camera location in this survey.

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
- observed mean encounter frequency: 1.040 events per 100 camera-days;
- map target: effort-weighted annualized 2024 surface.

Weakly informative priors:

- intercept: Gaussian centered on crude observed daily rate, SD 2.5 on log
  scale;
- month log-rate ratios: Gaussian(0, SD 1);
- negative-binomial log(size): Gaussian(log(2), SD 2);
- spatial range: PC prior, `P(range < 1000 m) = 0.5`;
- spatial marginal SD: PC prior, `P(SD > 1.5) = 0.05`.

Fitted hyperparameters:

- negative-binomial size posterior mean: 1.836;
- spatial range posterior mean: 584.1 m;
- spatial SD posterior mean: 0.816.

Month effects (rate ratio vs. reference month 2024-08):

- 2024-03: 0.34 (95% CrI 0.08 to 1.47);
- 2024-04: 0.65 (95% CrI 0.22 to 1.92);
- 2024-05: 0.91 (95% CrI 0.36 to 2.31);
- 2024-06: 1.69 (95% CrI 0.74 to 3.79);
- 2024-07: 1.28 (95% CrI 0.53 to 3.07);
- 2024-09: 0.67 (95% CrI 0.23 to 1.99).

Main diagnostics:

- posterior predictive row and camera total events / zero fraction /
  maximum count: all pass;
- row Pearson dispersion: 0.584; camera Pearson dispersion: 0.614;
- residual Moran's I: -0.035, two-sided p = 0.691;
- row PIT KS p-value: 0.7938; camera PIT KS p-value: 0.5125;
- required diagnostics pass: TRUE;
- temporal residual autocorrelation: within-camera lag-1 r = -0.046,
  p = 0.4211 (n = 303 pairs); no evidence of residual autocorrelation. This
  survey additionally has enough sampled months (seven, 2024-03 to 2024-09)
  to compute a month-level lag-1 ACF as a second, low-power supporting
  check: 0.144. The road-camera surveys sample only 3-4 months, too few for
  that check to carry any power, so it is not reported for them;
- spatial block cross-validation: row 90 percent coverage = 0.978, camera 90
  percent coverage = 0.925. This survey's cross-validation builds each
  fold's SPDE mesh from the training-fold cameras only, and simulates
  held-out counts from full joint posterior samples;
- prior sensitivity: all 12 variants pass required diagnostics (WAIC 269.49
  to 276.75; best variant `month_sd_0_5`, WAIC 269.49; final-current variant
  WAIC 270.21). This survey's sensitivity loop recomputes the full PPC/
  Moran's I diagnostic gate independently for each variant (see "Sensitivity
  checks" above);
- mesh sensitivity: final, finer, and coarser mesh variants all pass
  required diagnostics; WAIC range = 269.49 to 270.42.

Main limitation: the forest-camera dataset contains only 46 independent wolf
events. Weakly informative priors do not add information the data lack, so
with this few events the posterior for month and spatial effects stays wide
— the wide credible intervals above reflect data sparsity, not the priors
being informative. The model remains valid for relative encounter-frequency
mapping, but month and spatial effects should be read with that wide
uncertainty in mind.

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

`ZeroInflatedNegativeBinomial1` is INLA's type-1 zero-inflated
negative-binomial parameterization (see Martin et al. 2005 for the ecological
rationale for modelling structural zeros separately from count-process
zeros). `pi` is the additional structural-zero probability, while `mu_i` and
`size` describe the negative-binomial count component.

Key settings:

- 344 camera-month rows;
- 60 cameras;
- 479 independent wolf events;
- 3574.0 camera-days;
- observed mean encounter frequency: 13.402 events per 100 camera-days;
- map target: effort-weighted annualized 2024 surface (annualization factor
  1.195).

Weakly informative priors:

- intercept: Gaussian(mean = -2.010, SD 2.5 on log scale), centered on the
  crude observed daily rate;
- month log-rate ratios: Gaussian(0, SD 1);
- zero-inflation logit probability: Gaussian(mean = -2.94, SD approximately
  1.5 on logit scale);
- negative-binomial log(size): Gaussian(log(2), SD 2);
- spatial range: PC prior, `P(range < 5000 m) = 0.5`;
- spatial marginal SD: PC prior, `P(SD > 2.5) = 0.05`.

Model comparison:

| Model | WAIC | Delta WAIC |
| --- | ---: | ---: |
| ZINB spatial-month | 933.64 | 0.00 |
| NB spatial-month | 937.32 | 3.68 |
| Poisson spatial-month | 997.30 | 63.66 |

Main diagnostics:

- posterior predictive camera total events / zero fraction / maximum count:
  all pass;
- row Pearson dispersion: 0.589; camera Pearson dispersion: 0.232;
- residual Moran's I: -0.034 (expected -0.017), two-sided p = 0.384;
- row PIT KS p-value: 0.08442; camera PIT KS p-value: 0.0008276;
- zero-inflation probability posterior mean: 0.075;
- negative-binomial size posterior mean: 3.622;
- required diagnostics pass: TRUE;
- temporal residual autocorrelation: within-camera lag-1 r = -0.179,
  p = 0.002523 (n = 284 pairs); date-ordered mean-residual lag-1 ACF: 0.252.
  Residual deployment-order temporal structure remains detectable here,
  unlike the other two surveys. The originally hypothesized mechanism
  (staggered deployment timing correlated with camera location) was tested
  directly: a Spearman correlation of deployment start day-of-year against
  UTM northing (`results/road_2024/wolf_2024_exploratory_timing_vs_northing.csv`)
  gives rho = 0.058, p = 0.280 (n = 344) — not a meaningful correlation, so
  it does not explain the residual autocorrelation, and no specific
  mechanism is established. This is treated as an open temporal caution
  rather than corrected, on the grounds that it does not appear to distort
  the mapped spatial surface: the spatial field is fit jointly with, and net
  of, the month effect, and both spatial block cross-validation coverage and
  mesh sensitivity remain stable (below);
- spatial block cross-validation: row 90 percent coverage = 0.96, camera 90
  percent coverage = 0.93;
- prior sensitivity: WAIC, DIC, and posterior hyperparameters are stable
  across the retained prior variants (WAIC 933.40 to 933.89; delta WAIC 0.00
  to 0.50; stability checked, gate not recomputed per variant — see
  "Sensitivity checks" above);
- mesh sensitivity: WAIC and hyperparameters are stable across the final,
  finer, and coarser mesh variants (WAIC 933.32 to 933.64; delta WAIC 0.00
  to 0.32; same basis as above).

## Final Interpretation

All three models are final for relative encounter-frequency mapping and pass
the required diagnostic gate (camera-level PPC plus residual Moran's I).

The road-camera 2023 model shows no evidence of residual temporal
autocorrelation and is retained as a parsimonious NB model over the
marginally-better-fitting ZINB alternative. The forest-camera 2024 model's
prior and mesh sensitivity variants independently re-verify the same
diagnostic gate at every variant, the most thorough check of the three; its
main limitation is the small number of independent events (46), which widens
posterior uncertainty on month and spatial effects without indicating a
model problem. The road-camera 2024 model passes the required
posterior-predictive and spatial diagnostics and is supported over NB/Poisson
by WAIC; its one open issue is a residual within-camera temporal correlation
whose mechanism is not established, but which spatial block cross-validation
and mesh sensitivity both indicate does not distort the mapped spatial
surface.

## References

- Czado, C., Gneiting, T. & Held, L. (2009). Predictive model assessment for
  count data. *Biometrics*, 65(4), 1254–1261.
- Dormann, C. F. et al. (2007). Methods to account for spatial
  autocorrelation in the analysis of species distributional data: a review.
  *Ecography*, 30(5), 609–628.
- Fuglstad, G.-A., Simpson, D., Lindgren, F. & Rue, H. (2019). Constructing
  priors that penalize the complexity of Gaussian random fields. *Journal of
  the American Statistical Association*, 114(525), 445–452.
- Gelman, A., Meng, X.-L. & Stern, H. (1996). Posterior predictive assessment
  of model fitness via realized discrepancies. *Statistica Sinica*, 6(4),
  733–760.
- Hilbe, J. M. (2011). *Negative Binomial Regression* (2nd ed.). Cambridge
  University Press.
- Lindgren, F., Rue, H. & Lindström, J. (2011). An explicit link between
  Gaussian fields and Gaussian Markov random fields: the stochastic partial
  differential equation approach. *Journal of the Royal Statistical Society:
  Series B*, 73(4), 423–498.
- Martin, T. G. et al. (2005). Zero tolerance ecology: improving ecological
  inference by modelling the source of zero observations. *Ecology Letters*,
  8(11), 1235–1246.
- Moran, P. A. P. (1950). Notes on continuous stochastic phenomena.
  *Biometrika*, 37(1/2), 17–23.
- O'Brien, T. G. (2011). Abundance, density and relative abundance: a
  conceptual framework. In *Camera Traps in Animal Ecology* (pp. 71–96).
  Springer.
- Roberts, D. R. et al. (2017). Cross-validation strategies for data with
  temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*,
  40(8), 913–929.
- Rowcliffe, J. M., Field, J., Turvey, S. T. & Carbone, C. (2008). Estimating
  animal density using camera traps without the need for individual
  recognition. *Journal of Applied Ecology*, 45(4), 1228–1236.
- Rue, H., Martino, S. & Chopin, N. (2009). Approximate Bayesian inference
  for latent Gaussian models by using integrated nested Laplace
  approximations. *Journal of the Royal Statistical Society: Series B*,
  71(2), 319–392.
- Simpson, D., Rue, H., Riebler, A., Sørbye, S. H. & Fuglstad, G.-A. (2017).
  Penalising model component complexity: a principled, practical approach to
  constructing priors. *Statistical Science*, 32(1), 1–28.
- Spiegelhalter, D. J., Best, N. G., Carlin, B. P. & van der Linde, A.
  (2002). Bayesian measures of model complexity and fit. *Journal of the
  Royal Statistical Society: Series B*, 64(4), 583–639.
- Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and
  widely applicable information criterion in singular learning theory.
  *Journal of Machine Learning Research*, 11, 3571–3594.
