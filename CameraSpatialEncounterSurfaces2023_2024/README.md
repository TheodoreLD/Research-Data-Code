# Bayesian Spatial Encounter-Surface Models For 2023-2024 Wolf Camera-Trap Detections

This project contains the final 2023-2024 wolf relative encounter-frequency models
from camera-trap data. It is organized so a reader can understand the ecological
question, the input data structure, the three final statistical models, the
validation checks, and the outputs needed to reproduce or audit the analysis.

The analyses model the number of independent wolf event IDs recorded in
camera-month rows. The statistical approach is Bayesian count modelling with
INLA-SPDE spatial random fields. Active camera-days are used as an exposure
term, calendar month is included as a temporal adjustment, and outputs are
relative encounter-frequency surfaces expressed as expected wolf events per 100
camera-days across the sampled survey-year period. The maps should not be
interpreted as abundance, density, occupancy, or population size.

## Final Models

Three camera-specific analyses are included:

| Survey | Final model | Rows | Cameras | Events | Effort | Final output |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Road-camera 2023 | Negative-binomial spatial-month INLA-SPDE model | 490 camera-month rows | 60 | 586 | 5222.2 camera-days | `results/road_2023/` |
| Forest-camera 2024 | Negative-binomial spatial-month INLA-SPDE model | 356 camera-month rows | 53 | 46 | 4423.0 camera-days | `results/forest_2024/` |
| Road-camera 2024 | Zero-inflated negative-binomial spatial-month INLA-SPDE model | 344 camera-month rows | 60 | 479 | 3574.0 camera-days | `results/road_2024/` |

Rows and events are different quantities. A row is one camera location observed
during one calendar month. The event count in a row can be zero, one, or many
independent wolf events. It is therefore normal for the total number of events
to be larger than the number of rows.

## Model Structure And Camera-Month Rows

The analysis is fitted at the camera-month scale. Each row represents one camera
location during one calendar month, with:

- `y_i`: the number of independent wolf events detected by camera `i` during
  that month;
- `E_i`: the number of active camera-days for that camera during that month;
- `s_i`: the spatial location of the camera;
- `m_i`: the calendar month assigned to the row.

If a deployment crosses a month boundary, the deployment is split before
modelling. For example, a camera active from late August into September
contributes one August row and one September row. August camera-days and August
events are assigned to the August row; September camera-days and September
events are assigned to the September row. This keeps event counts and exposure
aligned before fitting month effects.

The shared linear predictor is:

```text
log(mu_i) = log(E_i) + beta_0 + gamma[m_i] + u(s_i)
```

where:

- `mu_i` is the expected number of wolf events in camera-month row `i`;
- `log(E_i)` is an offset for camera effort, so cameras active for more days are
  expected to record more events;
- `beta_0` is the baseline log encounter rate for the reference month;
- `gamma[m_i]` is a fixed calendar-month effect used as a temporal control;
- `u(s_i)` is the spatial INLA-SPDE random field, which estimates a smooth
  spatial surface while allowing nearby camera locations to be correlated.

The forest-camera 2024 model and road-camera 2023 model use a
negative-binomial likelihood:

```text
y_i ~ NegativeBinomial(mu_i, size)
```

The road-camera 2024 model uses INLA's `zeroinflatednbinomial1` likelihood:

```text
y_i ~ ZeroInflatedNegativeBinomial(mu_i, size, pi)
```

The negative-binomial component allows event counts to be more variable than a
Poisson model. The zero-inflation component in the road-camera 2024 model
allows for additional zero counts beyond those expected from the
negative-binomial count process. In this type-1 parameterization, `pi` is the
probability of an additional structural zero, while `mu_i` and `size` describe
the negative-binomial count component.

The final map target is an effort-weighted annualized survey-year encounter-frequency
surface. Month remains in the model as a temporal control because residual
diagnostics indicated temporal structure when time was not handled explicitly.
The public maps are not intended to represent one selected calendar month. For
each prediction cell, the mapped daily encounter rate is:

```text
lambda_year(s) = sum_m w_m * 100 * exp(beta_0 + gamma[m] + u(s))
```

where `w_m` is the proportion of total sampled camera-days in month `m`. This
keeps the annual spatial pattern aligned with the months that were actually
sampled in each camera dataset. The reference month is retained only as the
baseline for coding the month coefficients.

Month is treated as a fixed effect in the final models. A random month effect
is possible in INLA, for example with `f(month, model = "iid")`, but it was not
used here because each dataset contains only a small number of sampled months.
With so few levels, a random-effect variance would be weakly identified and
would shrink month contrasts toward zero. Fixed month effects are a simpler
temporal control for these final models.

The reference month for coefficients is August in each year: August 2023 for
the 2023 model and August 2024 for both 2024 models. Changing the reference
month changes coefficient coding, but it does not change the annualized mapped
surface.

The mapped central estimate is the posterior mean because the reported target
is an expected encounter rate. Posterior-SD maps are included as the matching
uncertainty maps. Posterior medians can be useful for skewed posteriors, but
they answer a different question from the expected number of events.

## Repository Layout

```text
CameraSpatialEncounterSurfaces2023_2024/
  README.md
  data/
    README.md
  docs/
    final-model-details.md
  results/
    README.md
    road_2023/
    forest_2024/
    road_2024/
  scripts/
    wolf_2023_nb_month_split_workflow.R
    wolf_forest_month_refit.R
    wolf_2024_zinb_month_split_workflow.R
    wolf_relative_frequency_inla_helpers.R
```

Raw camera-trap CSV files are not committed here. The `data/README.md` file
lists the expected input files and where to place them for reproduction.

## Main Scripts

### Road-Camera 2023 Final Model

```sh
Rscript scripts/wolf_2023_nb_month_split_workflow.R
```

Final model:

- likelihood: negative binomial;
- spatial component: INLA-SPDE spatial random field;
- temporal component: calendar-month fixed effects;
- reference month for coefficients: August 2023;
- map target: effort-weighted annualized 2023 surface;
- effort is split into camera-month rows before fitting.

### Forest-Camera 2024 Final Model

```sh
Rscript scripts/wolf_forest_month_refit.R
```

Final model:

- likelihood: negative binomial;
- spatial component: INLA-SPDE spatial random field;
- temporal component: calendar-month fixed effects;
- reference month for coefficients: August 2024;
- map target: effort-weighted annualized 2024 surface;
- spatial range is estimated with a weakly informative PC prior.

### Road-Camera 2024 Final Model

```sh
Rscript scripts/wolf_2024_zinb_month_split_workflow.R
```

Final model:

- likelihood: zero-inflated negative binomial type 1;
- spatial component: INLA-SPDE spatial random field;
- temporal component: calendar-month fixed effects;
- reference month for coefficients: August 2024;
- map target: effort-weighted annualized 2024 surface;
- effort is split into camera-month rows before fitting.

## Runtime Profiles

All workflows use the `WOLF_RUN_PROFILE` environment variable:

```sh
# Fast development run
WOLF_RUN_PROFILE=quick

# Recommended reproducible analysis
WOLF_RUN_PROFILE=balanced

# Heavier final run with more posterior simulations
WOLF_RUN_PROFILE=final
```

On Windows PowerShell, for example:

```powershell
$env:WOLF_RUN_PROFILE = "balanced"
& "C:\Program Files\R\R-4.5.2\bin\Rscript.exe" scripts\wolf_2024_zinb_month_split_workflow.R
```

Path overrides are available:

```powershell
$env:WOLF_PROJECT_DIR = "C:\path\to\CameraSpatialEncounterSurfaces2023_2024"
$env:WOLF_DATA_DIR = "C:\path\to\CameraSpatialEncounterSurfaces2023_2024\data"
$env:WOLF_OUTPUT_DIR = "C:\path\to\outputs"
```

## Key Results

### Road-Camera 2023

The corrected negative-binomial spatial-month model passes the required
diagnostics:

- posterior predictive camera total events: pass;
- posterior predictive camera zero fraction: pass;
- posterior predictive camera maximum count: pass;
- residual Moran's I: `I = -0.008`, `p = 0.638`;
- row PIT KS p-value: `0.268`;
- temporal lag-1 residual correlation: `r = 0.018`, `p = 0.711`;
- formal 7-day lag-1 check: `r = -0.072`, `p = 0.573`;
- formal 14-day lag-1 check: `r = -0.046`, `p = 0.412`;
- spatial block CV row mean log predictive density: `-1.407`;
- spatial block CV row 90 percent coverage: `0.933`;
- spatial block CV camera 90 percent coverage: `0.900`;
- mesh sensitivity: final, finer, and coarser mesh variants pass required diagnostics;
- prior sensitivity: retained variants pass required diagnostics.

Model comparison shows ZINB is only marginally lower by WAIC, with low estimated
zero inflation, so the NB model is retained for parsimony:

| Model | WAIC | Delta WAIC |
| --- | ---: | ---: |
| ZINB spatial-month | 1162.72 | 0.00 |
| NB spatial-month | 1162.97 | 0.25 |
| Poisson spatial-month | 1303.16 | 140.44 |

The non-spatial Poisson and NB controls are retained in
`results/road_2023/wolf_2023_model_comparison.csv`, but the README table is
limited to the spatial candidate models used to choose the mapped likelihood.

### Forest-Camera 2024

The final negative-binomial spatial-month model passes the required
diagnostics:

- posterior predictive total events: pass;
- posterior predictive zero fraction: pass;
- posterior predictive maximum camera count: pass;
- residual Moran's I: `I = -0.041`, `p = 0.656`;
- PIT KS p-value: `0.952`;
- month-level residual lag-1 ACF: `0.050`;
- spatial block CV mean log predictive density: `-0.412`;
- spatial block CV 90 percent coverage: `0.978`;
- prior sensitivity: all retained variants pass required diagnostics;
- spatial prior/range sensitivity: fixed and estimated spatial-range variants
  are reported in `results/forest_2024/wolf_forest_month_prior_sensitivity.csv`.

The main caveat is low information content: 46 independent wolf events, so
month and spatial effects have wider uncertainty.

### Road-Camera 2024

The corrected zero-inflated negative-binomial spatial-month model passes the
required diagnostics:

- posterior predictive camera total events: pass;
- posterior predictive camera zero fraction: pass;
- posterior predictive camera maximum count: pass;
- residual Moran's I: `I = -0.033`, `p = 0.370`;
- row PIT KS p-value: `0.118`;
- temporal lag-1 residual correlation: `r = -0.178`, `p = 0.00267`;
- formal 7-day lag-1 check: `r = -0.153`, `p = 0.118`;
- formal 14-day lag-1 check: `r = -0.136`, `p = 0.094`;
- spatial block CV row mean log predictive density: `-1.524`;
- spatial block CV row 90 percent coverage: `0.962`;
- spatial block CV camera 90 percent coverage: `0.933`;
- mesh sensitivity: final, finer, and coarser mesh variants pass required diagnostics;
- prior sensitivity: retained variants are stable and pass required diagnostics.

Model comparison supports the ZINB model:

| Model | WAIC | Delta WAIC |
| --- | ---: | ---: |
| ZINB spatial-month | 933.67 | 0.00 |
| NB spatial-month | 937.32 | 3.65 |
| Poisson spatial-month | 997.19 | 63.52 |

The non-spatial Poisson and NB controls are retained in
`results/road_2024/wolf_2024_model_comparison.csv`, but the README table is
limited to the spatial candidate models used to choose the mapped likelihood.
The deployment-order temporal correlation is retained as a warning-style
diagnostic because the gaps between observations are irregular. The equal-time
7-day and 14-day checks are the formal lag checks and are not significant.

## Outputs Included Here

Each final-results folder contains:

- final validation report;
- posterior predictive checks;
- hyperparameter summaries;
- month-effect summaries;
- prior sensitivity report and table;
- mesh sensitivity report and table for the road-camera models;
- model-comparison report and table for the road-camera models;
- spatial block cross-validation summary;
- temporal diagnostics;
- posterior mean encounter-frequency map as PNG and GeoTIFF;
- posterior-SD uncertainty map as PNG and GeoTIFF.

The full generated output folders also contain exploratory plots, full
prediction grids, and additional intermediate diagnostics. Those full scratch
archives are not part of the curated GitHub result set.

## Required R Packages

- `readr`
- `dplyr`
- `tidyr`
- `sf`
- `terra`
- `ggplot2`
- `viridis`
- `scales`
- `INLA`

The scripts do not install packages automatically unless explicitly changed by
the user. INLA usually requires installing from the INLA repository.
