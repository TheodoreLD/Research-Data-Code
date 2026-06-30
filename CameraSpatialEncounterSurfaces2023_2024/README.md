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
| Forest-camera 2024 | Negative-binomial spatial-month INLA-SPDE model | 356 camera-month rows | 53 | 46 | 4423.0 camera-days | `results/forest/` |
| Road-camera 2024 | Zero-inflated negative-binomial spatial-month INLA-SPDE model | 344 camera-month rows | 60 | 479 | 3574.0 camera-days | `results/road/` |

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
- `gamma[m_i]` is a fixed month effect, estimating how the encounter rate in
  month `m_i` differs from the reference month;
- `u(s_i)` is the spatial INLA-SPDE random field, which estimates a smooth
  spatial surface while allowing nearby camera locations to be correlated.

The forest-camera 2024 model and road-camera 2023 model use a
negative-binomial likelihood:

```text
y_i ~ NegativeBinomial(mu_i, size)
```

The road-camera 2024 model uses a zero-inflated negative-binomial likelihood:

```text
y_i ~ ZeroInflatedNegativeBinomial(mu_i, size, pi)
```

The negative-binomial component allows event counts to be more variable than a
Poisson model. The zero-inflation component in the road-camera 2024 model
allows for additional zero counts beyond those expected from the
negative-binomial count process.

The final map target is an effort-weighted annualized survey-year encounter-frequency
surface. Month remains in the model as an adjustment for seasonal differences in
encounter rate and sampling effort, but the public maps are not intended to
represent one selected calendar month. For each prediction cell, the mapped
daily encounter rate is:

```text
lambda_year(s) = sum_m w_m * 100 * exp(beta_0 + gamma[m] + u(s))
```

where `w_m` is the proportion of total sampled camera-days in month `m`. This
keeps the annual spatial pattern aligned with the months that were actually
sampled in each camera dataset. The reference month is retained only as the
baseline for coding the month coefficients.

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
    forest/
    road/
    road_2023/
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
- reference month for coefficients: June 2024;
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
- reference month for coefficients: September 2024;
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
$env:WOLF_OUTPUT_DIR = "C:\path\to\custom\outputs"
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
- spatial block CV row mean log predictive density: `-1.407`;
- spatial block CV camera 90 percent coverage: `0.900`.

Model comparison shows ZINB is only marginally lower by WAIC, with low estimated
zero inflation, so the NB model is retained for parsimony:

| Model | WAIC | Delta WAIC |
| --- | ---: | ---: |
| ZINB spatial-month | 1162.72 | 0.00 |
| NB spatial-month | 1162.97 | 0.25 |
| Poisson spatial-month | 1303.16 | 140.44 |
| NB non-spatial month | 1348.70 | 185.98 |
| Poisson non-spatial month | 1861.53 | 698.81 |

The formal equal-time lag-1 temporal checks are not significant:

- 7-day lag: `r = -0.072`, `p = 0.573`;
- 14-day lag: `r = -0.046`, `p = 0.412`.

### Forest-Camera 2024

The final negative-binomial spatial-month model passes the required
diagnostics:

- posterior predictive total events: pass;
- posterior predictive zero fraction: pass;
- posterior predictive maximum camera count: pass;
- residual Moran's I: `I = -0.041`, `p = 0.656`;
- PIT KS p-value: `0.952`;
- spatial block CV mean log predictive density: `-0.415`;
- month-level residual lag-1 ACF: `0.050`.

All prior sensitivity variants passed required diagnostics. The main caveat is
low information content: 46 independent wolf events, so month and spatial
effects have wide uncertainty.

### Road-Camera 2024

The corrected zero-inflated negative-binomial spatial-month model passes the
required diagnostics:

- posterior predictive camera total events: pass;
- posterior predictive camera zero fraction: pass;
- posterior predictive camera maximum count: pass;
- residual Moran's I: `I = -0.034`, `p = 0.336`;
- row PIT KS p-value: `0.0905`;
- spatial block CV row mean log predictive density: `-1.513`;
- spatial block CV camera 90 percent coverage: `0.933`.

Model comparison supports the ZINB model:

| Model | WAIC | Delta WAIC |
| --- | ---: | ---: |
| ZINB spatial-month | 933.67 | 0.00 |
| NB spatial-month | 937.27 | 3.60 |
| Poisson spatial-month | 997.23 | 63.56 |
| NB non-spatial month | 1042.52 | 108.85 |
| Poisson non-spatial month | 1327.08 | 393.41 |

The formal equal-time lag-1 temporal checks are not significant:

- 7-day lag: `r = -0.153`, `p = 0.118`;
- 14-day lag: `r = -0.136`, `p = 0.094`.

Weekly and biweekly time-bin models were fitted as sensitivity checks, but the
month model remained best by WAIC and was retained.

## Outputs Included Here

Each final-results folder contains:

- final validation report;
- posterior predictive checks;
- hyperparameter summaries;
- month-effect summaries;
- prior sensitivity outputs;
- spatial block cross-validation summary;
- temporal diagnostics;
- final map PNGs for mean encounter frequency, coefficient of variation, and
  exceedance probability.

Large GeoTIFF rasters and full prediction grids are not committed here by
default. They can be regenerated from the scripts.

## Public Release And Publication Notes

Raw camera-trap files and exact coordinate tables are not included in this
repository. The committed result tables contain model summaries and diagnostics,
not per-camera coordinates or raw event records.

The committed map PNGs show the relative spatial surface and camera positions,
but they do not include numeric coordinate ticks or raw coordinate tables. If
camera locations are considered sensitive, use anonymized publication figures or
restrict the raw spatial data through the journal data-availability statement.

For a manuscript-linked release, archive a tagged repository version and cite
that version in the paper. Add the final repository license and citation metadata
once authorship, journal requirements, and data-sharing constraints are fixed.

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

## Interpretation Limits

These analyses estimate relative encounter frequency. High predicted values may
reflect wolf activity, camera placement, trail or road use, detection
conditions, or camera effort. The models do not estimate abundance or density,
and they do not include habitat, prey, road, human disturbance, or detection
covariates.

Prediction maps should be interpreted within the camera sampling domain as
month-adjusted, effort-weighted annualized survey-year encounter-frequency
surfaces.
