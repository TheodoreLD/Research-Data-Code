# Wolf Relative Encounter Frequency, 2024

This project contains the final 2024 wolf relative encounter-frequency models
from camera-trap data. It is organized so a reader can understand the ecological
question, the input data structure, the two final statistical models, the
validation checks, and the outputs needed to reproduce or audit the analysis.

The response variable is the number of independent wolf event IDs recorded by
camera traps. Camera effort is used as an exposure term, so all predictions are
reported as expected wolf events per 100 camera-days. The maps are relative
encounter-frequency surfaces. They should not be interpreted as abundance,
density, occupancy, or population size.

## Final Models

Two 2024 analyses are included because they answer related but different
sampling questions.

| Survey | Final model | Rows | Cameras | Events | Effort | Final output |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Small/local 2024 | Negative-binomial spatial-month INLA-SPDE model | 356 camera-month rows | 53 | 46 | 4423.0 camera-days | `results/small_2024/` |
| Large 2024 | Zero-inflated negative-binomial spatial-month INLA-SPDE model | 344 camera-month rows | 60 | 479 | 3574.0 camera-days | `results/large_2024/` |

The small/local model is deliberately more conservative because the dataset
contains only 46 independent wolf events. The large 2024 model has enough data
to support a zero-inflated negative-binomial likelihood, and model comparison
favours this over negative-binomial and Poisson alternatives.

## Why Camera-Month Rows Matter

The final models use calendar-month temporal adjustment. That means effort and
events must be aligned to the month in which they occurred, not only to the
month in which a deployment started.

For both final 2024 models:

- deployment effort is split across calendar months when a deployment crosses a
  month boundary;
- wolf events are assigned to months using their `eventStart` timestamp;
- month fixed effects adjust for seasonal sampling differences;
- prediction maps are conditional on the selected prediction month.

This camera-month construction is especially important for the large 2024
survey, where many deployments crossed from September into October. After the
correction, formal equal-time temporal residual checks are no longer significant
at lag 1.

## Repository Layout

```text
WolfRelativeEncounterFrequency2024/
  README.md
  data/
    README.md
  docs/
    final-model-details.md
  results/
    README.md
    small_2024/
    large_2024/
  scripts/
    wolf_small_2024_month_refit.R
    wolf_2024_zinb_month_split_workflow.R
    wolf_relative_frequency_inla_helpers.R
```

Raw camera-trap CSV files are not committed here. The `data/README.md` file
lists the expected input files and where to place them for reproduction.

## Main Scripts

### Small/local 2024 final model

```sh
Rscript scripts/wolf_small_2024_month_refit.R
```

Final model:

- likelihood: negative binomial;
- spatial component: INLA-SPDE spatial random field;
- temporal component: calendar-month fixed effects;
- reference and prediction month: June 2024;
- spatial range is estimated with a weakly informative PC prior.

### Large 2024 final model

```sh
Rscript scripts/wolf_2024_zinb_month_split_workflow.R
```

Final model:

- likelihood: zero-inflated negative binomial type 1;
- spatial component: INLA-SPDE spatial random field;
- temporal component: calendar-month fixed effects;
- reference and prediction month: September 2024;
- effort is split into camera-month rows before fitting.

## Runtime Profiles

Both workflows use the `WOLF_RUN_PROFILE` environment variable:

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
$env:WOLF_PROJECT_DIR = "C:\path\to\WolfRelativeEncounterFrequency2024"
$env:WOLF_DATA_DIR = "C:\path\to\WolfRelativeEncounterFrequency2024\data"
$env:WOLF_OUTPUT_DIR = "C:\path\to\custom\outputs"
```

## Key Results

### Small/local 2024

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

### Large 2024

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

Prediction maps should be interpreted within the camera sampling domain and
conditional on the selected prediction month.

