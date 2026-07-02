# Results Included In This Repository

The `results/` directory contains compact final outputs for the three final
models. It is meant for review and interpretation, not as a complete
generated-output archive. For methodology, priors, and full diagnostic
numbers, see [`docs/final-model-details.md`](../docs/final-model-details.md).

## Road-Camera 2023

Folder:

```text
results/road_2023/
```

Key files:

- `wolf_2023_validation_report.txt` — model, diagnostic status, priors, CV
- `wolf_2023_run_manifest.csv`
- `wolf_2023_science_checks_summary.txt`
- `wolf_2023_model_choice_report.txt`
- `wolf_2023_model_comparison_report.txt` / `wolf_2023_model_comparison.csv`
- `wolf_2023_prior_sensitivity_report.txt` / `wolf_2023_prior_sensitivity.csv`
- `wolf_2023_mesh_sensitivity_report.txt` / `wolf_2023_mesh_sensitivity.csv`
- `wolf_2023_nb_spatial_month_temporal_autocorrelation_report.txt`
- `wolf_2023_exploratory_timing_vs_northing.csv` — Spearman correlation of
  deployment start day-of-year against UTM northing
- `wolf_2023_final_spatial_block_cv_summary.csv`
- `wolf_2023_annualization_weights.csv`
- `wolf_2023_hyperparameters.csv`
- `wolf_2023_month_coefficients.csv`
- `wolf_2023_month_observed_summary.csv`
- `wolf_2023_nb_spatial_month_posterior_predictive_check.csv`
- `wolf_2023_final_event_frequency_mean.png` / `.tif`
- `wolf_2023_final_event_frequency_sd.png` / `wolf_2023_final_predicted_events_per_100_days_sd.tif`

## Forest-Camera 2024

Folder:

```text
results/forest_2024/
```

Key files:

- `wolf_forest_2024_full_final_model_report.txt` — consolidated model,
  diagnostic, and prior/mesh sensitivity summary
- `wolf_forest_2024_validation_report.txt`
- `wolf_forest_2024_run_manifest.csv`
- `wolf_forest_2024_prior_sensitivity_report.txt` / `wolf_forest_2024_prior_sensitivity.csv`
- `wolf_forest_2024_mesh_sensitivity_report.txt` / `wolf_forest_2024_mesh_sensitivity.csv`
- `wolf_forest_2024_final_spatial_block_cv_summary.csv`
- `wolf_forest_2024_nb_spatial_month_temporal_autocorrelation_report.txt`
- `wolf_forest_2024_temporal_residual_diagnostics.csv`
- `wolf_forest_2024_temporal_within_camera_lag_correlation.csv`
- `wolf_forest_2024_annualization_weights.csv`
- `wolf_forest_2024_hyperparameters.csv`
- `wolf_forest_2024_month_coefficients.csv`
- `wolf_forest_2024_month_observed_summary.csv`
- `wolf_forest_2024_nb_spatial_month_posterior_predictive_check.csv`
- `wolf_forest_2024_final_event_frequency_mean.png` / `.tif`
- `wolf_forest_2024_final_event_frequency_sd.png` / `wolf_forest_2024_final_predicted_events_per_100_days_sd.tif`

## Road-Camera 2024

Folder:

```text
results/road_2024/
```

Key files:

- `wolf_2024_validation_report.txt`
- `wolf_2024_run_manifest.csv`
- `wolf_2024_science_checks_summary.txt`
- `wolf_2024_model_choice_report.txt`
- `wolf_2024_model_comparison_report.txt` / `wolf_2024_model_comparison.csv`
- `wolf_2024_prior_sensitivity_report.txt` / `wolf_2024_prior_sensitivity.csv`
- `wolf_2024_mesh_sensitivity_report.txt` / `wolf_2024_mesh_sensitivity.csv`
- `wolf_2024_zinb_spatial_month_temporal_autocorrelation_report.txt`
- `wolf_2024_exploratory_timing_vs_northing.csv` — Spearman correlation of
  deployment start day-of-year against UTM northing; used to test the
  hypothesized mechanism for this survey's residual temporal autocorrelation
- `wolf_2024_final_spatial_block_cv_summary.csv`
- `wolf_2024_annualization_weights.csv`
- `wolf_2024_hyperparameters.csv`
- `wolf_2024_month_coefficients.csv`
- `wolf_2024_month_observed_summary.csv`
- `wolf_2024_zinb_spatial_month_posterior_predictive_check.csv`
- `wolf_2024_final_event_frequency_mean.png` / `.tif`
- `wolf_2024_final_event_frequency_sd.png` / `wolf_2024_final_predicted_events_per_100_days_sd.tif`

## Not Included By Default

The full generated output folders contain full prediction-grid CSV files,
exploratory figures, and many intermediate diagnostic files. These are omitted
here to keep the GitHub project focused on the final diagnostics and map
products.
