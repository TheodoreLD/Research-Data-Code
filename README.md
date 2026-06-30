# Bayesian Spatial Encounter-Surface Modelling Of 2023-2024 Camera-Trap Detections

This repository contains a 2023-2024 camera-trap modelling project focused on
the spatial pattern of wolf relative encounter frequency. The analyses model
independent wolf detections recorded by camera traps using Bayesian count models
with INLA-SPDE spatial random fields. Active camera-days are used as exposure,
calendar month is included as a temporal adjustment, and outputs are relative
encounter-frequency surfaces expressed as expected independent wolf events per
100 camera-days across the sampled survey-year period.

## Project

- [CameraSpatialEncounterSurfaces2023_2024](CameraSpatialEncounterSurfaces2023_2024/):
  three Bayesian spatial encounter-surface models from camera-trap data: a
  2023 road-camera negative-binomial spatial-month model, a 2024 forest-camera
  negative-binomial spatial-month model, and a 2024 road-camera
  zero-inflated negative-binomial spatial-month model.
