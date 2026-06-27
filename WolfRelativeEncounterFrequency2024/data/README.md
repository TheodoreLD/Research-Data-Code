# Data Inputs

Raw camera-trap data are not committed to this repository. To rerun the final
models, place the required CSV files in this `data/` folder, or set
`WOLF_DATA_DIR` to the folder that contains them.

## Required Files

For the small/local 2024 model:

```text
theodata_1.0(1).csv
```

For the large 2024 model:

```text
deployments_2024.csv
observations_2024.csv
```

The helper workflow also knows about 2023 files, but the final project described
here is focused on the two final 2024 models.

## Required Fields

The small/local flat file must contain:

```text
deploymentID
eventID
eventStart
scientificName
plotID
deploymentEffort
latitude
longitude
startDate
endDate
```

The large 2024 deployment file must contain:

```text
deploymentID
locationID
latitude
longitude
deploymentStart
deploymentEnd
```

The large 2024 observation file must contain:

```text
deploymentID
eventID
scientificName
eventStart
```

## Event Definition

The scripts count distinct wolf `eventID` values where `scientificName` is
`Canis_lupus` or `Canis lupus`. The analysis assumes each event ID represents an
independent wolf event according to the camera-trap processing workflow.

## Temporal Alignment

Both final models split deployment effort by calendar month and assign wolf
events by `eventStart` month. This avoids assigning all effort and events from a
long deployment to only the deployment-start month.

