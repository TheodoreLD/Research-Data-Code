# Data Inputs

Raw camera-trap data are not committed to this repository. To rerun the final
models, place the required CSV files in this `data/` folder, or set
`WOLF_DATA_DIR` to the folder that contains them.

## Required Files

For the forest-camera 2024 model:

```text
forest_camera_trap_events.csv
```

A custom path or filename can be supplied with `WOLF_FOREST_FILE`.

For the road-camera 2024 model:

```text
deployments_2024.csv
observations_2024.csv
```

For the road-camera 2023 model:

```text
deployments_2023.csv
observations_2023.csv
```

## Required Fields

The forest-camera flat file must contain:

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

The road-camera deployment file must contain:

```text
deploymentID
locationID
latitude
longitude
deploymentStart
deploymentEnd
```

The road-camera observation file must contain:

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

## How Camera Effort And Events Are Assigned To Months

All final models use camera-month rows. If a camera deployment spans more than
one month, the active camera-days are divided among the months in which the
camera was active. Wolf events are then counted in the month indicated by their
`eventStart` timestamp.

This is required because the statistical model includes month effects: the wolf
count for a month must be paired with the camera effort from that same month.

This keeps the count data and the exposure data aligned before fitting month
effects.
