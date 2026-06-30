###############################################################################
# Wolf relative encounter frequency: 2024 INLA-SPDE ordered science workflow
# -----------------------------------------------------------------------------
# Purpose
#   Standalone 2024 workflow for modelling relative wolf encounter frequency
#   from camera-trap data, producing prediction maps and validation diagnostics.
#
# Model
#   Response: independent wolf eventID count per camera-month row.
#   Exposure: camera effort split across calendar months, passed to INLA through E.
#   Final 2024 model: zero-inflated negative-binomial type 1 likelihood,
#   calendar-month fixed effects, and a spatial SPDE field. Wolf events are
#   assigned to months using eventStart timestamps.
#
# Interpretation
#   Predictions are relative encounter frequency: wolf events per 100 camera-days.
#   They are not abundance, density, occupancy, or population size.
#
# Required inputs
#   deployments_2024.csv
#   observations_2024.csv
#
# Key corrections used in this final 2024 workflow
#   * Uses joint posterior samples for PPC, PIT, fitted counts, and CV.
#   * Uses camera-level and deployment-row-level validation summaries.
#   * Uses a two-sided Moran permutation test around the expected value.
#   * Computes map exceedance from the marginal posterior linear predictor, avoiding unstable full-grid latent sampling.
#   * Adds fitted-scale sanity checks against INLA fitted means when available.
#
# Prior sensitivity and science-check version
#   This v2 2024 script uses the configured final ZINB spatial-month model.
#   It still runs model comparison against Poisson/NB candidates, a prior-influence
#   screen, and targeted prior-sensitivity checks before final interpretation.
#   The intercept prior is centered after loading the data at the crude observed
#   daily event rate, with a broad SD.
#
# Additional science checks in v2
#   * Model comparison among spatial Poisson/NB/ZINB variants.
#   * Mesh sensitivity checks for the final spatial model.
#   * Full convex-hull prediction map only; disk/domain-sensitivity mapping has been removed.
#   * A short scientific limitations report for interpretation.
###############################################################################


## 01. User settings ----------------------------------------------------------

input_files_required <- c("deployments_2024.csv", "observations_2024.csv")

script_file <- tryCatch({
  ofile <- sys.frames()[[1]]$ofile
  if (is.null(ofile)) NA_character_ else normalizePath(ofile, winslash = "/", mustWork = FALSE)
}, error = function(e) NA_character_)

cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (is.na(script_file) && length(cmd_file)) {
  script_file <- normalizePath(sub("^--file=", "", cmd_file[[1]]),
                               winslash = "/", mustWork = FALSE)
}

script_dir <- if (is.na(script_file)) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  dirname(script_file)
}

has_input_layout <- function(dir) {
  all(file.exists(file.path(dir, input_files_required))) ||
    all(file.exists(file.path(dir, "data", input_files_required)))
}

project_candidates <- unique(normalizePath(c(script_dir, getwd()),
                                           winslash = "/", mustWork = FALSE))
if (basename(script_dir) %in% c("scripts", "R")) {
  project_candidates <- unique(c(
    normalizePath(dirname(script_dir), winslash = "/", mustWork = FALSE),
    project_candidates
  ))
}
detected_project <- project_candidates[
  vapply(project_candidates, has_input_layout, logical(1))
]
default_project_dir <- if (length(detected_project)) detected_project[[1]] else getwd()

PROJECT_DIR <- Sys.getenv("WOLF_PROJECT_DIR", unset = default_project_dir)
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)

DATA_DIR <- Sys.getenv("WOLF_DATA_DIR", unset = "")
if (!nzchar(DATA_DIR)) {
  data_subdir <- file.path(PROJECT_DIR, "data")
  DATA_DIR <- if (dir.exists(data_subdir)) data_subdir else PROJECT_DIR
}
DATA_DIR <- normalizePath(DATA_DIR, winslash = "/", mustWork = FALSE)

OUTPUT_DIR <- Sys.getenv(
  "WOLF_OUTPUT_DIR",
  unset = file.path(PROJECT_DIR, "outputs", "2024", "wolf_2024_ZINB_month_split_final_v1")
)
OUTPUT_DIR <- normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Spatial CRS for modelling and maps.
EPSG_UTM <- 32634L

# Species labels in the source files.
WOLF_NAMES <- c("Canis_lupus", "Canis lupus")

# Set TRUE only if you explicitly want this script to install missing packages.
INSTALL_MISSING_PACKAGES <- FALSE

# Runtime profile:
#   quick    : fast testing; fewer posterior samples and no spatial CV
#   balanced : recommended default
#   final    : heavier publication rerun
RUN_PROFILE <- tolower(Sys.getenv("WOLF_RUN_PROFILE", unset = "balanced"))
if (!RUN_PROFILE %in% c("quick", "balanced", "final")) {
  stop("WOLF_RUN_PROFILE must be one of: quick, balanced, final.")
}

PPC_NSIM <- switch(RUN_PROFILE, quick = 200L, balanced = 750L, final = 1500L)
# Retained for compatibility; maps use marginal posterior summaries to avoid unstable full-grid latent sampling.
PRED_NSIM <- switch(RUN_PROFILE, quick = 100L, balanced = 300L, final = 600L)
CV_NSIM <- switch(RUN_PROFILE, quick = 150L, balanced = 300L, final = 600L)
RUN_SPATIAL_CV <- RUN_PROFILE != "quick"
CV_K <- switch(RUN_PROFILE, quick = 3L, balanced = 4L, final = 5L)

# Diagnostics.
MORAN_ALPHA <- 0.05
MORAN_NPERM <- switch(RUN_PROFILE, quick = 199L, balanced = 499L, final = 999L)

# Prediction domain.
# Full map: buffered convex hull around all cameras. Disk-based maps are not produced.
PRED_DOMAIN <- "hull"
MAP_EXCEEDANCE <- FALSE
EXCEED_MULT <- 1.5

# Final 2024 model settings.
SURVEY_LABEL <- "Road-camera 2024 survey"
SURVEY_PREFIX <- "wolf_2024"
FINAL_FAMILY <- "zeroinflatednbinomial1"
FINAL_MODEL_NAME <- "zinb_spatial_month"
# Configured final model for 2024 after model comparison: ZINB spatial-month.
MONTH_REFERENCE <- "2024-08"
MONTH_PREDICTION <- "2024-08"

settings <- list(
  cell_size_m = 150,
  pred_buffer_m = 1500,
  max_dist_m = 2500,
  mesh_cutoff_m = 350,
  mesh_max_edge = c(700, 5000),
  mesh_offset = c(5000, 15000),
  fix_range_m = NULL,
  prior_range_m = c(5000, 0.5),       # P(range < 5000 m) = 0.5
  prior_sigma = c(2.50, 0.05),        # final 2024 prior: P(sigma > 2.00) = 0.05; widened after prior-influence screen
  include_grid_in_mesh = FALSE,
  use_month_effect = TRUE,
  month_reference = MONTH_REFERENCE,
  month_prediction = MONTH_PREDICTION
)

# Fixed-effect and likelihood priors.
# Intercept prior mean is filled after loading the 2024 data, using the crude
# observed daily event rate. The prior SD is deliberately broad but finite.
PRIOR_INTERCEPT_MEAN <- NA_real_
PRIOR_INTERCEPT_PREC <- 1 / 2.5^2

# Month effects remain weakly informative: SD = 1 on log-rate-ratio scale.
PRIOR_MONTH_LOG_RATE_RATIO_PREC <- 1

# Zero inflation: skeptical prior, but not hard-constrained.
# Center: 5% structural-zero probability; SD = 1.5 on logit scale.
PRIOR_ZI_LOGIT_MEAN <- qlogis(0.05)  # skeptical but flexible ZINB prior
PRIOR_ZI_LOGIT_PREC <- 1 / 1.5^2

# Negative-binomial size prior for the final ZINB likelihood. Smaller size = stronger overdispersion.
PRIOR_NB_LOGSIZE_MEAN <- log(2)
PRIOR_NB_LOGSIZE_PREC <- 1 / 2^2

set.seed(1)


## 02. Packages ---------------------------------------------------------------

cat("Project directory: ", PROJECT_DIR, "\n", sep = "")
cat("Data directory:    ", DATA_DIR, "\n", sep = "")
cat("Output directory:  ", OUTPUT_DIR, "\n", sep = "")
cat(sprintf(
  "Run profile:       %s | PPC_NSIM=%d | PRED_NSIM=%d | spatial_CV=%s | CV_K=%d | CV_NSIM=%d\n",
  RUN_PROFILE, PPC_NSIM, PRED_NSIM, RUN_SPATIAL_CV, CV_K, CV_NSIM
))

required_packages <- c(
  "readr", "dplyr", "tidyr", "sf", "terra",
  "ggplot2", "viridis", "scales", "INLA"
)

for (p in required_packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (!INSTALL_MISSING_PACKAGES) {
      stop("Package '", p, "' is missing. Install it or set INSTALL_MISSING_PACKAGES <- TRUE.")
    }
    if (p == "INLA") {
      install.packages(
        "INLA",
        dep = TRUE,
        repos = c(getOption("repos"),
                  INLA = "https://inla.r-inla-download.org/R/stable")
      )
    } else {
      install.packages(p)
    }
  }
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(sf)
  library(terra)
  library(ggplot2)
  library(viridis)
  library(scales)
  library(INLA)
})

sf::sf_use_s2(FALSE)
try(INLA::inla.setOption(fmesher.evolution.warn = FALSE), silent = TRUE)


## 03. General helpers --------------------------------------------------------

path_in <- function(...) file.path(DATA_DIR, ...)
path_out <- function(...) file.path(OUTPUT_DIR, ...)

stop_missing_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(label, " missing required column(s): ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

first_finite <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) x[[1]] else NA_real_
}

parse_time <- function(x) {
  if (inherits(x, "POSIXt")) return(x)
  if (inherits(x, "Date")) return(as.POSIXct(x, tz = "UTC"))

  s <- as.character(x)
  s <- gsub("T", " ", s, fixed = TRUE)
  s <- sub("([+-][0-9]{2}:?[0-9]{2}|Z)$", "", s)

  out <- as.POSIXct(substr(s, 1, 19),
                    format = "%Y-%m-%d %H:%M:%S",
                    tz = "UTC")
  bad <- is.na(out) & !is.na(s) & nzchar(s)
  if (any(bad)) {
    out[bad] <- as.POSIXct(substr(s[bad], 1, 10),
                           format = "%Y-%m-%d",
                           tz = "UTC")
  }
  out
}

log_mean_exp <- function(x) {
  m <- max(x, na.rm = TRUE)
  if (!is.finite(m)) return(NA_real_)
  m + log(mean(exp(x - m), na.rm = TRUE))
}

safe_cor <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3 || sd(a[ok]) == 0 || sd(b[ok]) == 0) return(NA_real_)
  cor(a[ok], b[ok])
}

safe_cor_p_value <- function(a, b, min_pairs = 5L) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < min_pairs || sd(a[ok]) == 0 || sd(b[ok]) == 0) return(NA_real_)
  suppressWarnings(tryCatch(
    stats::cor.test(a[ok], b[ok])$p.value,
    error = function(e) NA_real_
  ))
}

as_utc_date <- function(x) {
  if (inherits(x, "POSIXt")) return(as.Date(x, tz = "UTC"))
  as.Date(x)
}


## 04. Family, prior, and hyperparameter helpers ------------------------------

fit_family <- function(family) {
  key <- gsub("[^a-z0-9]", "", tolower(family))
  if (key %in% c("poisson")) {
    "poisson"
  } else if (key %in% c("nbinomial", "negativebinomial", "nb")) {
    "nbinomial"
  } else if (key %in% c("zeroinflatedpoisson1", "zip", "zip1")) {
    "zeroinflatedpoisson1"
  } else if (key %in% c("zeroinflatednbinomial1", "zinb", "zinb1")) {
    "zeroinflatednbinomial1"
  } else {
    family
  }
}

is_zi <- function(family) grepl("zeroinflated", family)
is_nb <- function(family) grepl("nbinomial", family)

fam_nb_size <- function(size) {
  ifelse(is.finite(size) & size > 0, size, 1e6)
}

fam_mean <- function(mu, pi, family, size = NA_real_) {
  p <- ifelse(is_zi(family) & is.finite(pi), pi, 0)
  (1 - p) * mu
}

fam_var <- function(mu, pi, family, size = NA_real_) {
  p <- ifelse(is_zi(family) & is.finite(pi), pi, 0)
  base_var <- if (is_nb(family)) {
    s <- fam_nb_size(size)
    mu + mu^2 / s
  } else {
    mu
  }
  (1 - p) * base_var + p * (1 - p) * mu^2
}

fam_logpmf <- function(y, mu, pi, family, size = NA_real_) {
  p <- ifelse(is_zi(family) & is.finite(pi), pi, 0)
  base_prob <- if (is_nb(family)) {
    dnbinom(y, mu = mu, size = fam_nb_size(size))
  } else {
    dpois(y, mu)
  }
  log(pmax(p * (y == 0) + (1 - p) * base_prob,
           .Machine$double.xmin))
}

fam_sim <- function(mu, pi, family, size = NA_real_) {
  y <- if (is_nb(family)) {
    rnbinom(length(mu), mu = mu, size = fam_nb_size(size))
  } else {
    rpois(length(mu), mu)
  }
  if (is_zi(family) && is.finite(pi) && pi > 0) {
    y * (1 - rbinom(length(mu), 1, pi))
  } else {
    y
  }
}

PAT_ZPROB <- "zero.*prob|probability.*zero|zero-probability"
PAT_NB_SIZE <- "size.*nbinomial|size for|size parameter"
PAT_RANGE <- "range.*spatial|range for spatial|spatial.*range"
PAT_SIGMA <- "stdev.*spatial|stdev for spatial|standard deviation.*spatial|sigma.*spatial"

hyp_point <- function(fit, pattern) {
  if (is.null(fit$summary.hyperpar)) return(NA_real_)
  i <- grep(pattern, rownames(fit$summary.hyperpar), ignore.case = TRUE)
  if (length(i)) fit$summary.hyperpar[i[[1]], "mean"] else NA_real_
}

hyp_marg <- function(fit, pattern) {
  if (is.null(fit$marginals.hyperpar)) return(NULL)
  i <- grep(pattern, names(fit$marginals.hyperpar), ignore.case = TRUE)
  if (length(i)) fit$marginals.hyperpar[[i[[1]]]] else NULL
}

nb_size_point <- function(fit) {
  x <- hyp_point(fit, PAT_NB_SIZE)
  if (is.finite(x) && x > 0) x else NA_real_
}

make_control_family <- function(family) {
  hyper <- list()

  if (is_nb(family)) {
    hyper$size <- list(
      prior = "gaussian",
      param = c(PRIOR_NB_LOGSIZE_MEAN, PRIOR_NB_LOGSIZE_PREC)
    )
  }

  if (is_zi(family)) {
    hyper$prob <- list(
      prior = "gaussian",
      param = c(PRIOR_ZI_LOGIT_MEAN, PRIOR_ZI_LOGIT_PREC)
    )
  }

  if (length(hyper)) list(hyper = hyper) else list()
}

make_control_fixed <- function(fixed_terms = "intercept") {
  fixed_terms <- unique(fixed_terms)

  mean_values <- rep(0, length(fixed_terms))
  names(mean_values) <- fixed_terms
  if ("intercept" %in% names(mean_values)) {
    if (!is.finite(PRIOR_INTERCEPT_MEAN)) {
      stop("PRIOR_INTERCEPT_MEAN must be set before fitting.")
    }
    mean_values[["intercept"]] <- PRIOR_INTERCEPT_MEAN
  }
  mean_prior <- as.list(mean_values)
  mean_prior$default <- 0

  prec_values <- rep(PRIOR_MONTH_LOG_RATE_RATIO_PREC, length(fixed_terms))
  names(prec_values) <- fixed_terms
  if ("intercept" %in% names(prec_values)) {
    prec_values[["intercept"]] <- PRIOR_INTERCEPT_PREC
  }

  prec_prior <- as.list(prec_values)
  prec_prior$default <- PRIOR_MONTH_LOG_RATE_RATIO_PREC

  list(mean = mean_prior, prec = prec_prior)
}

month_term_name <- function(month) {
  paste0("month_", gsub("[^A-Za-z0-9]", "_", month))
}

month_from_term <- function(term) {
  gsub("_", "-", sub("^month_", "", term), fixed = TRUE)
}

temporal_month_terms <- function(data) {
  grep("^month_[0-9]{4}_[0-9]{2}$", names(data), value = TRUE)
}

time_bin_term_name <- function(bin_label) {
  paste0("timebin_", gsub("[^A-Za-z0-9]", "_", bin_label))
}

time_bin_from_term <- function(term) {
  gsub("_", "-", sub("^timebin_", "", term), fixed = TRUE)
}

temporal_time_bin_terms <- function(data) {
  grep("^timebin_[0-9]{4}_[0-9]{2}_[0-9]{2}$", names(data), value = TRUE)
}

fixed_effect_terms <- function(data) {
  c("intercept", temporal_month_terms(data), temporal_time_bin_terms(data))
}


## 05. Data preparation -------------------------------------------------------

validate_inputs <- function() {
  missing <- input_files_required[!file.exists(path_in(input_files_required))]
  if (length(missing)) {
    stop("Missing input file(s): ", paste(missing, collapse = ", "),
         " in ", DATA_DIR)
  }
  invisible(TRUE)
}

month_reference_from_settings <- function(months, settings) {
  ref <- settings$month_reference
  if (is.null(ref) || !nzchar(ref)) ref <- months[[1]]
  if (!ref %in% months) {
    stop("month_reference '", ref, "' is not present in deployment months: ",
         paste(months, collapse = ", "))
  }
  ref
}

month_prediction_from_settings <- function(months, settings) {
  pred <- settings$month_prediction
  if (is.null(pred) || !nzchar(pred)) pred <- month_reference_from_settings(months, settings)
  if (!pred %in% months) {
    stop("month_prediction '", pred, "' is not present in deployment months: ",
         paste(months, collapse = ", "))
  }
  pred
}

add_month_design <- function(model_dat, settings) {
  months <- sort(unique(model_dat$month))
  if (length(months) < 2) {
    stop("Month effect requires at least two deployment months.")
  }

  reference_month <- month_reference_from_settings(months, settings)
  prediction_month <- month_prediction_from_settings(months, settings)

  for (m in setdiff(months, reference_month)) {
    term <- month_term_name(m)
    model_dat[[term]] <- as.integer(model_dat$month == m)
  }

  model_dat$month_reference <- reference_month
  model_dat$month_prediction <- prediction_month
  model_dat$model_row_type <- "deployment_month"
  model_dat
}

drop_temporal_fixed_effect_design <- function(model_dat) {
  drop_terms <- c(temporal_month_terms(model_dat), temporal_time_bin_terms(model_dat))
  if (length(drop_terms)) {
    model_dat <- model_dat[, setdiff(names(model_dat), drop_terms), drop = FALSE]
  }
  model_dat
}

add_time_bin_design <- function(model_dat, bin_days = 14L, prediction_date = NULL) {
  bin_days <- as.integer(bin_days)
  if (!is.finite(bin_days) || bin_days <= 0) {
    stop("bin_days must be a positive integer.")
  }

  start_date <- as_utc_date(model_dat$start)
  if (all(is.na(start_date))) stop("Cannot build time-bin design: start dates are missing.")

  origin <- min(start_date, na.rm = TRUE)
  bin_index <- as.integer(floor(as.numeric(start_date - origin) / bin_days) + 1L)
  bin_start <- origin + (bin_index - 1L) * bin_days
  bin_label <- format(bin_start, "%Y-%m-%d")

  reference_label <- NULL
  if (!is.null(prediction_date) && length(prediction_date)) {
    prediction_date <- as.Date(prediction_date[[1]])
    if (!is.na(prediction_date)) {
      pred_index <- as.integer(floor(as.numeric(prediction_date - origin) / bin_days) + 1L)
      pred_start <- origin + (pred_index - 1L) * bin_days
      candidate <- format(pred_start, "%Y-%m-%d")
      if (candidate %in% bin_label) reference_label <- candidate
    }
  }
  if (is.null(reference_label)) {
    bin_effort <- tapply(model_dat$total_effort_days, bin_label, sum, na.rm = TRUE)
    reference_label <- names(sort(bin_effort, decreasing = TRUE))[[1]]
  }

  model_dat <- drop_temporal_fixed_effect_design(model_dat)
  model_dat$time_bin_days <- bin_days
  model_dat$time_bin_index <- bin_index
  model_dat$time_bin_start <- bin_start
  model_dat$time_bin <- bin_label
  model_dat$time_bin_reference <- reference_label
  model_dat$time_bin_prediction <- reference_label
  model_dat$model_row_type <- paste0("deployment_", bin_days, "day_timebin")

  for (b in setdiff(sort(unique(bin_label)), reference_label)) {
    model_dat[[time_bin_term_name(b)]] <- as.integer(bin_label == b)
  }

  model_dat
}

camera_summary_from_model <- function(model_dat) {
  has_deployment_id <- "deploymentID" %in% names(model_dat)
  has_n_deployments <- "n_deployments" %in% names(model_dat)

  model_dat %>%
    group_by(plotID) %>%
    summarise(
      longitude = mean(longitude, na.rm = TRUE),
      latitude = mean(latitude, na.rm = TRUE),
      total_effort_days = sum(total_effort_days, na.rm = TRUE),
      wolf_events = sum(wolf_events, na.rm = TRUE),
      n_model_rows = n(),
      n_deployments = if (has_deployment_id) {
        n_distinct(deploymentID)
      } else if (has_n_deployments) {
        sum(n_deployments, na.rm = TRUE)
      } else {
        n()
      },
      wolf_events_per_100_days = 100 * wolf_events / total_effort_days,
      .groups = "drop"
    ) %>%
    filter(is.finite(longitude),
           is.finite(latitude),
           is.finite(total_effort_days),
           total_effort_days > 0) %>%
    arrange(plotID)
}

month_period_start <- function(x) {
  as.POSIXct(
    paste0(format(x, "%Y-%m", tz = "UTC"), "-01 00:00:00"),
    tz = "UTC"
  )
}

next_month_start <- function(x) {
  lt <- as.POSIXlt(month_period_start(x), tz = "UTC")
  lt$mon <- lt$mon + 1L
  as.POSIXct(lt, tz = "UTC")
}

split_one_deployment_by_month <- function(deployment_row) {
  start <- deployment_row$start[[1]]
  end <- deployment_row$end[[1]]
  if (is.na(start) || is.na(end) || end <= start) return(NULL)

  cuts <- start
  boundary <- next_month_start(start)
  while (!is.na(boundary) && boundary < end) {
    cuts <- c(cuts, boundary)
    boundary <- next_month_start(boundary)
  }
  cuts <- c(cuts, end)

  if (length(cuts) < 2L) return(NULL)

  out <- deployment_row[rep(1L, length(cuts) - 1L), , drop = FALSE]
  out$original_deployment_start <- deployment_row$start[[1]]
  out$original_deployment_end <- deployment_row$end[[1]]
  out$start <- cuts[-length(cuts)]
  out$end <- cuts[-1L]
  out$deploymentEffort <- as.numeric(difftime(out$end, out$start, units = "days"))
  out$month <- format(out$start, "%Y-%m", tz = "UTC")
  out$model_row_id <- paste0(out$deploymentID, "__", out$month)
  out
}

split_deployments_by_month <- function(deployments) {
  segments <- dplyr::bind_rows(lapply(seq_len(nrow(deployments)), function(i) {
    split_one_deployment_by_month(deployments[i, , drop = FALSE])
  }))

  segments %>%
    filter(is.finite(deploymentEffort), deploymentEffort > 0) %>%
    arrange(plotID, start, end)
}

count_segment_wolf_events <- function(segments, wolf_events) {
  if (!nrow(segments) || !nrow(wolf_events)) return(rep(0L, nrow(segments)))

  vapply(seq_len(nrow(segments)), function(i) {
    ev <- wolf_events[wolf_events$deploymentID == segments$deploymentID[[i]], , drop = FALSE]
    if (!nrow(ev)) return(0L)

    is_last_segment <- abs(as.numeric(difftime(
      segments$end[[i]],
      segments$original_deployment_end[[i]],
      units = "secs"
    ))) < 1e-6
    in_segment <- ev$event_start >= segments$start[[i]] &
      (ev$event_start < segments$end[[i]] |
         (is_last_segment & ev$event_start <= segments$end[[i]]))
    as.integer(sum(in_segment, na.rm = TRUE))
  }, integer(1))
}

load_2024_survey <- function(settings) {
  dep <- readr::read_csv(path_in("deployments_2024.csv"), show_col_types = FALSE)
  obs <- readr::read_csv(path_in("observations_2024.csv"), show_col_types = FALSE)

  required_dep <- c("deploymentID", "locationID", "latitude", "longitude",
                    "deploymentStart", "deploymentEnd")
  required_obs <- c("deploymentID", "eventID", "scientificName", "eventStart")
  stop_missing_columns(dep, required_dep, "[wolf_2024] deployments")
  stop_missing_columns(obs, required_obs, "[wolf_2024] observations")

  deployments <- dep %>%
    transmute(
      deploymentID = na_if(as.character(deploymentID), ""),
      plotID = na_if(as.character(locationID), ""),
      latitude = as.numeric(latitude),
      longitude = as.numeric(longitude),
      start = parse_time(deploymentStart),
      end = parse_time(deploymentEnd),
      deploymentEffort = as.numeric(difftime(end, start, units = "days")),
      month = format(start, "%Y-%m", tz = "UTC")
    ) %>%
    filter(!is.na(deploymentID),
           !is.na(plotID),
           is.finite(latitude),
           is.finite(longitude),
           !is.na(start),
           !is.na(end),
           is.finite(deploymentEffort),
           deploymentEffort > 0,
           !is.na(month),
           nzchar(month))

  if (!nrow(deployments)) stop("[wolf_2024] no valid dated deployments.")

  wolf_events <- obs %>%
    transmute(
      deploymentID = na_if(as.character(deploymentID), ""),
      eventID = na_if(as.character(eventID), ""),
      scientificName = as.character(scientificName),
      event_start = parse_time(eventStart)
    ) %>%
    filter(scientificName %in% WOLF_NAMES,
           !is.na(deploymentID),
           !is.na(eventID)) %>%
    distinct(deploymentID, eventID, .keep_all = TRUE)

  missing_event_time <- sum(is.na(wolf_events$event_start))
  if (missing_event_time > 0) {
    stop("[wolf_2024] cannot split events by month: ",
         missing_event_time, " wolf event(s) have missing eventStart.")
  }

  model_dat <- split_deployments_by_month(deployments)
  model_dat$wolf_events <- count_segment_wolf_events(model_dat, wolf_events)
  model_dat <- model_dat %>%
    mutate(
      total_effort_days = deploymentEffort,
      wolf_events_per_100_days = 100 * wolf_events / total_effort_days
    ) %>%
    arrange(plotID, start)

  model_dat <- add_month_design(model_dat, settings)

  camera_rate <- camera_summary_from_model(model_dat)

  month_summary <- model_dat %>%
    group_by(month) %>%
    summarise(
      model_rows = n(),
      source_deployments = n_distinct(deploymentID),
      cameras = n_distinct(plotID),
      events = sum(wolf_events),
      effort_days = sum(total_effort_days),
      rate_per_100 = 100 * events / effort_days,
      .groups = "drop"
    )

  readr::write_csv(model_dat,
                   path_out(paste0(SURVEY_PREFIX, "_deployment_month_effort_rates.csv")))
  readr::write_csv(camera_rate,
                   path_out(paste0(SURVEY_PREFIX, "_camera_effort_rates.csv")))
  readr::write_csv(month_summary,
                   path_out(paste0(SURVEY_PREFIX, "_month_observed_summary.csv")))

  cat(sprintf(
    "[wolf_2024] cameras %d | model rows %d | positive rows %d | events %d | effort %.1f camera-days | observed %.2f /100\n",
    nrow(camera_rate),
    nrow(model_dat),
    sum(model_dat$wolf_events > 0),
    sum(model_dat$wolf_events),
    sum(model_dat$total_effort_days),
    100 * sum(model_dat$wolf_events) / sum(model_dat$total_effort_days)
  ))
  cat(sprintf(
    "[wolf_2024] month effect: reference=%s | prediction=%s | months=%s\n",
    unique(model_dat$month_reference),
    unique(model_dat$month_prediction),
    paste(sort(unique(model_dat$month)), collapse = ", ")
  ))

  model_dat
}


## 06. Spatial domain, mesh, and grid ----------------------------------------

camera_to_utm <- function(camera_rate) {
  camera_rate %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(EPSG_UTM)
}

build_spatial <- function(coords, settings) {
  range_arg <- if (!is.null(settings$fix_range_m)) {
    c(settings$fix_range_m, NA)
  } else {
    settings$prior_range_m
  }

  mesh <- INLA::inla.mesh.2d(
    loc = coords,
    cutoff = settings$mesh_cutoff_m,
    max.edge = settings$mesh_max_edge,
    offset = settings$mesh_offset
  )

  spde <- INLA::inla.spde2.pcmatern(
    mesh,
    alpha = 2,
    prior.range = range_arg,
    prior.sigma = settings$prior_sigma
  )

  list(
    mesh = mesh,
    spde = spde,
    s_index = INLA::inla.spde.make.index("spatial", n.spde = spde$n.spde),
    range_fixed = !is.null(settings$fix_range_m)
  )
}

prediction_grid <- function(camera_sf, settings) {
  pts <- st_union(st_geometry(camera_sf))

  # Full prediction map: buffered convex hull around all cameras.
  area <- st_buffer(st_convex_hull(pts), settings$pred_buffer_m)

  grid <- st_make_grid(area, cellsize = settings$cell_size_m, what = "centers")
  pred_sf <- st_sf(grid_id = seq_along(grid), geometry = grid, crs = st_crs(camera_sf))
  pred_sf <- pred_sf[lengths(st_intersects(pred_sf, area)) > 0, ]


  if (!nrow(pred_sf)) stop("Prediction grid is empty.")
  pred_sf
}

prediction_fixed_effects <- function(model_dat, fixed_terms, settings, n_pred) {
  fixed_pred <- as.data.frame(matrix(0, nrow = n_pred, ncol = length(fixed_terms)))
  names(fixed_pred) <- fixed_terms

  if ("intercept" %in% fixed_terms) fixed_pred$intercept <- 1

  month_terms <- intersect(temporal_month_terms(model_dat), fixed_terms)
  if (length(month_terms)) {
    prediction_month <- settings$month_prediction
    prediction_term <- month_term_name(prediction_month)
    if (prediction_term %in% month_terms) fixed_pred[[prediction_term]] <- 1
  }

  time_bin_terms <- intersect(temporal_time_bin_terms(model_dat), fixed_terms)
  if (length(time_bin_terms)) {
    prediction_bin <- if ("time_bin_prediction" %in% names(model_dat)) {
      unique(model_dat$time_bin_prediction)[[1]]
    } else if ("time_bin_reference" %in% names(model_dat)) {
      unique(model_dat$time_bin_reference)[[1]]
    } else {
      NA_character_
    }
    prediction_term <- time_bin_term_name(prediction_bin)
    if (!is.na(prediction_bin) && prediction_term %in% time_bin_terms) {
      fixed_pred[[prediction_term]] <- 1
    }
  }

  fixed_pred
}


## 07. Posterior sample extraction -------------------------------------------

posterior_samples_safe <- function(fit, nsim) {
  nsim <- as.integer(nsim)
  samples <- tryCatch(
    INLA::inla.posterior.sample(nsim, fit, intern = FALSE),
    error = function(e) {
      tryCatch(INLA::inla.posterior.sample(nsim, fit), error = function(e2) NULL)
    }
  )
  if (is.null(samples) || !length(samples)) {
    stop("Posterior sampling failed. Make sure control.compute$config = TRUE.")
  }
  samples
}

predictor_rows_in_sample <- function(sample) {
  rn <- rownames(sample$latent)
  rows <- grep("^APredictor", rn)
  if (!length(rows)) rows <- grep("^Predictor", rn)
  rows
}

extract_eta_matrix <- function(samples, stack_index, expected_n_stack = NULL) {
  first_rows <- predictor_rows_in_sample(samples[[1]])
  if (!length(first_rows)) stop("Could not find APredictor/Predictor rows in posterior samples.")

  if (!is.null(expected_n_stack) && length(first_rows) < max(stack_index)) {
    stop("Posterior sample predictor length is shorter than requested stack index.")
  }

  out <- vapply(samples, function(s) {
    rows <- predictor_rows_in_sample(s)
    as.numeric(s$latent[rows[stack_index], 1])
  }, numeric(length(stack_index)))

  if (is.null(dim(out))) out <- matrix(out, ncol = 1)
  out
}

sample_hyper_vector <- function(samples, pattern, fallback, transform = identity) {
  vapply(samples, function(s) {
    hp <- s$hyperpar
    if (!is.null(hp) && length(hp)) {
      i <- grep(pattern, names(hp), ignore.case = TRUE)
      if (length(i)) {
        val <- as.numeric(hp[i[[1]]])
        val <- transform(val)
        if (is.finite(val)) return(val)
      }
    }
    fallback
  }, numeric(1))
}

extract_pi_samples <- function(fit, samples, family) {
  if (!is_zi(family)) return(rep(0, length(samples)))
  fallback <- hyp_point(fit, PAT_ZPROB)
  if (!is.finite(fallback)) fallback <- plogis(PRIOR_ZI_LOGIT_MEAN)
  p <- sample_hyper_vector(samples, PAT_ZPROB, fallback)

  # If samples accidentally came back on logit scale, transform them. This branch
  # is intentionally conservative and only triggers for values outside [0, 1].
  if (any(is.finite(p) & (p < 0 | p > 1))) {
    p <- plogis(p)
  }
  pmin(pmax(p, 0), 1)
}

extract_size_samples <- function(fit, samples, family) {
  if (!is_nb(family)) return(rep(NA_real_, length(samples)))
  fallback <- nb_size_point(fit)
  if (!is.finite(fallback) || fallback <= 0) fallback <- 1e6
  s <- sample_hyper_vector(samples, PAT_NB_SIZE, fallback)

  # If samples accidentally came back on log-size scale, exponentiate values that
  # look like internal-scale log(size). This is a fallback only.
  if (median(s, na.rm = TRUE) < 0 || any(s <= 0, na.rm = TRUE)) {
    s <- exp(s)
  }
  fam_nb_size(s)
}

build_posterior_draws <- function(fit, samples, index, effort, family,
                                  expected_n_stack = NULL) {
  eta <- extract_eta_matrix(samples, index, expected_n_stack)
  pi <- extract_pi_samples(fit, samples, family)
  size <- extract_size_samples(fit, samples, family)

  nsim <- ncol(eta)
  mu <- eta
  for (j in seq_len(nsim)) {
    mu[, j] <- effort * exp(eta[, j])
  }

  fitted <- mu
  fit_var <- mu
  for (j in seq_len(nsim)) {
    fitted[, j] <- fam_mean(mu[, j], pi[j], family, size[j])
    fit_var[, j] <- fam_var(mu[, j], pi[j], family, size[j])
  }

  list(eta = eta, mu = mu, fitted = fitted, fit_var = fit_var,
       pi = pi, size = size)
}

simulate_from_draws <- function(draws, family) {
  nsim <- ncol(draws$mu)
  sim <- draws$mu
  for (j in seq_len(nsim)) {
    sim[, j] <- fam_sim(draws$mu[, j], draws$pi[j], family, draws$size[j])
  }
  storage.mode(sim) <- "integer"
  sim
}


## 08. PPC, PIT, residual diagnostics ----------------------------------------

aggregate_matrix_by_group <- function(mat, group) {
  group <- as.factor(group)
  apply(mat, 2, function(x) as.numeric(rowsum(x, group)[, 1]))
}

summarise_ppc_simulations <- function(sim, model_dat, method) {
  yobs <- model_dat$y
  camera_group <- as.factor(model_dat$plotID)
  yobs_camera <- as.numeric(rowsum(yobs, camera_group)[, 1])

  if (is.null(dim(sim))) sim <- matrix(sim, ncol = 1)
  nsim <- ncol(sim)

  sim_camera <- aggregate_matrix_by_group(sim, camera_group)
  if (is.null(dim(sim_camera))) sim_camera <- matrix(sim_camera, ncol = 1)

  row_total <- colSums(sim)
  row_zero_fraction <- colMeans(sim == 0)
  row_max <- apply(sim, 2, max)

  camera_total <- colSums(sim_camera)
  camera_zero_fraction <- colMeans(sim_camera == 0)
  camera_max <- apply(sim_camera, 2, max)

  ppc_stat <- function(stat, observed, values, level) {
    data.frame(
      level = level,
      stat = stat,
      observed = observed,
      sim_median = median(values),
      sim_q025 = unname(quantile(values, 0.025)),
      sim_q975 = unname(quantile(values, 0.975)),
      pass = observed >= unname(quantile(values, 0.025)) &
        observed <= unname(quantile(values, 0.975)),
      method = method
    )
  }

  summary <- bind_rows(
    ppc_stat("total_events", sum(yobs), row_total, "model_row"),
    ppc_stat("zero_fraction", mean(yobs == 0), row_zero_fraction, "model_row"),
    ppc_stat("max_count", max(yobs), row_max, "model_row"),
    ppc_stat("total_events", sum(yobs_camera), camera_total, "camera"),
    ppc_stat("zero_fraction", mean(yobs_camera == 0), camera_zero_fraction, "camera"),
    ppc_stat("max_count", max(yobs_camera), camera_max, "camera")
  )

  row_pit <- (rowSums(sim < yobs) + runif(length(yobs)) * rowSums(sim == yobs)) / nsim
  cam_pit <- (rowSums(sim_camera < yobs_camera) +
                runif(length(yobs_camera)) * rowSums(sim_camera == yobs_camera)) / nsim

  list(summary = summary, row_pit = row_pit, camera_pit = cam_pit,
       sim = sim, sim_camera = sim_camera, nsim = nsim)
}

ks_uniform_p_value <- function(pit) {
  pit <- pit[is.finite(pit)]
  if (length(pit) < 5) return(NA_real_)
  suppressWarnings(tryCatch(
    ks.test(pit, "punif")$p.value,
    error = function(e) NA_real_
  ))
}

moran_perm <- function(coords, x, nperm = MORAN_NPERM, two_sided = TRUE) {
  ok <- is.finite(x)
  coords <- coords[ok, , drop = FALSE]
  x <- x[ok]
  n <- length(x)

  if (n < 5) {
    return(list(I = NA_real_, expected = NA_real_, p_value = NA_real_,
                alternative = if (two_sided) "two_sided" else "greater"))
  }

  D <- as.matrix(dist(coords))
  W <- 1 / D
  diag(W) <- 0
  W[!is.finite(W)] <- 0

  rs <- rowSums(W)
  rs[rs == 0] <- 1
  W <- W / rs

  z <- x - mean(x)
  z2 <- sum(z^2)
  S0 <- sum(W)
  expected <- -1 / (n - 1)
  if (!is.finite(z2) || z2 <= 0 || !is.finite(S0) || S0 <= 0) {
    return(list(I = NA_real_, expected = expected, p_value = NA_real_,
                alternative = if (two_sided) "two_sided" else "greater"))
  }

  I_stat <- function(v) {
    (n / S0) * sum(W * outer(v, v)) / sum(v^2)
  }

  I0 <- I_stat(z)
  perm <- replicate(nperm, I_stat(sample(z)))

  p_value <- if (two_sided) {
    obs_dev <- abs(I0 - expected)
    perm_dev <- abs(perm - expected)
    (1 + sum(perm_dev >= obs_dev)) / (nperm + 1)
  } else {
    (1 + sum(perm >= I0)) / (nperm + 1)
  }

  list(I = I0, expected = expected, p_value = p_value,
       alternative = if (two_sided) "two_sided" else "greater")
}

resid_variogram <- function(coords, residual, nbins = 12L) {
  ok <- is.finite(residual) & is.finite(coords[, 1]) & is.finite(coords[, 2])
  coords <- coords[ok, , drop = FALSE]
  residual <- residual[ok]
  if (length(residual) < 3) {
    return(data.frame(dist = numeric(), gamma = numeric(), n = integer()))
  }

  D <- as.matrix(dist(coords))
  G <- as.matrix(dist(residual))^2 / 2

  d <- D[lower.tri(D)]
  g <- G[lower.tri(G)]
  if (!length(d)) return(data.frame(dist = numeric(), gamma = numeric(), n = integer()))

  br <- seq(0, quantile(d, 0.9, na.rm = TRUE), length.out = nbins + 1)
  bin <- cut(d, br, include.lowest = TRUE)

  data.frame(
    dist = as.numeric(tapply(d, bin, mean)),
    gamma = as.numeric(tapply(g, bin, mean)),
    n = as.integer(table(bin))
  ) %>% filter(is.finite(dist), is.finite(gamma), n > 0)
}

camera_residual_diagnostics <- function(model_dat) {
  model_dat %>%
    group_by(plotID) %>%
    summarise(
      longitude = mean(longitude, na.rm = TRUE),
      latitude = mean(latitude, na.rm = TRUE),
      y = sum(y, na.rm = TRUE),
      fitted_count = sum(fitted_count, na.rm = TRUE),
      fit_var = sum(fit_var, na.rm = TRUE),
      total_effort_days = sum(total_effort_days, na.rm = TRUE),
      wolf_events_per_100_days = 100 * y / total_effort_days,
      fitted_rate_per_100 = 100 * fitted_count / total_effort_days,
      n_model_rows = n(),
      month_first = min(month, na.rm = TRUE),
      month_last = max(month, na.rm = TRUE),
      pearson = (y - fitted_count) / sqrt(pmax(fit_var, 1e-9)),
      .groups = "drop"
    ) %>%
    arrange(plotID)
}

ppc_pass_lookup <- function(ppc_summary, level, stat) {
  x <- ppc_summary$pass[ppc_summary$level == level & ppc_summary$stat == stat]
  if (length(x)) isTRUE(x[[1]]) else FALSE
}

compute_diagnostics <- function(fit, samples, model_dat, obs_index, camera_sf,
                                family, write_files = TRUE) {
  obs_draws <- build_posterior_draws(
    fit = fit,
    samples = samples,
    index = obs_index,
    effort = model_dat$total_effort_days,
    family = family,
    expected_n_stack = length(fit$summary.linear.predictor$mean)
  )

  sim <- simulate_from_draws(obs_draws, family)
  ppc <- summarise_ppc_simulations(sim, model_dat, "joint_posterior")

  model_dat$fitted_count <- rowMeans(obs_draws$fitted)
  model_dat$fit_var <- rowMeans(obs_draws$fit_var) +
    apply(obs_draws$fitted, 1, var, na.rm = TRUE)
  model_dat$eta_mean <- rowMeans(obs_draws$eta)
  model_dat$eta_sd <- apply(obs_draws$eta, 1, sd, na.rm = TRUE)
  model_dat$pearson <- (model_dat$y - model_dat$fitted_count) /
    sqrt(pmax(model_dat$fit_var, 1e-9))

  if (!is.null(fit$summary.fitted.values) &&
      nrow(fit$summary.fitted.values) >= max(obs_index)) {
    model_dat$inla_fitted_mean <- fit$summary.fitted.values$mean[obs_index]
  } else {
    model_dat$inla_fitted_mean <- NA_real_
  }

  camera_diag <- camera_residual_diagnostics(model_dat)
  coords_camera <- st_coordinates(camera_sf)
  moran <- moran_perm(coords_camera, camera_diag$pearson,
                      nperm = MORAN_NPERM, two_sided = TRUE)

  fitted_check <- data.frame(
    check = c("observed_total", "posterior_sample_fitted_total", "inla_fitted_total",
              "cor_posterior_fitted_vs_inla_fitted"),
    value = c(
      sum(model_dat$y),
      sum(model_dat$fitted_count),
      sum(model_dat$inla_fitted_mean, na.rm = TRUE),
      safe_cor(model_dat$fitted_count, model_dat$inla_fitted_mean)
    )
  )

  row_disp <- mean(model_dat$pearson^2, na.rm = TRUE)
  cam_disp <- mean(camera_diag$pearson^2, na.rm = TRUE)

  ppc_total_pass <- ppc_pass_lookup(ppc$summary, "camera", "total_events")
  ppc_zero_pass <- ppc_pass_lookup(ppc$summary, "camera", "zero_fraction")
  ppc_max_pass <- ppc_pass_lookup(ppc$summary, "camera", "max_count")
  moran_pass <- is.finite(moran$p_value) && moran$p_value >= MORAN_ALPHA

  diagnostics <- list(
    model_dat = model_dat,
    camera_diag = camera_diag,
    pi_hat = mean(obs_draws$pi, na.rm = TRUE),
    size_hat = mean(obs_draws$size, na.rm = TRUE),
    pearson_disp = row_disp,
    pearson_disp_camera = cam_disp,
    moran_I = moran$I,
    moran_expected = moran$expected,
    moran_p = moran$p_value,
    moran_alternative = moran$alternative,
    ppc_pit_ks_row = ks_uniform_p_value(ppc$row_pit),
    ppc_pit_ks_camera = ks_uniform_p_value(ppc$camera_pit),
    pit_mean_row = mean(ppc$row_pit, na.rm = TRUE),
    pit_mean_camera = mean(ppc$camera_pit, na.rm = TRUE),
    ppc_method = "joint_posterior",
    ppc_nsim = ppc$nsim,
    ppc = ppc$summary,
    row_pit = ppc$row_pit,
    camera_pit = ppc$camera_pit,
    ppc_total_pass = ppc_total_pass,
    ppc_zero_pass = ppc_zero_pass,
    ppc_max_pass = ppc_max_pass,
    moran_pass = moran_pass,
    diagnostics_ok = ppc_total_pass && ppc_zero_pass && ppc_max_pass && moran_pass,
    fitted_check = fitted_check,
    obs_draws = obs_draws
  )

  if (write_files) {
    readr::write_csv(model_dat,
                     path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                                     "_model_row_diagnostics.csv")))
    readr::write_csv(camera_diag,
                     path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                                     "_camera_residual_diagnostics.csv")))
    readr::write_csv(ppc$summary,
                     path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                                     "_posterior_predictive_check.csv")))
    readr::write_csv(fitted_check,
                     path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                                     "_fitted_scale_sanity_check.csv")))
  }

  diagnostics
}


## 09. Temporal autocorrelation diagnostics ----------------------------------

# These diagnostics test whether residuals remain temporally structured after
# the model has accounted for effort, spatial structure, and deployment-start
# month. They are diagnostic checks only; the fitted model remains the same.
#
# Outputs:
#   * residuals versus deployment start date
#   * pooled within-camera residual lag correlations, including lag 1
#   * ACF of mean residuals by deployment date
#   * a text note explaining limitations when there are too few repeated rows

temporal_autocorrelation_diagnostics <- function(model_dat,
                                                 prefix = SURVEY_PREFIX,
                                                 final_model = FINAL_MODEL_NAME,
                                                 max_lag = 5L,
                                                 min_pairs = 5L) {
  required <- c("plotID", "start", "month", "pearson")
  missing <- setdiff(required, names(model_dat))
  note_file <- path_out(paste0(prefix, "_", final_model,
                               "_temporal_autocorrelation_NOTE.txt"))

  if (length(missing)) {
    writeLines(
      c(
        "Temporal autocorrelation diagnostics skipped.",
        sprintf("Missing required column(s): %s", paste(missing, collapse = ", "))
      ),
      note_file
    )
    return(NULL)
  }

  dat <- model_dat %>%
    mutate(
      start_date = as_utc_date(start),
      start_doy = as.integer(format(start_date, "%j"))
    ) %>%
    filter(!is.na(plotID), !is.na(start_date), is.finite(pearson)) %>%
    arrange(plotID, start)

  if (!nrow(dat)) {
    writeLines(
      c("Temporal autocorrelation diagnostics skipped.",
        "No rows with finite Pearson residuals and valid deployment dates."),
      note_file
    )
    return(NULL)
  }

  readr::write_csv(
    dat,
    path_out(paste0(prefix, "_", final_model,
                    "_temporal_residual_input.csv"))
  )

  # --------------------------------------------------------------------------
  # 1. Residuals through calendar time.
  # --------------------------------------------------------------------------

  p_time <- ggplot(dat, aes(start_date, pearson)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    geom_point(alpha = 0.75) +
    geom_smooth(se = TRUE, method = "loess", formula = y ~ x) +
    labs(
      title = paste0("Temporal residual trend: ", prefix),
      subtitle = paste0(final_model, " | residuals by deployment start date"),
      x = "deployment start date",
      y = "row Pearson residual"
    ) +
    theme_minimal(base_size = 12)

  ggsave(
    path_out(paste0(prefix, "_", final_model,
                    "_diag_temporal_residuals_vs_date.png")),
    p_time,
    width = 7,
    height = 5,
    dpi = 220
  )

  # --------------------------------------------------------------------------
  # 2. Pooled within-camera lag residual correlations.
  #    Lag is deployment sequence at the same camera, not a fixed time interval.
  # --------------------------------------------------------------------------

  make_lag_pairs <- function(lag_k) {
    pieces <- lapply(split(dat, dat$plotID), function(d) {
      d <- d[order(d$start), , drop = FALSE]
      n <- nrow(d)
      if (n <= lag_k) return(NULL)

      previous <- seq_len(n - lag_k)
      current <- previous + lag_k

      data.frame(
        plotID = d$plotID[current],
        lag = lag_k,
        start_previous = d$start[previous],
        start_current = d$start[current],
        residual_previous = d$pearson[previous],
        residual_current = d$pearson[current],
        days_between = as.numeric(difftime(d$start[current],
                                           d$start[previous],
                                           units = "days")),
        stringsAsFactors = FALSE
      )
    })
    dplyr::bind_rows(pieces)
  }

  lag_pairs <- dplyr::bind_rows(lapply(seq_len(max_lag), make_lag_pairs))

  if (nrow(lag_pairs)) {
    readr::write_csv(
      lag_pairs,
      path_out(paste0(prefix, "_", final_model,
                      "_temporal_within_camera_lag_pairs.csv"))
    )
  }

  safe_cor_p <- function(a, b) {
    safe_cor_p_value(a, b, min_pairs)
  }

  lag_summary <- if (nrow(lag_pairs)) {
    lag_pairs %>%
      group_by(lag) %>%
      summarise(
        n_pairs = n(),
        n_cameras = n_distinct(plotID),
        correlation = safe_cor(residual_previous, residual_current),
        p_value = safe_cor_p(residual_previous, residual_current),
        mean_days_between = mean(days_between, na.rm = TRUE),
        median_days_between = median(days_between, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    data.frame(
      lag = integer(),
      n_pairs = integer(),
      n_cameras = integer(),
      correlation = numeric(),
      p_value = numeric(),
      mean_days_between = numeric(),
      median_days_between = numeric()
    )
  }

  readr::write_csv(
    lag_summary,
    path_out(paste0(prefix, "_", final_model,
                    "_temporal_within_camera_lag_correlation.csv"))
  )

  lag1 <- lag_pairs %>% filter(lag == 1L)
  if (nrow(lag1) >= min_pairs &&
      is.finite(safe_cor(lag1$residual_previous, lag1$residual_current))) {
    lag1_row <- lag_summary %>% filter(lag == 1L)

    p_lag1 <- ggplot(lag1, aes(residual_previous, residual_current)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey70") +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey70") +
      geom_point(alpha = 0.75) +
      geom_smooth(se = TRUE, method = "lm", formula = y ~ x) +
      labs(
        title = paste0("Within-camera lag-1 residual correlation: ", prefix),
        subtitle = sprintf(
          "%s | r = %.3f, p = %.3g, n = %d pairs, median gap = %.1f days",
          final_model,
          lag1_row$correlation[[1]],
          lag1_row$p_value[[1]],
          lag1_row$n_pairs[[1]],
          lag1_row$median_days_between[[1]]
        ),
        x = "previous residual at same camera",
        y = "current residual"
      ) +
      theme_minimal(base_size = 12)

    ggsave(
      path_out(paste0(prefix, "_", final_model,
                      "_diag_temporal_lag1_residual_correlation.png")),
      p_lag1,
      width = 6.5,
      height = 5.5,
      dpi = 220
    )
  } else {
    writeLines(
      c(
        "Within-camera lag-1 residual plot skipped.",
        sprintf("Usable lag-1 pairs: %d", nrow(lag1)),
        sprintf("Minimum required pairs: %d", min_pairs),
        "This usually means there are too few repeated camera-month rows per camera."
      ),
      path_out(paste0(prefix, "_", final_model,
                      "_temporal_lag1_residual_correlation_NOTE.txt"))
    )
  }

  if (nrow(lag_summary)) {
    p_lag_summary <- ggplot(lag_summary, aes(lag, correlation)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
      geom_point(aes(size = n_pairs)) +
      geom_line() +
      scale_x_continuous(breaks = seq_len(max_lag)) +
      labs(
        title = paste0("Within-camera residual autocorrelation by lag: ", prefix),
        subtitle = paste0(final_model, " | lag is deployment sequence within camera"),
        x = "within-camera deployment lag",
        y = "pooled residual correlation",
        size = "pairs"
      ) +
      theme_minimal(base_size = 12)

    ggsave(
      path_out(paste0(prefix, "_", final_model,
                      "_diag_temporal_lag_correlation_summary.png")),
      p_lag_summary,
      width = 6.5,
      height = 5,
      dpi = 220
    )
  }

  # --------------------------------------------------------------------------
  # 3. ACF of mean residuals by deployment start date.
  #    This is a broad date-ordered diagnostic. It mixes cameras, so it should
  #    be interpreted as evidence of survey-wide temporal residual structure,
  #    not camera-specific serial dependence.
  # --------------------------------------------------------------------------

  daily_resid <- dat %>%
    group_by(start_date) %>%
    summarise(
      mean_pearson = mean(pearson, na.rm = TRUE),
      median_pearson = median(pearson, na.rm = TRUE),
      n_rows = n(),
      .groups = "drop"
    ) %>%
    arrange(start_date)

  readr::write_csv(
    daily_resid,
    path_out(paste0(prefix, "_", final_model,
                    "_temporal_mean_residual_by_date.csv"))
  )

  acf_summary <- data.frame(lag = integer(), acf = numeric())
  if (nrow(daily_resid) >= 6 && sd(daily_resid$mean_pearson, na.rm = TRUE) > 0) {
    acf_obj <- stats::acf(daily_resid$mean_pearson,
                          plot = FALSE,
                          na.action = stats::na.pass)
    acf_summary <- data.frame(
      lag = as.integer(round(as.numeric(acf_obj$lag))),
      acf = as.numeric(acf_obj$acf)
    )
    readr::write_csv(
      acf_summary,
      path_out(paste0(prefix, "_", final_model,
                      "_temporal_acf_mean_residual_by_date.csv"))
    )

    png(
      filename = path_out(paste0(prefix, "_", final_model,
                                 "_diag_temporal_acf_mean_residual_by_date.png")),
      width = 1400,
      height = 1000,
      res = 180
    )
    stats::acf(
      daily_resid$mean_pearson,
      main = paste0("ACF of mean residuals by deployment date: ", prefix),
      xlab = "lag in ordered deployment dates"
    )
    dev.off()
  } else {
    writeLines(
      c(
        "Date-ordered residual ACF skipped.",
        sprintf("Unique deployment dates: %d", nrow(daily_resid)),
        "At least 6 unique dates and non-constant mean residuals are required."
      ),
      path_out(paste0(prefix, "_", final_model,
                      "_temporal_acf_mean_residual_by_date_NOTE.txt"))
    )
  }

  # --------------------------------------------------------------------------
  # 4. Compact text report.
  # --------------------------------------------------------------------------

  lag1_summary <- lag_summary %>% filter(lag == 1L)
  lag1_line <- if (nrow(lag1_summary)) {
    sprintf(
      "Within-camera lag-1 residual correlation: r = %.3f, p = %.4g, n pairs = %d, cameras = %d, median gap = %.1f days.",
      lag1_summary$correlation[[1]],
      lag1_summary$p_value[[1]],
      lag1_summary$n_pairs[[1]],
      lag1_summary$n_cameras[[1]],
      lag1_summary$median_days_between[[1]]
    )
  } else {
    "Within-camera lag-1 residual correlation was not evaluable."
  }

  acf1 <- if (nrow(acf_summary) >= 2 && any(acf_summary$lag == 1L)) {
    acf_summary$acf[acf_summary$lag == 1L][[1]]
  } else {
    NA_real_
  }
  acf_line <- if (is.finite(acf1)) {
    sprintf("Date-ordered ACF of mean residuals: lag-1 ACF = %.3f.", acf1)
  } else {
    "Date-ordered ACF of mean residuals was not evaluable."
  }

  report <- c(
    sprintf("Temporal autocorrelation diagnostics: %s", prefix),
    sprintf("Model: %s", final_model),
    sprintf("Rows used: %d", nrow(dat)),
    sprintf("Cameras used: %d", dplyr::n_distinct(dat$plotID)),
    sprintf("Unique deployment start dates: %d", dplyr::n_distinct(dat$start_date)),
    "",
    "Diagnostics are based on Pearson residuals from the fitted model.",
    "The month fixed effect already adjusts the mean by calendar camera-month; these diagnostics check residual temporal structure remaining after that adjustment.",
    "",
    lag1_line,
    acf_line
  )

  writeLines(
    report,
    path_out(paste0(prefix, "_", final_model,
                    "_TEMPORAL_AUTOCORRELATION_REPORT.txt"))
  )

  list(
    data = dat,
    lag_pairs = lag_pairs,
    lag_summary = lag_summary,
    daily_resid = daily_resid,
    acf_summary = acf_summary,
    report = report
  )
}


## 09. Diagnostic plots -------------------------------------------------------

write_diagnostic_plots <- function(diag, camera_sf) {
  model_dat <- diag$model_dat
  camera_diag <- diag$camera_diag
  coords_camera <- st_coordinates(camera_sf)

  obs_fit_plot <- ggplot(camera_diag, aes(fitted_count, y)) +
    geom_point(alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    scale_x_continuous(trans = "sqrt") +
    scale_y_continuous(trans = "sqrt") +
    labs(
      title = "Observed vs fitted camera counts: wolf_2024",
      subtitle = FINAL_MODEL_NAME,
      x = "posterior mean fitted wolf events",
      y = "observed wolf events"
    ) +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                         "_diag_obs_vs_fitted.png")),
         obs_fit_plot, width = 6, height = 5.5, dpi = 220)

  residual_sf <- camera_sf
  residual_sf$pearson <- camera_diag$pearson
  residual_plot <- ggplot(residual_sf) +
    geom_sf(aes(colour = pearson, size = abs(pearson)), alpha = 0.9) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_size_continuous(range = c(2, 8)) +
    coord_sf(datum = NA) +
    labs(title = "Spatial Pearson residuals: wolf_2024",
         subtitle = FINAL_MODEL_NAME,
         colour = "Pearson", size = "|Pearson|") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())
  ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                         "_diag_spatial_residuals.png")),
         residual_plot, width = 7, height = 7, dpi = 220)

  pit_row_plot <- ggplot(data.frame(pit = diag$row_pit), aes(pit)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20,
                   fill = "grey55", colour = "white", boundary = 0) +
    geom_hline(yintercept = 1, linetype = 2, colour = "red") +
    labs(
      title = "Posterior predictive PIT, model rows: wolf_2024",
      subtitle = sprintf("mean %.3f | KS p %.3g",
                         diag$pit_mean_row, diag$ppc_pit_ks_row),
      x = "PIT", y = "density"
    ) +
    xlim(0, 1) +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                         "_diag_pit_hist_model_rows.png")),
         pit_row_plot, width = 6, height = 5, dpi = 220)

  pit_camera_plot <- ggplot(data.frame(pit = diag$camera_pit), aes(pit)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20,
                   fill = "grey55", colour = "white", boundary = 0) +
    geom_hline(yintercept = 1, linetype = 2, colour = "red") +
    labs(
      title = "Posterior predictive PIT, camera aggregates: wolf_2024",
      subtitle = sprintf("mean %.3f | KS p %.3g",
                         diag$pit_mean_camera, diag$ppc_pit_ks_camera),
      x = "PIT", y = "density"
    ) +
    xlim(0, 1) +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                         "_diag_pit_hist_camera.png")),
         pit_camera_plot, width = 6, height = 5, dpi = 220)

  variogram <- resid_variogram(coords_camera, camera_diag$pearson)
  readr::write_csv(variogram,
                   path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                                   "_diag_resid_variogram.csv")))
  if (nrow(variogram)) {
    variogram_plot <- ggplot(variogram, aes(dist, gamma)) +
      geom_point(aes(size = n)) +
      geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
      labs(
        title = "Residual semivariogram: wolf_2024",
        subtitle = FINAL_MODEL_NAME,
        x = "distance (m)",
        y = "semivariance"
      ) +
      theme_minimal(base_size = 12)
    ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                           "_diag_resid_variogram.png")),
           variogram_plot, width = 6.5, height = 5, dpi = 220)
  }

  residual_covariate <- model_dat %>%
    mutate(
      fitted_rate_per_100 = 100 * fitted_count / total_effort_days,
      northing_proxy = latitude
    )

  effort_plot <- ggplot(residual_covariate,
                        aes(total_effort_days, pearson)) +
    geom_point(alpha = 0.75) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
    labs(title = "Residuals vs effort: wolf_2024",
         x = "deployment effort days", y = "row Pearson residual") +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                         "_diag_resid_vs_effort.png")),
         effort_plot, width = 6, height = 5, dpi = 220)

  month_plot <- ggplot(residual_covariate,
                       aes(month, pearson)) +
    geom_boxplot(outlier.alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = "Residuals by calendar month: wolf_2024",
         x = "month", y = "row Pearson residual") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                         "_diag_resid_by_month.png")),
         month_plot, width = 7, height = 5, dpi = 220)

  lat_plot <- ggplot(residual_covariate,
                     aes(latitude, pearson)) +
    geom_point(alpha = 0.75) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
    labs(title = "Residuals vs latitude: wolf_2024",
         x = "latitude", y = "row Pearson residual") +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_", FINAL_MODEL_NAME,
                         "_diag_resid_vs_latitude.png")),
         lat_plot, width = 6, height = 5, dpi = 220)

  invisible(TRUE)
}


## 10. Prior-posterior plots and summaries -----------------------------------

pc_range_density <- function(x, range0, prob_below_range0) {
  lambda <- -range0 * log(prob_below_range0)
  ifelse(x > 0, lambda * x^(-2) * exp(-lambda / x), 0)
}

pc_range_quantile <- function(p, range0, prob_below_range0) {
  lambda <- -range0 * log(prob_below_range0)
  lambda / (-log(p))
}

pc_sigma_density <- function(x, sigma0, prob_above_sigma0) {
  lambda <- -log(prob_above_sigma0) / sigma0
  ifelse(x >= 0, lambda * exp(-lambda * x), 0)
}

pc_sigma_quantile <- function(p, sigma0, prob_above_sigma0) {
  lambda <- -log(prob_above_sigma0) / sigma0
  -log(1 - p) / lambda
}

zip_prob_prior_density <- function(p) {
  sd_logit <- 1 / sqrt(PRIOR_ZI_LOGIT_PREC)
  dnorm(qlogis(p), mean = PRIOR_ZI_LOGIT_MEAN, sd = sd_logit) / (p * (1 - p))
}

nb_size_prior_density <- function(x) {
  sd_logsize <- 1 / sqrt(PRIOR_NB_LOGSIZE_PREC)
  dnorm(log(x), mean = PRIOR_NB_LOGSIZE_MEAN, sd = sd_logsize) / x
}

nb_size_prior_quantile <- function(p) {
  sd_logsize <- 1 / sqrt(PRIOR_NB_LOGSIZE_PREC)
  exp(qnorm(p, mean = PRIOR_NB_LOGSIZE_MEAN, sd = sd_logsize))
}

plot_prior_posterior_density <- function(parameter, prior_df,
                                         posterior_marginal, file_suffix,
                                         x_label, log_x = FALSE,
                                         reference_x = NULL,
                                         reference_label = NULL) {
  if (is.null(posterior_marginal) || nrow(posterior_marginal) < 2) {
    writeLines(
      paste("Posterior marginal not available for", parameter),
      path_out(paste0(SURVEY_PREFIX, "_prior_posterior_", file_suffix, "_NOTE.txt"))
    )
    return(invisible(NULL))
  }

  posterior_df <- data.frame(
    value = posterior_marginal[, 1],
    density = posterior_marginal[, 2],
    source = "posterior"
  ) %>% filter(is.finite(value), is.finite(density), density >= 0)

  prior_df <- prior_df %>%
    mutate(source = "prior") %>%
    filter(is.finite(value), is.finite(density), density >= 0)

  if (log_x) {
    posterior_df <- posterior_df %>% filter(value > 0)
    prior_df <- prior_df %>% filter(value > 0)
  }

  plot_df <- bind_rows(prior_df, posterior_df)
  if (!nrow(plot_df)) return(invisible(NULL))

  p <- ggplot(plot_df, aes(value, density, colour = source, linetype = source)) +
    geom_line(linewidth = 0.85) +
    labs(title = paste0("Prior vs posterior: ", parameter),
         subtitle = SURVEY_PREFIX,
         x = x_label, y = "density", colour = NULL, linetype = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")

  if (log_x) p <- p + scale_x_log10(labels = label_number())

  if (!is.null(reference_x) && is.finite(reference_x)) {
    p <- p + geom_vline(xintercept = reference_x, linetype = 2, colour = "grey35")
    if (!is.null(reference_label)) p <- p + labs(caption = reference_label)
  }

  ggsave(path_out(paste0(SURVEY_PREFIX, "_prior_posterior_", file_suffix, ".png")),
         p, width = 6.8, height = 4.8, dpi = 250)
  readr::write_csv(plot_df,
                   path_out(paste0(SURVEY_PREFIX, "_prior_posterior_", file_suffix, ".csv")))
  invisible(plot_df)
}

write_prior_posterior_plots <- function(fit, settings, family) {
  intercept_post <- fit$marginals.fixed[["intercept"]]
  intercept_sd_prior <- 1 / sqrt(PRIOR_INTERCEPT_PREC)
  intercept_grid <- seq(qnorm(0.001, PRIOR_INTERCEPT_MEAN, intercept_sd_prior),
                        qnorm(0.999, PRIOR_INTERCEPT_MEAN, intercept_sd_prior),
                        length.out = 700)
  plot_prior_posterior_density(
    parameter = "intercept",
    prior_df = data.frame(value = intercept_grid,
                          density = dnorm(intercept_grid, PRIOR_INTERCEPT_MEAN, intercept_sd_prior)),
    posterior_marginal = intercept_post,
    file_suffix = "intercept",
    x_label = "log-rate intercept",
    reference_x = PRIOR_INTERCEPT_MEAN,
    reference_label = sprintf("Prior center: crude observed log daily rate = %.2f",
                              PRIOR_INTERCEPT_MEAN)
  )

  month_terms <- grep("^month_[0-9]{4}_[0-9]{2}$", names(fit$marginals.fixed), value = TRUE)
  if (length(month_terms)) {
    month_sd_prior <- 1 / sqrt(PRIOR_MONTH_LOG_RATE_RATIO_PREC)
    month_grid <- seq(qnorm(0.001, 0, month_sd_prior),
                      qnorm(0.999, 0, month_sd_prior),
                      length.out = 700)
    month_prior <- data.frame(value = month_grid,
                              density = dnorm(month_grid, 0, month_sd_prior))
    for (term in month_terms) {
      plot_prior_posterior_density(
        parameter = paste0("month log-rate ratio: ", month_from_term(term),
                           " vs ", settings$month_reference),
        prior_df = month_prior,
        posterior_marginal = fit$marginals.fixed[[term]],
        file_suffix = paste0("fixed_", term),
        x_label = "log-rate ratio"
      )
    }
  }

  if (!is.null(settings$fix_range_m)) {
    writeLines(c("Spatial range was fixed.",
                 sprintf("fix_range_m = %s", settings$fix_range_m)),
               path_out(paste0(SURVEY_PREFIX, "_prior_posterior_range_NOTE.txt")))
  } else {
    range0 <- settings$prior_range_m[1]
    prob_range <- settings$prior_range_m[2]
    range_grid <- exp(seq(log(pc_range_quantile(0.001, range0, prob_range)),
                          log(pc_range_quantile(0.995, range0, prob_range)),
                          length.out = 700))
    range_prior <- data.frame(value = range_grid,
                              density = pc_range_density(range_grid, range0, prob_range))
    plot_prior_posterior_density(
      parameter = "spatial range",
      prior_df = range_prior,
      posterior_marginal = hyp_marg(fit, PAT_RANGE),
      file_suffix = "range",
      x_label = "range (m)",
      log_x = TRUE,
      reference_x = range0,
      reference_label = sprintf("Prior statement: P(range < %.0f m) = %.2f",
                                range0, prob_range)
    )
  }

  sigma0 <- settings$prior_sigma[1]
  prob_sigma <- settings$prior_sigma[2]
  sigma_grid <- seq(0, pc_sigma_quantile(0.995, sigma0, prob_sigma), length.out = 700)
  sigma_prior <- data.frame(value = sigma_grid,
                            density = pc_sigma_density(sigma_grid, sigma0, prob_sigma))
  plot_prior_posterior_density(
    parameter = "spatial marginal SD",
    prior_df = sigma_prior,
    posterior_marginal = hyp_marg(fit, PAT_SIGMA),
    file_suffix = "spatial_sd",
    x_label = "spatial marginal SD",
    reference_x = sigma0,
    reference_label = sprintf("Prior statement: P(SD > %.2f) = %.2f",
                              sigma0, prob_sigma)
  )

  if (is_zi(family)) {
    prob_grid <- seq(0.001, 0.999, length.out = 700)
    plot_prior_posterior_density(
      parameter = "zero-inflation probability",
      prior_df = data.frame(value = prob_grid,
                            density = zip_prob_prior_density(prob_grid)),
      posterior_marginal = hyp_marg(fit, PAT_ZPROB),
      file_suffix = "zero_inflation_probability",
      x_label = "zero-inflation probability",
      reference_x = plogis(PRIOR_ZI_LOGIT_MEAN),
      reference_label = sprintf("Prior center on probability scale: %.2f",
                                plogis(PRIOR_ZI_LOGIT_MEAN))
    )
  }

  if (is_nb(family)) {
    size_grid <- exp(seq(log(nb_size_prior_quantile(0.001)),
                         log(nb_size_prior_quantile(0.995)),
                         length.out = 700))
    prior_median_size <- exp(PRIOR_NB_LOGSIZE_MEAN)
    plot_prior_posterior_density(
      parameter = "negative-binomial size",
      prior_df = data.frame(value = size_grid,
                            density = nb_size_prior_density(size_grid)),
      posterior_marginal = hyp_marg(fit, PAT_NB_SIZE),
      file_suffix = "nb_size",
      x_label = "negative-binomial size",
      log_x = TRUE,
      reference_x = prior_median_size,
      reference_label = sprintf(
        "Prior on log(size): Normal(%.2f, SD %.2f); median size %.2f",
        PRIOR_NB_LOGSIZE_MEAN,
        1 / sqrt(PRIOR_NB_LOGSIZE_PREC),
        prior_median_size
      )
    )
  }

  invisible(TRUE)
}



# Quantitative prior-influence screen ---------------------------------------
# These diagnostics do not prove causality by themselves. They summarize how
# much the posterior narrowed relative to the prior and whether the posterior
# lies in a low-prior-density region. Parameters flagged here are the ones that
# should be prioritized in the prior-sensitivity reruns.

trapz_numeric <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2) return(NA_real_)
  o <- order(x)
  x <- x[o]
  y <- y[o]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

standardise_density <- function(df) {
  df <- df %>%
    filter(is.finite(value), is.finite(density), density >= 0) %>%
    arrange(value)
  if (nrow(df) < 2) return(df[FALSE, ])
  area <- trapz_numeric(df$value, df$density)
  if (!is.finite(area) || area <= 0) return(df[FALSE, ])
  df$density <- df$density / area
  df
}

transform_density_for_metric <- function(df, transform = "identity") {
  df <- df %>% filter(is.finite(value), is.finite(density), density >= 0)
  if (identical(transform, "log")) {
    df <- df %>% filter(value > 0)
    if (!nrow(df)) return(data.frame(value = numeric(), density = numeric()))
    z <- log(df$value)
    # density for z = log(x): f_Z(z) = f_X(exp(z)) * exp(z)
    out <- data.frame(value = z, density = df$density * exp(z))
  } else {
    out <- data.frame(value = df$value, density = df$density)
  }
  standardise_density(out)
}

density_summary_stats <- function(df) {
  df <- standardise_density(df)
  if (nrow(df) < 2) {
    return(list(mean = NA_real_, sd = NA_real_, q025 = NA_real_,
                median = NA_real_, q975 = NA_real_))
  }
  x <- df$value
  d <- df$density
  dx <- diff(x)
  increments <- dx * (head(d, -1) + tail(d, -1)) / 2
  cdf <- c(0, cumsum(increments))
  if (max(cdf, na.rm = TRUE) > 0) cdf <- cdf / max(cdf, na.rm = TRUE)
  qfun <- function(p) {
    keep <- !duplicated(cdf) & is.finite(cdf) & is.finite(x)
    if (sum(keep) < 2) return(NA_real_)
    approx(cdf[keep], x[keep], xout = p, rule = 2)$y
  }
  m <- trapz_numeric(x, x * d)
  v <- trapz_numeric(x, (x - m)^2 * d)
  list(
    mean = m,
    sd = sqrt(pmax(v, 0)),
    q025 = qfun(0.025),
    median = qfun(0.5),
    q975 = qfun(0.975)
  )
}

density_overlap <- function(prior_df, post_df) {
  prior_df <- standardise_density(prior_df)
  post_df <- standardise_density(post_df)
  if (nrow(prior_df) < 2 || nrow(post_df) < 2) return(NA_real_)

  lo <- max(min(prior_df$value), min(post_df$value))
  hi <- min(max(prior_df$value), max(post_df$value))
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) return(0)
  grid <- seq(lo, hi, length.out = 1500)
  p <- approx(prior_df$value, prior_df$density, xout = grid, rule = 1, yleft = 0, yright = 0)$y
  q <- approx(post_df$value, post_df$density, xout = grid, rule = 1, yleft = 0, yright = 0)$y
  ov <- trapz_numeric(grid, pmin(p, q))
  if (is.finite(ov)) pmin(pmax(ov, 0), 1) else NA_real_
}

density_mass_between <- function(df, lower, upper) {
  df <- standardise_density(df)
  if (nrow(df) < 2 || !is.finite(lower) || !is.finite(upper) || lower >= upper) return(NA_real_)
  lo <- max(min(df$value), lower)
  hi <- min(max(df$value), upper)
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) return(0)
  grid <- seq(lo, hi, length.out = 1000)
  dens <- approx(df$value, df$density, xout = grid, rule = 1, yleft = 0, yright = 0)$y
  mass <- trapz_numeric(grid, dens)
  if (is.finite(mass)) pmin(pmax(mass, 0), 1) else NA_real_
}

prior_influence_row <- function(parameter, prior_df, posterior_marginal,
                                metric_scale = "identity",
                                class = "model_parameter") {
  if (is.null(posterior_marginal) || nrow(posterior_marginal) < 2) {
    return(data.frame(
      parameter = parameter,
      class = class,
      metric_scale = metric_scale,
      prior_mean = NA_real_, prior_sd = NA_real_, prior_q025 = NA_real_,
      prior_median = NA_real_, prior_q975 = NA_real_,
      posterior_mean = NA_real_, posterior_sd = NA_real_, posterior_q025 = NA_real_,
      posterior_median = NA_real_, posterior_q975 = NA_real_,
      posterior_prior_sd_ratio = NA_real_,
      posterior_shift_prior_sd = NA_real_,
      prior_mass_in_posterior_95 = NA_real_,
      density_overlap = NA_real_,
      sensitivity_priority_score = NA_real_,
      sensitivity_priority = NA,
      reason = "posterior marginal not available",
      row.names = NULL
    ))
  }

  post_df <- data.frame(value = posterior_marginal[, 1], density = posterior_marginal[, 2])
  p0 <- transform_density_for_metric(prior_df, metric_scale)
  p1 <- transform_density_for_metric(post_df, metric_scale)
  s0 <- density_summary_stats(p0)
  s1 <- density_summary_stats(p1)

  sd_ratio <- if (is.finite(s0$sd) && s0$sd > 0 && is.finite(s1$sd)) s1$sd / s0$sd else NA_real_
  shift <- if (is.finite(s0$sd) && s0$sd > 0 && is.finite(s0$median) && is.finite(s1$median)) {
    abs(s1$median - s0$median) / s0$sd
  } else NA_real_
  mass <- density_mass_between(p0, s1$q025, s1$q975)
  ov <- density_overlap(p0, p1)

  # Heuristic priority score for sensitivity. A high value can mean either:
  #   1) posterior remains wide relative to prior -> possible prior influence;
  #   2) posterior is in low-prior-density space -> prior-data tension.
  comp_width <- if (is.finite(sd_ratio)) pmin(1, sd_ratio / 0.60) else 0
  comp_mass <- if (is.finite(mass)) pmin(1, pmax(0, (0.10 - mass) / 0.10)) else 0
  comp_overlap <- if (is.finite(ov)) pmin(1, pmax(0, (0.20 - ov) / 0.20)) else 0
  score <- max(comp_width, comp_mass, comp_overlap, na.rm = TRUE)

  reasons <- character()
  if (is.finite(sd_ratio) && sd_ratio > 0.60) {
    reasons <- c(reasons, "posterior SD is >60% of prior SD; data may not dominate this prior")
  }
  if (is.finite(mass) && mass < 0.10) {
    reasons <- c(reasons, "posterior 95% interval had <10% prior mass; prior-data tension")
  }
  if (is.finite(ov) && ov < 0.20) {
    reasons <- c(reasons, "prior-posterior density overlap <20%; strong update/tension")
  }
  if (!length(reasons)) reasons <- "no obvious prior-sensitivity signal"

  data.frame(
    parameter = parameter,
    class = class,
    metric_scale = metric_scale,
    prior_mean = s0$mean,
    prior_sd = s0$sd,
    prior_q025 = s0$q025,
    prior_median = s0$median,
    prior_q975 = s0$q975,
    posterior_mean = s1$mean,
    posterior_sd = s1$sd,
    posterior_q025 = s1$q025,
    posterior_median = s1$median,
    posterior_q975 = s1$q975,
    posterior_prior_sd_ratio = sd_ratio,
    posterior_shift_prior_sd = shift,
    prior_mass_in_posterior_95 = mass,
    density_overlap = ov,
    sensitivity_priority_score = score,
    sensitivity_priority = is.finite(score) && score >= 0.5,
    reason = paste(reasons, collapse = "; "),
    row.names = NULL
  )
}

write_prior_influence_diagnostics <- function(fit, settings, family,
                                              label = "pre_sensitivity") {
  rows <- list()

  # Intercept on the log-rate scale.
  intercept_sd_prior <- 1 / sqrt(PRIOR_INTERCEPT_PREC)
  intercept_grid <- seq(qnorm(0.001, PRIOR_INTERCEPT_MEAN, intercept_sd_prior),
                        qnorm(0.999, PRIOR_INTERCEPT_MEAN, intercept_sd_prior),
                        length.out = 1200)
  rows[[length(rows) + 1L]] <- prior_influence_row(
    parameter = "intercept",
    prior_df = data.frame(value = intercept_grid,
                          density = dnorm(intercept_grid, PRIOR_INTERCEPT_MEAN, intercept_sd_prior)),
    posterior_marginal = fit$marginals.fixed[["intercept"]],
    metric_scale = "identity",
    class = "fixed_effect"
  )

  # Month fixed effects on the log-rate-ratio scale.
  month_terms <- grep("^month_[0-9]{4}_[0-9]{2}$", names(fit$marginals.fixed), value = TRUE)
  if (length(month_terms)) {
    month_sd_prior <- 1 / sqrt(PRIOR_MONTH_LOG_RATE_RATIO_PREC)
    month_grid <- seq(qnorm(0.001, 0, month_sd_prior),
                      qnorm(0.999, 0, month_sd_prior),
                      length.out = 1200)
    month_prior <- data.frame(value = month_grid,
                              density = dnorm(month_grid, 0, month_sd_prior))
    for (term in month_terms) {
      rows[[length(rows) + 1L]] <- prior_influence_row(
        parameter = paste0("month_", month_from_term(term), "_vs_", settings$month_reference),
        prior_df = month_prior,
        posterior_marginal = fit$marginals.fixed[[term]],
        metric_scale = "identity",
        class = "fixed_effect"
      )
    }
  }

  # Spatial range, evaluated on log(range) scale for interpretability.
  if (is.null(settings$fix_range_m)) {
    range0 <- settings$prior_range_m[1]
    prob_range <- settings$prior_range_m[2]
    range_grid <- exp(seq(log(pc_range_quantile(0.001, range0, prob_range)),
                          log(pc_range_quantile(0.995, range0, prob_range)),
                          length.out = 1200))
    rows[[length(rows) + 1L]] <- prior_influence_row(
      parameter = "spatial_range",
      prior_df = data.frame(value = range_grid,
                            density = pc_range_density(range_grid, range0, prob_range)),
      posterior_marginal = hyp_marg(fit, PAT_RANGE),
      metric_scale = "log",
      class = "spatial_hyperparameter"
    )
  }

  # Spatial marginal SD, native scale.
  sigma0 <- settings$prior_sigma[1]
  prob_sigma <- settings$prior_sigma[2]
  sigma_grid <- seq(0, pc_sigma_quantile(0.995, sigma0, prob_sigma), length.out = 1200)
  rows[[length(rows) + 1L]] <- prior_influence_row(
    parameter = "spatial_sd",
    prior_df = data.frame(value = sigma_grid,
                          density = pc_sigma_density(sigma_grid, sigma0, prob_sigma)),
    posterior_marginal = hyp_marg(fit, PAT_SIGMA),
    metric_scale = "identity",
    class = "spatial_hyperparameter"
  )

  # Zero inflation, only if relevant to the family.
  if (is_zi(family)) {
    prob_grid <- seq(0.001, 0.999, length.out = 1200)
    rows[[length(rows) + 1L]] <- prior_influence_row(
      parameter = "zero_inflation_probability",
      prior_df = data.frame(value = prob_grid,
                            density = zip_prob_prior_density(prob_grid)),
      posterior_marginal = hyp_marg(fit, PAT_ZPROB),
      metric_scale = "identity",
      class = "likelihood_hyperparameter"
    )
  }

  # NB size, evaluated on log(size) scale.
  if (is_nb(family)) {
    size_grid <- exp(seq(log(nb_size_prior_quantile(0.001)),
                         log(nb_size_prior_quantile(0.995)),
                         length.out = 1200))
    rows[[length(rows) + 1L]] <- prior_influence_row(
      parameter = "negative_binomial_size",
      prior_df = data.frame(value = size_grid,
                            density = nb_size_prior_density(size_grid)),
      posterior_marginal = hyp_marg(fit, PAT_NB_SIZE),
      metric_scale = "log",
      class = "likelihood_hyperparameter"
    )
  }

  out <- do.call(rbind, rows) %>%
    arrange(desc(sensitivity_priority_score), parameter)

  readr::write_csv(out, path_out(paste0(SURVEY_PREFIX, "_prior_influence_screen_", label, ".csv")))

  plot_df <- out %>%
    mutate(parameter = factor(parameter, levels = rev(parameter)))

  p <- ggplot(plot_df, aes(parameter, sensitivity_priority_score)) +
    geom_col() +
    coord_flip() +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey35") +
    scale_y_continuous(limits = c(0, 1), labels = label_number(accuracy = 0.01)) +
    labs(
      title = paste0("Prior-sensitivity priority screen: ", SURVEY_PREFIX),
      subtitle = "Higher values flag parameters where sensitivity runs are most important",
      x = NULL,
      y = "sensitivity priority score"
    ) +
    theme_minimal(base_size = 12)

  ggsave(path_out(paste0(SURVEY_PREFIX, "_prior_influence_priority_", label, ".png")),
         p, width = 8, height = 5.5, dpi = 250)

  detail_plot <- ggplot(out,
                        aes(posterior_prior_sd_ratio,
                            prior_mass_in_posterior_95,
                            label = parameter)) +
    geom_vline(xintercept = 0.60, linetype = 2, colour = "grey55") +
    geom_hline(yintercept = 0.10, linetype = 2, colour = "grey55") +
    geom_point(size = 2.5) +
    geom_text(size = 3, check_overlap = TRUE, nudge_y = 0.03) +
    scale_x_continuous(labels = label_number(accuracy = 0.01)) +
    scale_y_continuous(limits = c(0, 1), labels = label_number(accuracy = 0.01)) +
    labs(
      title = paste0("Prior influence components: ", SURVEY_PREFIX),
      subtitle = "Upper-right means posterior remains prior-width; lower-left means posterior lies in low-prior region",
      x = "posterior SD / prior SD",
      y = "prior mass inside posterior 95% interval"
    ) +
    theme_minimal(base_size = 12)

  ggsave(path_out(paste0(SURVEY_PREFIX, "_prior_influence_components_", label, ".png")),
         detail_plot, width = 8, height = 5.5, dpi = 250)

  flagged <- out %>% filter(isTRUE(sensitivity_priority) | sensitivity_priority == TRUE)
  report <- c(
    "Prior-influence screen:",
    "  Purpose: decide which priors deserve explicit sensitivity reruns before final interpretation.",
    "  This is a heuristic screen, not a formal proof that a prior caused the posterior.",
    "",
    "How to read the metrics:",
    "  posterior_prior_sd_ratio close to 1 means the posterior remains almost as wide as the prior, so data may be weak for that parameter.",
    "  prior_mass_in_posterior_95 below 0.10 means the posterior sits in a low-prior-probability region, suggesting prior-data tension.",
    "  density_overlap below 0.20 means the prior and posterior distributions are very different.",
    "",
    if (nrow(flagged)) {
      c("Parameters prioritized for sensitivity:",
        paste0("  - ", flagged$parameter, ": ", flagged$reason))
    } else {
      "No parameter crossed the automatic sensitivity-priority thresholds."
    },
    "",
    "Full numeric table written to:",
    paste0("  ", path_out(paste0(SURVEY_PREFIX, "_prior_influence_screen_", label, ".csv")))
  )
  writeLines(report, path_out(paste0(SURVEY_PREFIX, "_PRIOR_INFLUENCE_SCREEN_", label, ".txt")))

  invisible(out)
}

write_month_coefficients <- function(fit, model_dat, settings) {
  month_terms <- temporal_month_terms(model_dat)
  if (!length(month_terms)) return(invisible(NULL))
  month_summary <- fit$summary.fixed[month_terms, , drop = FALSE]
  out <- data.frame(
    term = rownames(month_summary),
    month = month_from_term(rownames(month_summary)),
    reference_month = settings$month_reference,
    prediction_month = settings$month_prediction,
    mean_log_rate_ratio = month_summary[, "mean"],
    q025_log_rate_ratio = month_summary[, "0.025quant"],
    q975_log_rate_ratio = month_summary[, "0.975quant"],
    mean_rate_ratio = exp(month_summary[, "mean"]),
    q025_rate_ratio = exp(month_summary[, "0.025quant"]),
    q975_rate_ratio = exp(month_summary[, "0.975quant"]),
    row.names = NULL
  )
  readr::write_csv(out, path_out(paste0(SURVEY_PREFIX, "_month_coefficients.csv")))
  invisible(out)
}

month_log_rate_ratio_mean <- function(fit, month, settings) {
  if (is.null(month) || is.na(month) || !nzchar(month)) return(0)
  if (!is.null(settings$month_reference) &&
      !is.na(settings$month_reference) &&
      month == settings$month_reference) {
    return(0)
  }

  term <- month_term_name(month)
  fixed_names <- rownames(fit$summary.fixed)
  if (term %in% fixed_names) {
    return(as.numeric(fit$summary.fixed[term, "mean"]))
  }

  0
}

write_annualization_weights <- function(fit, model_dat, settings) {
  if (!isTRUE(settings$use_month_effect) || !"month" %in% names(model_dat)) {
    return(list(
      factor = 1,
      label = "single-period encounter-frequency surface",
      weights = NULL
    ))
  }

  effort_col <- if ("total_effort_days" %in% names(model_dat)) {
    "total_effort_days"
  } else if ("effort_days" %in% names(model_dat)) {
    "effort_days"
  } else {
    stop("No effort column found for annualized prediction weighting.")
  }

  month_effort <- model_dat %>%
    dplyr::group_by(month) %>%
    dplyr::summarise(
      effort_days = sum(.data[[effort_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(effort_days), effort_days > 0) %>%
    dplyr::arrange(month)

  total_effort <- sum(month_effort$effort_days)
  if (!nrow(month_effort) || !is.finite(total_effort) || total_effort <= 0) {
    return(list(
      factor = 1,
      label = "single-period encounter-frequency surface",
      weights = month_effort
    ))
  }

  prediction_month <- settings$month_prediction
  prediction_log_ratio <- month_log_rate_ratio_mean(fit, prediction_month, settings)

  month_effort$weight <- month_effort$effort_days / total_effort
  month_effort$log_rate_ratio_to_prediction_month <- vapply(
    month_effort$month,
    function(m) month_log_rate_ratio_mean(fit, m, settings) - prediction_log_ratio,
    numeric(1)
  )
  month_effort$mean_rate_ratio_to_prediction_month <-
    exp(month_effort$log_rate_ratio_to_prediction_month)
  factor <- sum(month_effort$weight *
                  month_effort$mean_rate_ratio_to_prediction_month)
  month_effort$annualization_factor <- factor
  month_effort$prediction_baseline_month <- prediction_month

  readr::write_csv(month_effort,
                   path_out(paste0(SURVEY_PREFIX, "_annualization_weights.csv")))

  list(
    factor = factor,
    label = sprintf(
      "annualized 2024 surface; baseline %s; factor %.3f",
      prediction_month,
      factor
    ),
    weights = month_effort
  )
}

write_model_hyperparameters <- function(fit) {
  if (is.null(fit$summary.hyperpar)) return(invisible(NULL))
  out <- data.frame(parameter = rownames(fit$summary.hyperpar),
                    fit$summary.hyperpar,
                    row.names = NULL,
                    check.names = FALSE)
  readr::write_csv(out, path_out(paste0(SURVEY_PREFIX, "_hyperparameters.csv")))
  invisible(out)
}


## 11. Fit model and prediction maps -----------------------------------------

fit_2024_model <- function(model_dat, settings, family) {
  cat(sprintf("\n[wolf_2024] fitting %s: family=%s\n", FINAL_MODEL_NAME, family))

  model_dat <- model_dat %>% mutate(y = as.integer(wolf_events), intercept = 1)

  camera_summary <- camera_summary_from_model(model_dat)
  camera_sf <- camera_to_utm(camera_summary)
  coords_camera <- st_coordinates(camera_sf)
  colnames(coords_camera) <- c("x", "y")

  obs_sf <- model_dat %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(EPSG_UTM)
  coords_obs <- st_coordinates(obs_sf)
  colnames(coords_obs) <- c("x", "y")

  pred_sf <- prediction_grid(camera_sf, settings)
  coords_pred <- st_coordinates(pred_sf)
  colnames(coords_pred) <- c("x", "y")
  cat(sprintf("[wolf_2024] prediction cells: %d at %.0f m\n",
              nrow(pred_sf), settings$cell_size_m))

  mesh_loc <- if (isTRUE(settings$include_grid_in_mesh)) {
    rbind(coords_camera, coords_pred)
  } else {
    coords_camera
  }
  spde_obj <- build_spatial(mesh_loc, settings)

  fixed_terms <- fixed_effect_terms(model_dat)
  fixed_obs <- as.data.frame(model_dat[, fixed_terms, drop = FALSE])
  fixed_pred <- prediction_fixed_effects(model_dat, fixed_terms, settings, nrow(pred_sf))

  A_obs <- INLA::inla.spde.make.A(spde_obj$mesh, loc = coords_obs)
  A_pred <- INLA::inla.spde.make.A(spde_obj$mesh, loc = coords_pred)

  stack_obs <- INLA::inla.stack(
    tag = "obs",
    data = list(y = model_dat$y, e = model_dat$total_effort_days),
    A = list(A_obs, 1),
    effects = list(spatial = spde_obj$s_index, fixed = fixed_obs)
  )

  stack_pred <- INLA::inla.stack(
    tag = "pred",
    data = list(y = rep(NA_real_, nrow(pred_sf)), e = rep(100, nrow(pred_sf))),
    A = list(A_pred, 1),
    effects = list(spatial = spde_obj$s_index, fixed = fixed_pred)
  )

  stack_all <- INLA::inla.stack(stack_obs, stack_pred)
  stack_data <- INLA::inla.stack.data(stack_all)

  obs_index <- INLA::inla.stack.index(stack_all, tag = "obs")$data
  pred_index <- INLA::inla.stack.index(stack_all, tag = "pred")$data

  formula <- as.formula(
    paste("y ~ 0 +",
          paste(c(fixed_terms, "f(spatial, model = spde_obj$spde)"), collapse = " + "))
  )

  fit <- INLA::inla(
    formula,
    family = family,
    data = stack_data,
    E = stack_data$e,
    control.predictor = list(
      A = INLA::inla.stack.A(stack_all),
      compute = TRUE,
      link = 1
    ),
    control.compute = list(config = FALSE, dic = TRUE, waic = TRUE, cpo = TRUE),
    control.fixed = make_control_fixed(fixed_terms),
    control.family = make_control_family(family),
    verbose = FALSE
  )

  list(fit = fit, stack_all = stack_all, obs_index = obs_index,
       pred_index = pred_index, model_dat = model_dat,
       camera_sf = camera_sf, pred_sf = pred_sf,
       spde_obj = spde_obj, fixed_terms = fixed_terms)
}


fit_2024_diagnostic_model <- function(model_dat, settings, family) {
  cat(sprintf("\n[wolf_2024] fitting observed-data diagnostic model for posterior sampling\n"))

  model_dat <- model_dat %>% mutate(y = as.integer(wolf_events), intercept = 1)

  camera_summary <- camera_summary_from_model(model_dat)
  camera_sf <- camera_to_utm(camera_summary)
  coords_camera <- st_coordinates(camera_sf)
  colnames(coords_camera) <- c("x", "y")

  obs_sf <- model_dat %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(EPSG_UTM)
  coords_obs <- st_coordinates(obs_sf)
  colnames(coords_obs) <- c("x", "y")

  spde_obj <- build_spatial(coords_camera, settings)
  fixed_terms <- fixed_effect_terms(model_dat)
  fixed_obs <- as.data.frame(model_dat[, fixed_terms, drop = FALSE])
  A_obs <- INLA::inla.spde.make.A(spde_obj$mesh, loc = coords_obs)

  stack_obs <- INLA::inla.stack(
    tag = "obs",
    data = list(y = model_dat$y, e = model_dat$total_effort_days),
    A = list(A_obs, 1),
    effects = list(spatial = spde_obj$s_index, fixed = fixed_obs)
  )

  stack_data <- INLA::inla.stack.data(stack_obs)
  obs_index <- INLA::inla.stack.index(stack_obs, tag = "obs")$data

  formula <- as.formula(
    paste("y ~ 0 +",
          paste(c(fixed_terms, "f(spatial, model = spde_obj$spde)"), collapse = " + "))
  )

  fit <- INLA::inla(
    formula,
    family = family,
    data = stack_data,
    E = stack_data$e,
    control.predictor = list(
      A = INLA::inla.stack.A(stack_obs),
      compute = TRUE,
      link = 1
    ),
    control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE),
    control.fixed = make_control_fixed(fixed_terms),
    control.family = make_control_family(family),
    verbose = FALSE
  )

  list(fit = fit, stack_obs = stack_obs, obs_index = obs_index,
       model_dat = model_dat, camera_sf = camera_sf,
       spde_obj = spde_obj, fixed_terms = fixed_terms)
}

make_prediction_outputs <- function(fit_obj, diag, settings, family) {
  fit <- fit_obj$fit
  pred_sf <- fit_obj$pred_sf
  pred_index <- fit_obj$pred_index
  eta_mean <- fit$summary.linear.predictor$mean[pred_index]
  eta_sd <- fit$summary.linear.predictor$sd[pred_index]
  eta_sd <- pmax(eta_sd, 1e-9)
  annualization <- write_annualization_weights(fit, diag$model_dat, settings)
  annual_factor <- annualization$factor

  # Marginal posterior mean of the expected encounter frequency per 100 camera-days.
  # This is the latent mean surface, not a draw of future observation noise.
  rate100_latent_mean <- annual_factor * 100 * exp(eta_mean + 0.5 * eta_sd^2)
  pred_sf$mean <- fam_mean(rate100_latent_mean, diag$pi_hat, family, diag$size_hat)

  # Approximate posterior uncertainty in the latent mean surface from the marginal
  # linear-predictor posterior. This avoids sampling the whole 36k-cell latent grid.
  pred_sf$cv <- sqrt(expm1(eta_sd^2))
  pred_sf$sd <- pred_sf$mean * pred_sf$cv
  pred_sf$annualization_factor <- annual_factor

  model_dat <- diag$model_dat
  overall_rate <- 100 * sum(model_dat$y) / sum(model_dat$total_effort_days)
  threshold <- EXCEED_MULT * overall_rate

  if (MAP_EXCEEDANCE) {
    denom <- annual_factor * 100 * pmax(1 - diag$pi_hat, 1e-12)
    latent_threshold <- threshold / denom
    pred_sf$exceed <- if (is.finite(latent_threshold) && latent_threshold > 0) {
      1 - pnorm((log(latent_threshold) - eta_mean) / eta_sd)
    } else {
      rep(1, length(eta_mean))
    }
    pred_sf$exceed <- pmin(pmax(pred_sf$exceed, 0), 1)
  }

  coords_pred <- st_coordinates(pred_sf)
  pred_sf$x <- coords_pred[, 1]
  pred_sf$y <- coords_pred[, 2]

  wkt <- st_crs(pred_sf)$wkt
  pred_table <- st_drop_geometry(pred_sf)

  readr::write_csv(pred_table,
                   path_out(paste0(SURVEY_PREFIX, "_final_prediction_grid.csv")))

  make_raster <- function(col) {
    terra::rast(pred_table[, c("x", "y", col)], type = "xyz", crs = wkt)
  }

  r_mean <- make_raster("mean")
  r_sd <- make_raster("sd")
  names(r_mean) <- "wolf_events_per_100_camera_days"
  names(r_sd) <- "posterior_sd"

  terra::writeRaster(r_mean,
                     path_out(paste0(SURVEY_PREFIX,
                                     "_final_predicted_events_per_100_days_mean.tif")),
                     overwrite = TRUE)
  terra::writeRaster(r_sd,
                     path_out(paste0(SURVEY_PREFIX,
                                     "_final_predicted_events_per_100_days_sd.tif")),
                     overwrite = TRUE)

  rasters <- list(mean = r_mean, sd = r_sd, exceed = NULL)

  if (MAP_EXCEEDANCE) {
    r_exceed <- make_raster("exceed")
    names(r_exceed) <- "exceedance_probability"
    rasters$exceed <- r_exceed
    terra::writeRaster(r_exceed,
                       path_out(paste0(SURVEY_PREFIX, "_final_exceedance_prob.tif")),
                       overwrite = TRUE)
  }

  plot_map_outputs(fit_obj$camera_sf, diag$model_dat, rasters, overall_rate,
                   annualization)
  invisible(list(pred_sf = pred_sf, rasters = rasters,
                 overall_rate = overall_rate, threshold = threshold,
                 annualization = annualization))
}

plot_map_outputs <- function(camera_sf, model_dat, rasters, overall_rate,
                             annualization = NULL) {
  plot_label <- sub(" survey$", "", SURVEY_LABEL)

  raster_to_df <- function(r, name) {
    d <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
    names(d) <- c("x", "y", name)
    d
  }

  mean_df <- raster_to_df(rasters$mean, "rate")
  cap <- quantile(mean_df$rate, 0.98, na.rm = TRUE)

  camera_obs <- camera_summary_from_model(model_dat)
  positive_sf <- camera_sf
  positive_sf$wolf_events_per_100_days <- camera_obs$wolf_events_per_100_days
  positive_sf$wolf_events <- camera_obs$wolf_events
  positive_sf <- positive_sf[camera_obs$wolf_events > 0, ]

  mean_plot <- ggplot() +
    geom_raster(data = mean_df,
                aes(x, y, fill = pmin(rate, cap)),
                interpolate = TRUE) +
    geom_sf(data = camera_sf,
            shape = 21, size = 1.4, fill = "white",
            colour = "grey35", stroke = 0.25) +
    geom_sf(data = positive_sf,
            aes(size = wolf_events_per_100_days),
            shape = 21, fill = "black", colour = "white",
            stroke = 0.25, alpha = 0.9) +
    scale_fill_viridis_c(option = "magma", na.value = NA,
                         name = "predicted events\n/100 camera-days",
                         labels = label_number(accuracy = 0.01)) +
    scale_size_continuous(range = c(2, 7),
                          name = "observed events\n/100 camera-days",
                          labels = label_number(accuracy = 0.01)) +
    coord_sf(datum = NA) +
    labs(title = paste0("Wolf encounter-frequency surface: ", plot_label),
         subtitle = sprintf(
           "%s (%s)\n%s",
           FINAL_MODEL_NAME,
           FINAL_FAMILY,
           if (!is.null(annualization)) {
             annualization$label
           } else {
             "prediction surface"
           }
         ),
         x = "Easting, UTM 34N", y = "Northing, UTM 34N") +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), legend.position = "right")

  ggsave(path_out(paste0(SURVEY_PREFIX, "_final_event_frequency_mean.png")),
         mean_plot, width = 9.5, height = 9, dpi = 350)

  sd_df <- raster_to_df(rasters$sd, "sd")
  sd_cap <- quantile(sd_df$sd, 0.98, na.rm = TRUE)
  sd_plot <- ggplot() +
    geom_raster(data = sd_df, aes(x, y, fill = pmin(sd, sd_cap)),
                interpolate = TRUE) +
    geom_sf(data = camera_sf,
            shape = 21, size = 1.4, fill = "white",
            colour = "grey35", stroke = 0.25) +
    scale_fill_viridis_c(option = "viridis", na.value = NA,
                         name = "posterior SD\n(events /100 camera-days)",
                         labels = label_number(accuracy = 0.01)) +
    coord_sf(datum = NA) +
    labs(title = paste0("Uncertainty surface: ", plot_label),
         subtitle = "posterior standard deviation of annualized expected encounter frequency",
         x = "Easting, UTM 34N", y = "Northing, UTM 34N") +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), legend.position = "right")

  ggsave(path_out(paste0(SURVEY_PREFIX, "_final_event_frequency_sd.png")),
         sd_plot, width = 9.5, height = 9, dpi = 350)

  if (!is.null(rasters$exceed)) {
    exceed_df <- raster_to_df(rasters$exceed, "p")
    exceed_plot <- ggplot() +
      geom_raster(data = exceed_df, aes(x, y, fill = p), interpolate = TRUE) +
      geom_sf(data = camera_sf,
              shape = 21, size = 1.4, fill = "white",
              colour = "grey35", stroke = 0.25) +
      scale_fill_viridis_c(option = "inferno", limits = c(0, 1),
                           na.value = NA, name = "P(rate > threshold)") +
      coord_sf(datum = NA) +
      labs(title = paste0("Elevated encounter-frequency probability: ", plot_label),
           subtitle = sprintf("annualized surface; threshold = %.2f events / 100 camera-days (%.1fx observed mean)",
                              EXCEED_MULT * overall_rate, EXCEED_MULT),
           x = "Easting, UTM 34N", y = "Northing, UTM 34N") +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank(), legend.position = "right")

    ggsave(path_out(paste0(SURVEY_PREFIX, "_final_exceedance_prob.png")),
           exceed_plot, width = 9.5, height = 9, dpi = 350)
  }

  invisible(TRUE)
}


## 12. Spatial block cross-validation ----------------------------------------

spatial_block_cv <- function(model_dat, settings, family, K = CV_K) {
  cat(sprintf("\n[wolf_2024] spatial block CV for %s (K=%d)\n", FINAL_MODEL_NAME, K))

  model_dat <- model_dat %>% mutate(y = as.integer(wolf_events), intercept = 1)
  fixed_terms <- fixed_effect_terms(model_dat)

  camera_summary <- camera_summary_from_model(model_dat)
  camera_sf <- camera_to_utm(camera_summary)
  coords_camera <- st_coordinates(camera_sf)
  colnames(coords_camera) <- c("x", "y")

  row_sf <- model_dat %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(EPSG_UTM)
  coords_row <- st_coordinates(row_sf)
  colnames(coords_row) <- c("x", "y")

  y <- model_dat$y
  effort <- model_dat$total_effort_days

  K_final <- min(K, nrow(coords_camera) - 1L)
  if (K_final < 2) {
    cat("[wolf_2024] spatial CV skipped: too few cameras.\n")
    return(NULL)
  }

  fold <- tryCatch(
    kmeans(scale(coords_camera), centers = K_final, nstart = 20)$cluster,
    error = function(e) sample(rep_len(seq_len(K_final), nrow(coords_camera)))
  )
  row_fold <- fold[match(model_dat$plotID, camera_summary$plotID)]

  rows <- list()
  cam_rows <- list()
  failed_folds <- character()

  for (f in sort(unique(fold))) {
    cat(sprintf("[wolf_2024]   CV fold %s\n", f))
    test <- which(row_fold == f)
    train <- which(row_fold != f)

    result <- tryCatch({
      # Strict CV: train mesh only. Offset is large enough for edge-fold prediction.
      train_camera_ids <- unique(model_dat$plotID[train])
      train_camera_coords <- coords_camera[camera_summary$plotID %in% train_camera_ids, , drop = FALSE]
      spde_obj <- build_spatial(train_camera_coords, settings)

      A_train <- INLA::inla.spde.make.A(spde_obj$mesh,
                                        loc = coords_row[train, , drop = FALSE])
      A_test <- INLA::inla.spde.make.A(spde_obj$mesh,
                                       loc = coords_row[test, , drop = FALSE])

      fixed_train <- as.data.frame(model_dat[train, fixed_terms, drop = FALSE])
      fixed_test <- as.data.frame(model_dat[test, fixed_terms, drop = FALSE])

      stack_train <- INLA::inla.stack(
        tag = "train",
        data = list(y = y[train], e = effort[train]),
        A = list(A_train, 1),
        effects = list(spatial = spde_obj$s_index, fixed = fixed_train)
      )

      stack_test <- INLA::inla.stack(
        tag = "test",
        data = list(y = rep(NA_real_, length(test)), e = effort[test]),
        A = list(A_test, 1),
        effects = list(spatial = spde_obj$s_index, fixed = fixed_test)
      )

      stack_all <- INLA::inla.stack(stack_train, stack_test)
      stack_data <- INLA::inla.stack.data(stack_all)

      formula <- as.formula(
        paste("y ~ 0 +",
              paste(c(fixed_terms, "f(spatial, model = spde_obj$spde)"),
                    collapse = " + "))
      )

      fit_fold <- INLA::inla(
        formula,
        family = family,
        data = stack_data,
        E = stack_data$e,
        control.predictor = list(
          A = INLA::inla.stack.A(stack_all),
          compute = TRUE,
          link = 1
        ),
        control.compute = list(config = TRUE),
        control.fixed = make_control_fixed(fixed_terms),
        control.family = make_control_family(family),
        verbose = FALSE
      )

      test_index <- INLA::inla.stack.index(stack_all, tag = "test")$data
      samples <- posterior_samples_safe(fit_fold, CV_NSIM)
      test_draws <- build_posterior_draws(fit_fold, samples, test_index,
                                          effort[test], family,
                                          expected_n_stack = length(fit_fold$summary.linear.predictor$mean))
      sim <- simulate_from_draws(test_draws, family)
      Ey <- rowMeans(test_draws$fitted)
      lo <- apply(sim, 1, quantile, 0.05, na.rm = TRUE)
      hi <- apply(sim, 1, quantile, 0.95, na.rm = TRUE)

      lpd <- vapply(seq_along(test), function(j) {
        log_mean_exp(fam_logpmf(y[test[j]], test_draws$mu[j, ],
                                test_draws$pi, family, test_draws$size))
      }, numeric(1))

      row_out <- data.frame(
        fold = f,
        plotID = model_dat$plotID[test],
        deploymentID = model_dat$deploymentID[test],
        month = model_dat$month[test],
        y = y[test],
        Ey = Ey,
        effort_days = effort[test],
        rate_obs = 100 * y[test] / effort[test],
        rate_pred = 100 * Ey / effort[test],
        lpd = lpd,
        lo90 = lo,
        hi90 = hi,
        covered_90 = y[test] >= lo & y[test] <= hi
      )

      group <- as.factor(model_dat$plotID[test])
      sim_cam <- aggregate_matrix_by_group(sim, group)
      if (is.null(dim(sim_cam))) sim_cam <- matrix(sim_cam, nrow = 1)
      y_cam <- as.numeric(rowsum(y[test], group)[, 1])
      effort_cam <- as.numeric(rowsum(effort[test], group)[, 1])
      Ey_cam <- rowMeans(sim_cam)
      lo_cam <- apply(sim_cam, 1, quantile, 0.05, na.rm = TRUE)
      hi_cam <- apply(sim_cam, 1, quantile, 0.95, na.rm = TRUE)
      cam_ids <- rownames(rowsum(y[test], group))

      cam_out <- data.frame(
        fold = f,
        plotID = cam_ids,
        y = y_cam,
        Ey = Ey_cam,
        effort_days = effort_cam,
        rate_obs = 100 * y_cam / effort_cam,
        rate_pred = 100 * Ey_cam / effort_cam,
        lo90 = lo_cam,
        hi90 = hi_cam,
        covered_90 = y_cam >= lo_cam & y_cam <= hi_cam
      )

      list(row = row_out, camera = cam_out)
    }, error = function(e) {
      failed_folds <<- c(failed_folds, paste0(f, ": ", conditionMessage(e)))
      NULL
    })

    if (!is.null(result)) {
      rows[[length(rows) + 1L]] <- result$row
      cam_rows[[length(cam_rows) + 1L]] <- result$camera
    }
  }

  if (length(failed_folds)) {
    stop("[wolf_2024] spatial CV failed for fold(s): ",
         paste(failed_folds, collapse = "; "))
  }
  if (!length(rows)) stop("[wolf_2024] spatial CV produced no successful folds.")

  cv_row <- do.call(rbind, rows)
  cv_cam <- do.call(rbind, cam_rows)

  readr::write_csv(cv_row, path_out(paste0(SURVEY_PREFIX, "_final_spatial_block_cv_rows.csv")))
  readr::write_csv(cv_cam, path_out(paste0(SURVEY_PREFIX, "_final_spatial_block_cv_camera.csv")))

  summary <- bind_rows(
    data.frame(
      level = "model_row",
      metric = c("mean_log_predictive_density", "rmse_count",
                 "rmse_rate_per100", "coverage_90"),
      value = c(mean(cv_row$lpd),
                sqrt(mean((cv_row$y - cv_row$Ey)^2)),
                sqrt(mean((cv_row$rate_obs - cv_row$rate_pred)^2)),
                mean(cv_row$covered_90))
    ),
    data.frame(
      level = "camera",
      metric = c("rmse_count", "rmse_rate_per100", "coverage_90"),
      value = c(sqrt(mean((cv_cam$y - cv_cam$Ey)^2)),
                sqrt(mean((cv_cam$rate_obs - cv_cam$rate_pred)^2)),
                mean(cv_cam$covered_90))
    )
  )

  readr::write_csv(summary,
                   path_out(paste0(SURVEY_PREFIX, "_final_spatial_block_cv_summary.csv")))

  cat(sprintf(
    "[wolf_2024] spatial CV rows: mean LPD %.3f | RMSE count %.2f | RMSE rate %.2f | 90%% coverage %.2f\n",
    summary$value[summary$level == "model_row" & summary$metric == "mean_log_predictive_density"],
    summary$value[summary$level == "model_row" & summary$metric == "rmse_count"],
    summary$value[summary$level == "model_row" & summary$metric == "rmse_rate_per100"],
    summary$value[summary$level == "model_row" & summary$metric == "coverage_90"]
  ))
  cat(sprintf(
    "[wolf_2024] spatial CV cameras: RMSE count %.2f | RMSE rate %.2f | 90%% coverage %.2f\n",
    summary$value[summary$level == "camera" & summary$metric == "rmse_count"],
    summary$value[summary$level == "camera" & summary$metric == "rmse_rate_per100"],
    summary$value[summary$level == "camera" & summary$metric == "coverage_90"]
  ))

  list(row = cv_row, camera = cv_cam, summ = summary)
}


## 13. Reporting --------------------------------------------------------------

diagnostic_failures <- function(diag) {
  c(
    if (!isTRUE(diag$ppc_total_pass)) "camera-level PPC total events" else NULL,
    if (!isTRUE(diag$ppc_zero_pass)) "camera-level PPC zero fraction" else NULL,
    if (!isTRUE(diag$ppc_max_pass)) "camera-level PPC max count" else NULL,
    if (!isTRUE(diag$moran_pass)) {
      if (is.finite(diag$moran_p)) {
        sprintf("residual spatial autocorrelation (Moran's I = %.3f, p = %.3f)",
                diag$moran_I, diag$moran_p)
      } else {
        "residual Moran's I was not evaluable"
      }
    } else NULL
  )
}

prior_lines_for_report <- function(settings, family) {
  c(
    "",
    "Priors:",
    sprintf("  Field range: PC prior, P(range < %d m) = %.2f",
            as.integer(settings$prior_range_m[1]), settings$prior_range_m[2]),
    sprintf("  Field marginal SD: PC prior, P(SD > %.2f) = %.2f",
            settings$prior_sigma[1], settings$prior_sigma[2]),
    sprintf("  Intercept: Gaussian(mean = %.3f, prec = %.3f), SD %.1f on log scale",
            PRIOR_INTERCEPT_MEAN, PRIOR_INTERCEPT_PREC, 1 / sqrt(PRIOR_INTERCEPT_PREC)),
    sprintf("  Month log-rate ratios: Gaussian(0, prec = %g), SD %.1f",
            PRIOR_MONTH_LOG_RATE_RATIO_PREC,
            1 / sqrt(PRIOR_MONTH_LOG_RATE_RATIO_PREC)),
    if (is_zi(family)) {
      sprintf("  Zero-inflation logit(p): Gaussian(%.2f, prec = %.2f)",
              PRIOR_ZI_LOGIT_MEAN, PRIOR_ZI_LOGIT_PREC)
    } else NULL,
    if (is_nb(family)) {
      sprintf("  Negative-binomial log(size): Gaussian(mean = %.2f, prec = %.3f), SD %.1f",
              PRIOR_NB_LOGSIZE_MEAN,
              PRIOR_NB_LOGSIZE_PREC,
              1 / sqrt(PRIOR_NB_LOGSIZE_PREC))
    } else NULL
  )
}

write_validation_report <- function(model_dat, diag, cv, prediction,
                                    temporal_diag = NULL) {
  failures <- diagnostic_failures(diag)
  passes_required <- isTRUE(diag$diagnostics_ok)
  n_cameras <- dplyr::n_distinct(model_dat$plotID)
  annualization <- prediction$annualization

  cv_lines <- if (!is.null(cv)) {
    c(
      "",
      "Spatial block cross-validation:",
      sprintf("  Row mean LPD: %.3f",
              cv$summ$value[cv$summ$level == "model_row" &
                              cv$summ$metric == "mean_log_predictive_density"]),
      sprintf("  Row RMSE count: %.2f",
              cv$summ$value[cv$summ$level == "model_row" &
                              cv$summ$metric == "rmse_count"]),
      sprintf("  Row RMSE rate /100: %.2f",
              cv$summ$value[cv$summ$level == "model_row" &
                              cv$summ$metric == "rmse_rate_per100"]),
      sprintf("  Row 90%% coverage: %.2f",
              cv$summ$value[cv$summ$level == "model_row" &
                              cv$summ$metric == "coverage_90"]),
      sprintf("  Camera RMSE count: %.2f",
              cv$summ$value[cv$summ$level == "camera" &
                              cv$summ$metric == "rmse_count"]),
      sprintf("  Camera RMSE rate /100: %.2f",
              cv$summ$value[cv$summ$level == "camera" &
                              cv$summ$metric == "rmse_rate_per100"]),
      sprintf("  Camera 90%% coverage: %.2f",
              cv$summ$value[cv$summ$level == "camera" &
                              cv$summ$metric == "coverage_90"])
    )
  } else {
    c("", "Spatial block cross-validation: not run")
  }



  temporal_lines <- if (!is.null(temporal_diag) &&
                        !is.null(temporal_diag$lag_summary)) {
    lag1 <- temporal_diag$lag_summary %>% filter(lag == 1L)
    acf1 <- if (!is.null(temporal_diag$acf_summary) &&
                nrow(temporal_diag$acf_summary) >= 2 &&
                any(temporal_diag$acf_summary$lag == 1L)) {
      temporal_diag$acf_summary$acf[temporal_diag$acf_summary$lag == 1L][[1]]
    } else {
      NA_real_
    }

    c(
      "",
      "Temporal residual autocorrelation diagnostics:",
      "  Month is included as a fixed effect; these diagnostics check residual temporal structure after that adjustment.",
      if (nrow(lag1)) {
        sprintf(
          "  Within-camera lag-1 residual correlation: r = %.3f, p = %.4g, n pairs = %d, median gap = %.1f days",
          lag1$correlation[[1]],
          lag1$p_value[[1]],
          lag1$n_pairs[[1]],
          lag1$median_days_between[[1]]
        )
      } else {
        "  Within-camera lag-1 residual correlation: not evaluable"
      },
      if (is.finite(acf1)) {
        sprintf("  Date-ordered mean residual ACF, lag 1: %.3f", acf1)
      } else {
        "  Date-ordered mean residual ACF: not evaluable"
      }
    )
  } else {
    c(
      "",
      "Temporal residual autocorrelation diagnostics: not run"
    )
  }

  report <- c(
    sprintf("Survey: %s", SURVEY_LABEL),
    sprintf("Final model: %s (family = %s)", FINAL_MODEL_NAME, FINAL_FAMILY),
    sprintf("Run profile: %s", RUN_PROFILE),
    sprintf("Joint posterior PPC simulations: %d", diag$ppc_nsim),
    sprintf("Prediction posterior samples: %d", PRED_NSIM),
    sprintf(
      "Cameras: %d | model rows: %d | positive rows: %d | events: %d | effort: %.1f camera-days | observed mean %.3f /100",
      n_cameras,
      nrow(model_dat),
      sum(model_dat$wolf_events > 0),
      sum(model_dat$wolf_events),
      sum(model_dat$total_effort_days),
      100 * sum(model_dat$wolf_events) / sum(model_dat$total_effort_days)
    ),
    "",
    "Validation:",
    sprintf("  PPC method: %s", diag$ppc_method),
    sprintf("  Pearson dispersion, model rows: %.3f", diag$pearson_disp),
    sprintf("  Pearson dispersion, camera aggregates: %.3f", diag$pearson_disp_camera),
    sprintf("  Camera-level PPC total events pass: %s", isTRUE(diag$ppc_total_pass)),
    sprintf("  Camera-level PPC zero fraction pass: %s", isTRUE(diag$ppc_zero_pass)),
    sprintf("  Camera-level PPC max count pass: %s", isTRUE(diag$ppc_max_pass)),
    sprintf("  Residual Moran's I: %.3f; expected %.3f; two-sided p = %.3f",
            diag$moran_I, diag$moran_expected, diag$moran_p),
    sprintf("  Residual Moran pass: %s", isTRUE(diag$moran_pass)),
    sprintf("  Row PIT KS p: %.4g", diag$ppc_pit_ks_row),
    sprintf("  Camera PIT KS p: %.4g", diag$ppc_pit_ks_camera),
    if (is_zi(FINAL_FAMILY)) sprintf("  Zero-inflation probability posterior mean: %.3f", diag$pi_hat) else NULL,
    if (is_nb(FINAL_FAMILY)) sprintf("  Negative-binomial size posterior mean: %.3f", diag$size_hat) else NULL,
    sprintf("  Passes required checks: %s", passes_required),
    if (!passes_required) {
      sprintf("  Stated limitation: %s", paste(failures, collapse = "; "))
    } else NULL,
    cv_lines,
    temporal_lines,
    prior_lines_for_report(settings, FINAL_FAMILY),
    "",
    "Temporal structure:",
    "  Calendar camera-month is included as a fixed effect after splitting effort by month and assigning events by eventStart month.",
    sprintf("  Reference month for coefficients: %s", MONTH_REFERENCE),
    sprintf("  Prediction-stack baseline month: %s", MONTH_PREDICTION),
    if (!is.null(annualization)) {
      sprintf("  Map aggregation: effort-weighted annualized 2024 surface; scale factor %.3f.",
              annualization$factor)
    } else {
      "  Map aggregation: single-period prediction surface."
    },
    "",
    "Prediction:",
    sprintf("  Map units: expected wolf events per 100 camera-days."),
    sprintf("  Exceedance threshold: %.3f events / 100 camera-days (%.1fx observed mean).",
            prediction$threshold, EXCEED_MULT),
    "",
    "Interpretation:",
    "  Relative wolf encounter frequency only.",
    "  Not abundance, density, occupancy, or population size.",
    "  Spatial surface is month-adjusted and annualized over the sampled 2024 months."
  )

  writeLines(report, path_out(paste0(SURVEY_PREFIX, "_VALIDATION_REPORT.txt")))
  invisible(report)
}

write_manifest <- function(model_dat, diag, cv) {
  manifest <- data.frame(
    survey = "2024",
    prefix = SURVEY_PREFIX,
    model = FINAL_MODEL_NAME,
    family = FINAL_FAMILY,
    cameras = dplyr::n_distinct(model_dat$plotID),
    model_rows = nrow(model_dat),
    events = sum(model_dat$wolf_events),
    effort_days = sum(model_dat$total_effort_days),
    month_effect = TRUE,
    diagnostics_ok = isTRUE(diag$diagnostics_ok),
    spatial_cv_run = !is.null(cv),
    run_profile = RUN_PROFILE
  )
  readr::write_csv(manifest, path_out("wolf_2024_run_manifest.csv"))
  invisible(manifest)
}



## 14. Additional science checks: model comparison and mesh ----------------

# These checks are deliberately separated from the main fit. They are not used
# to draw the final map directly; they document whether the selected model is
# actually needed and whether the spatial result is sensitive to common modelling
# choices that reviewers usually ask about.

RUN_MODEL_COMPARISON <- tolower(Sys.getenv("WOLF_RUN_MODEL_COMPARISON", unset = ifelse(RUN_PROFILE == "quick", "false", "true"))) %in%
  c("true", "1", "yes", "y")
RUN_MESH_SENSITIVITY <- tolower(Sys.getenv("WOLF_RUN_MESH_SENSITIVITY", unset = ifelse(RUN_PROFILE == "quick", "false", "true"))) %in%
  c("true", "1", "yes", "y")

safe_value <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x)) default else as.numeric(x[[1]])
}

inla_criteria <- function(fit) {
  c(
    dic = safe_value(fit$dic$dic),
    p_dic = safe_value(fit$dic$p.eff),
    waic = safe_value(fit$waic$waic),
    p_waic = safe_value(fit$waic$p.eff),
    marginal_loglik = safe_value(fit$mlik[1])
  )
}

cpo_failure_rate <- function(fit) {
  if (is.null(fit$cpo)) return(NA_real_)
  vals <- fit$cpo$cpo
  fail <- fit$cpo$failure
  bad_vals <- !is.finite(vals) | vals <= 0
  bad_fail <- if (is.null(fail)) rep(FALSE, length(vals)) else fail > 0
  mean(bad_vals | bad_fail, na.rm = TRUE)
}

mesh_n_vertices <- function(fit_result) {
  if (!is.null(fit_result$spde_obj) && !is.null(fit_result$spde_obj$mesh$n)) {
    return(as.integer(fit_result$spde_obj$mesh$n))
  }
  NA_integer_
}

fit_observed_model_generic <- function(model_dat, settings, family,
                                       spatial = TRUE,
                                       model_label = "model") {
  family <- fit_family(family)
  dat <- model_dat %>% mutate(y = as.integer(wolf_events), intercept = 1)
  fixed_terms <- fixed_effect_terms(dat)
  formula_fixed <- paste(fixed_terms, collapse = " + ")

  if (!isTRUE(spatial)) {
    cat(sprintf("[wolf_2024] model comparison fit: %s | family=%s | spatial=false\n",
                model_label, family))
    fit <- INLA::inla(
      as.formula(paste("y ~ 0 +", formula_fixed)),
      family = family,
      data = dat,
      E = dat$total_effort_days,
      control.predictor = list(compute = TRUE, link = 1),
      control.compute = list(config = FALSE, dic = TRUE, waic = TRUE, cpo = TRUE),
      control.fixed = make_control_fixed(fixed_terms),
      control.family = make_control_family(family),
      verbose = FALSE
    )
    return(list(fit = fit, model_label = model_label, family = family,
                spatial = FALSE, spde_obj = NULL))
  }

  cat(sprintf("[wolf_2024] model comparison fit: %s | family=%s | spatial=true\n",
              model_label, family))
  camera_summary <- camera_summary_from_model(dat)
  camera_sf <- camera_to_utm(camera_summary)
  coords_camera <- st_coordinates(camera_sf)
  colnames(coords_camera) <- c("x", "y")

  obs_sf <- dat %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(EPSG_UTM)
  coords_obs <- st_coordinates(obs_sf)
  colnames(coords_obs) <- c("x", "y")

  spde_obj <- build_spatial(coords_camera, settings)
  A_obs <- INLA::inla.spde.make.A(spde_obj$mesh, loc = coords_obs)
  fixed_obs <- as.data.frame(dat[, fixed_terms, drop = FALSE])

  stack_obs <- INLA::inla.stack(
    tag = "obs",
    data = list(y = dat$y, e = dat$total_effort_days),
    A = list(A_obs, 1),
    effects = list(spatial = spde_obj$s_index, fixed = fixed_obs)
  )
  stack_data <- INLA::inla.stack.data(stack_obs)

  formula <- as.formula(
    paste("y ~ 0 +",
          paste(c(fixed_terms, "f(spatial, model = spde_obj$spde)"), collapse = " + "))
  )

  fit <- INLA::inla(
    formula,
    family = family,
    data = stack_data,
    E = stack_data$e,
    control.predictor = list(
      A = INLA::inla.stack.A(stack_obs),
      compute = TRUE,
      link = 1
    ),
    control.compute = list(config = FALSE, dic = TRUE, waic = TRUE, cpo = TRUE),
    control.fixed = make_control_fixed(fixed_terms),
    control.family = make_control_family(family),
    verbose = FALSE
  )

  list(fit = fit, model_label = model_label, family = family,
       spatial = TRUE, spde_obj = spde_obj)
}

summarise_fitted_model <- function(fit_result, note = NA_character_) {
  fit <- fit_result$fit
  crit <- inla_criteria(fit)
  data.frame(
    model = fit_result$model_label,
    family = fit_result$family,
    spatial = isTRUE(fit_result$spatial),
    note = note,
    dic = crit[["dic"]],
    p_dic = crit[["p_dic"]],
    waic = crit[["waic"]],
    p_waic = crit[["p_waic"]],
    marginal_loglik = crit[["marginal_loglik"]],
    cpo_failure_rate = cpo_failure_rate(fit),
    zi_prob_mean = if (is_zi(fit_result$family)) hyp_point(fit, PAT_ZPROB) else NA_real_,
    nb_size_mean = if (is_nb(fit_result$family)) nb_size_point(fit) else NA_real_,
    spatial_range_mean_m = if (isTRUE(fit_result$spatial)) hyp_point(fit, PAT_RANGE) else NA_real_,
    spatial_sd_mean = if (isTRUE(fit_result$spatial)) hyp_point(fit, PAT_SIGMA) else NA_real_,
    mesh_vertices = mesh_n_vertices(fit_result),
    row.names = NULL
  )
}

run_model_comparison <- function(model_dat, settings) {
  if (!isTRUE(RUN_MODEL_COMPARISON)) {
    writeLines("Model comparison skipped by WOLF_RUN_MODEL_COMPARISON.",
               path_out(paste0(SURVEY_PREFIX, "_MODEL_COMPARISON_SKIPPED.txt")))
    return(NULL)
  }

  cat("\n[wolf_2024] running model comparison\n")
  specs <- data.frame(
    model = c(
      "poisson_spatial_month",
      "nb_spatial_month",
      "zinb_spatial_month"
    ),
    family = c(
      "poisson",
      "nbinomial",
      "zeroinflatednbinomial1"
    ),
    spatial = c(TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  rows <- list()
  failures <- character()
  for (i in seq_len(nrow(specs))) {
    sp <- specs[i, ]
    res <- tryCatch(
      fit_observed_model_generic(model_dat, settings, sp$family,
                                 spatial = sp$spatial,
                                 model_label = sp$model),
      error = function(e) {
        failures <<- c(failures, paste0(sp$model, ": ", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(res)) rows[[length(rows) + 1L]] <- summarise_fitted_model(res)
  }

  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(out)) {
    out <- out %>%
      mutate(delta_waic = waic - min(waic, na.rm = TRUE),
             delta_dic = dic - min(dic, na.rm = TRUE)) %>%
      arrange(waic)
  }
  readr::write_csv(out, path_out(paste0(SURVEY_PREFIX, "_model_comparison.csv")))

  report <- c(
    "Model comparison:",
    "  Purpose: compare spatial Poisson, NB, and ZINB month models to check whether overdispersion or zero inflation is needed.",
    "  Lower WAIC/DIC is better, but spatial block CV and PPC diagnostics should be given more weight than information criteria alone.",
    if (nrow(out)) capture.output(print(out, row.names = FALSE)) else "  No comparison models completed.",
    if (length(failures)) c("", "Failures:", paste0("  ", failures)) else NULL
  )
  writeLines(report, path_out(paste0(SURVEY_PREFIX, "_MODEL_COMPARISON_REPORT.txt")))
  invisible(out)
}

mesh_settings_variants <- function(settings) {
  list(
    current = settings,
    finer = modifyList(settings, list(
      mesh_cutoff_m = max(100, settings$mesh_cutoff_m * 0.70),
      mesh_max_edge = pmax(100, settings$mesh_max_edge * 0.70),
      mesh_offset = pmax(500, settings$mesh_offset * 0.80)
    )),
    coarser = modifyList(settings, list(
      mesh_cutoff_m = settings$mesh_cutoff_m * 1.45,
      mesh_max_edge = settings$mesh_max_edge * 1.45,
      mesh_offset = settings$mesh_offset * 1.20
    ))
  )
}

run_mesh_sensitivity <- function(model_dat, settings, family) {
  if (!isTRUE(RUN_MESH_SENSITIVITY)) {
    writeLines("Mesh sensitivity skipped by profile or WOLF_RUN_MESH_SENSITIVITY.",
               path_out(paste0(SURVEY_PREFIX, "_MESH_SENSITIVITY_SKIPPED.txt")))
    return(NULL)
  }

  cat("\n[wolf_2024] running mesh sensitivity\n")
  variants <- mesh_settings_variants(settings)
  rows <- list()
  failures <- character()

  for (nm in names(variants)) {
    stg <- variants[[nm]]
    label <- paste0("final_", nm, "_mesh")
    res <- tryCatch(
      fit_observed_model_generic(model_dat, stg, family,
                                 spatial = TRUE,
                                 model_label = label),
      error = function(e) {
        failures <<- c(failures, paste0(label, ": ", conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(res)) {
      sm <- summarise_fitted_model(res, note = nm)
      sm$mesh_cutoff_m <- stg$mesh_cutoff_m
      sm$mesh_max_edge_inner_m <- stg$mesh_max_edge[1]
      sm$mesh_max_edge_outer_m <- stg$mesh_max_edge[2]
      rows[[length(rows) + 1L]] <- sm
    }
  }

  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(out)) {
    out <- out %>% mutate(delta_waic = waic - min(waic, na.rm = TRUE))
  }
  readr::write_csv(out, path_out(paste0(SURVEY_PREFIX, "_mesh_sensitivity.csv")))

  report <- c(
    "Mesh sensitivity:",
    "  Purpose: check whether spatial hyperparameters and information criteria are stable under finer/coarser SPDE meshes.",
    "  This is a diagnostic, not an automatic model-selection rule.",
    if (nrow(out)) capture.output(print(out, row.names = FALSE)) else "  No mesh variants completed.",
    if (length(failures)) c("", "Failures:", paste0("  ", failures)) else NULL
  )
  writeLines(report, path_out(paste0(SURVEY_PREFIX, "_MESH_SENSITIVITY_REPORT.txt")))
  invisible(out)
}

extract_first_value <- function(data, field, default = NA_real_) {
  if (is.null(data) || !nrow(data) || !field %in% names(data)) return(default)
  value <- data[[field]][[1]]
  if (is.numeric(value) || is.integer(value)) return(as.numeric(value))
  value
}

temporal_lag1_row <- function(temporal_diag) {
  if (is.null(temporal_diag) || is.null(temporal_diag$lag_summary)) return(data.frame())
  temporal_diag$lag_summary %>% filter(lag == 1L)
}

write_scientific_limitations_report <- function(model_dat) {
  observed_rate <- 100 * sum(model_dat$wolf_events, na.rm = TRUE) /
    sum(model_dat$total_effort_days, na.rm = TRUE)
  n_cameras <- dplyr::n_distinct(model_dat$plotID)
  n_rows <- nrow(model_dat)
  months <- paste(sort(unique(model_dat$month)), collapse = ", ")

  lines <- c(
    "Scientific interpretation and limitations:",
    "",
    sprintf("Survey rows: %d camera-month rows at %d cameras.", n_rows, n_cameras),
    sprintf("Observed mean encounter frequency: %.3f wolf events per 100 camera-days.", observed_rate),
    sprintf("Calendar months represented in the model: %s.", months),
    "",
    "1. Response and interpretation",
    "   The response is independent wolf event count with camera-days as exposure.",
    "   The map is a relative encounter-frequency index, not abundance, density, occupancy, or population size.",
    "",
    "2. Detection and camera-placement bias",
    "   Encounter frequency can reflect wolf activity, placement on trails/roads, camera visibility, camera settings, and detection probability.",
    "   Unless these factors are modelled explicitly, high predicted values should be interpreted as high relative encounter frequency, not necessarily more wolves.",
    "",
    "3. Temporal interpretation",
    "   Calendar camera-month is included as a fixed effect after splitting effort by month and assigning events by eventStart month.",
    "   The final map is an effort-weighted annualized 2024 encounter-frequency surface, not a single-month map.",
    sprintf("   Month %s is retained only as the prediction-stack baseline used to express month-rate ratios.", MONTH_PREDICTION),
    "   Because sampling time and location may be correlated, the spatial surface should be described as month-adjusted rather than purely spatial.",
    "",
    "4. Spatial prediction domain",
    "   Predictions are produced as a full buffered convex-hull map around the camera array.",
    "   Predictions should not be interpreted far outside the sampled camera domain; edge and unsampled-gap areas remain more uncertain.",
    "",
    "5. Missing ecological covariates",
    "   This is a spatial smoothing model. It does not estimate effects of habitat, roads, prey, elevation, human disturbance, or other covariates.",
    "   If ecological explanation is required, those covariates should be added and checked separately.",
    "",
    "6. Event independence",
    "   The analysis assumes eventID represents independent wolf events.",
    "   The manuscript should cite or describe the camera-trap processing rule used to assign independent eventID values, because this rule is part of the data-generating process rather than a parameter estimated by the model."
  )
  writeLines(lines, path_out(paste0(SURVEY_PREFIX, "_SCIENTIFIC_LIMITATIONS.txt")))
  invisible(lines)
}

## 15. Ordered workflow helpers: exploratory checks, prior sensitivity --------

RUN_PRIOR_SENSITIVITY <- tolower(Sys.getenv("WOLF_RUN_PRIOR_SENSITIVITY", unset = ifelse(RUN_PROFILE == "quick", "false", "true"))) %in%
  c("true", "1", "yes", "y")

write_workflow_order_report <- function() {
  lines <- c(
    "Ordered 2024 spatial-modelling workflow:",
    "",
    "1. Data preparation and quality control",
    "   Load deployments and observations, check coordinates, effort, months, and independent wolf-event counts.",
    "",
    "2. Exploratory checks",
    "   Summarise observed encounter rates, month structure, effort, raw spatial pattern, and deployment timing versus northing.",
    "",
    "3. Candidate model comparison",
    "   Fit spatial Poisson/NB/ZINB month models to check whether overdispersion and zero inflation are needed.",
    "",
    "4. Final-model decision",
    "   Keep the configured final model only if it is supported by diagnostics, information criteria, parsimony, and ecological interpretation.",
    "",
    "5. Prior-influence screen",
    "   Fit the selected model with the chosen priors and quantify which priors may be influential or in tension with the data.",
    "",
    "6. Prior sensitivity",
    "   Rerun reasonable alternatives for the priors flagged by the screen, plus core spatial and likelihood priors.",
    "",
    "7. Final model fit and mapping",
    "   Fit the configured final model with the chosen priors and predict relative encounter frequency per 100 camera-days.",
    "",
    "8. Final model diagnostics",
    "   Run posterior predictive checks, PIT, Pearson residual checks, Moran's I, semivariogram, residual covariate plots, and temporal lag/ACF-style checks.",
    "",
    "9. Robustness checks",
    "   Run spatial block cross-validation and mesh sensitivity.",
    "",
    "10. Interpretation and limitations",
    "   Report that the output is relative encounter frequency, not abundance, density, occupancy, or population size."
  )
  writeLines(lines, path_out(paste0(SURVEY_PREFIX, "_ORDERED_WORKFLOW.txt")))
  invisible(lines)
}

write_exploratory_checks <- function(model_dat) {
  cat("\n[wolf_2024] writing exploratory checks\n")

  overall <- data.frame(
    cameras = dplyr::n_distinct(model_dat$plotID),
    deployments = nrow(model_dat),
    positive_deployment_rows = sum(model_dat$wolf_events > 0, na.rm = TRUE),
    zero_deployment_rows = sum(model_dat$wolf_events == 0, na.rm = TRUE),
    events = sum(model_dat$wolf_events, na.rm = TRUE),
    effort_days = sum(model_dat$total_effort_days, na.rm = TRUE),
    observed_rate_per_100_days = 100 * sum(model_dat$wolf_events, na.rm = TRUE) /
      sum(model_dat$total_effort_days, na.rm = TRUE),
    first_start = as.character(min(model_dat$start, na.rm = TRUE)),
    last_start = as.character(max(model_dat$start, na.rm = TRUE)),
    months = paste(sort(unique(model_dat$month)), collapse = ", "),
    row.names = NULL
  )
  readr::write_csv(overall, path_out(paste0(SURVEY_PREFIX, "_exploratory_overall_summary.csv")))

  by_month <- model_dat %>%
    group_by(month) %>%
    summarise(
      deployment_rows = n(),
      cameras = n_distinct(plotID),
      positive_rows = sum(wolf_events > 0, na.rm = TRUE),
      events = sum(wolf_events, na.rm = TRUE),
      effort_days = sum(total_effort_days, na.rm = TRUE),
      rate_per_100 = 100 * events / effort_days,
      .groups = "drop"
    ) %>%
    arrange(month)
  readr::write_csv(by_month, path_out(paste0(SURVEY_PREFIX, "_exploratory_month_summary.csv")))

  camera_summary <- camera_summary_from_model(model_dat %>% mutate(y = wolf_events))
  camera_sf <- camera_to_utm(camera_summary)
  camera_sf$wolf_events <- camera_summary$wolf_events
  camera_sf$wolf_events_per_100_days <- camera_summary$wolf_events_per_100_days

  raw_map <- ggplot(camera_sf) +
    geom_sf(aes(size = wolf_events_per_100_days), shape = 21,
            fill = "black", colour = "white", alpha = 0.85, stroke = 0.25) +
    scale_size_continuous(range = c(1.5, 7), labels = label_number(accuracy = 0.01)) +
    coord_sf(datum = NA) +
    labs(title = "Observed wolf encounter frequency by camera: wolf_2024",
         subtitle = "Raw observed events per 100 camera-days before model smoothing",
         size = "observed events\n/100 days",
         x = "Easting, UTM 34N", y = "Northing, UTM 34N") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())
  ggsave(path_out(paste0(SURVEY_PREFIX, "_exploratory_raw_camera_rates_map.png")),
         raw_map, width = 7.5, height = 7, dpi = 250)

  model_dat$start_date <- as.Date(model_dat$start)
  date_plot <- ggplot(model_dat, aes(start_date, wolf_events_per_100_days)) +
    geom_point(alpha = 0.65) +
    geom_smooth(se = TRUE, method = "loess", formula = y ~ x) +
    scale_y_continuous(trans = "sqrt") +
    labs(title = "Observed deployment-row rates through time: wolf_2024",
         x = "deployment start date",
         y = "observed wolf events / 100 camera-days, square-root scale") +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_exploratory_observed_rate_vs_date.png")),
         date_plot, width = 7, height = 5, dpi = 220)

  month_plot <- ggplot(by_month, aes(month, rate_per_100)) +
    geom_col() +
    labs(title = "Observed monthly wolf encounter frequency: wolf_2024",
         x = "calendar month", y = "observed events / 100 camera-days") +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_exploratory_month_rates.png")),
         month_plot, width = 6.5, height = 4.8, dpi = 220)

  obs_sf <- model_dat %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(EPSG_UTM)
  xy <- st_coordinates(obs_sf)
  timing <- model_dat %>%
    mutate(
      utm_x = xy[, 1],
      utm_y = xy[, 2],
      start_date = as.Date(start),
      start_doy = as.numeric(format(start, "%j"))
    )

  timing_cor <- tryCatch(
    suppressWarnings(cor.test(timing$start_doy, timing$utm_y, method = "spearman")),
    error = function(e) NULL
  )
  timing_summary <- data.frame(
    diagnostic = "deployment_start_day_of_year_vs_northing",
    spearman_rho = if (!is.null(timing_cor)) unname(timing_cor$estimate) else NA_real_,
    p_value = if (!is.null(timing_cor)) timing_cor$p.value else NA_real_,
    n = nrow(timing),
    row.names = NULL
  )
  readr::write_csv(timing_summary, path_out(paste0(SURVEY_PREFIX, "_exploratory_timing_vs_northing.csv")))

  timing_plot <- ggplot(timing, aes(start_doy, utm_y)) +
    geom_point(alpha = 0.65) +
    geom_smooth(se = TRUE, method = "lm", formula = y ~ x) +
    labs(title = "Deployment timing versus northing: wolf_2024",
         subtitle = sprintf("Spearman rho = %.3f, p = %.3g",
                            timing_summary$spearman_rho, timing_summary$p_value),
         x = "deployment start day of year", y = "Northing, UTM 34N") +
    theme_minimal(base_size = 12)
  ggsave(path_out(paste0(SURVEY_PREFIX, "_exploratory_deployment_timing_vs_northing.png")),
         timing_plot, width = 6.8, height = 5, dpi = 220)

  report <- c(
    "Exploratory checks:",
    "",
    capture.output(print(overall, row.names = FALSE)),
    "",
    "Month summary:",
    capture.output(print(by_month, row.names = FALSE)),
    "",
    sprintf("Deployment timing vs northing: Spearman rho = %.3f, p = %.3g, n = %d.",
            timing_summary$spearman_rho,
            timing_summary$p_value,
            timing_summary$n),
    "",
    "Interpretation:",
    "  These checks are descriptive. They justify including month and a spatial field but do not replace posterior diagnostics."
  )
  writeLines(report, path_out(paste0(SURVEY_PREFIX, "_EXPLORATORY_REPORT.txt")))
  invisible(list(overall = overall, by_month = by_month, timing = timing_summary))
}

write_model_choice_report <- function(model_comparison) {
  lines <- c(
    "Final model-choice report:",
    "",
    sprintf("Configured final model: %s, family = %s.", FINAL_MODEL_NAME, fit_family(FINAL_FAMILY)),
    "The configured final model is used for mapping, diagnostics, mesh sensitivity, and spatial CV.",
    "Use model comparison together with posterior predictive checks, spatial CV, residual diagnostics, prior sensitivity, and parsimony."
  )

  if (!is.null(model_comparison) && nrow(model_comparison)) {
    best <- model_comparison[order(model_comparison$waic), ][1, ]
    lines <- c(lines, "", sprintf("Best WAIC model in the comparison table: %s, WAIC = %.2f.", best$model, best$waic))

    pois <- model_comparison[model_comparison$model == "poisson_spatial_month", ]
    nb <- model_comparison[model_comparison$model == "nb_spatial_month", ]
    zinb <- model_comparison[model_comparison$model == "zinb_spatial_month", ]

    if (nrow(zinb)) {
      lines <- c(lines, sprintf("Configured ZINB spatial-month WAIC = %.2f, DIC = %.2f.", zinb$waic[[1]], zinb$dic[[1]]))
      if (is.finite(zinb$zi_prob_mean[[1]])) {
        lines <- c(lines, sprintf("ZINB estimated zero-inflation probability mean in comparison fit: %.3f.",
                                  zinb$zi_prob_mean[[1]]))
      }
      if (is.finite(zinb$nb_size_mean[[1]])) {
        lines <- c(lines, sprintf("ZINB estimated negative-binomial size mean in comparison fit: %.3f.",
                                  zinb$nb_size_mean[[1]]))
      }
    }
    if (nrow(nb) && nrow(zinb)) {
      lines <- c(lines, sprintf("NB spatial-month WAIC = %.2f; delta(NB - ZINB) = %.2f.",
                                nb$waic[[1]], nb$waic[[1]] - zinb$waic[[1]]))
    }
    if (nrow(pois) && nrow(zinb)) {
      lines <- c(lines, sprintf("Poisson spatial-month WAIC = %.2f; delta(Poisson - ZINB) = %.2f.",
                                pois$waic[[1]], pois$waic[[1]] - zinb$waic[[1]]))
      if (is.finite(pois$cpo_failure_rate[[1]])) {
        lines <- c(lines, sprintf("Poisson spatial-month CPO failure rate = %.3f.",
                                  pois$cpo_failure_rate[[1]]))
      }
    }

    if (nrow(zinb) && identical(as.character(best$model[[1]]), "zinb_spatial_month")) {
      lines <- c(lines, "ZINB spatial-month is the best WAIC model in the comparison set; retain it unless PPC/CV diagnostics fail.")
    } else if (nrow(zinb)) {
      lines <- c(lines, "The configured ZINB model is not the best WAIC model; decide final use after PPC, residual diagnostics, spatial CV, and parsimony.")
    }
  } else {
    lines <- c(lines, "", "Model comparison was skipped or failed, so the configured final model remains a prior decision.")
  }

  writeLines(lines, path_out(paste0(SURVEY_PREFIX, "_MODEL_CHOICE_REPORT.txt")))
  invisible(lines)
}

save_prior_state <- function(settings_current) {
  list(
    settings = settings_current,
    PRIOR_INTERCEPT_MEAN = PRIOR_INTERCEPT_MEAN,
    PRIOR_INTERCEPT_PREC = PRIOR_INTERCEPT_PREC,
    PRIOR_MONTH_LOG_RATE_RATIO_PREC = PRIOR_MONTH_LOG_RATE_RATIO_PREC,
    PRIOR_ZI_LOGIT_MEAN = PRIOR_ZI_LOGIT_MEAN,
    PRIOR_ZI_LOGIT_PREC = PRIOR_ZI_LOGIT_PREC,
    PRIOR_NB_LOGSIZE_MEAN = PRIOR_NB_LOGSIZE_MEAN,
    PRIOR_NB_LOGSIZE_PREC = PRIOR_NB_LOGSIZE_PREC
  )
}

restore_prior_state <- function(state) {
  PRIOR_INTERCEPT_MEAN <<- state$PRIOR_INTERCEPT_MEAN
  PRIOR_INTERCEPT_PREC <<- state$PRIOR_INTERCEPT_PREC
  PRIOR_MONTH_LOG_RATE_RATIO_PREC <<- state$PRIOR_MONTH_LOG_RATE_RATIO_PREC
  PRIOR_ZI_LOGIT_MEAN <<- state$PRIOR_ZI_LOGIT_MEAN
  PRIOR_ZI_LOGIT_PREC <<- state$PRIOR_ZI_LOGIT_PREC
  PRIOR_NB_LOGSIZE_MEAN <<- state$PRIOR_NB_LOGSIZE_MEAN
  PRIOR_NB_LOGSIZE_PREC <<- state$PRIOR_NB_LOGSIZE_PREC
  invisible(TRUE)
}

run_prior_sensitivity <- function(model_dat, settings, family, observed_daily_rate) {
  if (!isTRUE(RUN_PRIOR_SENSITIVITY)) {
    writeLines("Prior sensitivity skipped by profile or WOLF_RUN_PRIOR_SENSITIVITY.",
               path_out(paste0(SURVEY_PREFIX, "_PRIOR_SENSITIVITY_SKIPPED.txt")))
    return(NULL)
  }

  cat("\n[wolf_2024] running prior sensitivity\n")
  old <- save_prior_state(settings)
  on.exit(restore_prior_state(old), add = TRUE)

  variants <- list(
    final_current = list(
      settings = settings,
      intercept_mean = log(observed_daily_rate),
      intercept_prec = 1 / 2.5^2,
      zi_mean = qlogis(0.05),
      zi_prec = 1 / 1.5^2,
      nb_mean = log(2),
      nb_prec = 1 / 2^2,
      note = "final selected ZINB priors: wider spatial SD, current range, skeptical zero-inflation, broad log-size prior"
    ),
    wider_spatial_sd_1p5 = list(
      settings = modifyList(settings, list(prior_sigma = c(1.50, 0.05))),
      intercept_mean = log(observed_daily_rate),
      intercept_prec = 1 / 2.5^2,
      zi_mean = qlogis(0.05),
      zi_prec = 1 / 1.5^2,
      nb_mean = log(2),
      nb_prec = 1 / 2^2,
      note = "wider spatial SD prior for comparison"
    ),
    tighter_original_like = list(
      settings = modifyList(settings, list(prior_sigma = c(0.85, 0.05))),
      intercept_mean = 0,
      intercept_prec = 0.01,
      zi_mean = -1,
      zi_prec = 0.2,
      nb_mean = log(2),
      nb_prec = 1 / 2^2,
      note = "tighter original-like spatial SD and broad intercept priors"
    ),
    extra_wide_spatial_sd = list(
      settings = modifyList(settings, list(prior_sigma = c(2.50, 0.05))),
      intercept_mean = log(observed_daily_rate),
      intercept_prec = 1 / 2.5^2,
      zi_mean = qlogis(0.05),
      zi_prec = 1 / 1.5^2,
      nb_mean = log(2),
      nb_prec = 1 / 2^2,
      note = "extra-wide spatial SD prior, checks robustness beyond final prior"
    ),
    shorter_spatial_range = list(
      settings = modifyList(settings, list(prior_range_m = c(2500, 0.5))),
      intercept_mean = log(observed_daily_rate),
      intercept_prec = 1 / 2.5^2,
      zi_mean = qlogis(0.05),
      zi_prec = 1 / 1.5^2,
      nb_mean = log(2),
      nb_prec = 1 / 2^2,
      note = "shorter spatial-range prior median"
    ),
    wider_spatial_range = list(
      settings = modifyList(settings, list(prior_range_m = c(10000, 0.5))),
      intercept_mean = log(observed_daily_rate),
      intercept_prec = 1 / 2.5^2,
      zi_mean = qlogis(0.05),
      zi_prec = 1 / 1.5^2,
      nb_mean = log(2),
      nb_prec = 1 / 2^2,
      note = "wider spatial-range prior median"
    )
  )

  rows <- list()
  failures <- character()
  for (nm in names(variants)) {
    v <- variants[[nm]]
    PRIOR_INTERCEPT_MEAN <<- v$intercept_mean
    PRIOR_INTERCEPT_PREC <<- v$intercept_prec
    PRIOR_ZI_LOGIT_MEAN <<- v$zi_mean
    PRIOR_ZI_LOGIT_PREC <<- v$zi_prec
    PRIOR_NB_LOGSIZE_MEAN <<- v$nb_mean
    PRIOR_NB_LOGSIZE_PREC <<- v$nb_prec

    res <- tryCatch(
      fit_observed_model_generic(model_dat, v$settings, family,
                                 spatial = TRUE,
                                 model_label = paste0("prior_", nm)),
      error = function(e) {
        failures <<- c(failures, paste0(nm, ": ", conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(res)) {
      sm <- summarise_fitted_model(res, note = v$note)
      sm$prior_variant <- nm
      sm$intercept_prior_mean <- v$intercept_mean
      sm$intercept_prior_sd <- 1 / sqrt(v$intercept_prec)
      sm$spatial_sd_prior_sigma0 <- v$settings$prior_sigma[1]
      sm$spatial_sd_prior_prob_above <- v$settings$prior_sigma[2]
      sm$spatial_range_prior_range0 <- v$settings$prior_range_m[1]
      sm$spatial_range_prior_prob_below <- v$settings$prior_range_m[2]
      sm$zi_prior_center <- plogis(v$zi_mean)
      sm$zi_prior_logit_sd <- 1 / sqrt(v$zi_prec)
      sm$nb_logsize_prior_mean <- v$nb_mean
      sm$nb_logsize_prior_sd <- 1 / sqrt(v$nb_prec)
      rows[[length(rows) + 1L]] <- sm
    }
  }

  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(out)) {
    out <- out %>%
      mutate(delta_waic = waic - min(waic, na.rm = TRUE),
             delta_dic = dic - min(dic, na.rm = TRUE)) %>%
      arrange(waic)
  }
  readr::write_csv(out, path_out(paste0(SURVEY_PREFIX, "_prior_sensitivity.csv")))

  report <- c(
    "Prior sensitivity:",
    "  Purpose: check whether the selected spatial model changes under reasonable prior choices.",
    "  This observed-data sensitivity is not the final map; it screens spatial/likelihood hyperparameters and information criteria.",
    if (nrow(out)) capture.output(print(out, row.names = FALSE)) else "  No prior-sensitivity variants completed.",
    if (length(failures)) c("", "Failures:", paste0("  ", failures)) else NULL
  )
  writeLines(report, path_out(paste0(SURVEY_PREFIX, "_PRIOR_SENSITIVITY_REPORT.txt")))
  invisible(out)
}

write_science_checks_summary_ordered <- function(model_comparison, prior_influence, prior_sensitivity, mesh_sensitivity) {
  best_model_line <- if (!is.null(model_comparison) && nrow(model_comparison)) {
    best <- model_comparison[order(model_comparison$waic), ][1, ]
    sprintf("  Best WAIC model among comparison set: %s (WAIC %.2f).", best$model, best$waic)
  } else {
    "  Model comparison was not available."
  }

  prior_influence_line <- if (!is.null(prior_influence) && nrow(prior_influence)) {
    flagged <- prior_influence %>% filter(isTRUE(sensitivity_priority) | sensitivity_priority == TRUE)
    if (nrow(flagged)) {
      paste0("  Prior-influence screen flagged: ", paste(flagged$parameter, collapse = ", "), ".")
    } else {
      "  Prior-influence screen did not flag any parameter at the automatic thresholds."
    }
  } else {
    "  Prior-influence screen was not available."
  }

  prior_line <- if (!is.null(prior_sensitivity) && nrow(prior_sensitivity)) {
    rng <- range(prior_sensitivity$waic, na.rm = TRUE)
    sprintf("  Prior-sensitivity WAIC range across variants: %.2f to %.2f.", rng[1], rng[2])
  } else {
    "  Prior sensitivity was not available."
  }

  mesh_line <- if (!is.null(mesh_sensitivity) && nrow(mesh_sensitivity)) {
    rng <- range(mesh_sensitivity$waic, na.rm = TRUE)
    sprintf("  Mesh-sensitivity WAIC range across variants: %.2f to %.2f.", rng[1], rng[2])
  } else {
    "  Mesh sensitivity was not available."
  }


  lines <- c(
    "Science-check summary, ordered workflow:",
    "",
    "1. Model comparison:",
    best_model_line,
    "  Use this together with PPC and spatial block CV; do not select solely by WAIC/DIC.",
    "",
    "2. Prior-influence screen:",
    prior_influence_line,
    "  Use this screen to decide which priors need explicit sensitivity reruns.",
    "",
    "3. Prior sensitivity:",
    prior_line,
    "  Stable hyperparameters and similar criteria across priors support prior robustness.",
    "",
    "4. Mesh sensitivity:",
    mesh_line,
    "  Stable spatial range/SD and similar WAIC across meshes support mesh robustness.",
    "",
    "5. Full-map prediction domain:",
    "  Final maps use the buffered convex-hull domain only. Disk-based prediction maps have been removed by design.",
    "",
    "Interpretation reminder:",
    "  These models estimate relative encounter frequency only."
  )
  writeLines(lines, path_out(paste0(SURVEY_PREFIX, "_SCIENCE_CHECKS_SUMMARY.txt")))
  invisible(lines)
}


## 16. Main run ---------------------------------------------------------------

validate_inputs()
write_workflow_order_report()

# 1. Data preparation and quality control.
family <- fit_family(FINAL_FAMILY)
model_dat <- load_2024_survey(settings)

observed_daily_rate <- sum(model_dat$wolf_events, na.rm = TRUE) /
  sum(model_dat$total_effort_days, na.rm = TRUE)
if (!is.finite(observed_daily_rate) || observed_daily_rate <= 0) {
  stop("Cannot set intercept prior: observed daily rate is not positive.")
}

# Intercept prior must be set before any model fit because all model-comparison
# and diagnostic fits use the same fixed-effect prior machinery.
PRIOR_INTERCEPT_MEAN <- log(observed_daily_rate)
cat(sprintf(
  "[%s] prior setup: intercept mean %.3f = log(observed daily rate %.4f); spatial SD prior P(SD > %.2f) = %.2f; ZI candidate prior center %.3f; NB candidate log-size prior mean %.3f SD %.2f\n",
  SURVEY_PREFIX,
  PRIOR_INTERCEPT_MEAN,
  observed_daily_rate,
  settings$prior_sigma[1],
  settings$prior_sigma[2],
  plogis(PRIOR_ZI_LOGIT_MEAN),
  PRIOR_NB_LOGSIZE_MEAN,
  1 / sqrt(PRIOR_NB_LOGSIZE_PREC)
))

# 2. Exploratory checks before model selection.
exploratory <- write_exploratory_checks(model_dat)

# 3. Candidate model comparison before deciding whether ZINB is supported relative to Poisson/NB alternatives.
model_comparison <- run_model_comparison(model_dat, settings)

# 4. Provisional final-model decision report. The configured family is fixed in this script; this report documents whether the configured ZINB model is supported.
write_model_choice_report(model_comparison)

# 5. Prior-influence screen after model comparison and before sensitivity runs.
# This observed-data fit is inexpensive relative to the final mapping fit and is
# used only to decide which priors deserve explicit sensitivity checks.
cat("\n[wolf_2024] running prior-influence screen before sensitivity reruns\n")
prior_screen_fit <- fit_observed_model_generic(
  model_dat,
  settings,
  family,
  spatial = TRUE,
  model_label = "prior_influence_screen"
)
prior_influence <- write_prior_influence_diagnostics(
  prior_screen_fit$fit,
  settings,
  family,
  label = "pre_sensitivity"
)

# 6. Prior sensitivity after the influence screen and before the final map.
prior_sensitivity <- run_prior_sensitivity(model_dat, settings, family, observed_daily_rate)

# Restore the chosen final relaxed priors after the sensitivity loop.
PRIOR_INTERCEPT_MEAN <- log(observed_daily_rate)
PRIOR_INTERCEPT_PREC <- 1 / 2.5^2
PRIOR_ZI_LOGIT_MEAN <- qlogis(0.05)  # skeptical but flexible ZINB prior
PRIOR_ZI_LOGIT_PREC <- 1 / 1.5^2
PRIOR_NB_LOGSIZE_MEAN <- log(2)
PRIOR_NB_LOGSIZE_PREC <- 1 / 2^2

# 7. Main final model fit and mapping fit.
fit_obj <- fit_2024_model(model_dat, settings, family)

# 8. Final diagnostics from observed-data diagnostic fit with posterior samples.
diag_fit_obj <- fit_2024_diagnostic_model(model_dat, settings, family)

cat(sprintf("\n[wolf_2024] drawing %d posterior samples for PPC and diagnostics from observed-data fit\n", PPC_NSIM))
ppc_samples <- posterior_samples_safe(diag_fit_obj$fit, PPC_NSIM)

diag <- compute_diagnostics(
  fit = diag_fit_obj$fit,
  samples = ppc_samples,
  model_dat = diag_fit_obj$model_dat,
  obs_index = diag_fit_obj$obs_index,
  camera_sf = diag_fit_obj$camera_sf,
  family = family,
  write_files = TRUE
)

write_diagnostic_plots(diag, diag_fit_obj$camera_sf)
temporal_diag <- temporal_autocorrelation_diagnostics(diag$model_dat)
write_prior_posterior_plots(diag_fit_obj$fit, settings, family)
write_month_coefficients(diag_fit_obj$fit, diag$model_dat, settings)
write_model_hyperparameters(diag_fit_obj$fit)

prediction <- make_prediction_outputs(fit_obj, diag, settings, family)

# 9. Robustness checks after the final model is fitted.
cv <- NULL
if (RUN_SPATIAL_CV) {
  cv <- spatial_block_cv(model_dat, settings, family, K = CV_K)
}
mesh_sensitivity <- run_mesh_sensitivity(model_dat, settings, family)
# 10. Final summaries.
write_science_checks_summary_ordered(model_comparison, prior_influence, prior_sensitivity, mesh_sensitivity)
write_validation_report(model_dat, diag, cv, prediction, temporal_diag)
write_manifest(model_dat, diag, cv)

cat(sprintf("\n[wolf_2024] VALIDATION\n"))
cat(sprintf("[wolf_2024]   PPC method: %s\n", diag$ppc_method))
cat(sprintf("[wolf_2024]   Pearson dispersion, model rows: %.3f\n", diag$pearson_disp))
cat(sprintf("[wolf_2024]   Pearson dispersion, camera aggregates: %.3f\n", diag$pearson_disp_camera))
cat(sprintf("[wolf_2024]   PPC pass, camera total / zero / max: %s / %s / %s\n",
            isTRUE(diag$ppc_total_pass),
            isTRUE(diag$ppc_zero_pass),
            isTRUE(diag$ppc_max_pass)))
cat(sprintf("[wolf_2024]   residual Moran's I: %.3f (expected %.3f, two-sided p = %.3f)\n",
            diag$moran_I, diag$moran_expected, diag$moran_p))
cat(sprintf("[wolf_2024]   row PIT KS p: %.3g | camera PIT KS p: %.3g\n",
            diag$ppc_pit_ks_row, diag$ppc_pit_ks_camera))
if (!is.null(temporal_diag) && nrow(temporal_diag$lag_summary)) {
  lag1_console <- temporal_diag$lag_summary %>% filter(lag == 1L)
  if (nrow(lag1_console)) {
    cat(sprintf("[wolf_2024]   temporal lag-1 residual correlation: r = %.3f (p = %.3g, n pairs = %d)\n",
                lag1_console$correlation[[1]],
                lag1_console$p_value[[1]],
                lag1_console$n_pairs[[1]]))
  }
}
if (is_zi(family)) {
  cat(sprintf("[wolf_2024]   zero-inflation probability mean: %.3f\n", diag$pi_hat))
}
if (is_nb(family)) {
  cat(sprintf("[wolf_2024]   negative-binomial size mean: %.3f\n", diag$size_hat))
}
cat(sprintf("[wolf_2024]   passes required checks: %s\n", isTRUE(diag$diagnostics_ok)))

if (!isTRUE(diag$diagnostics_ok)) {
  cat(sprintf("[wolf_2024]   NOTE: mapped despite failing %s.\n",
              paste(diagnostic_failures(diag), collapse = "; ")))
  cat("[wolf_2024]         Treat flagged aspects as stated limitations.\n")
}

if (!is.null(model_comparison) && nrow(model_comparison)) {
  cat(sprintf("[wolf_2024]   model comparison written: %s\n",
              path_out(paste0(SURVEY_PREFIX, "_model_comparison.csv"))))
}
if (!is.null(prior_sensitivity) && nrow(prior_sensitivity)) {
  cat(sprintf("[wolf_2024]   prior sensitivity written: %s\n",
              path_out(paste0(SURVEY_PREFIX, "_prior_sensitivity.csv"))))
}
if (!is.null(mesh_sensitivity) && nrow(mesh_sensitivity)) {
  cat(sprintf("[wolf_2024]   mesh sensitivity written: %s\n",
              path_out(paste0(SURVEY_PREFIX, "_mesh_sensitivity.csv"))))
}
cat("\nAll 2024 workflow steps completed in the ordered analysis sequence.\n")
cat("Final outputs are in:\n  ", OUTPUT_DIR, "\n", sep = "")
cat("Key files:\n")
cat("  wolf_2024_ORDERED_WORKFLOW.txt\n")
cat("  wolf_2024_EXPLORATORY_REPORT.txt\n")
cat("  wolf_2024_MODEL_CHOICE_REPORT.txt\n")
cat("  wolf_2024_VALIDATION_REPORT.txt\n")
cat("  wolf_2024_SCIENCE_CHECKS_SUMMARY.txt\n")
cat("  wolf_2024_model_comparison.csv / wolf_2024_MODEL_COMPARISON_REPORT.txt\n")
cat("  wolf_2024_prior_sensitivity.csv / wolf_2024_PRIOR_SENSITIVITY_REPORT.txt\n")
cat("  wolf_2024_mesh_sensitivity.csv / wolf_2024_MESH_SENSITIVITY_REPORT.txt\n")
cat("  wolf_2024_final_predicted_events_per_100_days_mean.tif / _sd.tif\n")
cat("  wolf_2024_final_event_frequency_mean.png / _sd.png\n")
cat("  wolf_2024_zinb_spatial_month_posterior_predictive_check.csv\n")
cat("  wolf_2024_zinb_spatial_month_model_row_diagnostics.csv\n")
cat("  wolf_2024_zinb_spatial_month_camera_residual_diagnostics.csv\n")
cat("  wolf_2024_zinb_spatial_month_fitted_scale_sanity_check.csv\n")
cat("  wolf_2024_prior_posterior_*.png / .csv\n")
cat("  wolf_2024_hyperparameters.csv\n")
cat("  wolf_2024_month_coefficients.csv\n")
cat("  wolf_2024_month_observed_summary.csv\n")
cat("  wolf_2024_zinb_spatial_month_TEMPORAL_AUTOCORRELATION_REPORT.txt\n")
cat("  wolf_2024_zinb_spatial_month_temporal_within_camera_lag_correlation.csv\n")
cat("  wolf_2024_final_spatial_block_cv_summary.csv, if CV was run\n")
cat("  wolf_2024_run_manifest.csv\n")

###############################################################################
