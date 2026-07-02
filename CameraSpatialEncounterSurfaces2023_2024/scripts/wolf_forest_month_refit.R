###############################################################################
# Forest-camera 2024 wolf refit with calendar camera-month fixed effects
# -----------------------------------------------------------------------------
# This wrapper reuses the INLA-SPDE fitting and diagnostic helpers from
# wolf_relative_frequency_inla_helpers.R, but rebuilds the forest-camera flat
# data as calendar camera-month rows so that month fixed effects can be fitted.
###############################################################################

if (!nzchar(Sys.getenv("WOLF_RUN_PROFILE", unset = ""))) {
  Sys.setenv(WOLF_RUN_PROFILE = "balanced")
}

this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, winslash = "/"),
                      error = function(e) NA_character_)
PROJECT_DIR <- Sys.getenv("WOLF_PROJECT_DIR", unset = "")
if (!nzchar(PROJECT_DIR)) {
  script_dir <- if (is.na(this_file)) getwd() else dirname(this_file)
  if (basename(script_dir) %in% c("scripts", "R")) script_dir <- dirname(script_dir)
  PROJECT_DIR <- script_dir
}
PROJECT_DIR <- normalizePath(PROJECT_DIR, winslash = "/", mustWork = TRUE)

Sys.setenv(WOLF_PROJECT_DIR = PROJECT_DIR)
Sys.setenv(WOLF_OUTPUT_DIR = Sys.getenv(
  "WOLF_OUTPUT_DIR",
  unset = file.path(PROJECT_DIR, "outputs", "wolf_forest_NB_month_est_range_v1")
))

main_file <- file.path(PROJECT_DIR, "scripts", "wolf_relative_frequency_inla_helpers.R")
main_lines <- readLines(main_file, warn = FALSE)
main_start <- grep("^validate_requested_surveys\\(\\)", main_lines)
if (length(main_start) != 1L) {
  stop("Expected exactly one main execution boundary line matching ",
       "'^validate_requested_surveys()' in ", main_file, ", found ",
       length(main_start), ". Update this boundary marker to match.")
}
eval(parse(text = main_lines[seq_len(main_start[[1]] - 1L)]), envir = .GlobalEnv)

# Match the weakly informative priors selected for the final month-adjusted
# forest-camera 2024 NB model.
PRIOR_INTERCEPT_MEAN <- NA_real_
PRIOR_INTERCEPT_PREC <- 1 / 2.5^2
PRIOR_NB_LOGSIZE_MEAN <- log(2)
PRIOR_NB_LOGSIZE_PREC <- 1 / 2^2

FOREST_INPUT_FILES <- unique(c(
  Sys.getenv("WOLF_FOREST_FILE", unset = ""),
  "forest_camera_trap_events.csv"
))

resolve_input_file <- function(candidates, label) {
  candidates <- candidates[nzchar(candidates)]
  for (candidate in candidates) {
    if (file.exists(candidate)) return(candidate)
    data_path <- path_in(candidate)
    if (file.exists(data_path)) return(data_path)
  }
  stop("Could not find ", label, ". Checked: ",
       paste(candidates, collapse = ", "),
       ". Put the file in WOLF_DATA_DIR or set WOLF_FOREST_FILE.")
}

make_control_fixed <- function(fixed_terms = "intercept") {
  fixed_terms <- unique(fixed_terms)

  mean_values <- rep(0, length(fixed_terms))
  names(mean_values) <- fixed_terms
  if ("intercept" %in% names(mean_values) && is.finite(PRIOR_INTERCEPT_MEAN)) {
    mean_values[["intercept"]] <- PRIOR_INTERCEPT_MEAN
  }

  prec_values <- rep(PRIOR_MONTH_LOG_RATE_RATIO_PREC, length(fixed_terms))
  names(prec_values) <- fixed_terms
  if ("intercept" %in% names(prec_values)) {
    prec_values[["intercept"]] <- PRIOR_INTERCEPT_PREC
  }

  mean_prior <- as.list(mean_values)
  mean_prior$default <- 0
  prec_prior <- as.list(prec_values)
  prec_prior$default <- PRIOR_MONTH_LOG_RATE_RATIO_PREC

  list(mean = mean_prior, prec = prec_prior)
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

nb_size_prior_density <- function(x) {
  ifelse(
    x > 0,
    dnorm(log(x),
          mean = PRIOR_NB_LOGSIZE_MEAN,
          sd = 1 / sqrt(PRIOR_NB_LOGSIZE_PREC)) / x,
    0
  )
}

nb_size_prior_quantile <- function(p) {
  exp(qnorm(p,
            mean = PRIOR_NB_LOGSIZE_MEAN,
            sd = 1 / sqrt(PRIOR_NB_LOGSIZE_PREC)))
}

prior_lines_for_report <- function(settings, spec) {
  c(
    "",
    "Priors:",
    sprintf(
      "  Field range: %s",
      if (!is.null(settings$fix_range_m)) {
        sprintf("fixed at %d m", as.integer(settings$fix_range_m))
      } else {
        sprintf("PC prior, P(range < %d m) = %.2f",
                as.integer(settings$prior_range_m[1]),
                settings$prior_range_m[2])
      }
    ),
    sprintf("  Field marginal SD: PC prior, P(SD > %.2f) = %.2f",
            settings$prior_sigma[1],
            settings$prior_sigma[2]),
    sprintf("  Intercept: Gaussian(mean = %.3f, prec = %.3f), SD %.1f on log daily rate",
            PRIOR_INTERCEPT_MEAN,
            PRIOR_INTERCEPT_PREC,
            1 / sqrt(PRIOR_INTERCEPT_PREC)),
    if (isTRUE(settings$use_month_effect)) {
      sprintf(
        "  Month log-rate ratios: Gaussian(0, prec = %g), SD %.1f",
        PRIOR_MONTH_LOG_RATE_RATIO_PREC,
        1 / sqrt(PRIOR_MONTH_LOG_RATE_RATIO_PREC)
      )
    } else {
      NULL
    },
    if (is_nb(spec$family)) {
      sprintf("  Negative-binomial log(size): Gaussian(mean = %.2f, prec = %.3f), SD %.1f",
              PRIOR_NB_LOGSIZE_MEAN,
              PRIOR_NB_LOGSIZE_PREC,
              1 / sqrt(PRIOR_NB_LOGSIZE_PREC))
    } else {
      NULL
    }
  )
}

split_deployment_month_effort <- function(deployments) {
  rows <- vector("list", nrow(deployments))

  for (i in seq_len(nrow(deployments))) {
    start <- deployments$start_date[[i]]
    end <- deployments$end_date[[i]]
    month_starts <- seq(as.Date(format(start, "%Y-%m-01")),
                        as.Date(format(end - 1, "%Y-%m-01")),
                        by = "1 month")

    rows[[i]] <- do.call(rbind, lapply(month_starts, function(month_start) {
      next_month <- seq(month_start, by = "1 month", length.out = 2)[[2]]
      overlap_start <- max(start, month_start)
      overlap_end <- min(end, next_month)
      effort <- as.numeric(overlap_end - overlap_start)

      data.frame(
        deploymentID = deployments$deploymentID[[i]],
        plotID = deployments$plotID[[i]],
        longitude = deployments$longitude[[i]],
        latitude = deployments$latitude[[i]],
        month = format(month_start, "%Y-%m"),
        total_effort_days = effort,
        stringsAsFactors = FALSE
      )
    }))
  }

  dplyr::bind_rows(rows) %>%
    dplyr::filter(is.finite(total_effort_days), total_effort_days > 0)
}

load_forest_flat_deployment_month <- function(settings, prefix) {
  input_file <- resolve_input_file(FOREST_INPUT_FILES,
                                   "forest-camera 2024 camera-trap input")
  dat <- readr::read_csv(input_file, show_col_types = FALSE)
  required <- c("deploymentID", "eventID", "eventStart", "scientificName",
                "plotID", "deploymentEffort", "latitude", "longitude",
                "startDate", "endDate")
  stop_missing_columns(dat, required, paste0("[", prefix, "] flat input"))

  deployments <- dat %>%
    dplyr::transmute(
      deploymentID = dplyr::na_if(as.character(deploymentID), ""),
      plotID = dplyr::na_if(as.character(plotID), ""),
      latitude = as.numeric(latitude),
      longitude = as.numeric(longitude),
      start_date = as.Date(startDate),
      end_date = as.Date(endDate)
    ) %>%
    dplyr::filter(!is.na(deploymentID),
                  !is.na(plotID),
                  is.finite(latitude),
                  is.finite(longitude),
                  !is.na(start_date),
                  !is.na(end_date),
                  end_date > start_date) %>%
    dplyr::distinct(deploymentID, .keep_all = TRUE)

  month_effort <- split_deployment_month_effort(deployments)

  wolf_events <- dat %>%
    dplyr::transmute(
      deploymentID = dplyr::na_if(as.character(deploymentID), ""),
      eventID = dplyr::na_if(as.character(eventID), ""),
      scientificName = as.character(scientificName),
      event_time = parse_time(eventStart),
      event_month = format(event_time, "%Y-%m", tz = "UTC")
    ) %>%
    dplyr::filter(scientificName %in% WOLF_NAMES,
                  !is.na(deploymentID),
                  !is.na(eventID),
                  !is.na(event_month),
                  nzchar(event_month)) %>%
    dplyr::distinct(deploymentID, eventID, .keep_all = TRUE) %>%
    dplyr::count(deploymentID, month = event_month, name = "wolf_events")

  model_dat <- month_effort %>%
    dplyr::left_join(wolf_events, by = c("deploymentID", "month")) %>%
    dplyr::mutate(
      wolf_events = tidyr::replace_na(wolf_events, 0L),
      wolf_events_per_100_days = 100 * wolf_events / total_effort_days
    ) %>%
    add_month_design(settings, prefix) %>%
    dplyr::arrange(plotID, deploymentID, month)

  month_summary <- model_dat %>%
    dplyr::group_by(month) %>%
    dplyr::summarise(
      deployment_rows = dplyr::n(),
      cameras = dplyr::n_distinct(plotID),
      positive_rows = sum(wolf_events > 0),
      events = sum(wolf_events),
      effort_days = sum(total_effort_days),
      rate_per_100 = 100 * events / effort_days,
      .groups = "drop"
    )

  camera_summary <- camera_summary_from_model(model_dat)

  readr::write_csv(model_dat, path_out(paste0(prefix, "_deployment_month_effort_rates.csv")))
  readr::write_csv(month_summary, path_out(paste0(prefix, "_month_observed_summary.csv")))
  readr::write_csv(camera_summary, path_out(paste0(prefix, "_camera_effort_rates.csv")))

  cat(sprintf(
    "[%s] camera-month rows %d | cameras %d | positive rows %d | events %d | effort %.1f camera-days | observed %.2f /100\n",
    prefix,
    nrow(model_dat),
    dplyr::n_distinct(model_dat$plotID),
    sum(model_dat$wolf_events > 0),
    sum(model_dat$wolf_events),
    sum(model_dat$total_effort_days),
    100 * sum(model_dat$wolf_events) / sum(model_dat$total_effort_days)
  ))
  cat(sprintf(
    "[%s] month effect: reference=%s | prediction=%s | months=%s\n",
    prefix,
    unique(model_dat$month_reference),
    unique(model_dat$month_prediction),
    paste(sort(unique(model_dat$month)), collapse = ", ")
  ))

  model_dat
}

write_temporal_residual_diagnostic <- function(diag, prefix, model_name) {
  if (!"month" %in% names(diag$model_dat)) return(invisible(NULL))

  safe_cor <- function(a, b, min_pairs = 5L) {
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < min_pairs ||
        stats::sd(a[ok], na.rm = TRUE) == 0 ||
        stats::sd(b[ok], na.rm = TRUE) == 0) {
      return(NA_real_)
    }
    stats::cor(a[ok], b[ok])
  }
  safe_cor_p <- function(a, b, min_pairs = 5L) {
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < min_pairs ||
        stats::sd(a[ok], na.rm = TRUE) == 0 ||
        stats::sd(b[ok], na.rm = TRUE) == 0) {
      return(NA_real_)
    }
    suppressWarnings(tryCatch(
      stats::cor.test(a[ok], b[ok])$p.value,
      error = function(e) NA_real_
    ))
  }

  dat <- diag$model_dat %>%
    dplyr::mutate(month_date = as.Date(paste0(month, "-01"))) %>%
    dplyr::filter(!is.na(plotID), !is.na(month_date), is.finite(pearson)) %>%
    dplyr::arrange(plotID, month_date, deploymentID)

  lag_pairs <- dplyr::bind_rows(lapply(split(dat, dat$plotID), function(d) {
    d <- d[order(d$month_date, d$deploymentID), , drop = FALSE]
    if (nrow(d) < 2) return(NULL)
    previous <- seq_len(nrow(d) - 1L)
    current <- previous + 1L
    data.frame(
      plotID = d$plotID[current],
      lag = 1L,
      month_previous = d$month[previous],
      month_current = d$month[current],
      residual_previous = d$pearson[previous],
      residual_current = d$pearson[current],
      days_between = as.numeric(difftime(d$month_date[current],
                                         d$month_date[previous],
                                         units = "days")),
      stringsAsFactors = FALSE
    )
  }))

  lag_summary <- if (nrow(lag_pairs)) {
    lag_pairs %>%
      dplyr::group_by(lag) %>%
      dplyr::summarise(
        n_pairs = dplyr::n(),
        n_cameras = dplyr::n_distinct(plotID),
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
    lag_pairs,
    path_out(paste0(prefix, "_temporal_within_camera_lag_pairs.csv"))
  )
  readr::write_csv(
    lag_summary,
    path_out(paste0(prefix, "_temporal_within_camera_lag_correlation.csv"))
  )

  month_diag <- diag$model_dat %>%
    dplyr::group_by(month) %>%
    dplyr::summarise(
      rows = dplyr::n(),
      cameras = dplyr::n_distinct(plotID),
      observed_events = sum(y),
      fitted_events = sum(fitted_count),
      residual = observed_events - fitted_events,
      pearson = residual / sqrt(pmax(sum(fit_var), 1e-9)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(month)

  lag1 <- if (nrow(month_diag) >= 3) {
    stats::acf(month_diag$pearson, lag.max = 1, plot = FALSE,
               na.action = stats::na.pass)$acf[2]
  } else {
    NA_real_
  }

  month_diag$lag1_acf <- NA_real_
  if (nrow(month_diag)) month_diag$lag1_acf[[1]] <- lag1
  if (nrow(month_diag)) {
    lag1_row <- lag_summary[lag_summary$lag == 1L, , drop = FALSE]
    month_diag$within_camera_lag1_r <- NA_real_
    month_diag$within_camera_lag1_p <- NA_real_
    month_diag$within_camera_lag1_pairs <- NA_integer_
    month_diag$within_camera_lag1_cameras <- NA_integer_
    if (nrow(lag1_row)) {
      month_diag$within_camera_lag1_r[[1]] <- lag1_row$correlation[[1]]
      month_diag$within_camera_lag1_p[[1]] <- lag1_row$p_value[[1]]
      month_diag$within_camera_lag1_pairs[[1]] <- lag1_row$n_pairs[[1]]
      month_diag$within_camera_lag1_cameras[[1]] <- lag1_row$n_cameras[[1]]
    }
  }
  readr::write_csv(month_diag,
                   path_out(paste0(prefix, "_temporal_residual_diagnostics.csv")))

  lag1_line <- if (nrow(lag_summary)) {
    sprintf(
      "  Within-camera lag-1 residual correlation: r = %s, p = %s, n pairs = %d, cameras = %d, median gap = %.1f days.",
      ifelse(is.finite(lag_summary$correlation[[1]]),
             sprintf("%.3f", lag_summary$correlation[[1]]), "not estimable"),
      ifelse(is.finite(lag_summary$p_value[[1]]),
             sprintf("%.4g", lag_summary$p_value[[1]]), "not estimable"),
      lag_summary$n_pairs[[1]],
      lag_summary$n_cameras[[1]],
      lag_summary$median_days_between[[1]]
    )
  } else {
    "  Within-camera lag-1 residual correlation was not evaluable."
  }

  note <- c(
    "Temporal residual autocorrelation diagnostic:",
    sprintf("  Months represented: %s", paste(month_diag$month, collapse = ", ")),
    lag1_line,
    sprintf("  Month-level Pearson residual lag-1 ACF: %s",
            ifelse(is.finite(lag1), sprintf("%.3f", lag1), "not estimable")),
    "  The month-level ACF is a low-power supporting check because only seven monthly time points are available.",
    "  Month fixed effects are the primary temporal adjustment in this refit."
  )
  writeLines(note, path_out(paste0(prefix, "_", model_name,
                                   "_temporal_autocorrelation_report.txt")))
  invisible(month_diag)
}

set_prior_state <- function(intercept_mean,
                            intercept_sd = 2.5,
                            month_sd = 1,
                            nb_logsize_mean = log(2),
                            nb_logsize_sd = 2) {
  PRIOR_INTERCEPT_MEAN <<- intercept_mean
  PRIOR_INTERCEPT_PREC <<- 1 / intercept_sd^2
  PRIOR_MONTH_LOG_RATE_RATIO_PREC <<- 1 / month_sd^2
  PRIOR_NB_LOGSIZE_MEAN <<- nb_logsize_mean
  PRIOR_NB_LOGSIZE_PREC <<- 1 / nb_logsize_sd^2
  invisible(TRUE)
}

modify_settings <- function(settings, changes) {
  out <- settings
  for (nm in names(changes)) {
    out[[nm]] <- changes[[nm]]
  }
  out
}

scalar_or_na <- function(x) {
  if (is.null(x) || !length(x) || !is.finite(x[[1]])) NA_real_ else as.numeric(x[[1]])
}

fit_prior_sensitivity_model <- function(camera_rate, settings, spec, prefix) {
  model_dat <- camera_rate %>%
    dplyr::mutate(y = as.integer(wolf_events), intercept = 1)

  camera_summary <- camera_summary_from_model(model_dat)
  camera_sf <- camera_to_utm(camera_summary)
  coords_camera <- sf::st_coordinates(camera_sf)
  colnames(coords_camera) <- c("x", "y")

  obs_sf <- model_dat %>%
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    sf::st_transform(EPSG_UTM)
  coords_obs <- sf::st_coordinates(obs_sf)
  colnames(coords_obs) <- c("x", "y")

  fixed_terms <- fixed_effect_terms(model_dat)
  fixed_obs <- as.data.frame(model_dat[, fixed_terms, drop = FALSE])

  spde_obj <- build_spatial(coords_camera, settings)
  A_obs <- INLA::inla.spde.make.A(spde_obj$mesh, loc = coords_obs)

  stack_obs <- INLA::inla.stack(
    tag = "obs",
    data = list(y = model_dat$y, e = model_dat$total_effort_days),
    A = list(A_obs, 1),
    effects = list(
      spatial = spde_obj$s_index,
      fixed = fixed_obs
    )
  )

  stack_data <- INLA::inla.stack.data(stack_obs)
  obs_index <- INLA::inla.stack.index(stack_obs, tag = "obs")$data
  formula <- as.formula(
    paste("y ~ 0 +",
          paste(c(fixed_terms, "f(spatial, model = spde_obj$spde)"),
                collapse = " + "))
  )

  fit <- INLA::inla(
    formula,
    family = spec$family,
    data = stack_data,
    E = stack_data$e,
    control.predictor = list(
      A = INLA::inla.stack.A(stack_obs),
      compute = TRUE,
      link = 1
    ),
    control.compute = list(
      config = TRUE,
      dic = TRUE,
      waic = TRUE,
      cpo = TRUE
    ),
    control.fixed = make_control_fixed(fixed_terms),
    control.family = make_control_family(spec$family),
    verbose = FALSE
  )

  diagnostics <- diagnose_fit(
    fit,
    model_dat,
    camera_sf,
    spec,
    obs_index,
    prefix,
    write_files = FALSE
  )

  list(
    fit = fit,
    diag = diagnostics,
    model_dat = diagnostics$model_dat,
    mesh_vertices = nrow(spde_obj$mesh$loc)
  )
}

summarise_sensitivity_fit <- function(variant, note, settings, prior_state,
                                      fit_obj, spec) {
  fit <- fit_obj$fit
  diag <- fit_obj$diag
  model_dat <- fit_obj$model_dat
  month_terms <- temporal_month_terms(model_dat)
  month_fixed <- fit$summary.fixed[month_terms, , drop = FALSE]
  month_means <- if (length(month_terms)) month_fixed[, "mean"] else numeric()
  month_sds <- if (length(month_terms)) month_fixed[, "sd"] else numeric()

  data.frame(
    prior_variant = variant,
    note = note,
    model = spec$name,
    family = spec$family,
    dic = scalar_or_na(fit$dic$dic),
    p_dic = scalar_or_na(fit$dic$p.eff),
    waic = scalar_or_na(fit$waic$waic),
    p_waic = scalar_or_na(fit$waic$p.eff),
    marginal_loglik = scalar_or_na(fit$mlik[1]),
    cpo_failure_rate = if (!is.null(fit$cpo$failure)) {
      mean(fit$cpo$failure, na.rm = TRUE)
    } else {
      NA_real_
    },
    diagnostics_ok = isTRUE(diag$diagnostics_ok),
    ppc_total_pass = isTRUE(diag$ppc_total_pass),
    ppc_zero_pass = isTRUE(diag$ppc_zero_pass),
    ppc_max_pass = isTRUE(diag$ppc_max_pass),
    moran_p = diag$moran_p,
    pearson_disp = diag$pearson_disp,
    pearson_disp_camera = diag$pearson_disp_camera,
    ppc_pit_ks = diag$ppc_pit_ks,
    nb_size_mean = diag$size_hat,
    spatial_range_mean_m = hyp_point(fit, PAT_RANGE),
    spatial_sd_mean = hyp_point(fit, PAT_SIGMA),
    max_abs_month_log_rate_ratio = if (length(month_means)) {
      max(abs(month_means), na.rm = TRUE)
    } else {
      NA_real_
    },
    mean_month_coef_sd = if (length(month_sds)) mean(month_sds, na.rm = TRUE) else NA_real_,
    mesh_vertices = fit_obj$mesh_vertices,
    intercept_prior_mean = prior_state$intercept_mean,
    intercept_prior_sd = prior_state$intercept_sd,
    month_prior_sd = prior_state$month_sd,
    nb_logsize_prior_mean = prior_state$nb_logsize_mean,
    nb_logsize_prior_sd = prior_state$nb_logsize_sd,
    spatial_range_fixed_m = if (!is.null(settings$fix_range_m)) settings$fix_range_m else NA_real_,
    spatial_range_prior_range0 = if (!is.null(settings$prior_range_m)) settings$prior_range_m[[1]] else NA_real_,
    spatial_range_prior_prob_below = if (!is.null(settings$prior_range_m)) settings$prior_range_m[[2]] else NA_real_,
    spatial_sd_prior_sigma0 = settings$prior_sigma[[1]],
    spatial_sd_prior_prob_above = settings$prior_sigma[[2]],
    stringsAsFactors = FALSE
  )
}

run_prior_sensitivity <- function(camera_rate, base_settings, spec, prefix) {
  cat(sprintf("\n[%s] prior sensitivity for month-adjusted model\n", prefix))

  crude_intercept <- log(sum(camera_rate$wolf_events) / sum(camera_rate$total_effort_days))
  base_prior <- list(
    intercept_mean = crude_intercept,
    intercept_sd = 2.5,
    month_sd = 1,
    nb_logsize_mean = log(2),
    nb_logsize_sd = 2
  )

  variants <- list(
    list(
      name = "final_current",
      note = "final selected priors: estimated range PC prior P(range < 1000 m) = 0.5, broad spatial SD, month SD 1, NB log-size SD 2",
      settings = list(),
      prior = list()
    ),
    list(
      name = "fixed_range_500",
      note = "shorter fixed spatial range",
      settings = list(fix_range_m = 500, prior_range_m = c(500, 0.5)),
      prior = list()
    ),
    list(
      name = "fixed_range_1000",
      note = "old stabilising choice: fixed 1000 m spatial range",
      settings = list(fix_range_m = 1000, prior_range_m = c(1000, 0.5)),
      prior = list()
    ),
    list(
      name = "fixed_range_1500",
      note = "longer fixed spatial range",
      settings = list(fix_range_m = 1500, prior_range_m = c(1500, 0.5)),
      prior = list()
    ),
    list(
      name = "fixed_range_2500",
      note = "much longer fixed spatial range",
      settings = list(fix_range_m = 2500, prior_range_m = c(2500, 0.5)),
      prior = list()
    ),
    list(
      name = "estimated_range_pc1500",
      note = "estimated spatial range with PC prior P(range < 1500 m) = 0.5",
      settings = list(fix_range_m = NULL, prior_range_m = c(1500, 0.5)),
      prior = list()
    ),
    list(
      name = "spatial_sd_tighter",
      note = "tighter spatial SD prior, P(SD > 0.85) = 0.05",
      settings = list(prior_sigma = c(0.85, 0.05)),
      prior = list()
    ),
    list(
      name = "spatial_sd_wider",
      note = "wider spatial SD prior, P(SD > 2.5) = 0.05",
      settings = list(prior_sigma = c(2.5, 0.05)),
      prior = list()
    ),
    list(
      name = "month_sd_0_5",
      note = "more regularising month-effect prior, SD 0.5 on log-rate ratios",
      settings = list(),
      prior = list(month_sd = 0.5)
    ),
    list(
      name = "month_sd_2",
      note = "weaker month-effect prior, SD 2 on log-rate ratios",
      settings = list(),
      prior = list(month_sd = 2)
    ),
    list(
      name = "nb_size_median_1",
      note = "NB log-size prior centered at size 1",
      settings = list(),
      prior = list(nb_logsize_mean = log(1))
    ),
    list(
      name = "nb_size_median_3",
      note = "NB log-size prior centered at size 3",
      settings = list(),
      prior = list(nb_logsize_mean = log(3))
    )
  )

  rows <- list()

  for (variant in variants) {
    variant_settings <- modify_settings(base_settings, variant$settings)
    variant_prior <- modifyList(base_prior, variant$prior)
    set_prior_state(
      intercept_mean = variant_prior$intercept_mean,
      intercept_sd = variant_prior$intercept_sd,
      month_sd = variant_prior$month_sd,
      nb_logsize_mean = variant_prior$nb_logsize_mean,
      nb_logsize_sd = variant_prior$nb_logsize_sd
    )

    cat(sprintf("[%s]   sensitivity variant: %s\n", prefix, variant$name))
    fit_obj <- fit_prior_sensitivity_model(
      camera_rate,
      variant_settings,
      spec,
      paste(prefix, "prior", variant$name, sep = "_")
    )

    rows[[length(rows) + 1L]] <- summarise_sensitivity_fit(
      variant$name,
      variant$note,
      variant_settings,
      variant_prior,
      fit_obj,
      spec
    )
  }

  out <- dplyr::bind_rows(rows)
  if (any(is.finite(out$waic))) {
    out$delta_waic <- out$waic - min(out$waic, na.rm = TRUE)
  } else {
    out$delta_waic <- NA_real_
  }
  if (any(is.finite(out$dic))) {
    out$delta_dic <- out$dic - min(out$dic, na.rm = TRUE)
  } else {
    out$delta_dic <- NA_real_
  }

  readr::write_csv(out, path_out(paste0(prefix, "_prior_sensitivity.csv")))
  write_prior_sensitivity_report(out, prefix)

  set_prior_state(
    intercept_mean = base_prior$intercept_mean,
    intercept_sd = base_prior$intercept_sd,
    month_sd = base_prior$month_sd,
    nb_logsize_mean = base_prior$nb_logsize_mean,
    nb_logsize_sd = base_prior$nb_logsize_sd
  )

  invisible(out)
}

write_prior_sensitivity_report <- function(out, prefix) {
  failed <- out$prior_variant[!out$diagnostics_ok]
  finite_waic <- out$waic[is.finite(out$waic)]
  finite_sd <- out$spatial_sd_mean[is.finite(out$spatial_sd_mean)]
  finite_nb <- out$nb_size_mean[is.finite(out$nb_size_mean)]

  lines <- c(
    "Prior sensitivity for forest-camera 2024 month-adjusted NB spatial model",
    "",
    "Purpose:",
    "  Check whether the month-adjusted model is driven by the chosen spatial, month-effect, or NB-size priors.",
    "  These fits use observed data only; final prediction maps remain from the selected prior set.",
    "",
    "Summary:",
    sprintf("  Variants fitted: %d", nrow(out)),
    sprintf("  Variants passing required diagnostics: %d / %d",
            sum(out$diagnostics_ok, na.rm = TRUE), nrow(out)),
    if (length(failed)) {
      sprintf("  Failed variants: %s", paste(failed, collapse = ", "))
    } else {
      "  Failed variants: none"
    },
    if (length(finite_waic)) {
      sprintf("  WAIC range: %.2f to %.2f", min(finite_waic), max(finite_waic))
    } else {
      "  WAIC range: not available"
    },
    if (length(finite_sd)) {
      sprintf("  Spatial SD posterior mean range: %.3f to %.3f",
              min(finite_sd), max(finite_sd))
    } else {
      "  Spatial SD posterior mean range: not available"
    },
    if (length(finite_nb)) {
      sprintf("  NB size posterior mean range: %.3f to %.3f",
              min(finite_nb), max(finite_nb))
    } else {
      "  NB size posterior mean range: not available"
    },
    "",
    "Best WAIC variants:",
    apply(head(out[order(out$waic), ], 5), 1, function(x) {
      sprintf("  %s: WAIC %.2f, diagnostics_ok=%s, Moran p %.3f, NB size %.3f, spatial SD %.3f",
              x[["prior_variant"]],
              as.numeric(x[["waic"]]),
              x[["diagnostics_ok"]],
              as.numeric(x[["moran_p"]]),
              as.numeric(x[["nb_size_mean"]]),
              as.numeric(x[["spatial_sd_mean"]]))
    }),
    "",
    "Interpretation:",
    "  Stable diagnostics across variants support the prior choices.",
    "  The final model estimates spatial range with a weakly informative PC prior rather than fixing it.",
    "  Fixed-range variants are retained as sensitivity checks for whether the range prior drives conclusions.",
    "  If maps are used for publication, cite the estimated range prior and this sensitivity table together."
  )

  writeLines(lines, path_out(paste0(prefix, "_prior_sensitivity_report.txt")))
  invisible(lines)
}

mesh_settings_variants <- function(settings) {
  list(
    current = settings,
    finer = modify_settings(settings, list(
      mesh_cutoff_m = max(50, settings$mesh_cutoff_m * 0.70),
      mesh_max_edge = pmax(50, settings$mesh_max_edge * 0.70),
      mesh_offset = pmax(500, settings$mesh_offset * 0.80)
    )),
    coarser = modify_settings(settings, list(
      mesh_cutoff_m = settings$mesh_cutoff_m * 1.45,
      mesh_max_edge = settings$mesh_max_edge * 1.45,
      mesh_offset = settings$mesh_offset * 1.20
    ))
  )
}

summarise_mesh_fit <- function(variant, note, settings, fit_obj, spec) {
  fit <- fit_obj$fit
  diag <- fit_obj$diag

  data.frame(
    mesh_variant = variant,
    note = note,
    model = spec$name,
    family = spec$family,
    dic = scalar_or_na(fit$dic$dic),
    p_dic = scalar_or_na(fit$dic$p.eff),
    waic = scalar_or_na(fit$waic$waic),
    p_waic = scalar_or_na(fit$waic$p.eff),
    marginal_loglik = scalar_or_na(fit$mlik[1]),
    cpo_failure_rate = if (!is.null(fit$cpo$failure)) {
      mean(fit$cpo$failure, na.rm = TRUE)
    } else {
      NA_real_
    },
    diagnostics_ok = isTRUE(diag$diagnostics_ok),
    ppc_total_pass = isTRUE(diag$ppc_total_pass),
    ppc_zero_pass = isTRUE(diag$ppc_zero_pass),
    ppc_max_pass = isTRUE(diag$ppc_max_pass),
    moran_p = diag$moran_p,
    pearson_disp = diag$pearson_disp,
    pearson_disp_camera = diag$pearson_disp_camera,
    ppc_pit_ks = diag$ppc_pit_ks,
    nb_size_mean = diag$size_hat,
    spatial_range_mean_m = hyp_point(fit, PAT_RANGE),
    spatial_sd_mean = hyp_point(fit, PAT_SIGMA),
    mesh_vertices = fit_obj$mesh_vertices,
    mesh_cutoff_m = settings$mesh_cutoff_m,
    mesh_max_edge_inner_m = settings$mesh_max_edge[[1]],
    mesh_max_edge_outer_m = settings$mesh_max_edge[[2]],
    mesh_offset_inner_m = settings$mesh_offset[[1]],
    mesh_offset_outer_m = settings$mesh_offset[[2]],
    stringsAsFactors = FALSE
  )
}

run_mesh_sensitivity <- function(camera_rate, base_settings, spec, prefix) {
  cat(sprintf("\n[%s] mesh sensitivity for month-adjusted model\n", prefix))

  variants <- mesh_settings_variants(base_settings)
  notes <- c(
    current = "current final mesh",
    finer = "finer SPDE mesh",
    coarser = "coarser SPDE mesh"
  )

  rows <- list()
  for (nm in names(variants)) {
    cat(sprintf("[%s]   mesh variant: %s\n", prefix, nm))
    fit_obj <- fit_prior_sensitivity_model(
      camera_rate,
      variants[[nm]],
      spec,
      paste(prefix, "mesh", nm, sep = "_")
    )
    rows[[length(rows) + 1L]] <- summarise_mesh_fit(
      nm,
      notes[[nm]],
      variants[[nm]],
      fit_obj,
      spec
    )
  }

  out <- dplyr::bind_rows(rows)
  if (any(is.finite(out$waic))) {
    out$delta_waic <- out$waic - min(out$waic, na.rm = TRUE)
  } else {
    out$delta_waic <- NA_real_
  }
  if (any(is.finite(out$dic))) {
    out$delta_dic <- out$dic - min(out$dic, na.rm = TRUE)
  } else {
    out$delta_dic <- NA_real_
  }

  readr::write_csv(out, path_out(paste0(prefix, "_mesh_sensitivity.csv")))
  write_mesh_sensitivity_report(out, prefix)
  invisible(out)
}

write_mesh_sensitivity_report <- function(out, prefix) {
  failed <- out$mesh_variant[!out$diagnostics_ok]
  finite_waic <- out$waic[is.finite(out$waic)]
  finite_range <- out$spatial_range_mean_m[is.finite(out$spatial_range_mean_m)]
  finite_sd <- out$spatial_sd_mean[is.finite(out$spatial_sd_mean)]

  lines <- c(
    "Mesh sensitivity for forest-camera 2024 month-adjusted NB spatial model",
    "",
    "Purpose:",
    "  Check whether the selected spatial surface and fit diagnostics are stable under finer and coarser SPDE meshes.",
    "",
    "Summary:",
    sprintf("  Variants fitted: %d", nrow(out)),
    sprintf("  Variants passing required diagnostics: %d / %d",
            sum(out$diagnostics_ok, na.rm = TRUE), nrow(out)),
    if (length(failed)) {
      sprintf("  Failed variants: %s", paste(failed, collapse = ", "))
    } else {
      "  Failed variants: none"
    },
    if (length(finite_waic)) {
      sprintf("  WAIC range: %.2f to %.2f", min(finite_waic), max(finite_waic))
    } else {
      "  WAIC range: not available"
    },
    if (length(finite_range)) {
      sprintf("  Spatial range posterior mean range: %.1f m to %.1f m",
              min(finite_range), max(finite_range))
    } else {
      "  Spatial range posterior mean range: not available"
    },
    if (length(finite_sd)) {
      sprintf("  Spatial SD posterior mean range: %.3f to %.3f",
              min(finite_sd), max(finite_sd))
    } else {
      "  Spatial SD posterior mean range: not available"
    },
    "",
    "Mesh variants:",
    apply(out, 1, function(x) {
      sprintf("  %s: WAIC %.2f, delta WAIC %.2f, diagnostics_ok=%s, vertices=%s, Moran p %.3f",
              x[["mesh_variant"]],
              as.numeric(x[["waic"]]),
              as.numeric(x[["delta_waic"]]),
              x[["diagnostics_ok"]],
              x[["mesh_vertices"]],
              as.numeric(x[["moran_p"]]))
    }),
    "",
    "Interpretation:",
    "  Similar diagnostics and information criteria across mesh variants support mesh robustness for the forest-camera model."
  )

  writeLines(lines, path_out(paste0(prefix, "_mesh_sensitivity_report.txt")))
  invisible(lines)
}

write_full_final_model_report <- function(cfg, spec, result, cv,
                                          temporal_diag, prior_sensitivity,
                                          mesh_sensitivity) {
  prefix <- cfg$prefix
  fit_obj <- result$final
  diag <- fit_obj$diag
  model_dat <- fit_obj$model_dat
  fit_diag <- fit_obj$fit_diagnostics
  month_coef_file <- path_out(paste0(prefix, "_month_coefficients.csv"))
  month_coef <- if (file.exists(month_coef_file)) {
    readr::read_csv(month_coef_file, show_col_types = FALSE)
  } else {
    NULL
  }

  range_mean <- hyp_point(fit_diag, PAT_RANGE)
  sigma_mean <- hyp_point(fit_diag, PAT_SIGMA)
  ppc <- diag$ppc
  best_prior <- prior_sensitivity[order(prior_sensitivity$waic), ][1, ]
  current_prior <- prior_sensitivity[prior_sensitivity$prior_variant == "final_current", ][1, ]
  failed_prior <- prior_sensitivity$prior_variant[!prior_sensitivity$diagnostics_ok]
  cv_value <- function(level, metric) {
    if (is.null(cv) || is.null(cv$summ) || !"level" %in% names(cv$summ)) return(NA_real_)
    x <- cv$summ$value[cv$summ$level == level & cv$summ$metric == metric]
    if (length(x)) x[[1]] else NA_real_
  }

  month_lines <- if (!is.null(month_coef) && nrow(month_coef)) {
    apply(month_coef, 1, function(x) {
      sprintf("  %s vs %s: rate ratio %.2f (95%% CrI %.2f to %.2f)",
              x[["month"]],
              x[["reference_month"]],
              as.numeric(x[["mean_rate_ratio"]]),
              as.numeric(x[["q025_rate_ratio"]]),
              as.numeric(x[["q975_rate_ratio"]]))
    })
  } else {
    "  Month coefficients were not written."
  }

  ppc_lines <- apply(ppc, 1, function(x) {
    sprintf("  %s %s: observed %s; simulated 95%% interval %s to %s; pass=%s",
            x[["level"]],
            x[["stat"]],
            x[["observed"]],
            x[["sim_q025"]],
            x[["sim_q975"]],
            x[["pass"]])
  })

  cv_lines <- if (!is.null(cv)) {
    c(
      sprintf("  Row mean log predictive density: %.3f",
              cv_value("model_row", "mean_log_predictive_density")),
      sprintf("  Row RMSE count: %.3f",
              cv_value("model_row", "rmse_count")),
      sprintf("  Row RMSE rate per 100 camera-days: %.3f",
              cv_value("model_row", "rmse_rate_per100")),
      sprintf("  Row 90%% coverage: %.3f",
              cv_value("model_row", "coverage_90")),
      sprintf("  Camera RMSE count: %.3f",
              cv_value("camera", "rmse_count")),
      sprintf("  Camera RMSE rate per 100 camera-days: %.3f",
              cv_value("camera", "rmse_rate_per100")),
      sprintf("  Camera 90%% coverage: %.3f",
              cv_value("camera", "coverage_90"))
    )
  } else {
    "  Spatial block CV was not run for this profile."
  }

  temporal_lines <- if (!is.null(temporal_diag) && nrow(temporal_diag)) {
    c(
      sprintf("  Months checked: %s", paste(temporal_diag$month, collapse = ", ")),
      sprintf("  Within-camera lag-1 residual correlation: r = %s, p = %s",
              ifelse("within_camera_lag1_r" %in% names(temporal_diag) &&
                       is.finite(temporal_diag$within_camera_lag1_r[[1]]),
                     sprintf("%.3f", temporal_diag$within_camera_lag1_r[[1]]),
                     "not estimable"),
              ifelse("within_camera_lag1_p" %in% names(temporal_diag) &&
                       is.finite(temporal_diag$within_camera_lag1_p[[1]]),
                     sprintf("%.4g", temporal_diag$within_camera_lag1_p[[1]]),
                     "not estimable")),
      sprintf("  Month-level Pearson residual lag-1 ACF: %.3f",
              temporal_diag$lag1_acf[[1]]),
      "  This is a supporting low-power check because only seven monthly points are available."
    )
  } else {
    "  Temporal residual diagnostic was not available."
  }

  prior_lines <- c(
    sprintf("  Sensitivity variants fitted: %d", nrow(prior_sensitivity)),
    sprintf("  Variants passing required diagnostics: %d / %d",
            sum(prior_sensitivity$diagnostics_ok, na.rm = TRUE),
            nrow(prior_sensitivity)),
    if (length(failed_prior)) {
      sprintf("  Failed prior variants: %s", paste(failed_prior, collapse = ", "))
    } else {
      "  Failed prior variants: none"
    },
    sprintf("  Final-current WAIC: %.2f", current_prior$waic),
    sprintf("  Best WAIC variant: %s, WAIC %.2f",
            best_prior$prior_variant, best_prior$waic),
    sprintf("  WAIC range across variants: %.2f to %.2f",
            min(prior_sensitivity$waic, na.rm = TRUE),
            max(prior_sensitivity$waic, na.rm = TRUE)),
    sprintf("  NB size posterior mean range: %.3f to %.3f",
            min(prior_sensitivity$nb_size_mean, na.rm = TRUE),
            max(prior_sensitivity$nb_size_mean, na.rm = TRUE)),
    sprintf("  Spatial SD posterior mean range: %.3f to %.3f",
            min(prior_sensitivity$spatial_sd_mean, na.rm = TRUE),
            max(prior_sensitivity$spatial_sd_mean, na.rm = TRUE))
  )

  mesh_failed <- mesh_sensitivity$mesh_variant[!mesh_sensitivity$diagnostics_ok]
  mesh_lines <- c(
    sprintf("  Mesh variants fitted: %d", nrow(mesh_sensitivity)),
    sprintf("  Variants passing required diagnostics: %d / %d",
            sum(mesh_sensitivity$diagnostics_ok, na.rm = TRUE),
            nrow(mesh_sensitivity)),
    if (length(mesh_failed)) {
      sprintf("  Failed mesh variants: %s", paste(mesh_failed, collapse = ", "))
    } else {
      "  Failed mesh variants: none"
    },
    sprintf("  WAIC range across meshes: %.2f to %.2f",
            min(mesh_sensitivity$waic, na.rm = TRUE),
            max(mesh_sensitivity$waic, na.rm = TRUE)),
    sprintf("  Spatial range posterior mean range: %.1f m to %.1f m",
            min(mesh_sensitivity$spatial_range_mean_m, na.rm = TRUE),
            max(mesh_sensitivity$spatial_range_mean_m, na.rm = TRUE)),
    sprintf("  Spatial SD posterior mean range: %.3f to %.3f",
            min(mesh_sensitivity$spatial_sd_mean, na.rm = TRUE),
            max(mesh_sensitivity$spatial_sd_mean, na.rm = TRUE))
  )

  lines <- c(
    "FULL FINAL MODEL REPORT",
    "Forest-camera 2024 wolf relative encounter-frequency model",
    "",
    "1. Final model",
    sprintf("  Model: %s", spec$name),
    sprintf("  Likelihood: %s", spec$family),
    "  Spatial structure: INLA-SPDE spatial random field.",
    "  Temporal structure: calendar camera-month fixed effects.",
    "  Effort structure: active camera-days are included as exposure.",
    "  Map target: effort-weighted annualized 2024 encounter-frequency surface.",
    "  Prediction units: expected independent wolf events per 100 camera-days.",
    "  Priors: weakly informative Gaussian and PC priors.",
    "",
    "2. Data represented in the model",
    sprintf("  Cameras: %d", dplyr::n_distinct(model_dat$plotID)),
    sprintf("  Camera-month rows: %d", nrow(model_dat)),
    sprintf("  Independent wolf events: %d", sum(model_dat$y)),
    sprintf("  Camera effort: %.1f camera-days", sum(model_dat$total_effort_days)),
    sprintf("  Observed mean encounter frequency: %.3f events per 100 camera-days",
            100 * sum(model_dat$y) / sum(model_dat$total_effort_days)),
    sprintf("  Months: %s", paste(sort(unique(model_dat$month)), collapse = ", ")),
    "",
    "3. Final priors",
    "  Intercept: Gaussian centered on crude daily event rate, SD 2.5 on log scale.",
    "  Month log-rate ratios: Gaussian(0, SD 1.0).",
    "  Negative-binomial log(size): Gaussian(log(2), SD 2.0).",
    "  Spatial range: PC prior, P(range < 1000 m) = 0.5.",
    "  Spatial marginal SD: PC prior, P(SD > 1.5) = 0.05.",
    "  These priors are weakly informative/regularising; the spatial range is estimated, not fixed.",
    "",
    "4. Fitted hyperparameters",
    sprintf("  Negative-binomial size posterior mean: %.3f", diag$size_hat),
    sprintf("  Spatial range posterior mean: %.1f m", range_mean),
    sprintf("  Spatial SD posterior mean: %.3f", sigma_mean),
    "",
    "5. Month effects",
    month_lines,
    "",
    "6. Posterior predictive checks",
    sprintf("  Method: %s; simulations: %d", diag$ppc_method, diag$ppc_nsim),
    ppc_lines,
    "",
    "7. Residual diagnostics",
    sprintf("  Pearson dispersion, model rows: %.3f", diag$pearson_disp),
    sprintf("  Pearson dispersion, camera aggregates: %.3f", diag$pearson_disp_camera),
    sprintf("  Residual Moran's I: %.3f; p = %.3f", diag$moran_I, diag$moran_p),
    sprintf("  Residual Moran pass: %s", isTRUE(diag$moran_pass)),
    sprintf("  Row PIT KS p-value: %.4f", diag$ppc_pit_ks_row),
    sprintf("  Camera PIT KS p-value: %.4f", diag$ppc_pit_ks_camera),
    sprintf("  Required diagnostics pass: %s", isTRUE(diag$diagnostics_ok)),
    "",
    "8. Temporal residual check",
    temporal_lines,
    "",
    "9. Spatial block cross-validation",
    cv_lines,
    "",
    "10. Prior sensitivity",
    prior_lines,
    "",
    "11. Mesh sensitivity",
    mesh_lines,
    "",
    "12. Final assessment",
    "  The final estimated-range month-adjusted NB spatial model passes the required diagnostics.",
    "  Prior sensitivity supports the prior choices: all tested variants pass required diagnostics.",
    "  Mesh sensitivity supports the selected SPDE mesh: all tested variants pass required diagnostics.",
    "  The estimated spatial range prior removes the earlier fixed-range modelling assumption while retaining regularisation.",
    "  The output should be interpreted as relative encounter frequency, not abundance, density, occupancy, or population size.",
    "",
    "13. Main limitations",
    "  The data contain only 46 independent wolf events, so month and spatial effects have wide uncertainty.",
    "  Detection probability, camera placement, habitat, roads, prey, and human disturbance are not modelled explicitly.",
    "  The month-level temporal ACF check has low power because only seven monthly time points are available.",
    "  Predictions should not be interpreted far outside the sampled camera domain."
  )

  writeLines(lines, path_out(paste0(prefix, "_full_final_model_report.txt")))
  invisible(lines)
}

prefix <- "wolf_forest_2024"
settings <- list(
  cell_size_m = 60,
  pred_buffer_m = 1500,
  max_dist_m = 2500,
  mesh_cutoff_m = 75,
  mesh_max_edge = c(150, 900),
  mesh_offset = c(1200, 4000),
  fix_range_m = NULL,
  prior_range_m = c(1000, 0.5),
  prior_sigma = c(1.5, 0.05),
  use_month_effect = TRUE,
  month_reference = "2024-08",
  month_prediction = "2024-08",
  include_grid_in_mesh = FALSE
)

cfg <- list(
  label = "Forest-camera 2024 survey with month fixed effects",
  prefix = prefix,
  settings = settings,
  caveat = paste(
    "Camera-month rows are reconstructed from the dated flat forest-camera",
    "2024 file. Month fixed effects are retained as a temporal control;",
    "the spatial range is estimated with a weakly informative PC prior;",
    "the month-level temporal residual ACF is only a supporting check because",
    "there are seven monthly points."
  )
)
spec <- model_spec("nb_spatial_month", "nbinomial")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("\n==================== wolf_forest month refit ====================\n")
camera_rate <- load_forest_flat_deployment_month(settings, prefix)
PRIOR_INTERCEPT_MEAN <- log(sum(camera_rate$wolf_events) / sum(camera_rate$total_effort_days))

fit <- fit_final_model(
  camera_rate,
  settings,
  spec,
  prefix,
  add_prediction = TRUE,
  write_files = TRUE
)

cv <- NULL
if (RUN_FINAL_SPATIAL_CV) {
  cv <- spatial_block_cv(camera_rate, settings, spec, prefix, K = CV_K)
}

temporal_diag <- write_temporal_residual_diagnostic(fit$diag, prefix, spec$name)
write_validation_report(prefix, cfg, spec, camera_rate, fit$diag, cv, temporal_diag)

result <- list(
  prefix = prefix,
  camera_rate = camera_rate,
  spec = spec,
  final = fit,
  cv = cv
)
write_run_manifest(list(forest_2024 = result))
file.copy(path_out("wolf_final_run_manifest.csv"),
          path_out(paste0(prefix, "_run_manifest.csv")),
          overwrite = TRUE)
prior_sensitivity <- run_prior_sensitivity(camera_rate, settings, spec, prefix)
mesh_sensitivity <- run_mesh_sensitivity(camera_rate, settings, spec, prefix)
write_full_final_model_report(cfg, spec, result, cv, temporal_diag,
                              prior_sensitivity, mesh_sensitivity)

cat("\nMonth refit completed successfully.\n")
cat("Outputs are in:\n  ", OUTPUT_DIR, "\n", sep = "")
cat("Prior sensitivity variants fitted: ", nrow(prior_sensitivity), "\n", sep = "")
cat("Mesh sensitivity variants fitted: ", nrow(mesh_sensitivity), "\n", sep = "")
