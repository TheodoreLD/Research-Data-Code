###############################################################################
# Wolf relative encounter frequency: helper functions for the small 2024 model
# -----------------------------------------------------------------------------
# Purpose
#   This file is retained as a helper dependency for
#   scripts/wolf_small_2024_month_refit.R. The small-model wrapper reads and
#   evaluates the shared function definitions above the main execution boundary.
#   Do not use this file as the final 2024 workflow entry point; use the two
#   scripts documented in the project README instead.
#
# Response variable
#   The final project scripts model independent wolf eventIDs with camera effort
#   as the exposure term. This helper file contains older shared functions used
#   by the small/local 2024 wrapper; see README.md for the final 2024 model
#   definitions.
#
# Interpretation
#   The output is a relative encounter-frequency index. It is not abundance,
#   density, occupancy, or population size.
#
# Final model note
#   This helper file is not the final model-selection statement. The final 2024
#   models are documented in README.md and docs/final-model-details.md.
#
# Main outputs
#   * GeoTIFF prediction rasters: mean, SD, CV, and exceedance probability
#   * PNG maps for mean encounter frequency, uncertainty, and exceedance
#   * posterior predictive checks for total events, zero fraction, and max count
#   * residual diagnostics: Pearson dispersion, Moran's I, PIT, semivariogram
#   * spatial block cross-validation summaries
#   * fixed-effect summaries for temporal month effects
#   * per-survey validation reports
###############################################################################


## 01. User Settings: File Paths, Surveys, And Runtime Profile -----------------

# The script is portable. It first looks for input files next to this script or
# in a data/ subfolder. Paths can be overridden with WOLF_* environment variables.

input_files_required <- c(
  "theodata_1.0(1).csv",
  "deployments_2024.csv",
  "observations_2024.csv",
  "deployments_2023.csv",
  "observations_2023.csv"
)

script_file <- tryCatch({
  ofile <- sys.frames()[[1]]$ofile
  if (is.null(ofile)) {
    NA_character_
  } else {
    normalizePath(ofile, winslash = "/", mustWork = FALSE)
  }
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
  unset = file.path(PROJECT_DIR, "outputs", "wolf_final_maps_diagnostics")
)
OUTPUT_DIR <- normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE)

# Coordinate reference system used for spatial modelling and maps.
EPSG_UTM <- 32634L

# Species labels observed in the input files.
WOLF_NAMES <- c("Canis_lupus", "Canis lupus")

# Set TRUE only if you explicitly want this script to install missing packages.
INSTALL_MISSING_PACKAGES <- FALSE

# Surveys to run. Each survey below must have a final_model entry.
RUN_SURVEYS <- c("2023", "2024", "small")

# Runtime profile:
#   quick    : fast code testing; skips spatial block CV
#   balanced : default; still runs PPC and spatial block CV, with fewer draws
#   final    : original heavier settings for final publication reruns
RUN_PROFILE <- tolower(Sys.getenv("WOLF_RUN_PROFILE", unset = "balanced"))
if (!nzchar(RUN_PROFILE)) RUN_PROFILE <- "balanced"
if (!RUN_PROFILE %in% c("quick", "balanced", "final")) {
  stop("WOLF_RUN_PROFILE must be one of: quick, balanced, final.")
}

PPC_NSIM <- switch(RUN_PROFILE,
                   quick = 400L,
                   balanced = 1000L,
                   final = 2000L)

JOINT_PPC_NSIM <- switch(RUN_PROFILE,
                         quick = 100L,
                         balanced = 250L,
                         final = 500L)

RUN_FINAL_SPATIAL_CV <- RUN_PROFILE != "quick"
CV_K <- switch(RUN_PROFILE,
               quick = 3L,
               balanced = 4L,
               final = 5L)
CV_NSIM <- switch(RUN_PROFILE,
                  quick = 200L,
                  balanced = 400L,
                  final = 800L)

# Required diagnostic threshold for residual spatial autocorrelation.
MORAN_ALPHA <- 0.05

# Prediction domain:
#   "hull"  : buffered convex hull around cameras
#   "disks" : buffered disks around cameras, trimmed by max_dist_m
PRED_DOMAIN <- "hull"
MAP_EXCEEDANCE <- TRUE
EXCEED_MULT <- 1.5

# Priors. Field priors are survey-specific in each settings list below.
PRIOR_INTERCEPT_PREC <- 0.01
PRIOR_MONTH_LOG_RATE_RATIO_PREC <- 1
PRIOR_ZI_LOGIT_MEAN <- -1
PRIOR_ZI_LOGIT_PREC <- 0.2
PRIOR_NB_SIZE_LOGGAMMA <- c(1, 0.01)

# Small helper used only while the survey list is created. The canonical
# model_spec() function is defined later after family-name standardisation.
.model_spec <- function(name, family = "poisson") {
  list(name = name, family = family)
}


## 02. Survey Definitions: Final Pinned Models And Mesh Settings --------------

surveys <- list(
  small = list(
    label = "Small/local survey",
    type = "flat",
    prefix = "wolf_small",
    file = "theodata_1.0(1).csv",
    final_model = .model_spec("pois_field", "poisson"),
    caveat = paste(
      "Sensitivity checks found no evidence that residual temporal",
      "autocorrelation remained after using the simple spatial model.",
      "Month-adjusted and zero-inflated alternatives passed diagnostics but did",
      "not materially improve the small-survey model, so the simpler Poisson",
      "spatial model is retained for parsimony."
    ),
    settings = list(
      cell_size_m = 60,
      pred_buffer_m = 1500,
      max_dist_m = 2500,
      mesh_cutoff_m = 150,
      mesh_max_edge = c(300, 1500),
      mesh_offset = c(1800, 6000),
      fix_range_m = 1500,
      prior_range_m = c(1500, 0.5),
      prior_sigma = c(0.85, 0.05),
      include_grid_in_mesh = FALSE
    )
  ),

  `2024` = list(
    label = "Large 2024 survey",
    type = "camtrap",
    prefix = "wolf_2024",
    deployments = "deployments_2024.csv",
    observations = "observations_2024.csv",
    final_model = .model_spec("pois_field_month", "poisson"),
    caveat = paste(
      "Deployment timing is strongly correlated with latitude (Spearman rho +0.53,",
      "p < 0.001). This legacy helper model adjusts for month and maps",
      "the spatial surface at the reference month 2024-09."
    ),
    settings = list(
      cell_size_m = 150,
      pred_buffer_m = 1500,
      max_dist_m = 2500,
      mesh_cutoff_m = 350,
      mesh_max_edge = c(700, 5000),
      mesh_offset = c(5000, 15000),
      fix_range_m = NULL,
      prior_range_m = c(5000, 0.5),
      prior_sigma = c(0.85, 0.05),
      use_month_effect = TRUE,
      month_reference = "2024-09",
      month_prediction = "2024-09",
      include_grid_in_mesh = FALSE
    )
  ),

  `2023` = list(
    label = "Large 2023 survey",
    type = "camtrap",
    prefix = "wolf_2023",
    deployments = "deployments_2023.csv",
    observations = "observations_2023.csv",
    final_model = .model_spec("zinb_field_month", "zeroinflatednbinomial1"),
    caveat = paste(
      "Deployment timing correlates with latitude (Spearman rho +0.35, p 0.007)",
      "and encounter rates rise through the season. This final zero-inflated",
      "negative-binomial helper model adjusts for month and maps the",
      "spatial surface at the reference month 2023-08."
    ),
    settings = list(
      cell_size_m = 150,
      pred_buffer_m = 1500,
      max_dist_m = 2500,
      mesh_cutoff_m = 350,
      mesh_max_edge = c(700, 5000),
      mesh_offset = c(5000, 15000),
      fix_range_m = NULL,
      prior_range_m = c(5000, 0.5),
      prior_sigma = c(0.85, 0.05),
      use_month_effect = TRUE,
      month_reference = "2023-08",
      month_prediction = "2023-08",
      include_grid_in_mesh = FALSE
    )
  )
)


## 03. Package Setup: Load Required Libraries ---------------------------------

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Project directory: ", PROJECT_DIR, "\n", sep = "")
cat("Data directory:    ", DATA_DIR, "\n", sep = "")
cat("Output directory:  ", OUTPUT_DIR, "\n", sep = "")
cat(sprintf(
  "Run profile:       %s | PPC_NSIM=%d | JOINT_PPC_NSIM=%d | spatial_CV=%s | CV_K=%d | CV_NSIM=%d\n",
  RUN_PROFILE, PPC_NSIM, JOINT_PPC_NSIM, RUN_FINAL_SPATIAL_CV, CV_K, CV_NSIM
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
set.seed(1)
try(INLA::inla.setOption(fmesher.evolution.warn = FALSE), silent = TRUE)


## 04. General Helpers: Paths, Validation, Time Parsing ------------------------

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


## 05. Model-Family Helpers: Poisson, NB, ZIP, And ZINB -----------------------

fit_family <- function(family) {
  key <- gsub("[^a-z0-9]", "", tolower(family))
  switch(key,
         poisson = "poisson",
         nbinomial = ,
         negativebinomial = ,
         nb = "nbinomial",
         zeroinflatedpoisson1 = ,
         zip = ,
         zip1 = "zeroinflatedpoisson1",
         zeroinflatednbinomial1 = ,
         zinb = ,
         zinb1 = "zeroinflatednbinomial1",
         family)
}

is_zi <- function(family) {
  grepl("zeroinflated", family)
}

is_nb <- function(family) {
  grepl("nbinomial", family)
}

fam_nb_size <- function(size) {
  if (is.finite(size) && size > 0) size else 1e6
}

fam_mean <- function(mu, pi, family, size = NA_real_) {
  p <- if (is_zi(family) && is.finite(pi)) pi else 0
  (1 - p) * mu
}

fam_var <- function(mu, pi, family, size = NA_real_) {
  p <- if (is_zi(family) && is.finite(pi)) pi else 0
  base_var <- if (is_nb(family)) {
    s <- fam_nb_size(size)
    mu + mu^2 / s
  } else {
    mu
  }
  (1 - p) * base_var + p * (1 - p) * mu^2
}

fam_logpmf <- function(y, mu, pi, family, size = NA_real_) {
  p <- if (is_zi(family) && is.finite(pi)) pi else 0
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

PAT_ZPROB <- "zero.*prob|probability.*zero|zero-probability"
PAT_NB_SIZE <- "size.*nbinomial|size for"
PAT_RANGE <- "range.*spatial|range for spatial|spatial.*range"
PAT_SIGMA <- "stdev.*spatial|stdev for spatial|standard deviation.*spatial|sigma.*spatial"

nb_size_point <- function(fit) {
  x <- hyp_point(fit, PAT_NB_SIZE)
  if (is.finite(x) && x > 0) x else NA_real_
}

use_joint_ppc <- function(family) {
  is_nb(family)
}

model_spec <- function(name, family = "poisson") {
  list(name = name, family = fit_family(family))
}

finalise_spec <- function(spec) {
  model_spec(spec$name, spec$family)
}


## 06. Prior Helpers: Fixed Effects And Likelihood Hyperparameters ------------

month_term_name <- function(month) {
  paste0("month_", gsub("[^A-Za-z0-9]", "_", month))
}

month_from_term <- function(term) {
  gsub("_", "-", sub("^month_", "", term), fixed = TRUE)
}

temporal_month_terms <- function(data) {
  grep("^month_[0-9]{4}_[0-9]{2}$", names(data), value = TRUE)
}

fixed_effect_terms <- function(data) {
  c("intercept", temporal_month_terms(data))
}

make_control_fixed <- function(fixed_terms = "intercept") {
  fixed_terms <- unique(fixed_terms)

  mean_prior <- as.list(setNames(rep(0, length(fixed_terms)), fixed_terms))
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

make_control_family <- function(family) {
  hyper <- list()

  if (is_nb(family)) {
    hyper$size <- list(
      prior = "loggamma",
      param = PRIOR_NB_SIZE_LOGGAMMA
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


## 07. Prior-Posterior Plots: How Strongly The Data Update The Priors ----------

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
  # INLA loggamma prior: log(size) ~ Gamma(shape, rate), so Jacobian 1/x.
  dgamma(log(x), shape = PRIOR_NB_SIZE_LOGGAMMA[1], rate = PRIOR_NB_SIZE_LOGGAMMA[2]) / x
}

nb_size_prior_quantile <- function(p) {
  # Returns quantiles on the size scale (exp of Gamma quantiles on log-size scale).
  exp(qgamma(p, shape = PRIOR_NB_SIZE_LOGGAMMA[1], rate = PRIOR_NB_SIZE_LOGGAMMA[2]))
}

plot_prior_posterior_density <- function(prefix, parameter, prior_df,
                                         posterior_marginal, file_suffix,
                                         x_label, log_x = FALSE,
                                         reference_x = NULL,
                                         reference_label = NULL) {
  if (is.null(posterior_marginal) || nrow(posterior_marginal) < 2) {
    writeLines(
      paste("Posterior marginal not available for", parameter),
      path_out(paste0(prefix, "_prior_posterior_", file_suffix, "_NOTE.txt"))
    )
    return(invisible(NULL))
  }

  posterior_df <- data.frame(
    value = posterior_marginal[, 1],
    density = posterior_marginal[, 2],
    source = "posterior"
  ) %>%
    filter(is.finite(value), is.finite(density), density >= 0)

  prior_df <- prior_df %>%
    mutate(source = "prior") %>%
    filter(is.finite(value), is.finite(density), density >= 0)

  if (log_x) {
    posterior_df <- posterior_df %>% filter(value > 0)
    prior_df <- prior_df %>% filter(value > 0)
  }

  plot_df <- bind_rows(prior_df, posterior_df)
  if (!nrow(plot_df)) {
    writeLines(
      paste("No finite prior/posterior density values for", parameter),
      path_out(paste0(prefix, "_prior_posterior_", file_suffix, "_NOTE.txt"))
    )
    return(invisible(NULL))
  }

  p <- ggplot(plot_df, aes(value, density, colour = source, linetype = source)) +
    geom_line(linewidth = 0.85) +
    labs(
      title = paste0("Prior vs posterior: ", parameter),
      subtitle = prefix,
      x = x_label,
      y = "density",
      colour = NULL,
      linetype = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")

  if (log_x) {
    p <- p + scale_x_log10(labels = label_number())
  }

  if (!is.null(reference_x) && is.finite(reference_x)) {
    p <- p + geom_vline(xintercept = reference_x, linetype = 2, colour = "grey35")
    if (!is.null(reference_label)) {
      p <- p + labs(caption = reference_label)
    }
  }

  ggsave(
    path_out(paste0(prefix, "_prior_posterior_", file_suffix, ".png")),
    p,
    width = 6.8,
    height = 4.8,
    dpi = 250
  )

  readr::write_csv(
    plot_df,
    path_out(paste0(prefix, "_prior_posterior_", file_suffix, ".csv"))
  )

  invisible(plot_df)
}

write_prior_posterior_plots <- function(fit, settings, spec, prefix) {
  # Intercept: broad Gaussian prior on the log-rate intercept.
  intercept_post <- fit$marginals.fixed[["intercept"]]
  intercept_sd_prior <- 1 / sqrt(PRIOR_INTERCEPT_PREC)
  intercept_grid <- seq(
    qnorm(0.001, 0, intercept_sd_prior),
    qnorm(0.999, 0, intercept_sd_prior),
    length.out = 700
  )
  intercept_prior <- data.frame(
    value = intercept_grid,
    density = dnorm(intercept_grid, mean = 0, sd = intercept_sd_prior)
  )
  plot_prior_posterior_density(
    prefix = prefix,
    parameter = "intercept",
    prior_df = intercept_prior,
    posterior_marginal = intercept_post,
    file_suffix = "intercept",
    x_label = "log-rate intercept"
  )

  # Month fixed effects: Gaussian priors on log-rate ratios relative to the
  # configured reference month.
  month_terms <- grep("^month_[0-9]{4}_[0-9]{2}$",
                      names(fit$marginals.fixed),
                      value = TRUE)
  if (length(month_terms)) {
    month_sd_prior <- 1 / sqrt(PRIOR_MONTH_LOG_RATE_RATIO_PREC)
    month_grid <- seq(
      qnorm(0.001, 0, month_sd_prior),
      qnorm(0.999, 0, month_sd_prior),
      length.out = 700
    )
    month_prior <- data.frame(
      value = month_grid,
      density = dnorm(month_grid, mean = 0, sd = month_sd_prior)
    )

    for (term in month_terms) {
      plot_prior_posterior_density(
        prefix = prefix,
        parameter = paste0(
          "month log-rate ratio: ",
          month_from_term(term),
          " vs ",
          settings$month_reference
        ),
        prior_df = month_prior,
        posterior_marginal = fit$marginals.fixed[[term]],
        file_suffix = paste0("fixed_", term),
        x_label = "log-rate ratio"
      )
    }
  }

  # Spatial range: only plotted when range is estimated. For the small survey it
  # is fixed by design, so no posterior range distribution exists.
  if (!is.null(settings$fix_range_m)) {
    writeLines(
      c(
        "Spatial range was fixed in this model.",
        sprintf("fix_range_m = %s", settings$fix_range_m),
        "No prior-posterior overlay is available for a fixed parameter."
      ),
      path_out(paste0(prefix, "_prior_posterior_range_NOTE.txt"))
    )
  } else {
    range0 <- settings$prior_range_m[1]
    prob_range <- settings$prior_range_m[2]
    range_grid <- exp(seq(
      log(pc_range_quantile(0.001, range0, prob_range)),
      log(pc_range_quantile(0.995, range0, prob_range)),
      length.out = 700
    ))
    range_prior <- data.frame(
      value = range_grid,
      density = pc_range_density(range_grid, range0, prob_range)
    )
    plot_prior_posterior_density(
      prefix = prefix,
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

  # Spatial marginal SD: PC prior P(SD > sigma0) = alpha.
  sigma0 <- settings$prior_sigma[1]
  prob_sigma <- settings$prior_sigma[2]
  sigma_grid <- seq(
    0,
    pc_sigma_quantile(0.995, sigma0, prob_sigma),
    length.out = 700
  )
  sigma_prior <- data.frame(
    value = sigma_grid,
    density = pc_sigma_density(sigma_grid, sigma0, prob_sigma)
  )
  plot_prior_posterior_density(
    prefix = prefix,
    parameter = "spatial marginal SD",
    prior_df = sigma_prior,
    posterior_marginal = hyp_marg(fit, PAT_SIGMA),
    file_suffix = "spatial_sd",
    x_label = "spatial marginal SD",
    reference_x = sigma0,
    reference_label = sprintf("Prior statement: P(SD > %.2f) = %.2f",
                              sigma0, prob_sigma)
  )

  # Zero-inflation probability: relevant for ZIP and ZINB models.
  if (is_zi(spec$family)) {
    prob_grid <- seq(0.001, 0.999, length.out = 700)
    zip_prior <- data.frame(
      value = prob_grid,
      density = zip_prob_prior_density(prob_grid)
    )
    plot_prior_posterior_density(
      prefix = prefix,
      parameter = "zero-inflation probability",
      prior_df = zip_prior,
      posterior_marginal = hyp_marg(fit, PAT_ZPROB),
      file_suffix = "zero_inflation_probability",
      x_label = "zero-inflation probability",
      reference_x = plogis(PRIOR_ZI_LOGIT_MEAN),
      reference_label = sprintf("Prior center on probability scale: %.2f",
                                plogis(PRIOR_ZI_LOGIT_MEAN))
    )
  }

  # Negative-binomial size: relevant for NB and ZINB models.
  if (is_nb(spec$family)) {
    size_grid <- seq(
      nb_size_prior_quantile(0.001),
      nb_size_prior_quantile(0.995),
      length.out = 700
    )
    size_prior <- data.frame(
      value = size_grid,
      density = nb_size_prior_density(size_grid)
    )
    prior_mode_size <- exp(
      (PRIOR_NB_SIZE_LOGGAMMA[1] - 1) / PRIOR_NB_SIZE_LOGGAMMA[2]
    )
    plot_prior_posterior_density(
      prefix = prefix,
      parameter = "negative-binomial size",
      prior_df = size_prior,
      posterior_marginal = hyp_marg(fit, PAT_NB_SIZE),
      file_suffix = "nb_size",
      x_label = "negative-binomial size",
      log_x = TRUE,
      reference_x = prior_mode_size,
      reference_label = sprintf(
        "Prior mode on size scale: %.1f (log(size) ~ Gamma(%.2f, %.2f))",
        prior_mode_size,
        PRIOR_NB_SIZE_LOGGAMMA[1],
        PRIOR_NB_SIZE_LOGGAMMA[2]
      )
    )
  }

  invisible(TRUE)
}


## 08. Data Preparation: Collapse Events To Camera-Level Counts ---------------

summarise_camera_rate <- function(deployments, wolf_events, prefix) {
  camera_rate <- deployments %>%
    left_join(wolf_events, by = "deploymentID") %>%
    mutate(wolf_events = tidyr::replace_na(wolf_events, 0L)) %>%
    group_by(plotID) %>%
    summarise(
      longitude = mean(longitude, na.rm = TRUE),
      latitude = mean(latitude, na.rm = TRUE),
      total_effort_days = sum(deploymentEffort, na.rm = TRUE),
      wolf_events = sum(wolf_events, na.rm = TRUE),
      n_deployments = n_distinct(deploymentID),
      wolf_events_per_100_days = 100 * wolf_events / total_effort_days,
      .groups = "drop"
    ) %>%
    filter(is.finite(longitude),
           is.finite(latitude),
           is.finite(total_effort_days),
           total_effort_days > 0) %>%
    arrange(desc(wolf_events_per_100_days)) %>%
    mutate(model_row_type = "camera")

  if (!nrow(camera_rate)) {
    stop("[", prefix, "] no valid camera-level records after filtering.")
  }

  readr::write_csv(camera_rate, path_out(paste0(prefix, "_camera_effort_rates.csv")))

  cat(sprintf(
    "[%s] cameras %d | positive %d | events %d | effort %.1f camera-days | observed %.2f /100\n",
    prefix,
    nrow(camera_rate),
    sum(camera_rate$wolf_events > 0),
    sum(camera_rate$wolf_events),
    sum(camera_rate$total_effort_days),
    100 * sum(camera_rate$wolf_events) / sum(camera_rate$total_effort_days)
  ))

  camera_rate
}

camera_summary_from_model <- function(model_dat) {
  has_deployment_id <- "deploymentID" %in% names(model_dat)

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
      } else {
        sum(n_deployments, na.rm = TRUE)
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

month_reference_from_settings <- function(months, settings, prefix) {
  ref <- settings$month_reference
  if (is.null(ref) || !nzchar(ref)) ref <- months[[1]]
  if (!ref %in% months) {
    stop("[", prefix, "] month_reference '", ref,
         "' is not present in deployment months: ",
         paste(months, collapse = ", "))
  }
  ref
}

month_prediction_from_settings <- function(months, settings, prefix) {
  pred <- settings$month_prediction
  if (is.null(pred) || !nzchar(pred)) {
    pred <- month_reference_from_settings(months, settings, prefix)
  }
  if (!pred %in% months) {
    stop("[", prefix, "] month_prediction '", pred,
         "' is not present in deployment months: ",
         paste(months, collapse = ", "))
  }
  pred
}

add_month_design <- function(model_dat, settings, prefix) {
  months <- sort(unique(model_dat$month))
  if (length(months) < 2) {
    stop("[", prefix, "] month effect requires at least two deployment months.")
  }

  reference_month <- month_reference_from_settings(months, settings, prefix)
  prediction_month <- month_prediction_from_settings(months, settings, prefix)
  month_terms <- character()

  for (m in setdiff(months, reference_month)) {
    term <- month_term_name(m)
    model_dat[[term]] <- as.integer(model_dat$month == m)
    month_terms <- c(month_terms, term)
  }

  model_dat$month_reference <- reference_month
  model_dat$month_prediction <- prediction_month
  model_dat$model_row_type <- "deployment_month"
  model_dat
}

summarise_deployment_month_rate <- function(deployments, wolf_events,
                                            settings, prefix) {
  model_dat <- deployments %>%
    left_join(wolf_events, by = "deploymentID") %>%
    mutate(
      wolf_events = tidyr::replace_na(wolf_events, 0L),
      total_effort_days = deploymentEffort,
      wolf_events_per_100_days = 100 * wolf_events / total_effort_days,
      month = format(start, "%Y-%m", tz = "UTC")
    ) %>%
    filter(!is.na(month), nzchar(month)) %>%
    arrange(plotID, start)

  if (!nrow(model_dat)) {
    stop("[", prefix, "] no valid deployment-month records after filtering.")
  }

  model_dat <- add_month_design(model_dat, settings, prefix)
  camera_rate <- camera_summary_from_model(model_dat)

  month_summary <- model_dat %>%
    group_by(month) %>%
    summarise(
      deployments = n_distinct(deploymentID),
      cameras = n_distinct(plotID),
      events = sum(wolf_events),
      effort_days = sum(total_effort_days),
      rate_per_100 = 100 * events / effort_days,
      .groups = "drop"
    )

  readr::write_csv(model_dat,
                   path_out(paste0(prefix, "_deployment_month_effort_rates.csv")))
  readr::write_csv(camera_rate,
                   path_out(paste0(prefix, "_camera_effort_rates.csv")))
  readr::write_csv(month_summary,
                   path_out(paste0(prefix, "_month_observed_summary.csv")))

  cat(sprintf(
    "[%s] cameras %d | deployments %d | positive rows %d | events %d | effort %.1f camera-days | observed %.2f /100\n",
    prefix,
    nrow(camera_rate),
    nrow(model_dat),
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

load_flat_survey <- function(cfg) {
  dat <- readr::read_csv(path_in(cfg$file), show_col_types = FALSE)

  required <- c("deploymentID", "eventID", "scientificName", "plotID",
                "deploymentEffort", "latitude", "longitude")
  stop_missing_columns(dat, required, paste0("[", cfg$prefix, "] flat file"))

  dat <- dat %>%
    mutate(
      deploymentID = na_if(as.character(deploymentID), ""),
      eventID = na_if(as.character(eventID), ""),
      plotID = na_if(as.character(plotID), ""),
      scientificName = as.character(scientificName),
      deploymentEffort = as.numeric(deploymentEffort),
      latitude = as.numeric(latitude),
      longitude = as.numeric(longitude)
    )

  deployments <- dat %>%
    filter(!is.na(deploymentID), !is.na(plotID)) %>%
    group_by(deploymentID, plotID) %>%
    summarise(
      latitude = mean(latitude, na.rm = TRUE),
      longitude = mean(longitude, na.rm = TRUE),
      deploymentEffort = first_finite(deploymentEffort),
      .groups = "drop"
    ) %>%
    filter(is.finite(latitude),
           is.finite(longitude),
           is.finite(deploymentEffort),
           deploymentEffort > 0)

  wolf_events <- dat %>%
    filter(scientificName %in% WOLF_NAMES,
           !is.na(deploymentID),
           !is.na(eventID)) %>%
    distinct(deploymentID, eventID) %>%
    count(deploymentID, name = "wolf_events")

  summarise_camera_rate(deployments, wolf_events, cfg$prefix)
}

load_camtrap_survey <- function(cfg) {
  dep <- readr::read_csv(path_in(cfg$deployments), show_col_types = FALSE)
  obs <- readr::read_csv(path_in(cfg$observations), show_col_types = FALSE)

  required_dep <- c("deploymentID", "locationID", "latitude", "longitude",
                    "deploymentStart", "deploymentEnd")
  required_obs <- c("deploymentID", "eventID", "scientificName")
  stop_missing_columns(dep, required_dep, paste0("[", cfg$prefix, "] deployments"))
  stop_missing_columns(obs, required_obs, paste0("[", cfg$prefix, "] observations"))

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
           deploymentEffort > 0)

  if (!nrow(deployments)) {
    stop("[", cfg$prefix, "] no valid dated deployments.")
  }

  wolf_events <- obs %>%
    transmute(
      deploymentID = na_if(as.character(deploymentID), ""),
      eventID = na_if(as.character(eventID), ""),
      scientificName = as.character(scientificName)
    ) %>%
    filter(scientificName %in% WOLF_NAMES,
           !is.na(deploymentID),
           !is.na(eventID)) %>%
    distinct(deploymentID, eventID) %>%
    count(deploymentID, name = "wolf_events")

  if (isTRUE(cfg$settings$use_month_effect)) {
    summarise_deployment_month_rate(deployments, wolf_events,
                                    cfg$settings, cfg$prefix)
  } else {
    summarise_camera_rate(deployments, wolf_events, cfg$prefix)
  }
}

load_survey <- function(cfg) {
  if (identical(cfg$type, "flat")) {
    load_flat_survey(cfg)
  } else {
    load_camtrap_survey(cfg)
  }
}


## 09. Spatial Domain: Camera Coordinates, Mesh, And Prediction Grid -----------

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

  area <- if (identical(PRED_DOMAIN, "hull")) {
    st_buffer(st_convex_hull(pts), settings$pred_buffer_m)
  } else {
    st_buffer(pts, settings$pred_buffer_m)
  }

  grid <- st_make_grid(area, cellsize = settings$cell_size_m, what = "centers")
  pred_sf <- st_sf(grid_id = seq_along(grid), geometry = grid, crs = st_crs(camera_sf))
  pred_sf <- pred_sf[lengths(st_intersects(pred_sf, area)) > 0, ]

  if (!identical(PRED_DOMAIN, "hull")) {
    dist_to_camera <- as.numeric(apply(st_distance(pred_sf, camera_sf), 1, min))
    pred_sf <- pred_sf[dist_to_camera <= settings$max_dist_m, ]
  }

  if (!nrow(pred_sf)) {
    stop("Prediction grid is empty. Check PRED_DOMAIN, buffers, and input coordinates.")
  }

  pred_sf
}


## 10. Residual Diagnostics: Moran's I, Semivariogram, And PPC -----------------

moran_perm <- function(coords, x, nperm = 499L) {
  ok <- is.finite(x)
  coords <- coords[ok, , drop = FALSE]
  x <- x[ok]
  n <- length(x)

  if (n < 5) {
    return(list(I = NA_real_, expected = NA_real_, p_value = NA_real_))
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
  if (!is.finite(z2) || z2 <= 0 || !is.finite(S0) || S0 <= 0) {
    return(list(I = NA_real_, expected = -1 / (n - 1), p_value = NA_real_))
  }

  I_stat <- function(v) {
    (n / S0) * sum(W * outer(v, v)) / sum(v^2)
  }

  I0 <- I_stat(z)
  perm <- replicate(nperm, I_stat(sample(z)))

  list(
    I = I0,
    expected = -1 / (n - 1),
    p_value = (1 + sum(perm >= I0)) / (nperm + 1)
  )
}

resid_variogram <- function(coords, residual, nbins = 12L) {
  D <- as.matrix(dist(coords))
  G <- as.matrix(dist(residual))^2 / 2

  d <- D[lower.tri(D)]
  g <- G[lower.tri(G)]
  if (!length(d)) {
    return(data.frame(dist = numeric(), gamma = numeric(), n = integer()))
  }

  br <- seq(0, quantile(d, 0.9, na.rm = TRUE), length.out = nbins + 1)
  bin <- cut(d, br, include.lowest = TRUE)

  data.frame(
    dist = tapply(d, bin, mean),
    gamma = tapply(g, bin, mean),
    n = as.integer(table(bin))
  )
}

summarise_ppc_simulations <- function(sim, model_dat, method) {
  yobs <- model_dat$y
  n <- length(yobs)
  camera_group <- as.factor(model_dat$plotID)
  yobs_camera <- as.numeric(rowsum(yobs, camera_group)[, 1])

  if (is.null(dim(sim))) {
    sim <- matrix(sim, ncol = 1)
  }
  nsim <- ncol(sim)

  sim_camera <- apply(sim, 2, function(x) as.numeric(rowsum(x, camera_group)[, 1]))
  if (is.null(dim(sim_camera))) {
    sim_camera <- matrix(sim_camera, ncol = 1)
  }

  total_events <- colSums(sim_camera)
  zero_fraction <- colMeans(sim_camera == 0)
  max_count <- apply(sim_camera, 2, max)

  pit <- (rowSums(sim < yobs) + runif(n) * rowSums(sim == yobs)) / nsim

  summary <- data.frame(
    stat = c("total_events", "zero_fraction", "max_camera_count"),
    observed = c(sum(yobs_camera), mean(yobs_camera == 0), max(yobs_camera)),
    sim_median = c(median(total_events), median(zero_fraction), median(max_count)),
    sim_q025 = c(quantile(total_events, 0.025),
                 quantile(zero_fraction, 0.025),
                 quantile(max_count, 0.025)),
    sim_q975 = c(quantile(total_events, 0.975),
                 quantile(zero_fraction, 0.975),
                 quantile(max_count, 0.975))
  )
  summary$pass <- summary$observed >= summary$sim_q025 &
    summary$observed <= summary$sim_q975
  summary$method <- method

  list(summary = summary, pit = pit, nsim = ncol(sim))
}

ppc_compute_marginal <- function(fit, model_dat, family, pi_hat, size_hat,
                                 nsim = PPC_NSIM) {
  effort <- model_dat$total_effort_days
  n <- nrow(model_dat)

  eta_mean <- model_dat$eta_mean
  eta_sd <- model_dat$eta_sd
  pi_marg <- if (is_zi(family)) hyp_marg(fit, PAT_ZPROB) else NULL
  size_marg <- if (is_nb(family)) hyp_marg(fit, PAT_NB_SIZE) else NULL

  draw_one_dataset <- function() {
    eta <- rnorm(n, eta_mean, eta_sd)
    pp <- if (!is.null(pi_marg)) {
      tryCatch(INLA::inla.rmarginal(1, pi_marg), error = function(e) pi_hat)
    } else {
      pi_hat
    }
    ss <- if (!is.null(size_marg)) {
      tryCatch(INLA::inla.rmarginal(1, size_marg), error = function(e) size_hat)
    } else {
      size_hat
    }
    fam_sim(effort * exp(eta), pp, family, ss)
  }

  sim <- replicate(nsim, draw_one_dataset())
  summarise_ppc_simulations(sim, model_dat, "marginal_predictor")
}

ppc_compute_joint <- function(fit, model_dat, family, pi_hat, size_hat,
                              obs_index, nsim = JOINT_PPC_NSIM) {
  samples <- tryCatch(
    INLA::inla.posterior.sample(nsim, fit),
    error = function(e) NULL
  )
  if (is.null(samples) || !length(samples)) return(NULL)

  first_latent <- samples[[1]]$latent
  pred_rows <- grep("^APredictor", rownames(first_latent))
  if (!length(pred_rows)) {
    pred_rows <- grep("^Predictor", rownames(first_latent))
  }
  if (!length(pred_rows)) return(NULL)

  if (length(pred_rows) == nrow(model_dat)) {
    pred_rows <- pred_rows
  } else if (length(pred_rows) >= max(obs_index)) {
    pred_rows <- pred_rows[obs_index]
  } else if (length(pred_rows) >= nrow(model_dat)) {
    pred_rows <- pred_rows[seq_len(nrow(model_dat))]
  } else {
    return(NULL)
  }

  effort <- model_dat$total_effort_days
  sim <- vapply(samples, function(s) {
    eta <- as.numeric(s$latent[pred_rows, 1])

    pp <- pi_hat
    if (is_zi(family)) {
      pi_idx <- grep(PAT_ZPROB, names(s$hyperpar), ignore.case = TRUE)
      if (length(pi_idx)) pp <- as.numeric(s$hyperpar[pi_idx[[1]]])
    }

    ss <- size_hat
    if (is_nb(family)) {
      size_idx <- grep(PAT_NB_SIZE, names(s$hyperpar), ignore.case = TRUE)
      if (length(size_idx)) ss <- as.numeric(s$hyperpar[size_idx[[1]]])
    }

    fam_sim(effort * exp(eta), pp, family, ss)
  }, numeric(nrow(model_dat)))

  summarise_ppc_simulations(sim, model_dat, "joint_posterior")
}

ppc_compute <- function(fit, model_dat, family, pi_hat, size_hat,
                        obs_index, nsim = PPC_NSIM) {
  if (use_joint_ppc(family)) {
    joint <- ppc_compute_joint(
      fit,
      model_dat,
      family,
      pi_hat,
      size_hat,
      obs_index,
      nsim = min(nsim, JOINT_PPC_NSIM)
    )
    if (!is.null(joint)) return(joint)
  }

  ppc_compute_marginal(fit, model_dat, family, pi_hat, size_hat, nsim = nsim)
}

ks_uniform_p_value <- function(pit) {
  pit <- pit[is.finite(pit)]
  if (length(pit) < 5) return(NA_real_)
  suppressWarnings(tryCatch(
    ks.test(pit, "punif")$p.value,
    error = function(e) NA_real_
  ))
}


## 11. Diagnostic Plots: Observed-Fitted, Residual Maps, PIT, Variogram --------

write_diagnostic_plots <- function(prefix, spec, model_dat,
                                   camera_diag, camera_sf, coords, pit) {
  readr::write_csv(
    model_dat,
    path_out(paste0(prefix, "_", spec$name, "_model_row_diagnostics.csv"))
  )

  readr::write_csv(
    camera_diag,
    path_out(paste0(prefix, "_", spec$name, "_camera_residual_diagnostics.csv"))
  )

  obs_fit_plot <- ggplot(camera_diag, aes(fitted_count, y)) +
    geom_point(alpha = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    scale_x_continuous(trans = "sqrt") +
    scale_y_continuous(trans = "sqrt") +
    labs(
      title = paste0("Observed vs fitted: ", prefix),
      subtitle = spec$name,
      x = "fitted wolf events",
      y = "observed wolf events"
    ) +
    theme_minimal(base_size = 12)
  ggsave(
    path_out(paste0(prefix, "_", spec$name, "_diag_obs_vs_fitted.png")),
    obs_fit_plot,
    width = 6,
    height = 5.5,
    dpi = 200
  )

  residual_sf <- camera_sf
  residual_sf$pearson <- camera_diag$pearson
  residual_plot <- ggplot(residual_sf) +
    geom_sf(aes(colour = pearson, size = abs(pearson)), alpha = 0.9) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_size_continuous(range = c(2, 8)) +
    coord_sf(datum = NA) +
    labs(title = paste0("Spatial residuals: ", prefix), subtitle = spec$name) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())
  ggsave(
    path_out(paste0(prefix, "_", spec$name, "_diag_spatial_residuals.png")),
    residual_plot,
    width = 7,
    height = 7,
    dpi = 200
  )

  pit_plot <- ggplot(data.frame(pit = pit), aes(pit)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 20,
                   fill = "grey55",
                   colour = "white",
                   boundary = 0) +
    geom_hline(yintercept = 1, linetype = 2, colour = "red") +
    labs(
      title = paste0("Posterior predictive PIT: ", prefix),
      subtitle = sprintf("%s | mean %.3f | KS p %.3f",
                         spec$name,
                         mean(pit, na.rm = TRUE),
                         ks_uniform_p_value(pit)),
      x = "PIT",
      y = "density"
    ) +
    xlim(0, 1) +
    theme_minimal(base_size = 12)
  ggsave(
    path_out(paste0(prefix, "_", spec$name, "_diag_pit_hist.png")),
    pit_plot,
    width = 6,
    height = 5,
    dpi = 200
  )

  variogram <- resid_variogram(coords, camera_diag$pearson)
  variogram_plot <- ggplot(variogram, aes(dist, gamma)) +
    geom_point(aes(size = n)) +
    geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
    labs(
      title = paste0("Residual semivariogram: ", prefix),
      subtitle = spec$name,
      x = "distance (m)",
      y = "semivariance"
    ) +
    theme_minimal(base_size = 12)
  ggsave(
    path_out(paste0(prefix, "_", spec$name, "_diag_resid_variogram.png")),
    variogram_plot,
    width = 6.5,
    height = 5,
    dpi = 200
  )
}


## 12. Fit Diagnostics: Compute Fitted Counts And Required Checks --------------

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
      n_model_rows = n(),
      pearson = (y - fitted_count) / sqrt(pmax(fit_var, 1e-9)),
      .groups = "drop"
    ) %>%
    arrange(plotID)
}

diagnose_fit <- function(fit, model_dat, camera_sf, spec, obs_index,
                         prefix, write_files = FALSE) {
  family <- spec$family
  pi_hat <- if (is_zi(family)) hyp_point(fit, PAT_ZPROB) else NA_real_
  size_hat <- if (is_nb(family)) nb_size_point(fit) else NA_real_

  eta_mean <- fit$summary.linear.predictor$mean[obs_index]
  eta_sd <- fit$summary.linear.predictor$sd[obs_index]

  model_dat$eta_mean <- eta_mean
  model_dat$eta_sd <- eta_sd
  model_dat$mu_count <- model_dat$total_effort_days * exp(eta_mean + 0.5 * eta_sd^2)
  model_dat$fitted_count <- fam_mean(model_dat$mu_count, pi_hat, family, size_hat)
  model_dat$fit_var <- fam_var(model_dat$mu_count, pi_hat, family, size_hat)
  model_dat$pearson <- (model_dat$y - model_dat$fitted_count) /
    sqrt(pmax(model_dat$fit_var, 1e-9))

  pearson_disp <- mean(model_dat$pearson^2, na.rm = TRUE)
  camera_diag <- camera_residual_diagnostics(model_dat)
  coords <- st_coordinates(camera_sf)
  moran <- moran_perm(coords, camera_diag$pearson)

  ppc <- ppc_compute(fit, model_dat, family, pi_hat, size_hat,
                     obs_index, nsim = PPC_NSIM)
  pit <- ppc$pit[is.finite(ppc$pit)]
  ppc_pit_ks <- ks_uniform_p_value(pit)

  ppc_row <- function(stat, col) {
    ppc$summary[ppc$summary$stat == stat, col][[1]]
  }

  total_pass <- isTRUE(ppc_row("total_events", "pass"))
  zero_pass <- isTRUE(ppc_row("zero_fraction", "pass"))
  max_pass <- isTRUE(ppc_row("max_camera_count", "pass"))
  moran_pass <- is.finite(moran$p_value) && moran$p_value >= MORAN_ALPHA

  if (write_files) {
    readr::write_csv(
      ppc$summary,
      path_out(paste0(prefix, "_", spec$name, "_posterior_predictive_check.csv"))
    )
    write_diagnostic_plots(prefix, spec, model_dat,
                           camera_diag, camera_sf, coords, pit)
  }

  list(
    model_dat = model_dat,
    camera_diag = camera_diag,
    pi_hat = pi_hat,
    size_hat = size_hat,
    pearson_disp = pearson_disp,
    pearson_disp_camera = mean(camera_diag$pearson^2, na.rm = TRUE),
    moran_I = moran$I,
    moran_p = moran$p_value,
    ppc_pit_ks = ppc_pit_ks,
    pit_mean = mean(pit, na.rm = TRUE),
    ppc_method = unique(ppc$summary$method)[[1]],
    ppc_nsim = ppc$nsim,
    ppc = ppc$summary,
    ppc_total_pass = total_pass,
    ppc_zero_pass = zero_pass,
    ppc_max_pass = max_pass,
    moran_pass = moran_pass,
    diagnostics_ok = total_pass && zero_pass && max_pass && moran_pass
  )
}


## 13. Model Fitting: INLA Stack, SPDE Field, Prediction Surface ---------------

prediction_fixed_effects <- function(model_dat, fixed_terms, settings, n_pred) {
  fixed_pred <- as.data.frame(matrix(0, nrow = n_pred, ncol = length(fixed_terms)))
  names(fixed_pred) <- fixed_terms

  if ("intercept" %in% fixed_terms) {
    fixed_pred$intercept <- 1
  }

  month_terms <- intersect(temporal_month_terms(model_dat), fixed_terms)
  if (length(month_terms)) {
    prediction_month <- settings$month_prediction
    if (is.null(prediction_month) || !nzchar(prediction_month)) {
      prediction_month <- unique(model_dat$month_prediction)
      prediction_month <- prediction_month[!is.na(prediction_month)]
      prediction_month <- if (length(prediction_month)) prediction_month[[1]] else NA_character_
    }

    prediction_term <- month_term_name(prediction_month)
    if (prediction_term %in% month_terms) {
      fixed_pred[[prediction_term]] <- 1
    }
  }

  fixed_pred
}

write_month_coefficients <- function(fit, model_dat, settings, prefix) {
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

  readr::write_csv(out, path_out(paste0(prefix, "_month_coefficients.csv")))
  invisible(out)
}

write_model_hyperparameters <- function(fit, prefix) {
  if (is.null(fit$summary.hyperpar)) return(invisible(NULL))

  out <- data.frame(
    parameter = rownames(fit$summary.hyperpar),
    fit$summary.hyperpar,
    row.names = NULL,
    check.names = FALSE
  )
  readr::write_csv(out, path_out(paste0(prefix, "_hyperparameters.csv")))
  invisible(out)
}

fit_final_model <- function(camera_rate, settings, spec, survey_prefix,
                            add_prediction = TRUE, write_files = TRUE) {
  prefix <- survey_prefix
  cat(sprintf(
    "\n[%s] fitting %s: family=%s, prediction_grid=%s\n",
    prefix, spec$name, spec$family, add_prediction
  ))

  model_dat <- camera_rate %>%
    mutate(y = as.integer(wolf_events), intercept = 1)

  camera_summary <- camera_summary_from_model(model_dat)
  camera_sf <- camera_to_utm(camera_summary)
  coords_camera <- st_coordinates(camera_sf)
  colnames(coords_camera) <- c("x", "y")

  obs_sf <- model_dat %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(EPSG_UTM)
  coords_obs <- st_coordinates(obs_sf)
  colnames(coords_obs) <- c("x", "y")
  n_obs <- nrow(model_dat)
  fixed_terms <- fixed_effect_terms(model_dat)
  fixed_obs <- as.data.frame(model_dat[, fixed_terms, drop = FALSE])

  pred_sf <- NULL
  coords_pred <- NULL
  n_pred <- 0L
  if (add_prediction) {
    pred_sf <- prediction_grid(camera_sf, settings)
    coords_pred <- st_coordinates(pred_sf)
    colnames(coords_pred) <- c("x", "y")
    n_pred <- nrow(pred_sf)
    cat(sprintf("[%s] prediction cells: %d at %.0f m\n",
                prefix, n_pred, settings$cell_size_m))
  }

  mesh_loc <- if (add_prediction && isTRUE(settings$include_grid_in_mesh)) {
    rbind(coords_camera, coords_pred)
  } else {
    coords_camera
  }

  spde_obj <- build_spatial(mesh_loc, settings)
  A_obs <- INLA::inla.spde.make.A(spde_obj$mesh, loc = coords_obs)

  A_obs_list <- list(A_obs, 1)
  effects_obs <- list(
    spatial = spde_obj$s_index,
    fixed = fixed_obs
  )

  formula_terms <- c(fixed_terms, "f(spatial, model = spde_obj$spde)")

  if (add_prediction) {
    A_pred <- INLA::inla.spde.make.A(spde_obj$mesh, loc = coords_pred)
    A_pred_list <- list(A_pred, 1)
    fixed_pred <- prediction_fixed_effects(model_dat, fixed_terms, settings, n_pred)
    effects_pred <- list(
      spatial = spde_obj$s_index,
      fixed = fixed_pred
    )
  }

  stack_obs <- INLA::inla.stack(
    tag = "obs",
    data = list(y = model_dat$y, e = model_dat$total_effort_days),
    A = A_obs_list,
    effects = effects_obs
  )

  if (add_prediction) {
    stack_pred <- INLA::inla.stack(
      tag = "pred",
      data = list(y = rep(NA_real_, n_pred), e = rep(100, n_pred)),
      A = A_pred_list,
      effects = effects_pred
    )
    stack_all <- INLA::inla.stack(stack_obs, stack_pred)
  } else {
    stack_all <- stack_obs
  }

  stack_data <- INLA::inla.stack.data(stack_all)
  obs_index <- INLA::inla.stack.index(stack_all, tag = "obs")$data
  formula <- as.formula(paste("y ~ 0 +", paste(formula_terms, collapse = " + ")))

  fit <- INLA::inla(
    formula,
    family = spec$family,
    data = stack_data,
    E = stack_data$e,
    control.predictor = list(
      A = INLA::inla.stack.A(stack_all),
      compute = TRUE,
      link = 1
    ),
    control.compute = list(config = use_joint_ppc(spec$family) && !add_prediction),
    control.fixed = make_control_fixed(fixed_terms),
    control.family = make_control_family(spec$family),
    verbose = FALSE
  )

  fit_diag <- fit
  obs_index_diag <- obs_index
  if (use_joint_ppc(spec$family) && add_prediction) {
    cat(sprintf(
      "[%s] fitting observed-data diagnostic refit for joint posterior PPC\n",
      prefix
    ))

    stack_obs_data <- INLA::inla.stack.data(stack_obs)
    fit_diag <- INLA::inla(
      formula,
      family = spec$family,
      data = stack_obs_data,
      E = stack_obs_data$e,
      control.predictor = list(
        A = INLA::inla.stack.A(stack_obs),
        compute = TRUE,
        link = 1
      ),
      control.compute = list(config = TRUE),
      control.fixed = make_control_fixed(fixed_terms),
      control.family = make_control_family(spec$family),
      verbose = FALSE
    )
    obs_index_diag <- INLA::inla.stack.index(stack_obs, tag = "obs")$data
  }

  diagnostics <- diagnose_fit(
    fit_diag,
    model_dat,
    camera_sf,
    spec,
    obs_index_diag,
    prefix,
    write_files = write_files
  )

  if (write_files) {
    write_prior_posterior_plots(fit_diag, settings, spec, prefix)
    write_month_coefficients(fit_diag, model_dat, settings, prefix)
    write_model_hyperparameters(fit_diag, prefix)
  }

  rasters <- NULL
  if (add_prediction) {
    pred_index <- INLA::inla.stack.index(stack_all, tag = "pred")$data
    eta_mean <- fit$summary.linear.predictor$mean[pred_index]
    eta_sd <- fit$summary.linear.predictor$sd[pred_index]

    rate100_count <- 100 * exp(eta_mean + 0.5 * eta_sd^2)
    pred_sf$mean <- fam_mean(rate100_count, diagnostics$pi_hat,
                             spec$family, diagnostics$size_hat)
    pred_sf$cv <- sqrt(expm1(eta_sd^2))
    pred_sf$sd <- pred_sf$mean * pred_sf$cv
    pred_sf$x <- coords_pred[, 1]
    pred_sf$y <- coords_pred[, 2]

    overall_rate <- 100 * sum(model_dat$y) / sum(model_dat$total_effort_days)
    if (MAP_EXCEEDANCE) {
      pi0 <- if (is_zi(spec$family) && is.finite(diagnostics$pi_hat)) {
        diagnostics$pi_hat
      } else {
        0
      }
      threshold_latent_rate <- (EXCEED_MULT * overall_rate) / (100 * (1 - pi0))
      pred_sf$exceed <- 1 - pnorm((log(threshold_latent_rate) - eta_mean) / eta_sd)
    }

    wkt <- st_crs(pred_sf)$wkt
    pred_table <- st_drop_geometry(pred_sf)
    make_raster <- function(col) {
      terra::rast(pred_table[, c("x", "y", col)], type = "xyz", crs = wkt)
    }

    r_mean <- make_raster("mean")
    r_sd <- make_raster("sd")
    r_cv <- make_raster("cv")
    names(r_mean) <- "wolf_events_per_100_camera_days"
    names(r_sd) <- "posterior_sd"
    names(r_cv) <- "posterior_cv"

    rasters <- list(mean = r_mean, sd = r_sd, cv = r_cv, exceed = NULL)

    terra::writeRaster(
      r_mean,
      path_out(paste0(prefix, "_final_predicted_events_per_100_days_mean.tif")),
      overwrite = TRUE
    )
    terra::writeRaster(
      r_sd,
      path_out(paste0(prefix, "_final_predicted_events_per_100_days_sd.tif")),
      overwrite = TRUE
    )
    terra::writeRaster(
      r_cv,
      path_out(paste0(prefix, "_final_predicted_events_per_100_days_cv.tif")),
      overwrite = TRUE
    )

    if (MAP_EXCEEDANCE) {
      r_exceed <- make_raster("exceed")
      names(r_exceed) <- "exceedance_probability"
      rasters$exceed <- r_exceed
      terra::writeRaster(
        r_exceed,
        path_out(paste0(prefix, "_final_exceedance_prob.tif")),
        overwrite = TRUE
      )
    }

    plot_map_outputs(prefix, spec, camera_sf, model_dat, rasters, overall_rate)
  }

  list(
    fit = fit,
    fit_diagnostics = fit_diag,
    spec = spec,
    diag = diagnostics,
    camera_sf = camera_sf,
    model_dat = diagnostics$model_dat,
    rasters = rasters
  )
}


## 14. Map Outputs: Mean Surface, Uncertainty, And Exceedance ------------------

plot_map_outputs <- function(prefix, spec, camera_sf, model_dat,
                             rasters, overall_rate) {
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
            shape = 21,
            size = 1.4,
            fill = "white",
            colour = "grey35",
            stroke = 0.25) +
    geom_sf(data = positive_sf,
            aes(size = wolf_events_per_100_days),
            shape = 21,
            fill = "black",
            colour = "white",
            stroke = 0.25,
            alpha = 0.9) +
    scale_fill_viridis_c(
      option = "magma",
      na.value = NA,
      name = "predicted events\n/100 camera-days",
      labels = label_number(accuracy = 0.01)
    ) +
    scale_size_continuous(
      range = c(2, 7),
      name = "observed events\n/100 camera-days",
      labels = label_number(accuracy = 0.01)
    ) +
    coord_sf(datum = NA) +
    labs(
      title = paste0("Wolf encounter-frequency surface: ", prefix),
      subtitle = sprintf("final model: %s (%s)", spec$name, spec$family),
      x = "Easting, UTM 34N",
      y = "Northing, UTM 34N"
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), legend.position = "right")

  ggsave(
    path_out(paste0(prefix, "_final_event_frequency_mean.png")),
    mean_plot,
    width = 9.5,
    height = 9,
    dpi = 350
  )

  cv_df <- raster_to_df(rasters$cv, "cv")
  cv_plot <- ggplot() +
    geom_raster(data = cv_df, aes(x, y, fill = cv), interpolate = TRUE) +
    geom_sf(data = camera_sf,
            shape = 21,
            size = 1.4,
            fill = "white",
            colour = "grey35",
            stroke = 0.25) +
    scale_fill_viridis_c(
      option = "viridis",
      na.value = NA,
      name = "prediction CV",
      labels = label_number(accuracy = 0.01)
    ) +
    coord_sf(datum = NA) +
    labs(
      title = paste0("Uncertainty surface: ", prefix),
      subtitle = "posterior coefficient of variation of latent rate",
      x = "Easting, UTM 34N",
      y = "Northing, UTM 34N"
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), legend.position = "right")

  ggsave(
    path_out(paste0(prefix, "_final_event_frequency_cv.png")),
    cv_plot,
    width = 9.5,
    height = 9,
    dpi = 350
  )

  if (!is.null(rasters$exceed)) {
    exceed_df <- raster_to_df(rasters$exceed, "p")
    exceed_plot <- ggplot() +
      geom_raster(data = exceed_df, aes(x, y, fill = p), interpolate = TRUE) +
      geom_sf(data = camera_sf,
              shape = 21,
              size = 1.4,
              fill = "white",
              colour = "grey35",
              stroke = 0.25) +
      scale_fill_viridis_c(
        option = "inferno",
        limits = c(0, 1),
        na.value = NA,
        name = "P(rate > threshold)"
      ) +
      coord_sf(datum = NA) +
      labs(
        title = paste0("Robustly elevated encounter rate: ", prefix),
        subtitle = sprintf(
          "threshold = %.2f events / 100 camera-days (%.1fx observed mean)",
          EXCEED_MULT * overall_rate,
          EXCEED_MULT
        ),
        x = "Easting, UTM 34N",
        y = "Northing, UTM 34N"
      ) +
      theme_minimal(base_size = 13) +
      theme(panel.grid = element_blank(), legend.position = "right")

    ggsave(
      path_out(paste0(prefix, "_final_exceedance_prob.png")),
      exceed_plot,
      width = 9.5,
      height = 9,
      dpi = 350
    )
  }
}


## 15. Spatial Block Cross-Validation: Held-Out Spatial Folds ------------------

spatial_block_cv <- function(camera_rate, settings, spec, prefix, K = 5L) {
  cat(sprintf("\n[%s] spatial block CV for final model %s (K=%d)\n",
              prefix, spec$name, K))

  model_dat <- camera_rate %>%
    mutate(y = as.integer(wolf_events), intercept = 1)
  fixed_terms <- fixed_effect_terms(model_dat)

  camera_summary <- camera_summary_from_model(model_dat)
  camera_sf <- camera_to_utm(camera_summary)
  coords_camera <- st_coordinates(camera_sf)
  row_sf <- model_dat %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>%
    st_transform(EPSG_UTM)
  coords_row <- st_coordinates(row_sf)

  y <- model_dat$y
  effort <- model_dat$total_effort_days

  K_final <- min(K, nrow(coords_camera) - 1L)
  if (K_final < 2) {
    cat(sprintf("[%s] spatial CV skipped: too few cameras.\n", prefix))
    return(invisible(NULL))
  }

  fold <- tryCatch(
    kmeans(scale(coords_camera), centers = K_final, nstart = 10)$cluster,
    error = function(e) sample(rep_len(seq_len(K_final), nrow(coords_camera)))
  )
  row_fold <- fold[match(model_dat$plotID, camera_summary$plotID)]

  rows <- list()
  failed_folds <- character()

  for (f in sort(unique(fold))) {
    test <- which(row_fold == f)
    train <- which(row_fold != f)

    result <- tryCatch({
      spde_obj <- build_spatial(coords_camera, settings)
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
        effects = list(
          spatial = spde_obj$s_index,
          fixed = fixed_train
        )
      )

      stack_test <- INLA::inla.stack(
        tag = "test",
        data = list(y = rep(NA_real_, length(test)), e = effort[test]),
        A = list(A_test, 1),
        effects = list(
          spatial = spde_obj$s_index,
          fixed = fixed_test
        )
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
        family = spec$family,
        data = stack_data,
        E = stack_data$e,
        control.predictor = list(
          A = INLA::inla.stack.A(stack_all),
          compute = TRUE,
          link = 1
        ),
        control.compute = list(config = FALSE),
        control.fixed = make_control_fixed(fixed_terms),
        control.family = make_control_family(spec$family),
        verbose = FALSE
      )

      test_index <- INLA::inla.stack.index(stack_all, tag = "test")$data
      eta_mean <- fit_fold$summary.linear.predictor$mean[test_index]
      eta_sd <- fit_fold$summary.linear.predictor$sd[test_index]
      pi_fold <- if (is_zi(spec$family)) hyp_point(fit_fold, PAT_ZPROB) else NA_real_
      size_fold <- if (is_nb(spec$family)) nb_size_point(fit_fold) else NA_real_

      mu_mean <- effort[test] * exp(eta_mean + 0.5 * eta_sd^2)
      expected_y <- fam_mean(mu_mean, pi_fold, spec$family, size_fold)

      sim_list <- lapply(seq_along(test), function(j) {
        eta_draw <- rnorm(CV_NSIM, eta_mean[j], eta_sd[j])
        mu_draw <- effort[test[j]] * exp(eta_draw)
        fam_sim(mu_draw, pi_fold, spec$family, size_fold)
      })

      lo <- vapply(sim_list, quantile, numeric(1), 0.05)
      hi <- vapply(sim_list, quantile, numeric(1), 0.95)

      lpd <- vapply(seq_along(test), function(j) {
        eta_draw <- rnorm(CV_NSIM, eta_mean[j], eta_sd[j])
        mu_draw <- effort[test[j]] * exp(eta_draw)
        log_mean_exp(fam_logpmf(y[test[j]], mu_draw, pi_fold,
                                spec$family, size_fold))
      }, numeric(1))

      data.frame(
        fold = f,
        plotID = model_dat$plotID[test],
        model_row_type = model_dat$model_row_type[test],
        n = length(test),
        y = y[test],
        Ey = expected_y,
        rate_obs = 100 * y[test] / effort[test],
        rate_pred = 100 * expected_y / effort[test],
        lpd = lpd,
        covered_90 = y[test] >= lo & y[test] <= hi
      )
    }, error = function(e) {
      failed_folds <<- c(failed_folds, paste0(f, ": ", conditionMessage(e)))
      NULL
    })

    if (!is.null(result)) {
      rows[[length(rows) + 1L]] <- result
    }
  }

  if (length(failed_folds)) {
    stop("[", prefix, "] spatial CV failed for fold(s): ",
         paste(failed_folds, collapse = "; "))
  }
  if (!length(rows)) {
    stop("[", prefix, "] spatial CV produced no successful folds.")
  }

  cv <- do.call(rbind, rows)
  readr::write_csv(cv, path_out(paste0(prefix, "_final_spatial_block_cv.csv")))

  summary <- data.frame(
    metric = c("mean_log_predictive_density", "rmse_count",
               "rmse_rate_per100", "coverage_90"),
    value = c(
      mean(cv$lpd),
      sqrt(mean((cv$y - cv$Ey)^2)),
      sqrt(mean((cv$rate_obs - cv$rate_pred)^2)),
      mean(cv$covered_90)
    )
  )
  readr::write_csv(summary,
                   path_out(paste0(prefix, "_final_spatial_block_cv_summary.csv")))

  cat(sprintf(
    "[%s] spatial CV: mean LPD %.3f | RMSE count %.2f | RMSE rate %.2f | 90%% coverage %.2f\n",
    prefix,
    summary$value[1],
    summary$value[2],
    summary$value[3],
    summary$value[4]
  ))

  invisible(list(cv = cv, summ = summary))
}


## 16. Reporting: Validation Text, Priors, Caveats, And Failure Flags ----------

diagnostic_failures <- function(diag) {
  c(
    if (!isTRUE(diag$ppc_total_pass)) "PPC total events" else NULL,
    if (!isTRUE(diag$ppc_zero_pass)) "PPC zero fraction" else NULL,
    if (!isTRUE(diag$ppc_max_pass)) "PPC max camera count" else NULL,
    if (!isTRUE(diag$moran_pass)) {
      if (is.finite(diag$moran_p)) {
        sprintf("residual spatial autocorrelation (Moran's I = %.3f, p = %.3f)",
                diag$moran_I, diag$moran_p)
      } else {
        "residual Moran's I was not evaluable"
      }
    } else {
      NULL
    }
  )
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
    sprintf("  Intercept: Gaussian(0, prec = %g), SD %.1f on the log scale",
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
    if (is_zi(spec$family)) {
      sprintf("  Zero-inflation logit(p): Gaussian(%.2f, prec = %.2f)",
              PRIOR_ZI_LOGIT_MEAN,
              PRIOR_ZI_LOGIT_PREC)
    } else {
      NULL
    },
    if (is_nb(spec$family)) {
      sprintf("  Negative-binomial size: loggamma(shape = %.2f, rate = %.2f)",
              PRIOR_NB_SIZE_LOGGAMMA[1],
              PRIOR_NB_SIZE_LOGGAMMA[2])
    } else {
      NULL
    }
  )
}

hyper_lines_for_report <- function(spec, diag) {
  lines <- c()

  if (is_zi(spec$family)) {
    lines <- c(lines, sprintf("  Zero-inflation probability mean: %.3f",
                              diag$pi_hat))
  }

  if (is_nb(spec$family)) {
    lines <- c(lines, sprintf("  Negative-binomial size mean: %.3f",
                              diag$size_hat))
  }

  if (!length(lines)) return(NULL)
  c("", "Estimated likelihood hyperparameters:", lines)
}

temporal_lines_for_report <- function(settings) {
  if (!isTRUE(settings$use_month_effect)) return(NULL)

  c(
    "",
    "Temporal structure:",
    "  Month is included as a fixed effect.",
    sprintf("  Reference month for coefficients: %s", settings$month_reference),
    sprintf("  Month used for prediction maps: %s", settings$month_prediction)
  )
}

write_validation_report <- function(prefix, cfg, spec, camera_rate, diag, cv) {
  failures <- diagnostic_failures(diag)
  passes_required <- isTRUE(diag$diagnostics_ok)
  n_cameras <- dplyr::n_distinct(camera_rate$plotID)
  n_model_rows <- nrow(camera_rate)

  report <- c(
    sprintf("Survey: %s", cfg$label),
    sprintf("Final model: %s (family = %s)", spec$name, spec$family),
    sprintf("Run profile: %s", RUN_PROFILE),
    sprintf("PPC simulations: %d", diag$ppc_nsim),
    sprintf("Spatial CV: %s", if (RUN_FINAL_SPATIAL_CV) "run" else "not run"),
    sprintf(
      "Cameras: %d | model rows: %d | positive rows: %d | events: %d | effort: %.1f camera-days",
      n_cameras,
      n_model_rows,
      sum(camera_rate$wolf_events > 0),
      sum(camera_rate$wolf_events),
      sum(camera_rate$total_effort_days)
    ),
    "",
    "Validation:",
    sprintf("  PPC method: %s", diag$ppc_method),
    sprintf("  Pearson dispersion, model rows: %.3f", diag$pearson_disp),
    sprintf("  Pearson dispersion, camera aggregates: %.3f", diag$pearson_disp_camera),
    sprintf("  PPC total events pass: %s", isTRUE(diag$ppc_total_pass)),
    sprintf("  PPC zero fraction pass: %s", isTRUE(diag$ppc_zero_pass)),
    sprintf("  PPC max count pass: %s", isTRUE(diag$ppc_max_pass)),
    sprintf("  Residual Moran's I: %.3f (p = %.3f)", diag$moran_I, diag$moran_p),
    sprintf("  Residual Moran pass: %s", isTRUE(diag$moran_pass)),
    sprintf("  PPC PIT KS p: %.4g", diag$ppc_pit_ks),
    sprintf("  Passes required checks: %s", passes_required),
    if (!passes_required) {
      sprintf("  Stated limitation: %s", paste(failures, collapse = "; "))
    } else {
      NULL
    },
    if (!is.null(cv)) "" else NULL,
    if (!is.null(cv)) "Spatial block cross-validation:" else NULL,
    if (!is.null(cv)) sprintf("  mean LPD: %.3f", cv$summ$value[1]) else NULL,
    if (!is.null(cv)) sprintf("  RMSE count: %.2f", cv$summ$value[2]) else NULL,
    if (!is.null(cv)) sprintf("  RMSE rate /100: %.2f", cv$summ$value[3]) else NULL,
    if (!is.null(cv)) sprintf("  90%% coverage: %.2f", cv$summ$value[4]) else NULL,
    prior_lines_for_report(cfg$settings, spec),
    hyper_lines_for_report(spec, diag),
    temporal_lines_for_report(cfg$settings),
    "",
    "Interpretation:",
    "  Relative wolf encounter frequency in events per 100 camera-days.",
    "  Not abundance, density, occupancy, or population size.",
    if (!is.null(cfg$caveat)) "" else NULL,
    if (!is.null(cfg$caveat)) sprintf("Temporal note: %s", cfg$caveat) else NULL
  )

  writeLines(report, path_out(paste0(prefix, "_VALIDATION_REPORT.txt")))
  invisible(report)
}


## 17. Per-Survey Workflow: Load, Fit, Map, Validate, Report -------------------

run_final_survey <- function(name, cfg) {
  prefix <- cfg$prefix
  cat(sprintf("\n\n==================== %s: %s ====================\n",
              prefix, cfg$label))

  if (is.null(cfg$final_model)) {
    stop("[", prefix, "] no final_model set.")
  }

  spec <- finalise_spec(cfg$final_model)
  settings <- cfg$settings
  camera_rate <- load_survey(cfg)

  cat(sprintf("[%s] FINAL model: %s (family = %s)\n",
              prefix, spec$name, spec$family))

  fit <- fit_final_model(
    camera_rate,
    settings,
    spec,
    prefix,
    add_prediction = TRUE,
    write_files = TRUE
  )

  diag <- fit$diag
  failures <- diagnostic_failures(diag)
  passes_required <- isTRUE(diag$diagnostics_ok)

  cat(sprintf("\n[%s] VALIDATION\n", prefix))
  cat(sprintf("[%s]   PPC method: %s\n", prefix, diag$ppc_method))
  cat(sprintf("[%s]   Pearson dispersion, model rows: %.3f\n",
              prefix, diag$pearson_disp))
  cat(sprintf("[%s]   Pearson dispersion, camera aggregates: %.3f\n",
              prefix, diag$pearson_disp_camera))
  cat(sprintf(
    "[%s]   PPC pass (total / zero / max): %s / %s / %s\n",
    prefix,
    isTRUE(diag$ppc_total_pass),
    isTRUE(diag$ppc_zero_pass),
    isTRUE(diag$ppc_max_pass)
  ))
  cat(sprintf("[%s]   residual Moran's I: %.3f (p = %.3f)\n",
              prefix, diag$moran_I, diag$moran_p))
  cat(sprintf("[%s]   PPC PIT KS p: %.3g\n", prefix, diag$ppc_pit_ks))
  if (is_zi(spec$family)) {
    cat(sprintf("[%s]   zero-inflation probability mean: %.3f\n",
                prefix, diag$pi_hat))
  }
  if (is_nb(spec$family)) {
    cat(sprintf("[%s]   negative-binomial size mean: %.3f\n",
                prefix, diag$size_hat))
  }
  cat(sprintf("[%s]   passes required checks: %s\n", prefix, passes_required))

  if (!passes_required) {
    cat(sprintf("[%s]   NOTE: mapped despite failing %s.\n",
                prefix, paste(failures, collapse = "; ")))
    cat(sprintf("[%s]         Treat the flagged aspect as a stated limitation.\n",
                prefix))
  }
  if (!is.null(cfg$caveat)) {
    cat(sprintf("[%s]   temporal note: %s\n", prefix, cfg$caveat))
  }
  if (isTRUE(settings$use_month_effect)) {
    cat(sprintf("[%s]   month effect: reference=%s | prediction=%s\n",
                prefix, settings$month_reference, settings$month_prediction))
  }

  cv <- NULL
  if (RUN_FINAL_SPATIAL_CV) {
    cv <- spatial_block_cv(camera_rate, settings, spec, prefix, K = CV_K)
  }

  write_validation_report(prefix, cfg, spec, camera_rate, diag, cv)

  invisible(list(
    prefix = prefix,
    camera_rate = camera_rate,
    spec = spec,
    final = fit,
    cv = cv
  ))
}


## 18. Main Run: Validate Inputs, Run All Surveys, Stop On Failures ------------

validate_requested_surveys <- function() {
  unknown <- setdiff(RUN_SURVEYS, names(surveys))
  if (length(unknown)) {
    stop("Unknown survey name(s) in RUN_SURVEYS: ", paste(unknown, collapse = ", "))
  }

  for (nm in RUN_SURVEYS) {
    cfg <- surveys[[nm]]
    needed <- if (identical(cfg$type, "camtrap")) {
      c(cfg$deployments, cfg$observations)
    } else {
      cfg$file
    }

    missing <- needed[!file.exists(path_in(needed))]
    if (length(missing)) {
      stop("[", cfg$prefix, "] missing input file(s): ",
           paste(missing, collapse = ", "))
    }
  }

  invisible(TRUE)
}

write_run_manifest <- function(results) {
  manifest <- data.frame(
    survey = names(results),
    prefix = vapply(results, function(x) x$prefix, character(1)),
    model = vapply(results, function(x) x$spec$name, character(1)),
    family = vapply(results, function(x) x$spec$family, character(1)),
    cameras = vapply(results, function(x) dplyr::n_distinct(x$camera_rate$plotID), integer(1)),
    model_rows = vapply(results, function(x) nrow(x$camera_rate), integer(1)),
    events = vapply(results, function(x) sum(x$camera_rate$wolf_events), numeric(1)),
    effort_days = vapply(results, function(x) sum(x$camera_rate$total_effort_days), numeric(1)),
    month_effect = vapply(results, function(x) {
      length(temporal_month_terms(x$final$model_dat)) > 0
    }, logical(1)),
    diagnostics_ok = vapply(results, function(x) isTRUE(x$final$diag$diagnostics_ok), logical(1)),
    run_profile = RUN_PROFILE
  )

  readr::write_csv(manifest, path_out("wolf_final_run_manifest.csv"))
  invisible(manifest)
}

validate_requested_surveys()

results <- list()
failures <- character()

for (nm in RUN_SURVEYS) {
  cfg <- surveys[[nm]]
  result <- tryCatch(
    run_final_survey(nm, cfg),
    error = function(e) {
      failures <<- c(failures, sprintf("[%s] %s", cfg$prefix, conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(result)) {
    results[[nm]] <- result
  }
}

if (length(failures)) {
  stop("One or more surveys failed:\n  ", paste(failures, collapse = "\n  "))
}

manifest <- write_run_manifest(results)

cat("\nAll requested surveys completed successfully.\n")
cat("Final outputs are in:\n  ", OUTPUT_DIR, "\n", sep = "")
cat("Key files per survey:\n")
cat("  wolf_*_VALIDATION_REPORT.txt              (model, diagnostic status, priors, CV)\n")
cat("  wolf_*_final_predicted_events_per_100_days_mean.tif / _sd.tif / _cv.tif\n")
cat("  wolf_*_final_event_frequency_mean.png / _cv.png\n")
cat("  wolf_*_final_exceedance_prob.tif / .png    (if MAP_EXCEEDANCE)\n")
cat("  wolf_*_*_posterior_predictive_check.csv, _diag_*.png\n")
cat("  wolf_*_prior_posterior_*.png / .csv        (prior-posterior overlays)\n")
cat("  wolf_*_hyperparameters.csv                 (likelihood and spatial hyperparameters)\n")
cat("  wolf_*_month_coefficients.csv              (for month-adjusted models)\n")
cat("  wolf_*_month_observed_summary.csv          (for month-adjusted models)\n")
cat("  wolf_*_final_spatial_block_cv_summary.csv\n")
cat("  wolf_final_run_manifest.csv\n")

###############################################################################
