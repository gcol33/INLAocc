# =============================================================================
# occu_data.R — Data formatting for INLA occupancy models
# =============================================================================

#' Format data for INLA occupancy models
#'
#' Accepts data in spOccupancy-compatible format and validates it.
#'
#' @param y N x J detection history matrix (0/1/NA). Rows = sites, cols = visits.
#' @param occ.covs data.frame with N rows of site-level occupancy covariates.
#'   Can also be a named list of vectors, each length N.
#' @param det.covs named list of detection covariates. Each element is either:
#'   - a vector of length N (constant across visits)
#'   - an N x J matrix (varies by visit)
#' @param coords optional N x 2 matrix of site coordinates (for spatial models)
#' @param species optional character or integer species identifier (for multi-species)
#'
#' @return An object of class \code{"occu_data"} with validated components
#'
#' @examples
#' y <- matrix(rbinom(200, 1, 0.4), nrow = 50, ncol = 4)
#' occ_covs <- data.frame(elev = rnorm(50), forest = runif(50))
#' det_covs <- list(
#'   effort = matrix(runif(200, 1, 8), 50, 4),
#'   date   = matrix(rnorm(200), 50, 4)
#' )
#' dat <- occu_format(y, occ_covs, det_covs)
#' @export
occu_format <- function(y, occ.covs = NULL, det.covs = NULL,
                        coords = NULL, species = NULL) {

  # --- Validate y ---
  if (!is.matrix(y)) {
    if (is.data.frame(y)) {
      y <- as.matrix(y)
    } else {
      stop("y must be an N x J matrix of detection histories (0/1/NA)")
    }
  }
  if (!all(y[!is.na(y)] %in% c(0L, 1L))) {
    stop("y must contain only 0, 1, or NA")
  }

  N <- nrow(y)
  J <- ncol(y)

  # --- Validate occupancy covariates ---
  if (is.null(occ.covs)) {
    occ.covs <- data.frame(.intercept = rep(1, N))
  } else if (is.list(occ.covs) && !is.data.frame(occ.covs)) {
    # Convert named list to data.frame
    lens <- vapply(occ.covs, length, integer(1))
    if (!all(lens == N)) {
      stop(sprintf(
        "All occ.covs must have length N = %d. Got lengths: %s",
        N, paste(lens, collapse = ", ")
      ))
    }
    occ.covs <- as.data.frame(occ.covs)
  }
  if (nrow(occ.covs) != N) {
    stop(sprintf("occ.covs has %d rows but y has %d sites", nrow(occ.covs), N))
  }

  # --- Validate detection covariates ---
  if (!is.null(det.covs)) {
    if (!is.list(det.covs)) {
      stop("det.covs must be a named list of vectors (length N) or matrices (N x J)")
    }
    for (nm in names(det.covs)) {
      cov <- det.covs[[nm]]
      if (is.vector(cov) && length(cov) == N) {
        # Expand site-level covariate to N x J matrix
        det.covs[[nm]] <- matrix(cov, nrow = N, ncol = J)
      } else if (is.matrix(cov)) {
        if (nrow(cov) != N || ncol(cov) != J) {
          stop(sprintf(
            "det.covs$%s is %d x %d but expected %d x %d",
            nm, nrow(cov), ncol(cov), N, J
          ))
        }
      } else {
        stop(sprintf(
          "det.covs$%s must be a vector of length %d or a %d x %d matrix",
          nm, N, N, J
        ))
      }
    }
  }

  # --- Validate coordinates ---
  if (!is.null(coords)) {
    if (!is.matrix(coords)) coords <- as.matrix(coords)
    if (nrow(coords) != N) {
      stop(sprintf("coords has %d rows but y has %d sites", nrow(coords), N))
    }
    if (ncol(coords) != 2) {
      stop("coords must have exactly 2 columns (x, y or lon, lat)")
    }
  }

  # --- Check for NAs in covariates ---
  if (!is.null(occ.covs)) {
    na_cols <- names(occ.covs)[vapply(occ.covs, function(x) any(is.na(x)), logical(1))]
    if (length(na_cols) > 0) {
      warning(sprintf(
        "NAs found in occupancy covariates: %s. Sites with NA covariates will be dropped from the model.",
        paste(na_cols, collapse = ", ")
      ))
    }
  }
  if (!is.null(det.covs)) {
    na_det <- names(det.covs)[vapply(det.covs, function(x) any(is.na(x)), logical(1))]
    if (length(na_det) > 0) {
      warning(sprintf(
        "NAs found in detection covariates: %s. Visits with NA covariates will be dropped from the model.",
        paste(na_det, collapse = ", ")
      ))
    }
  }

  # --- Derived quantities ---
  n_visits <- rowSums(!is.na(y))
  n_det    <- rowSums(y, na.rm = TRUE)
  detected <- n_det > 0

  out <- list(
    y         = y,
    occ.covs  = occ.covs,
    det.covs  = det.covs,
    coords    = coords,
    species   = species,
    N         = N,
    J         = J,
    n_visits  = n_visits,
    n_det     = n_det,
    detected  = detected,
    n_occupied = sum(detected),
    naive_occ  = mean(detected),
    naive_det  = if (sum(n_visits[detected]) > 0) {
      sum(n_det) / sum(n_visits[detected])
    } else NA_real_
  )
  class(out) <- "occu_data"
  out
}


#' Convert a data.frame to occupancy data
#'
#' Converts a long-format data.frame (one row per site-visit) into the
#' structured list format required by \code{\link{occu}}. Pipe-friendly:
#' the data.frame is the first argument.
#'
#' @param df a data.frame in long format (one row per site-visit).
#' @param y character: name of the detection column (0/1/NA).
#' @param site character: name of the site identifier column.
#' @param visit character: name of the visit/replicate column (integer
#'   or factor indicating visit number within each site).
#' @param occ.covs character vector of column names for site-level
#'   occupancy covariates. If \code{NULL} (default), all columns that
#'   are constant within each site (excluding \code{y}, \code{site},
#'   \code{visit}, and \code{det.covs}) are used.
#' @param det.covs character vector of column names for visit-level
#'   detection covariates. If \code{NULL} (default), all columns that
#'   vary within at least one site are used.
#' @param coords character vector of length 2 giving coordinate column
#'   names, or \code{NULL}.
#' @param impute how to handle NAs in site-level covariates after
#'   collapsing visit-rows.  \code{"none"} (default) leaves NAs as-is.
#'   \code{"mean"} or \code{"median"} fills remaining NAs with the
#'   column mean or median.  Only affects site-level (occupancy)
#'   covariates, not detection covariates or the detection matrix.
#'
#' @return An object of class \code{"occu_data"} ready for
#'   \code{\link{occu}}.
#'
#' @examples
#' # Long-format data
#' df <- data.frame(
#'   site = rep(1:50, each = 4),
#'   visit = rep(1:4, 50),
#'   detected = rbinom(200, 1, 0.4),
#'   elev = rep(rnorm(50), each = 4),
#'   effort = runif(200, 1, 8)
#' )
#' dat <- occu_data(df, y = "detected", site = "site", visit = "visit")
#'
#' # With pipe
#' # dat <- df |> occu_data(y = "detected", site = "site", visit = "visit")
#'
#' @export
occu_data <- function(df, y, site, visit,
                      occ.covs = NULL, det.covs = NULL,
                      coords = NULL,
                      impute = c("none", "mean", "median")) {
  impute <- match.arg(impute)

  if (!is.data.frame(df)) stop("df must be a data.frame")

  # Validate required columns exist
  required <- c(y, site, visit)
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(sprintf("Column(s) not found in df: %s", paste(missing, collapse = ", ")))
  }

  sites <- df[[site]]
  visits <- df[[visit]]
  y_vec <- df[[y]]

  # Determine unique sites (preserve original order)
  site_ids <- unique(sites)
  N <- length(site_ids)
  site_map <- match(sites, site_ids)

  # Determine visit indices
  if (is.factor(visits)) visits <- as.integer(visits)
  if (!is.numeric(visits)) {
    # Convert character/factor visit IDs to integers
    visit_ids <- unique(visits)
    visits <- match(visits, visit_ids)
  }
  J <- max(visits, na.rm = TRUE)

  # Columns to exclude from auto-detection
  meta_cols <- c(y, site, visit)
  if (!is.null(coords)) meta_cols <- c(meta_cols, coords)
  other_cols <- setdiff(names(df), meta_cols)

  # Auto-detect occ vs det covariates if not specified
  if (is.null(occ.covs) && is.null(det.covs) && length(other_cols) > 0) {
    is_site_level <- vapply(other_cols, function(nm) {
      # Check if values are constant within each site
      !any(vapply(split(df[[nm]], sites), function(v) {
        length(unique(v[!is.na(v)])) > 1L
      }, logical(1)))
    }, logical(1))
    occ.covs <- other_cols[is_site_level]
    det.covs <- other_cols[!is_site_level]
    if (length(occ.covs) == 0) occ.covs <- NULL
    if (length(det.covs) == 0) det.covs <- NULL
    message(sprintf(
      "No manual occ.covs/det.covs set. Attempted to auto-detect:\n  occ.covs (site-level):  %s\n  det.covs (visit-level): %s",
      if (!is.null(occ.covs)) paste(occ.covs, collapse = ", ") else "(none)",
      if (!is.null(det.covs)) paste(det.covs, collapse = ", ") else "(none)"
    ))
  }

  # Validate user-specified columns exist
  if (!is.null(occ.covs)) {
    bad <- setdiff(occ.covs, names(df))
    if (length(bad) > 0)
      stop(sprintf("occ.covs column(s) not found: %s", paste(bad, collapse = ", ")))
  }
  if (!is.null(det.covs)) {
    bad <- setdiff(det.covs, names(df))
    if (length(bad) > 0)
      stop(sprintf("det.covs column(s) not found: %s", paste(bad, collapse = ", ")))
  }

  # Build y matrix (N x J)
  y_mat <- matrix(NA_integer_, N, J)
  y_mat[cbind(site_map, visits)] <- as.integer(y_vec)

  # Build occ.covs data.frame (site-level: first non-NA value per site)
  occ_df <- NULL
  if (!is.null(occ.covs) && length(occ.covs) > 0) {
    occ_df <- data.frame(row.names = seq_len(N))
    for (nm in occ.covs) {
      vals <- df[[nm]]
      # Check that values are actually constant within each site
      n_varying <- sum(vapply(split(vals, sites), function(v) {
        v <- v[!is.na(v)]
        length(unique(v)) > 1L
      }, logical(1)))
      if (n_varying > 0) {
        warning(sprintf(
          "'%s' is listed as a site-level covariate but varies within %d site(s). Using first non-NA value per site.",
          nm, n_varying
        ))
      }
      occ_df[[nm]] <- vapply(site_ids, function(sid) {
        v <- vals[sites == sid]
        v <- v[!is.na(v)]
        if (length(v) > 0) v[1] else NA_real_
      }, numeric(1))
    }

    # Optionally impute remaining NAs
    if (impute != "none") {
      fn <- if (impute == "mean") mean else median
      for (nm in occ.covs) {
        col <- occ_df[[nm]]
        na_idx <- is.na(col)
        if (any(na_idx) && !all(na_idx)) {
          fill <- fn(col, na.rm = TRUE)
          occ_df[[nm]][na_idx] <- fill
          message(sprintf("Imputed %d NA(s) in '%s' with %s = %.4f",
                          sum(na_idx), nm, impute, fill))
        }
      }
    }
  }

  # Build det.covs list (each element is N x J matrix)
  det_list <- NULL
  if (!is.null(det.covs) && length(det.covs) > 0) {
    det_list <- lapply(setNames(det.covs, det.covs), function(nm) {
      m <- matrix(NA_real_, N, J)
      m[cbind(site_map, visits)] <- as.numeric(df[[nm]])
      m
    })
  }

  # Build coords matrix
  coords_mat <- NULL
  if (!is.null(coords)) {
    bad <- setdiff(coords, names(df))
    if (length(bad) > 0)
      stop(sprintf("coords column(s) not found: %s", paste(bad, collapse = ", ")))
    first_row <- !duplicated(sites)
    coords_mat <- as.matrix(df[first_row, coords, drop = FALSE])
    rownames(coords_mat) <- NULL
  }

  occu_format(y_mat, occ_df, det_list, coords_mat)
}


#' Print method for occu_data
#' @param x an \code{occu_data} object
#' @param ... additional arguments (ignored)
#' @return The \code{occu_data} object \code{x}, returned invisibly.
#' @export
print.occu_data <- function(x, ...) {
  cat("Occupancy data (occu_data)\n")
  cat(sprintf("  Sites:    %d\n", x$N))
  cat(sprintf("  Visits:   %d (max)\n", x$J))
  cat(sprintf("  Detected: %d / %d (naive psi = %.3f)\n",
              x$n_occupied, x$N, x$naive_occ))
  if (!is.na(x$naive_det)) {
    cat(sprintf("  Naive p:  %.3f\n", x$naive_det))
  }
  cat(sprintf("  Occ covs: %s\n",
              paste(names(x$occ.covs), collapse = ", ")))
  if (!is.null(x$det.covs)) {
    cat(sprintf("  Det covs: %s\n",
                paste(names(x$det.covs), collapse = ", ")))
  }
  if (!is.null(x$coords)) {
    cat("  Spatial:  coordinates provided\n")
  }
  invisible(x)
}


#' Summary statistics for occupancy detection histories
#'
#' Computes descriptive statistics for the detection matrix before model
#' fitting: detection frequencies, visit completeness, per-visit detection
#' rates, and spatial patterns (if coordinates are available).
#'
#' @param object an \code{occu_data} object or a raw data list with
#'   component \code{y}
#' @param ... ignored
#'
#' @return A list of class \code{"occu_data_summary"} with:
#'   \describe{
#'     \item{N}{number of sites}
#'     \item{J}{max visits per site}
#'     \item{n_obs}{total non-NA observations}
#'     \item{n_missing}{total NAs in detection matrix}
#'     \item{naive_psi}{proportion of sites with at least one detection}
#'     \item{naive_p}{detection rate conditional on occupied sites}
#'     \item{det_freq}{table of detection counts per site (0, 1, 2, ...)}
#'     \item{visit_rate}{proportion of non-NA cells per visit column}
#'     \item{det_per_visit}{detection rate per visit column}
#'     \item{has_coords}{logical}
#'   }
#' @export
summary.occu_data <- function(object, ...) {
  y <- object$y
  N <- nrow(y)
  J <- ncol(y)

  not_na <- !is.na(y)
  n_obs <- sum(not_na)
  n_missing <- sum(is.na(y))

  # Per-site detection count
  det_count <- rowSums(y == 1, na.rm = TRUE)
  n_visits  <- rowSums(not_na)
  detected  <- det_count > 0

  naive_psi <- mean(detected)
  naive_p   <- if (sum(n_visits[detected]) > 0) {
    sum(det_count[detected]) / sum(n_visits[detected])
  } else NA_real_

  # Detection frequency table
  det_freq <- table(factor(det_count, levels = 0:J))

  # Per-visit stats
  visit_rate    <- colMeans(not_na)
  det_per_visit <- colMeans(y == 1, na.rm = TRUE)

  out <- list(
    N             = N,
    J             = J,
    n_obs         = n_obs,
    n_missing     = n_missing,
    naive_psi     = naive_psi,
    naive_p       = naive_p,
    det_freq      = det_freq,
    det_count     = det_count,
    n_visits      = n_visits,
    visit_rate    = visit_rate,
    det_per_visit = det_per_visit,
    has_coords    = !is.null(object$coords)
  )
  class(out) <- "occu_data_summary"
  out
}


#' @export
print.occu_data_summary <- function(x, ...) {
  cat("Occupancy data summary\n")
  cat(sprintf("  Sites: %d | Max visits: %d\n", x$N, x$J))
  cat(sprintf("  Observations: %d | Missing: %d (%.1f%%)\n",
              x$n_obs, x$n_missing,
              100 * x$n_missing / (x$n_obs + x$n_missing)))
  cat(sprintf("  Naive occupancy: %.3f (%d / %d sites)\n",
              x$naive_psi, sum(x$det_freq[-1]), x$N))
  if (!is.na(x$naive_p)) {
    cat(sprintf("  Naive detection: %.3f\n", x$naive_p))
  }
  cat("\n  Detection frequency (detections per site):\n")
  df <- as.data.frame(x$det_freq)
  names(df) <- c("detections", "sites")
  print(df, row.names = FALSE)
  cat(sprintf("\n  Per-visit detection rate: %s\n",
              paste(sprintf("V%d=%.2f", seq_along(x$det_per_visit),
                            x$det_per_visit), collapse = "  ")))
  if (x$has_coords) cat("  Coordinates: available\n")
  invisible(x)
}


#' Plot detection history patterns
#'
#' Produces a 2x2 panel: detection frequency histogram, per-visit detection
#' rates, visit completeness, and a spatial detection map (if coordinates
#' are available).
#'
#' @param x an \code{occu_data} object
#' @param ... ignored
#'
#' @return Invisible \code{NULL}.
#' @export
plot.occu_data <- function(x, ...) {
  y <- x$y
  N <- nrow(y)
  J <- ncol(y)

  det_count <- rowSums(y == 1, na.rm = TRUE)
  n_visits  <- rowSums(!is.na(y))
  detected  <- det_count > 0

  has_coords <- !is.null(x$coords)
  old_par <- par(mfrow = if (has_coords) c(2, 2) else c(1, 3),
                 mar = c(4, 4, 2.5, 1))
  on.exit(par(old_par))

  # Panel 1: detection frequency histogram
  hist(det_count, breaks = seq(-0.5, max(det_count) + 0.5, by = 1),
       main = "Detections per site", xlab = "Number of detections",
       col = "grey80", border = "grey50")

  # Panel 2: per-visit detection rate
  det_per_visit <- colMeans(y == 1, na.rm = TRUE)
  barplot(det_per_visit, names.arg = paste0("V", seq_len(J)),
          main = "Detection rate by visit",
          ylab = "P(detect)", col = "steelblue", ylim = c(0, max(det_per_visit) * 1.2))

  # Panel 3: visit completeness
  completeness <- colMeans(!is.na(y))
  barplot(completeness, names.arg = paste0("V", seq_len(J)),
          main = "Visit completeness",
          ylab = "Proportion surveyed", col = "darkseagreen",
          ylim = c(0, 1))

  # Panel 4: spatial map (if coords)
  if (has_coords) {
    cols <- ifelse(detected, "tomato", "grey70")
    plot(x$coords[, 1], x$coords[, 2],
         col = cols, pch = 19,
         cex = 0.5 + det_count / max(max(det_count), 1),
         xlab = "X", ylab = "Y",
         main = "Detection map")
    legend("topright",
           legend = c("Detected", "Not detected"),
           col = c("tomato", "grey70"), pch = 19, cex = 0.8)
  }

  invisible(NULL)
}


#' Create multi-species occupancy data
#'
#' @param y_list named list of N x J detection matrices (one per species)
#' @param occ.covs site-level covariates (shared across species)
#' @param det.covs detection covariates (shared or species-specific)
#' @param coords optional coordinates
#'
#' @return object of class \code{"occu_data_ms"}
#' @export
occu_format_ms <- function(y_list, occ.covs = NULL, det.covs = NULL,
                           coords = NULL) {
  # Accept 3D array (species x sites x visits) — spOccupancy format
  if (is.array(y_list) && length(dim(y_list)) == 3) {
    arr <- y_list
    n_sp <- dim(arr)[1]
    sp_names <- dimnames(arr)[[1]] %||% paste0("sp", seq_len(n_sp))
    y_list <- lapply(seq_len(n_sp), function(s) arr[s, , ])
    names(y_list) <- sp_names
  }

  if (!is.list(y_list)) {
    stop("y_list must be a named list of detection matrices (one per species)")
  }

  species_names <- names(y_list)
  if (is.null(species_names)) {
    species_names <- paste0("sp", seq_along(y_list))
  }

  # Validate all species have same dimensions
  dims <- lapply(y_list, dim)
  ref <- dims[[1]]
  for (i in seq_along(dims)) {
    if (!identical(dims[[i]], ref)) {
      stop(sprintf(
        "All species must have same dimensions. Species 1 is %d x %d but species %d is %d x %d",
        ref[1], ref[2], i, dims[[i]][1], dims[[i]][2]
      ))
    }
  }

  # Format each species
  species_data <- lapply(seq_along(y_list), function(i) {
    occu_format(
      y_list[[i]], occ.covs, det.covs, coords,
      species = species_names[i]
    )
  })
  names(species_data) <- species_names

  out <- list(
    species_data  = species_data,
    species_names = species_names,
    n_species     = length(y_list),
    N             = ref[1],
    J             = ref[2],
    occ.covs      = occ.covs,
    det.covs      = det.covs,
    coords        = coords
  )
  class(out) <- "occu_data_ms"
  out
}


#' Simulate occupancy data for testing
#'
#' @param N number of sites
#' @param J number of visits per site
#' @param n_occ_covs number of occupancy covariates
#' @param n_det_covs number of detection covariates
#' @param beta_occ occupancy coefficients (intercept + slopes)
#' @param beta_det detection coefficients (intercept + slopes)
#' @param random_occ_sd SD of site-level random intercept on occupancy (0 = none)
#' @param random_det_sd SD of site-level random intercept on detection (0 = none)
#' @param spatial_range spatial range parameter (NULL = no spatial effect)
#' @param spatial_var spatial variance
#' @param seed random seed
#'
#' @return list with occu_data object and true parameter values
#' @export
simulate_occu <- function(N = 100, J = 4,
                          n_occ_covs = 2, n_det_covs = 1,
                          beta_occ = NULL, beta_det = NULL,
                          random_occ_sd = 0, random_det_sd = 0,
                          spatial_range = NULL, spatial_var = 1,
                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Default coefficients
  if (is.null(beta_occ)) {
    beta_occ <- c(0.5, seq(-0.8, 0.8, length.out = n_occ_covs))
  }
  if (is.null(beta_det)) {
    beta_det <- c(0, seq(-0.5, 0.5, length.out = n_det_covs))
  }

  # Covariates
  X_occ <- matrix(rnorm(N * n_occ_covs), N, n_occ_covs)
  colnames(X_occ) <- paste0("occ_x", seq_len(n_occ_covs))

  X_det_list <- lapply(seq_len(n_det_covs), function(k) {
    matrix(rnorm(N * J), N, J)
  })
  names(X_det_list) <- paste0("det_x", seq_len(n_det_covs))

  # Coordinates (for spatial models)
  coords <- cbind(x = runif(N), y = runif(N))

  # Occupancy linear predictor
  eta_occ <- beta_occ[1] + X_occ %*% beta_occ[-1]

  # Random intercept on occupancy
  re_occ <- rep(0, N)
  if (random_occ_sd > 0) {
    re_occ <- rnorm(N, 0, random_occ_sd)
    eta_occ <- eta_occ + re_occ
  }

  # Spatial effect on occupancy
  spatial_w <- rep(0, N)
  if (!is.null(spatial_range)) {
    D <- as.matrix(dist(coords))
    Sigma <- spatial_var * exp(-D / spatial_range)
    spatial_w <- MASS::mvrnorm(1, rep(0, N), Sigma)
    eta_occ <- eta_occ + spatial_w
  }

  psi <- as.vector(1 / (1 + exp(-eta_occ)))
  z   <- rbinom(N, 1, psi)

  # Detection linear predictor
  p_mat <- matrix(NA, N, J)
  re_det <- rep(0, N)
  if (random_det_sd > 0) {
    re_det <- rnorm(N, 0, random_det_sd)
  }

  for (i in seq_len(N)) {
    for (j in seq_len(J)) {
      eta_det <- beta_det[1] + re_det[i]
      for (k in seq_len(n_det_covs)) {
        eta_det <- eta_det + beta_det[k + 1] * X_det_list[[k]][i, j]
      }
      p_mat[i, j] <- 1 / (1 + exp(-eta_det))
    }
  }

  # Generate detections
  y <- matrix(NA, N, J)
  for (i in seq_len(N)) {
    for (j in seq_len(J)) {
      y[i, j] <- rbinom(1, 1, z[i] * p_mat[i, j])
    }
  }

  # Format
  dat <- occu_format(
    y        = y,
    occ.covs = as.data.frame(X_occ),
    det.covs = X_det_list,
    coords   = coords
  )

  list(
    data = dat,
    truth = list(
      z          = z,
      psi        = psi,
      p          = p_mat,
      beta_occ   = beta_occ,
      beta_det   = beta_det,
      re_occ     = re_occ,
      re_det     = re_det,
      spatial_w  = spatial_w
    )
  )
}


#' Simulate multi-species occupancy data (cf. simMsOcc)
#'
#' @param N number of sites
#' @param J number of visits
#' @param n_species number of species
#' @param n_occ_covs number of occupancy covariates
#' @param n_det_covs number of detection covariates
#' @param beta_comm_mean community mean for occupancy coefficients
#' @param beta_comm_sd community SD for occupancy coefficients
#' @param alpha_comm_mean community mean for detection coefficients
#' @param alpha_comm_sd community SD for detection coefficients
#' @param spatial_range spatial range (NULL = no spatial effect)
#' @param spatial_var spatial variance
#' @param seed random seed
#'
#' @return list with occu_data_ms object and true parameters
#' @export
simMsOcc <- function(N = 100, J = 4, n_species = 10,
                     n_occ_covs = 1, n_det_covs = 1,
                     beta_comm_mean = NULL, beta_comm_sd = NULL,
                     alpha_comm_mean = NULL, alpha_comm_sd = NULL,
                     spatial_range = NULL, spatial_var = 1,
                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_beta <- 1 + n_occ_covs
  n_alpha <- 1 + n_det_covs

  if (is.null(beta_comm_mean))  beta_comm_mean  <- c(0.3, rep(0, n_occ_covs))
  if (is.null(beta_comm_sd))    beta_comm_sd    <- c(1.0, rep(0.5, n_occ_covs))
  if (is.null(alpha_comm_mean)) alpha_comm_mean <- c(0.0, rep(0, n_det_covs))
  if (is.null(alpha_comm_sd))   alpha_comm_sd   <- c(0.8, rep(0.3, n_det_covs))

  # Covariates (shared across species)
  X_occ <- matrix(rnorm(N * n_occ_covs), N, n_occ_covs)
  colnames(X_occ) <- paste0("occ_x", seq_len(n_occ_covs))

  X_det <- lapply(seq_len(n_det_covs), function(k) matrix(rnorm(N * J), N, J))
  names(X_det) <- paste0("det_x", seq_len(n_det_covs))

  coords <- cbind(x = runif(N), y = runif(N))

  # Spatial field (shared)
  spatial_w <- rep(0, N)
  if (!is.null(spatial_range)) {
    D <- as.matrix(dist(coords))
    Sigma <- spatial_var * exp(-D / spatial_range)
    spatial_w <- MASS::mvrnorm(1, rep(0, N), Sigma)
  }

  # Species-specific parameters
  betas  <- matrix(NA, n_species, n_beta)
  alphas <- matrix(NA, n_species, n_alpha)
  for (s in seq_len(n_species)) {
    betas[s, ]  <- rnorm(n_beta,  beta_comm_mean,  beta_comm_sd)
    alphas[s, ] <- rnorm(n_alpha, alpha_comm_mean, alpha_comm_sd)
  }

  # Generate data per species
  y_array <- array(NA, dim = c(n_species, N, J))
  z_mat   <- matrix(NA, n_species, N)
  psi_mat <- matrix(NA, n_species, N)
  p_array <- array(NA, dim = c(n_species, N, J))

  sp_names <- paste0("sp", seq_len(n_species))

  for (s in seq_len(n_species)) {
    eta_occ <- betas[s, 1] + X_occ %*% betas[s, -1] + spatial_w
    psi_s <- as.vector(1 / (1 + exp(-eta_occ)))
    z_s   <- rbinom(N, 1, psi_s)

    psi_mat[s, ] <- psi_s
    z_mat[s, ]   <- z_s

    for (i in seq_len(N)) {
      for (j in seq_len(J)) {
        eta_det <- alphas[s, 1]
        for (k in seq_len(n_det_covs)) {
          eta_det <- eta_det + alphas[s, k + 1] * X_det[[k]][i, j]
        }
        p_ij <- 1 / (1 + exp(-eta_det))
        p_array[s, i, j] <- p_ij
        y_array[s, i, j] <- rbinom(1, 1, z_s[i] * p_ij)
      }
    }
  }

  dimnames(y_array)[[1]] <- sp_names

  # Format as occu_data_ms
  y_list <- lapply(seq_len(n_species), function(s) y_array[s, , ])
  names(y_list) <- sp_names

  ms_data <- occu_format_ms(y_list, as.data.frame(X_occ), X_det, coords)

  list(
    data  = ms_data,
    truth = list(
      z              = z_mat,
      psi            = psi_mat,
      p              = p_array,
      betas          = betas,
      alphas         = alphas,
      beta_comm_mean = beta_comm_mean,
      beta_comm_sd   = beta_comm_sd,
      alpha_comm_mean = alpha_comm_mean,
      alpha_comm_sd   = alpha_comm_sd,
      spatial_w      = spatial_w
    )
  )
}


#' Simulate multi-season occupancy data (cf. simTOcc)
#'
#' @param N number of sites
#' @param J number of visits per season
#' @param n_seasons number of primary periods
#' @param beta_occ occupancy coefficients
#' @param beta_det detection coefficients
#' @param ar1 logical: use AR(1) temporal correlation on occupancy
#' @param rho AR(1) correlation parameter
#' @param sigma_t temporal innovation SD
#' @param seed random seed
#'
#' @return list with data suitable for temporal_occu_inla and true values
#' @export
simTOcc <- function(N = 100, J = 4, n_seasons = 5,
                    beta_occ = c(0.5, 0.3), beta_det = c(0, -0.3),
                    ar1 = TRUE, rho = 0.7, sigma_t = 0.5,
                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_occ_covs <- length(beta_occ) - 1
  n_det_covs <- length(beta_det) - 1

  coords <- cbind(x = runif(N), y = runif(N))

  # Covariates (time-invariant and time-varying)
  X_occ_static <- if (n_occ_covs > 0) {
    m <- matrix(rnorm(N * n_occ_covs), N, n_occ_covs)
    colnames(m) <- paste0("occ_x", seq_len(n_occ_covs))
    m
  } else NULL

  X_det_list <- if (n_det_covs > 0) {
    lapply(seq_len(n_det_covs), function(k) {
      array(rnorm(N * n_seasons * J), dim = c(N, n_seasons, J))
    })
  } else NULL
  if (!is.null(X_det_list)) names(X_det_list) <- paste0("det_x", seq_len(n_det_covs))

  # Temporal random effect
  temporal_re <- matrix(0, N, n_seasons)
  if (ar1) {
    for (i in seq_len(N)) {
      temporal_re[i, 1] <- rnorm(1, 0, sigma_t / sqrt(1 - rho^2))
      for (t in 2:n_seasons) {
        temporal_re[i, t] <- rho * temporal_re[i, t - 1] + rnorm(1, 0, sigma_t)
      }
    }
  }

  # Generate data
  y_array <- array(NA, dim = c(N, n_seasons, J))
  z_mat   <- matrix(NA, N, n_seasons)
  psi_mat <- matrix(NA, N, n_seasons)
  p_array <- array(NA, dim = c(N, n_seasons, J))

  for (t in seq_len(n_seasons)) {
    eta_occ <- beta_occ[1] + temporal_re[, t]
    if (!is.null(X_occ_static)) {
      eta_occ <- eta_occ + X_occ_static %*% beta_occ[-1]
    }

    psi_t <- as.vector(1 / (1 + exp(-eta_occ)))
    z_t   <- rbinom(N, 1, psi_t)

    psi_mat[, t] <- psi_t
    z_mat[, t]   <- z_t

    for (i in seq_len(N)) {
      for (j in seq_len(J)) {
        eta_det <- beta_det[1]
        if (!is.null(X_det_list)) {
          for (k in seq_len(n_det_covs)) {
            eta_det <- eta_det + beta_det[k + 1] * X_det_list[[k]][i, t, j]
          }
        }
        p_ij <- 1 / (1 + exp(-eta_det))
        p_array[i, t, j] <- p_ij
        y_array[i, t, j] <- rbinom(1, 1, z_t[i] * p_ij)
      }
    }
  }

  # Format covariates
  occ_covs <- if (!is.null(X_occ_static)) {
    as.list(as.data.frame(X_occ_static))
  } else list()

  det_covs <- if (!is.null(X_det_list)) {
    lapply(X_det_list, identity)  # already arrays
  } else list()

  list(
    data = list(
      y        = y_array,
      occ.covs = occ_covs,
      det.covs = det_covs,
      coords   = coords
    ),
    truth = list(
      z          = z_mat,
      psi        = psi_mat,
      p          = p_array,
      beta_occ   = beta_occ,
      beta_det   = beta_det,
      temporal_re = temporal_re,
      rho        = rho,
      sigma_t    = sigma_t
    )
  )
}


#' Simulate integrated (multi-source) occupancy data (cf. simIntOcc)
#'
#' @param N_total total unique sites
#' @param n_data number of data sources
#' @param J vector of visits per source
#' @param n_shared number of sites shared across sources
#' @param beta_occ occupancy coefficients
#' @param beta_det list of detection coefficients per source
#' @param seed random seed
#'
#' @return list with data for intOccu_inla and true values
#' @export
simIntOcc <- function(N_total = 150, n_data = 2, J = c(4, 3),
                      n_shared = 20,
                      beta_occ = c(0.5, 0.3),
                      beta_det = list(c(0.2, -0.4), c(-0.1, 0.3)),
                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_occ_covs <- length(beta_occ) - 1
  coords <- cbind(x = runif(N_total), y = runif(N_total))

  # Covariates
  X_occ <- matrix(rnorm(N_total * n_occ_covs), N_total, n_occ_covs)
  colnames(X_occ) <- paste0("occ_x", seq_len(n_occ_covs))

  # True occupancy
  eta_occ <- beta_occ[1] + X_occ %*% beta_occ[-1]
  psi <- as.vector(1 / (1 + exp(-eta_occ)))
  z   <- rbinom(N_total, 1, psi)

  # Assign sites to data sources
  all_sites <- seq_len(N_total)
  shared <- sample(all_sites, n_shared)
  remaining <- setdiff(all_sites, shared)
  per_source <- (N_total - n_shared) %/% n_data

  sites_list <- list()
  y_list     <- list()
  det_covs_list <- list()

  start <- 1
  for (d in seq_len(n_data)) {
    if (d < n_data) {
      unique_d <- remaining[start:(start + per_source - 1)]
      start <- start + per_source
    } else {
      unique_d <- remaining[start:length(remaining)]
    }
    sites_d <- sort(c(shared, unique_d))
    sites_list[[d]] <- sites_d

    N_d <- length(sites_d)
    J_d <- J[d]
    n_det_covs_d <- length(beta_det[[d]]) - 1

    det_covs_d <- list()
    if (n_det_covs_d > 0) {
      for (k in seq_len(n_det_covs_d)) {
        det_covs_d[[paste0("det_x", k)]] <- matrix(rnorm(N_d * J_d), N_d, J_d)
      }
    }

    y_d <- matrix(NA, N_d, J_d)
    for (i in seq_len(N_d)) {
      si <- sites_d[i]
      for (j in seq_len(J_d)) {
        eta_det <- beta_det[[d]][1]
        for (k in seq_len(n_det_covs_d)) {
          eta_det <- eta_det + beta_det[[d]][k + 1] * det_covs_d[[k]][i, j]
        }
        p_ij <- 1 / (1 + exp(-eta_det))
        y_d[i, j] <- rbinom(1, 1, z[si] * p_ij)
      }
    }

    y_list[[d]] <- y_d
    det_covs_list[[d]] <- det_covs_d
  }

  list(
    data = list(
      y        = y_list,
      occ.covs = as.data.frame(X_occ),
      det.covs = det_covs_list,
      sites    = sites_list,
      coords   = coords
    ),
    truth = list(
      z        = z,
      psi      = psi,
      beta_occ = beta_occ,
      beta_det = beta_det
    )
  )
}


#' Simulate temporal multi-species occupancy data (cf. simTMsOcc)
#'
#' @param N number of sites
#' @param J number of visits per season
#' @param n_species number of species
#' @param n_seasons number of primary periods
#' @param n_occ_covs number of occupancy covariates
#' @param n_det_covs number of detection covariates
#' @param beta_comm_mean community mean for occupancy coefficients
#' @param beta_comm_sd community SD for occupancy coefficients
#' @param alpha_comm_mean community mean for detection coefficients
#' @param alpha_comm_sd community SD for detection coefficients
#' @param ar1 logical: use AR(1) temporal correlation
#' @param rho AR(1) correlation parameter
#' @param sigma_t temporal innovation SD
#' @param seed random seed
#'
#' @return list with 4D data and true parameters
#' @export
simTMsOcc <- function(N = 100, J = 4, n_species = 10, n_seasons = 5,
                       n_occ_covs = 1, n_det_covs = 1,
                       beta_comm_mean = NULL, beta_comm_sd = NULL,
                       alpha_comm_mean = NULL, alpha_comm_sd = NULL,
                       ar1 = TRUE, rho = 0.7, sigma_t = 0.5,
                       seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_beta <- 1 + n_occ_covs
  n_alpha <- 1 + n_det_covs

  if (is.null(beta_comm_mean))  beta_comm_mean  <- c(0.3, rep(0, n_occ_covs))
  if (is.null(beta_comm_sd))    beta_comm_sd    <- c(1.0, rep(0.5, n_occ_covs))
  if (is.null(alpha_comm_mean)) alpha_comm_mean <- c(0.0, rep(0, n_det_covs))
  if (is.null(alpha_comm_sd))   alpha_comm_sd   <- c(0.8, rep(0.3, n_det_covs))

  X_occ <- matrix(rnorm(N * n_occ_covs), N, n_occ_covs)
  colnames(X_occ) <- paste0("occ_x", seq_len(n_occ_covs))
  coords <- cbind(x = runif(N), y = runif(N))

  betas  <- matrix(rnorm(n_species * n_beta, beta_comm_mean, beta_comm_sd),
                    n_species, n_beta, byrow = TRUE)
  alphas <- matrix(rnorm(n_species * n_alpha, alpha_comm_mean, alpha_comm_sd),
                    n_species, n_alpha, byrow = TRUE)
  sp_names <- paste0("sp", seq_len(n_species))

  # Temporal RE per species
  temporal_re <- array(0, dim = c(n_species, N, n_seasons))
  if (ar1) {
    for (s in seq_len(n_species)) {
      for (i in seq_len(N)) {
        temporal_re[s, i, 1] <- rnorm(1, 0, sigma_t / sqrt(1 - rho^2))
        for (t in 2:n_seasons) {
          temporal_re[s, i, t] <- rho * temporal_re[s, i, t - 1] + rnorm(1, 0, sigma_t)
        }
      }
    }
  }

  # Generate data: 4D array (species x sites x seasons x visits)
  y_array <- array(NA, dim = c(n_species, N, n_seasons, J))
  z_array <- array(NA, dim = c(n_species, N, n_seasons))
  psi_array <- array(NA, dim = c(n_species, N, n_seasons))

  det_covs <- lapply(seq_len(n_det_covs), function(k) {
    array(rnorm(N * n_seasons * J), dim = c(N, n_seasons, J))
  })
  names(det_covs) <- paste0("det_x", seq_len(n_det_covs))

  for (s in seq_len(n_species)) {
    for (t in seq_len(n_seasons)) {
      eta_occ <- betas[s, 1] + X_occ %*% betas[s, -1] + temporal_re[s, , t]
      psi_st <- as.vector(1 / (1 + exp(-eta_occ)))
      z_st <- rbinom(N, 1, psi_st)
      psi_array[s, , t] <- psi_st
      z_array[s, , t] <- z_st

      for (i in seq_len(N)) {
        for (j in seq_len(J)) {
          eta_det <- alphas[s, 1]
          for (k in seq_len(n_det_covs)) {
            eta_det <- eta_det + alphas[s, k + 1] * det_covs[[k]][i, t, j]
          }
          p_ij <- 1 / (1 + exp(-eta_det))
          y_array[s, i, t, j] <- rbinom(1, 1, z_st[i] * p_ij)
        }
      }
    }
  }
  dimnames(y_array)[[1]] <- sp_names

  list(
    data = list(
      y        = y_array,
      occ.covs = as.data.frame(X_occ),
      det.covs = det_covs,
      coords   = coords
    ),
    truth = list(
      z     = z_array,
      psi   = psi_array,
      betas = betas, alphas = alphas,
      beta_comm_mean = beta_comm_mean, beta_comm_sd = beta_comm_sd,
      alpha_comm_mean = alpha_comm_mean, alpha_comm_sd = alpha_comm_sd,
      temporal_re = temporal_re, rho = rho, sigma_t = sigma_t
    )
  )
}


#' Simulate integrated multi-species occupancy data (cf. simIntMsOcc)
#'
#' @param N_total total unique sites
#' @param n_species number of species
#' @param n_data number of data sources
#' @param J vector of visits per source
#' @param n_shared number of sites shared across sources
#' @param beta_comm_mean community mean for occupancy coefficients
#' @param beta_comm_sd community SD for occupancy coefficients
#' @param beta_det list of detection coefficient vectors (one per source)
#' @param seed random seed
#'
#' @return list with data for integrated multi-species models and true values
#' @export
simIntMsOcc <- function(N_total = 150, n_species = 10, n_data = 2,
                         J = c(4, 3), n_shared = 20,
                         beta_comm_mean = NULL, beta_comm_sd = NULL,
                         beta_det = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n_occ_covs <- 1
  n_beta <- 1 + n_occ_covs

  if (is.null(beta_comm_mean)) beta_comm_mean <- c(0.3, 0.2)
  if (is.null(beta_comm_sd))   beta_comm_sd   <- c(1.0, 0.5)
  if (is.null(beta_det)) {
    beta_det <- lapply(seq_len(n_data), function(d) c(0.2 * d, -0.3))
  }

  X_occ <- matrix(rnorm(N_total * n_occ_covs), N_total, n_occ_covs)
  colnames(X_occ) <- paste0("occ_x", seq_len(n_occ_covs))
  coords <- cbind(x = runif(N_total), y = runif(N_total))

  betas <- matrix(rnorm(n_species * n_beta, beta_comm_mean, beta_comm_sd),
                   n_species, n_beta, byrow = TRUE)
  sp_names <- paste0("sp", seq_len(n_species))

  # Site assignment
  all_sites <- seq_len(N_total)
  shared <- sample(all_sites, n_shared)
  remaining <- setdiff(all_sites, shared)
  per_source <- (N_total - n_shared) %/% n_data

  sites_list <- list()
  start <- 1
  for (d in seq_len(n_data)) {
    if (d < n_data) {
      unique_d <- remaining[start:(start + per_source - 1)]
      start <- start + per_source
    } else {
      unique_d <- remaining[start:length(remaining)]
    }
    sites_list[[d]] <- sort(c(shared, unique_d))
  }

  # True occupancy per species
  z_mat <- matrix(NA, n_species, N_total)
  psi_mat <- matrix(NA, n_species, N_total)
  for (s in seq_len(n_species)) {
    eta <- betas[s, 1] + X_occ %*% betas[s, -1]
    psi_mat[s, ] <- as.vector(1 / (1 + exp(-eta)))
    z_mat[s, ] <- rbinom(N_total, 1, psi_mat[s, ])
  }

  # Generate per-source detection data as 3D arrays (species x sites_d x J_d)
  y_list <- list()
  det_covs_list <- list()
  for (d in seq_len(n_data)) {
    sites_d <- sites_list[[d]]
    N_d <- length(sites_d)
    J_d <- J[d]
    n_det_covs_d <- length(beta_det[[d]]) - 1

    det_covs_d <- list()
    if (n_det_covs_d > 0) {
      for (k in seq_len(n_det_covs_d)) {
        det_covs_d[[paste0("det_x", k)]] <- matrix(rnorm(N_d * J_d), N_d, J_d)
      }
    }

    y_d <- array(NA, dim = c(n_species, N_d, J_d))
    for (s in seq_len(n_species)) {
      for (i in seq_len(N_d)) {
        si <- sites_d[i]
        for (j in seq_len(J_d)) {
          eta_det <- beta_det[[d]][1]
          for (k in seq_len(n_det_covs_d)) {
            eta_det <- eta_det + beta_det[[d]][k + 1] * det_covs_d[[k]][i, j]
          }
          p_ij <- 1 / (1 + exp(-eta_det))
          y_d[s, i, j] <- rbinom(1, 1, z_mat[s, si] * p_ij)
        }
      }
    }
    dimnames(y_d)[[1]] <- sp_names
    y_list[[d]] <- y_d
    det_covs_list[[d]] <- det_covs_d
  }

  list(
    data = list(
      y        = y_list,
      occ.covs = as.data.frame(X_occ),
      det.covs = det_covs_list,
      sites    = sites_list,
      coords   = coords
    ),
    truth = list(
      z = z_mat, psi = psi_mat, betas = betas,
      beta_comm_mean = beta_comm_mean, beta_comm_sd = beta_comm_sd,
      beta_det = beta_det
    )
  )
}
