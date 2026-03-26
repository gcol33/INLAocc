# =============================================================================
# occu_random.R — Random effects specification for INLA occupancy models
# =============================================================================

#' @noRd
validate_re_args <- function(type, covariate, prior, label = "Random") {
  if (type == "slope" && is.null(covariate)) {
    stop(sprintf("%s slopes require a 'covariate' argument", label))
  }
  if (is.null(prior)) {
    prior <- list(
      prec = list(prior = "pc.prec", param = c(1, 0.05))
    )
  }
  prior
}

#' Specify random effects for occupancy or detection process
#'
#' Creates a random effects specification that gets translated to INLA f() terms.
#'
#' @param type "intercept", "slope", or "iid"
#' @param group character: grouping variable name (column in covariates)
#' @param covariate character: for random slopes, the covariate name
#' @param model INLA model for the random effect ("iid", "ar1", "rw1", "rw2", "besag")
#' @param prior list with prior specification (passed to INLA hyper)
#' @param constr logical: sum-to-zero constraint
#' @param n_groups optional: number of groups (auto-detected if NULL)
#'
#' @return object of class \code{"occu_re"}
#' @export
occu_re <- function(type = c("intercept", "slope", "iid"),
                    group,
                    covariate = NULL,
                    model = "iid",
                    prior = NULL,
                    constr = FALSE,
                    n_groups = NULL) {

  type <- match.arg(type)
  prior <- validate_re_args(type, covariate, prior)

  out <- list(
    type      = type,
    group     = group,
    covariate = covariate,
    model     = model,
    prior     = prior,
    constr    = constr,
    n_groups  = n_groups
  )
  class(out) <- "occu_re"
  out
}


#' @noRd
re_to_inla_term <- function(re, prefix = "occ") {
  idx_name <- paste0(prefix, "_re_", re$group)

  if (re$type == "intercept" || re$type == "iid") {
    # f(idx, model = "iid", hyper = ...)
    sprintf(
      'f(%s, model = "%s", constr = %s, hyper = list(prec = list(prior = "%s", param = c(%s))))',
      idx_name,
      re$model,
      as.character(re$constr),
      re$prior$prec$prior,
      paste(re$prior$prec$param, collapse = ", ")
    )
  } else if (re$type == "slope") {
    # f(idx, covariate, model = "iid", ...)
    cov_name <- paste0(prefix, "_re_cov_", re$covariate)
    sprintf(
      'f(%s, %s, model = "%s", constr = %s, hyper = list(prec = list(prior = "%s", param = c(%s))))',
      idx_name,
      cov_name,
      re$model,
      as.character(re$constr),
      re$prior$prec$prior,
      paste(re$prior$prec$param, collapse = ", ")
    )
  }
}


#' @noRd
attach_re_columns <- function(df, re, group_values, covariate_values = NULL,
                              prefix = "occ") {
  idx_name <- paste0(prefix, "_re_", re$group)

  # Convert group labels to integer indices
  if (is.character(group_values) || is.factor(group_values)) {
    group_values <- as.integer(as.factor(group_values))
  }

  df[[idx_name]] <- group_values

  if (re$type == "slope" && !is.null(covariate_values)) {
    cov_name <- paste0(prefix, "_re_cov_", re$covariate)
    df[[cov_name]] <- covariate_values
  }

  df
}


#' Create a community-level (multi-species) random effect
#'
#' For multi-species models: species-specific intercepts/slopes drawn
#' from a community distribution.
#'
#' @param type "intercept" or "slope"
#' @param covariate for slopes, the covariate name
#' @param model "iid" for exchangeable species effects, "iid" with group for correlated
#' @param prior prior specification
#'
#' @return object of class \code{"occu_community_re"}
#' @export
occu_community_re <- function(type = c("intercept", "slope"),
                              covariate = NULL,
                              model = "iid",
                              prior = NULL) {
  type <- match.arg(type)
  prior <- validate_re_args(type, covariate, prior, label = "Community random")

  out <- list(
    type      = type,
    covariate = covariate,
    model     = model,
    prior     = prior
  )
  class(out) <- "occu_community_re"
  out
}


#' Resolve group variable values, handling interaction terms (e.g., "site:year")
#' @noRd
resolve_group_vals <- function(group_name, data1, data2 = NULL) {
  # Direct lookup first
  vals <- data1[[group_name]]
  if (!is.null(vals)) return(vals)
  if (!is.null(data2)) {
    vals <- data2[[group_name]]
    if (!is.null(vals)) return(vals)
  }

  # Interaction group: "a:b" → build from separate columns
  if (grepl(":", group_name, fixed = TRUE)) {
    parts <- strsplit(group_name, ":", fixed = TRUE)[[1]]
    cols <- list()
    for (p in parts) {
      v <- data1[[p]]
      if (is.null(v) && !is.null(data2)) v <- data2[[p]]
      if (is.null(v)) return(NULL)  # can't resolve a component
      cols[[p]] <- v
    }
    return(as.integer(interaction(cols, drop = TRUE)))
  }

  NULL
}


#' @noRd
build_re_components <- function(re_list, df, source_data, prefix = "occ",
                                process = c("occ", "det")) {
  process <- match.arg(process)

  if (is.null(re_list) || length(re_list) == 0) {
    return(list(formula_terms = character(0), df = df))
  }

  terms <- character(length(re_list))

  for (k in seq_along(re_list)) {
    re <- re_list[[k]]
    terms[k] <- re_to_inla_term(re, prefix = prefix)

    # Extract group values from appropriate source
    if (process == "occ") {
      group_vals <- resolve_group_vals(re$group, source_data$occ.covs, df)
      if (is.null(group_vals)) {
        group_vals <- df[[re$group]]
      }

      cov_vals <- NULL
      if (re$type == "slope") {
        cov_vals <- source_data$occ.covs[[re$covariate]]
        if (is.null(cov_vals)) cov_vals <- df[[re$covariate]]
      }

      # For augmented occ data, need to map back to original site indices
      site_indices <- df$site
      group_vals_expanded <- group_vals[site_indices]
      cov_vals_expanded <- if (!is.null(cov_vals)) cov_vals[site_indices] else NULL

      df <- attach_re_columns(df, re, group_vals_expanded, cov_vals_expanded,
                              prefix = prefix)

    } else {
      # Detection data: group values might be site-level or visit-level
      group_vals <- resolve_group_vals(re$group, df, NULL)
      if (is.null(group_vals)) {
        # Try to get from source_data at site level
        site_vals <- resolve_group_vals(re$group, source_data$occ.covs, NULL)
        if (!is.null(site_vals)) {
          group_vals <- site_vals[df$site]
        }
      }

      cov_vals <- NULL
      if (re$type == "slope") {
        cov_vals <- df[[re$covariate]]
        if (is.null(cov_vals)) {
          # Try det.covs — already in long format in df
          if (re$covariate %in% names(df)) {
            cov_vals <- df[[re$covariate]]
          }
        }
      }

      df <- attach_re_columns(df, re, group_vals, cov_vals, prefix = prefix)
    }
  }

  list(formula_terms = terms, df = df)
}
