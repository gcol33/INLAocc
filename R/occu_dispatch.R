# =============================================================================
# occu_dispatch.R — Model type resolution, registry, and data preparation
# =============================================================================

#' @noRd
resolve_model_type <- function(spatial, temporal, multispecies,
                                integrated, n.factors, svc, det.formula) {
  is_spatial    <- !is.null(spatial)
  is_temporal   <- !is.null(temporal)
  is_ms         <- isTRUE(multispecies) || identical(multispecies, "jsdm")
  is_jsdm       <- identical(multispecies, "jsdm")
  is_integrated <- isTRUE(integrated)
  is_lf         <- !is.null(n.factors)
  is_svc        <- !is.null(svc)


  # --- Validation ---
  if (is_svc && !is_spatial)
    stop("svc requires spatial (provide coordinates or an occu_spatial object)")
  if (is_jsdm && !is.null(det.formula))
    warning("det.formula ignored for JSDM models (no detection process)")
  if (is_jsdm && !is_lf)
    stop("JSDM models require n.factors")
  if (is_jsdm && is_integrated)
    stop("JSDM and integrated models cannot be combined")
  if (is_jsdm && is_temporal)
    stop("Temporal JSDM is not yet supported")
  if (is_lf && !is_ms)
    stop("n.factors requires multispecies = TRUE or multispecies = 'jsdm'")
  if (!is_jsdm && is.null(det.formula))
    stop("det.formula is required (use multispecies = 'jsdm' for models without detection)")

  # --- Build dispatch key ---
  key <- paste0(
    if (is_jsdm) "jsdm" else if (is_ms) "ms" else "ss",
    if (is_integrated) "_int" else "",
    if (is_spatial)    "_sp"  else "",
    if (is_temporal)   "_t"   else "",
    if (is_lf) "_lf" else "",
    if (is_svc)        "_svc" else ""
  )
  key
}


# --- Model registry: key → engine function name + S3 class hierarchy ---
MODEL_REGISTRY <- list(
  # Single-species
  ss           = list(engine = "engine_ss",           class = c("occu_inla", "occu_em")),
  ss_sp        = list(engine = "engine_ss_spatial",   class = c("occu_inla_spatial", "occu_inla", "occu_em")),
  ss_t         = list(engine = "engine_temporal",     class = c("occu_inla_temporal", "occu_em")),
  ss_sp_svc    = list(engine = "engine_svc",          class = c("occu_inla_svc", "occu_inla_spatial", "occu_inla", "occu_em")),
  ss_sp_t_svc  = list(engine = "engine_svc_temporal", class = c("occu_inla_svct", "occu_inla_temporal", "occu_em")),

  # Integrated
  ss_int       = list(engine = "engine_int",          class = c("occu_inla_int", "occu_inla", "occu_em")),
  ss_int_sp    = list(engine = "engine_int_spatial",   class = c("occu_inla_spint", "occu_inla_int", "occu_inla", "occu_em")),
  ss_int_sp_t  = list(engine = "engine_stint",         class = c("occu_inla_stint", "occu_inla_int", "occu_inla", "occu_em")),

  # Multi-species
  ms           = list(engine = "engine_ms",           class = c("occu_inla_ms", "occu_em")),
  ms_lf        = list(engine = "engine_ms_lf",        class = c("occu_inla_lfms", "occu_inla_ms", "occu_em")),
  ms_sp_lf     = list(engine = "engine_ms_sf",        class = c("occu_inla_sfms", "occu_inla_lfms", "occu_inla_ms", "occu_em")),
  ms_t         = list(engine = "engine_ms_temporal",   class = c("occu_inla_tms", "occu_inla_ms", "occu_em")),
  ms_sp_t      = list(engine = "engine_ms_st",         class = c("occu_inla_stms", "occu_inla_ms", "occu_em")),
  ms_int       = list(engine = "engine_ms_int",        class = c("occu_inla_intms", "occu_inla_ms", "occu_em")),
  ms_sp_t_svc  = list(engine = "engine_ms_svc_t",     class = c("occu_inla_svctms", "occu_inla_ms", "occu_em")),

  # JSDM
  jsdm_lf      = list(engine = "engine_jsdm_lf",     class = c("occu_inla_jsdm", "occu_em")),
  jsdm_sp_lf   = list(engine = "engine_jsdm_sf",     class = c("occu_inla_sfjsdm", "occu_inla_jsdm", "occu_em"))
)


#' @noRd
prepare_spatial <- function(spatial, spde.args) {
  if (is.null(spatial)) return(NULL)
  if (inherits(spatial, "occu_spatial")) return(spatial)
  if (inherits(spatial, "occu_areal")) return(spatial)
  # Treat as coordinate matrix
  if (is.matrix(spatial) || is.data.frame(spatial)) {
    return(do.call(occu_spatial, c(list(coords = as.matrix(spatial)), spde.args)))
  }
  stop("spatial must be an N x 2 coordinate matrix, an occu_spatial(), or an occu_areal() object")
}


#' @noRd
prepare_data <- function(data, multispecies, integrated, temporal) {
  is_jsdm <- identical(multispecies, "jsdm")

  if (is_jsdm) return(data)

  if (isTRUE(multispecies)) {
    # Multi-species + temporal: 4D array (species x sites x seasons x visits)
    # Pass through directly — engine_ms_temporal handles the 4D structure
    if (!is.null(temporal)) return(data)
    if (inherits(data, "occu_data_ms")) return(data)
    y_input <- data$y
    if (is.array(y_input) && length(dim(y_input)) == 3) {
      n_sp <- dim(y_input)[1]
      sp_names <- dimnames(y_input)[[1]] %||% paste0("sp", seq_len(n_sp))
      y_list <- lapply(seq_len(n_sp), function(s) y_input[s, , ])
      names(y_list) <- sp_names
    } else if (is.list(y_input)) {
      y_list <- y_input
    } else {
      stop("Multi-species data$y must be a 3D array or a named list of matrices")
    }
    return(occu_format_ms(y_list, data$occ.covs, data$det.covs, data$coords))
  }

  if (isTRUE(integrated)) {
    if (!is.list(data$y))
      stop("Integrated models require data$y to be a list of detection matrices")
    return(data)
  }

  if (!is.null(temporal)) return(data)

  # Single-species: coerce to occu_data
  if (!inherits(data, "occu_data")) return(coerce_data(data))
  data
}
