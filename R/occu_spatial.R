# =============================================================================
# occu_spatial.R — Spatial SPDE components for INLA occupancy models
# =============================================================================

#' Create SPDE mesh and spatial effect for occupancy models
#'
#' Wraps INLA's mesh construction and SPDE model setup for use in
#' occupancy models. The spatial effect is added to the occupancy
#' linear predictor.
#'
#' @param coords N x 2 matrix of coordinates
#' @param max.edge numeric vector of length 2: max triangle edge length
#'   (inner domain, outer extension). Default auto-scaled from data extent.
#' @param cutoff minimum distance between mesh nodes. Default: `max.edge[1]/5`.
#' @param offset numeric vector of length 2: inner and outer extension distances.
#' @param prior.range numeric(2): PC prior for range. c(r0, p) means P(range < r0) = p.
#' @param prior.sigma numeric(2): PC prior for marginal SD. c(s0, p) means P(sigma > s0) = p.
#' @param alpha SPDE smoothness parameter (default 2, Matern nu = 1).
#'
#' @return object of class \code{"occu_spatial"} with mesh, spde, and A matrix components
#' @export
occu_spatial <- function(coords,
                         max.edge = NULL,
                         cutoff = NULL,
                         offset = NULL,
                         prior.range = c(0.3, 0.5),
                         prior.sigma = c(1, 0.05),
                         alpha = 2) {
  check_inla()

  if (!is.matrix(coords)) coords <- as.matrix(coords)
  if (ncol(coords) != 2) stop("coords must have exactly 2 columns")

  N <- nrow(coords)

  # Auto-scale mesh parameters from data extent
  extent <- diff(range(coords[, 1])) + diff(range(coords[, 2]))
  if (is.null(max.edge)) {
    max.edge <- c(extent / 20, extent / 5)
  }
  if (is.null(cutoff)) {
    cutoff <- max.edge[1] / 5
  }
  if (is.null(offset)) {
    offset <- c(extent / 10, extent / 3)
  }

  # Build mesh
  mesh <- INLA::inla.mesh.2d(
    loc      = coords,
    max.edge = max.edge,
    cutoff   = cutoff,
    offset   = offset
  )

  # SPDE model with PC priors
  spde <- INLA::inla.spde2.pcmatern(
    mesh        = mesh,
    alpha       = alpha,
    prior.range = prior.range,
    prior.sigma = prior.sigma
  )

  # Projection matrix: maps mesh nodes -> observation locations
  A <- INLA::inla.spde.make.A(mesh = mesh, loc = coords)

  out <- list(
    mesh        = mesh,
    spde        = spde,
    A           = A,
    coords      = coords,
    n_mesh      = mesh$n,
    n_sites     = N,
    max.edge    = max.edge,
    cutoff      = cutoff,
    offset      = offset,
    prior.range = prior.range,
    prior.sigma = prior.sigma
  )
  class(out) <- "occu_spatial"
  out
}


#' Print method for occu_spatial
#' @param x an \code{occu_spatial} object
#' @param ... additional arguments (ignored)
#' @return The \code{occu_spatial} object \code{x}, returned invisibly.
#' @export
print.occu_spatial <- function(x, ...) {
  cat("SPDE spatial component (occu_spatial)\n")
  cat(sprintf("  Sites:      %d\n", x$n_sites))
  cat(sprintf("  Mesh nodes: %d\n", x$n_mesh))
  cat(sprintf("  Max edge:   [%.3f, %.3f]\n", x$max.edge[1], x$max.edge[2]))
  cat(sprintf("  PC prior range: P(range < %.2f) = %.2f\n",
              x$prior.range[1], x$prior.range[2]))
  cat(sprintf("  PC prior sigma: P(sigma > %.2f) = %.2f\n",
              x$prior.sigma[1], x$prior.sigma[2]))
  invisible(x)
}


#' @noRd
build_spatial_stack <- function(df, spatial, site_col = "site",
                                tag = "occ_spatial") {


  sites <- df[[site_col]]
  A_obs <- spatial$A[sites, , drop = FALSE]

  # For augmented occupancy data, sites may repeat — the A matrix handles this
  # because each row of A_obs maps an observation to the mesh

  INLA::inla.stack(
    data    = list(y = df$z),
    A       = list(A_obs, 1),
    effects = list(
      spatial = seq_len(ncol(A_obs)),
      df[, setdiff(names(df), c("z", site_col)), drop = FALSE]
    ),
    tag     = tag
  )
}


#' @noRd
predict_spatial_A <- function(spatial, newcoords) {

  if (!is.matrix(newcoords)) newcoords <- as.matrix(newcoords)
  INLA::inla.spde.make.A(mesh = spatial$mesh, loc = newcoords)
}


#' @noRd
extract_spatial_field <- function(fit, spatial, name = "spatial") {


  idx <- fit$summary.random[[name]]

  list(
    mean = idx$mean,
    sd   = idx$sd,
    q025 = idx$`0.025quant`,
    q500 = idx$`0.5quant`,
    q975 = idx$`0.975quant`,
    mesh = spatial$mesh,
    spde_params = if (!is.null(fit$summary.hyperpar)) {
      rows <- grep(name, rownames(fit$summary.hyperpar))
      fit$summary.hyperpar[rows, ]
    }
  )
}


#' @noRd
project_spatial_grid <- function(spatial_field, spatial,
                                 n_grid = 100, xlim = NULL, ylim = NULL) {


  coords <- spatial$coords
  if (is.null(xlim)) xlim <- range(coords[, 1])
  if (is.null(ylim)) ylim <- range(coords[, 2])

  proj <- INLA::inla.mesh.projector(
    spatial$mesh,
    xlim = xlim, ylim = ylim,
    dims = c(n_grid, n_grid)
  )

  z_mean <- INLA::inla.mesh.project(proj, spatial_field$mean)
  z_sd   <- INLA::inla.mesh.project(proj, spatial_field$sd)

  list(
    x      = proj$x,
    y      = proj$y,
    z_mean = z_mean,
    z_sd   = z_sd,
    proj   = proj
  )
}
