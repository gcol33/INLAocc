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

  # Build mesh (suppressWarnings: inla.mesh.2d deprecated in favour of fmesher)
  mesh <- suppressWarnings(INLA::inla.mesh.2d(
    loc      = coords,
    max.edge = max.edge,
    cutoff   = cutoff,
    offset   = offset
  ))

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


#' Create an areal spatial effect for occupancy models (CAR/BYM2)
#'
#' For data on grids, administrative units, or any areal structure. Uses
#' INLA's \code{besag} or \code{bym2} model with an adjacency graph.
#'
#' @param adj Adjacency structure: a symmetric matrix, a \code{nb} object
#'   (from \code{spdep}), or a path to an INLA graph file.
#' @param model \code{"bym2"} (default, recommended) or \code{"besag"}.
#'   BYM2 separates structured (spatial) and unstructured (iid) components
#'   with a mixing parameter.
#' @param prior.sigma numeric(2): PC prior on marginal SD. c(s0, p) means
#'   P(sigma > s0) = p. Default c(1, 0.05).
#' @param scale.model Logical: scale the precision matrix so the generalized
#'   variance is 1 (recommended for bym2). Default TRUE.
#'
#' @return object of class \code{"occu_areal"} with graph and model components
#' @export
occu_areal <- function(adj,
                        model = c("bym2", "besag"),
                        prior.sigma = c(1, 0.05),
                        scale.model = TRUE) {
  check_inla()
  model <- match.arg(model)

  # Convert adjacency input to INLA graph
  if (is.character(adj) && file.exists(adj)) {
    # Path to INLA graph file
    graph <- adj
    # Read to determine N
    g <- INLA::inla.read.graph(adj)
    N <- g$n
  } else if (is.matrix(adj)) {
    # Adjacency matrix → INLA graph via temp file
    if (nrow(adj) != ncol(adj))
      stop("Adjacency matrix must be square")
    if (!isSymmetric(adj))
      stop("Adjacency matrix must be symmetric")
    N <- nrow(adj)
    graph_file <- tempfile(fileext = ".graph")
    # Write INLA graph format (vectorized neighbor extraction)
    lines <- character(N + 1)
    lines[1] <- as.character(N)
    seq_N <- seq_len(N)
    for (i in seq_N) {
      neighbors <- seq_N[adj[i, ] > 0 & seq_N != i]
      lines[i + 1] <- paste(i, length(neighbors), paste(neighbors, collapse = " "))
    }
    writeLines(lines, graph_file)
    graph <- graph_file
  } else if (inherits(adj, "nb")) {
    # spdep nb object → INLA graph
    graph_file <- tempfile(fileext = ".graph")
    N <- length(adj)
    lines <- character(N + 1)
    lines[1] <- as.character(N)
    for (i in seq_len(N)) {
      nb_i <- adj[[i]]
      nb_i <- nb_i[nb_i > 0]  # remove 0 (no-neighbor indicator)
      lines[i + 1] <- paste(i, length(nb_i), paste(nb_i, collapse = " "))
    }
    writeLines(lines, graph_file)
    graph <- graph_file
  } else {
    stop("adj must be an adjacency matrix, an spdep 'nb' object, or a path to an INLA graph file")
  }

  # Build INLA model
  hyper <- list(
    prec = list(prior = "pc.prec", param = prior.sigma)
  )

  out <- list(
    graph       = graph,
    model       = model,
    n_regions   = N,
    hyper       = hyper,
    scale.model = scale.model,
    prior.sigma = prior.sigma,
    type        = "areal"
  )
  class(out) <- "occu_areal"
  out
}


#' @export
print.occu_areal <- function(x, ...) {
  cat(sprintf("Areal spatial component (occu_areal, %s)\n", x$model))
  cat(sprintf("  Regions: %d\n", x$n_regions))
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
