# =============================================================================
# occu_spatial.R — Spatial SPDE components, prediction grids, and maps
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
  mesh <- fmesher::fm_mesh_2d_inla(
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
  A <- fmesher::fm_basis(mesh, loc = coords)

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
  fmesher::fm_basis(spatial$mesh, loc = newcoords)
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

  proj <- fmesher::fm_evaluator(
    spatial$mesh,
    xlim = xlim, ylim = ylim,
    dims = c(n_grid, n_grid)
  )

  z_mean <- fmesher::fm_evaluate(proj, spatial_field$mean)
  z_sd   <- fmesher::fm_evaluate(proj, spatial_field$sd)

  list(
    x      = proj$x,
    y      = proj$y,
    z_mean = z_mean,
    z_sd   = z_sd,
    proj   = proj
  )
}


# =============================================================================
# occuMap — publication-quality spatial occupancy maps
# =============================================================================

#' @noRd
.validate_spatial_fit <- function(fit, species = NULL) {
  # Multi-species spatial
  if (inherits(fit, "occu_inla_ms") && !is.null(fit$species_fits)) {
    if (is.null(species)) {
      stop("Multi-species model: specify species = 'name' or species = index")
    }
    sp <- if (is.numeric(species)) fit$species_names[species] else species
    if (!sp %in% fit$species_names) {
      stop(sprintf("Species '%s' not found. Available: %s",
                   sp, paste(fit$species_names, collapse = ", ")))
    }
    sp_fit <- fit$species_fits[[sp]]
    if (is.null(sp_fit)) stop(sprintf("No fit for species '%s'", sp))
    return(list(
      occ_fit = sp_fit$occ_fit,
      spatial = fit$spatial,
      psi_hat = sp_fit$psi_hat,
      data    = sp_fit$data %||% fit$data,
      species = sp
    ))
  }

  # Single-species spatial
  if (is.null(fit$spatial)) {
    stop("occuMap requires a spatial model. Fit with spatial = coords.")
  }
  if (inherits(fit$spatial, "occu_areal")) {
    stop("occuMap does not support areal models. Use plot() for region-level maps.")
  }

  list(
    occ_fit = fit$occ_fit,
    spatial = fit$spatial,
    psi_hat = fit$psi_hat,
    data    = fit$data,
    species = NULL
  )
}


#' @noRd
.draw_colorbar <- function(zlim, col, label = NULL) {
  usr <- par("usr")
  n <- length(col)

  # Bar position: right margin
  bar_w <- (usr[2] - usr[1]) * 0.025
  gap   <- (usr[2] - usr[1]) * 0.015
  x0 <- usr[2] + gap
  x1 <- x0 + bar_w

  y_seq <- seq(usr[3], usr[4], length.out = n + 1L)

  old_xpd <- par(xpd = NA)
  on.exit(par(xpd = old_xpd$xpd), add = TRUE)

  for (i in seq_len(n)) {
    rect(x0, y_seq[i], x1, y_seq[i + 1L], col = col[i], border = NA)
  }
  rect(x0, usr[3], x1, usr[4], border = "black", lwd = 0.5)

  # Tick labels
  at_vals <- pretty(zlim, n = 5)
  at_vals <- at_vals[at_vals >= zlim[1] & at_vals <= zlim[2]]
  at_y <- usr[3] + (at_vals - zlim[1]) / diff(zlim) * (usr[4] - usr[3])
  segments(x1, at_y, x1 + bar_w * 0.4, at_y, lwd = 0.5)
  text(x1 + bar_w * 0.6, at_y, labels = format(at_vals, digits = 2),
       adj = 0, cex = 0.65)

  if (!is.null(label)) {
    mtext(label, side = 4, line = 4, cex = 0.75)
  }
}


#' Spatial occupancy map
#'
#' Produces a publication-quality map of predicted occupancy probability,
#' uncertainty, or the spatial random effect from a fitted spatial model.
#' Survey sites are overlaid with detection status.
#'
#' @param fit fitted spatial occupancy model (class \code{occu_inla_spatial}
#'   or multi-species with spatial component)
#' @param type \code{"psi"} (default): occupancy probability.
#'   \code{"sd"}: posterior SD.
#'   \code{"spatial"}: spatial random effect only.
#'   \code{"all"}: 2-panel (psi + sd side by side).
#' @param species for multi-species models, the species name or index
#' @param n_grid grid resolution per axis (default 100)
#' @param xlim,ylim coordinate limits (default: data extent)
#' @param col color palette vector (default chosen per \code{type})
#' @param sites logical: overlay survey sites? (default TRUE)
#' @param main plot title (default chosen per \code{type})
#' @param ... passed to \code{image()}
#'
#' @return Invisibly, a list with \code{x}, \code{y}, \code{z} (the
#'   gridded values), \code{type}, \code{zlim}, \code{coords},
#'   \code{detected}, and \code{col}.
#'
#' @export
occuMap <- function(fit,
                    type = c("psi", "sd", "spatial", "all"),
                    species = NULL,
                    n_grid = 100L,
                    xlim = NULL, ylim = NULL,
                    col = NULL,
                    sites = TRUE,
                    main = NULL,
                    ...) {

  type <- match.arg(type)

  # --- "all" mode: 2-panel ---
  if (type == "all") {
    old_par <- par(mfrow = c(1, 2))
    on.exit(par(old_par))
    r1 <- occuMap(fit, type = "psi", species = species, n_grid = n_grid,
                  xlim = xlim, ylim = ylim, col = col, sites = sites,
                  main = main, ...)
    r2 <- occuMap(fit, type = "sd", species = species, n_grid = n_grid,
                  xlim = xlim, ylim = ylim, sites = sites, ...)
    return(invisible(list(psi = r1, sd = r2)))
  }

  # --- Validate and extract ---
  v <- .validate_spatial_fit(fit, species)

  # --- Compute grid ---
  sf <- extract_spatial_field(v$occ_fit, v$spatial)
  grid <- project_spatial_grid(sf, v$spatial, n_grid = n_grid,
                                xlim = xlim, ylim = ylim)

  # --- Build the plotted quantity ---
  if (type == "psi") {
    # Spatial field + fixed-effects baseline on logit scale
    beta <- v$occ_fit$summary.fixed$mean
    z_logit <- grid$z_mean  # spatial RE on log scale
    # Add intercept + mean covariate effects
    if (length(beta) > 0) {
      z_logit <- z_logit + sum(beta)  # approximate: covariates at mean
    }
    z <- expit(z_logit)
    z_label <- expression(psi)
    default_col <- hcl.colors(100, "YlGnBu")
    default_zlim <- c(0, 1)
    default_main <- if (!is.null(v$species)) {
      sprintf("Occupancy: %s", v$species)
    } else "Occupancy probability"

  } else if (type == "sd") {
    z <- grid$z_sd
    z_label <- "SD"
    default_col <- hcl.colors(100, "YlOrRd")
    default_zlim <- range(z, na.rm = TRUE)
    default_main <- if (!is.null(v$species)) {
      sprintf("Uncertainty: %s", v$species)
    } else "Posterior SD"

  } else {
    z <- grid$z_mean
    z_label <- "RE"
    default_col <- hcl.colors(100, "Blue-Red 3")
    z_range <- max(abs(range(z, na.rm = TRUE)))
    default_zlim <- c(-z_range, z_range)  # symmetric around 0
    default_main <- if (!is.null(v$species)) {
      sprintf("Spatial RE: %s", v$species)
    } else "Spatial random effect"
  }

  col  <- col %||% default_col
  zlim <- default_zlim
  main <- main %||% default_main

  # --- Plot ---
  old_mar <- par(mar = c(4, 4, 3, 6))
  on.exit(par(mar = old_mar$mar), add = TRUE)

  image(grid$x, grid$y, z,
        col = col, zlim = zlim, asp = 1,
        xlab = "X", ylab = "Y", main = main, ...)

  # Site overlay
  coords <- v$spatial$coords
  detected <- NULL
  if (sites && !is.null(coords) && !is.null(v$data$y)) {
    detected <- rowSums(v$data$y == 1, na.rm = TRUE) > 0
    points(coords[!detected, , drop = FALSE],
           pch = 1, col = "grey40", cex = 0.5)
    points(coords[detected, , drop = FALSE],
           pch = 16, col = "black", cex = 0.5)
    legend("bottomright",
           legend = c("Detected", "Not detected"),
           pch = c(16, 1), col = c("black", "grey40"),
           cex = 0.65, bg = "white", box.lwd = 0.5)
  }

  # Color bar
  .draw_colorbar(zlim, col, label = z_label)

  invisible(list(
    x        = grid$x,
    y        = grid$y,
    z        = z,
    type     = type,
    zlim     = zlim,
    coords   = coords,
    detected = detected,
    col      = col
  ))
}
