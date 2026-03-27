# =============================================================================
# Create the scaling comparison figure for INLAocc README
# Top row: computation time.  Bottom row: accuracy (correlation with truth).
# =============================================================================

library(ggplot2)
library(scales)
library(svglite)
library(patchwork)

root <- "C:/Users/Gilles Colling/Documents/dev/INLAocc"

# --- Load and merge results ---
df1 <- read.csv(file.path(root, "benchmark", "scaling_results.csv"),
                stringsAsFactors = FALSE)
ms_file <- file.path(root, "benchmark", "scaling_results_ms.csv")
if (file.exists(ms_file)) {
  df2 <- read.csv(ms_file, stringsAsFactors = FALSE)
  df  <- rbind(df1, df2)
} else {
  df <- df1
}

# --- Map parallel multi-species into the same panel with a distinct method ---
par_rows <- df$type == "Multi-species (parallel)"
df$method[par_rows] <- "INLAocc (parallel)"
df$type[par_rows]   <- "Multi-species"

# --- Colours (Wong palette, colourblind-safe) ---
method_cols  <- c(INLAocc = "#0072B2", `INLAocc (parallel)` = "#56B4E9",
                  spOccupancy = "#D55E00", Stan = "#CC79A7")
method_shape <- c(INLAocc = 16, `INLAocc (parallel)` = 1,
                  spOccupancy = 15, Stan = 17)
method_lty   <- c(INLAocc = "solid", `INLAocc (parallel)` = "dashed",
                  spOccupancy = "solid", Stan = "solid")

method_levels <- c("INLAocc", "INLAocc (parallel)", "spOccupancy", "Stan")
df$method <- factor(df$method, levels = method_levels)

# --- Facet labels ---
panel_map <- c(
  "Non-spatial"   = "Non-spatial",
  "Spatial"       = "Spatial (SPDE / NNGP)",
  "Multi-species" = "Multi-species (10 spp.)"
)
df$panel <- panel_map[df$type]
panel_order <- intersect(panel_map, unique(df$panel))
df$panel <- factor(df$panel, levels = panel_order)

# --- Speedup annotations (bracket at rightmost shared N per panel) ---
anno_list <- lapply(split(df, df$panel), function(sub) {
  inla_method <- if (any(sub$method == "INLAocc (parallel)"))
    "INLAocc (parallel)" else "INLAocc"
  shared <- intersect(
    sub$N[sub$method == inla_method],
    sub$N[sub$method == "spOccupancy"]
  )
  if (length(shared) == 0) return(NULL)
  max_N   <- max(shared)
  t_inla  <- sub$time[sub$method == inla_method    & sub$N == max_N]
  t_spocc <- sub$time[sub$method == "spOccupancy"  & sub$N == max_N]
  if (length(t_inla) == 0 || length(t_spocc) == 0) return(NULL)
  speedup <- t_spocc / t_inla
  data.frame(
    panel = sub$panel[1], N = max_N,
    y = sqrt(t_inla * t_spocc), ymin = t_inla, ymax = t_spocc,
    label = sprintf("%.0f\u00d7", speedup), stringsAsFactors = FALSE
  )
})
anno <- do.call(rbind, anno_list)
anno$panel <- factor(anno$panel, levels = levels(df$panel))

# --- Top row: computation time ---
p_time <- ggplot(df, aes(x = N, y = time, colour = method, shape = method,
                         linetype = method)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 3) +
  geom_segment(data = anno, inherit.aes = FALSE,
               aes(x = N * 1.15, xend = N * 1.15, y = ymin, yend = ymax),
               colour = "grey50", linewidth = 0.4,
               arrow = arrow(ends = "both", length = unit(0.06, "inches"),
                             type = "open")) +
  geom_label(data = anno, inherit.aes = FALSE,
             aes(x = N * 1.15, y = y, label = label),
             size = 3.1, fontface = "bold",
             colour = "grey25", fill = "white",
             linewidth = 0, label.padding = unit(0.15, "lines")) +
  facet_wrap(~ panel, scales = "free_x", nrow = 1) +
  scale_x_log10(labels = label_comma(),
                expand = expansion(mult = c(0.05, 0.12))) +
  scale_y_log10(labels = label_comma(suffix = "s"),
                breaks = c(1, 3, 10, 30, 100, 300, 1000, 3000)) +
  scale_colour_manual(values = method_cols) +
  scale_shape_manual(values = method_shape) +
  scale_linetype_manual(values = method_lty) +
  labs(x = NULL, y = "Computation time",
       colour = NULL, shape = NULL, linetype = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major  = element_line(colour = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold", size = 11),
    legend.position   = "top",
    legend.text       = element_text(size = 11),
    legend.key.width  = unit(1.8, "lines"),
    axis.title        = element_text(size = 11),
    axis.text.x       = element_blank(),
    plot.margin       = margin(10, 14, 0, 6)
  )

# --- Bottom row: accuracy (correlation with truth) ---
# Drop parallel from accuracy (same model, same accuracy as sequential)
df_acc <- df[df$method != "INLAocc (parallel)" & !is.na(df$cor_psi), ]

p_acc <- ggplot(df_acc, aes(x = N, y = cor_psi, colour = method, shape = method)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 3) +
  facet_wrap(~ panel, scales = "free_x", nrow = 1) +
  scale_x_log10(labels = label_comma(),
                expand = expansion(mult = c(0.05, 0.12))) +
  scale_y_continuous(limits = c(0.5, 1.02), breaks = seq(0.5, 1.0, 0.1),
                     labels = function(x) sprintf("%.1f", x)) +
  scale_colour_manual(values = method_cols) +
  scale_shape_manual(values = method_shape) +
  labs(x = "Number of sites", y = "Cor(\u03C8, truth)",
       colour = NULL, shape = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major  = element_line(colour = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    strip.text        = element_blank(),
    legend.position   = "none",
    axis.title        = element_text(size = 11),
    plot.margin       = margin(2, 14, 6, 6)
  )

# --- Combine ---
p <- p_time / p_acc + plot_layout(heights = c(3, 1.5))

# --- Save ---
dir.create(file.path(root, "man", "figures"), showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(root, "man", "figures", "benchmark.svg"),
       p, width = 12, height = 6.5, device = svglite)
ggsave(file.path(root, "man", "figures", "benchmark.png"),
       p, width = 12, height = 6.5, dpi = 200)

cat("Saved:\n")
cat("  man/figures/benchmark.svg\n")
cat("  man/figures/benchmark.png\n")
