# =============================================================================
# Parameter recovery figure: estimated vs true coefficients
# Reads benchmark/param_recovery.csv, writes man/figures/parameters.png
# =============================================================================

library(ggplot2)
library(svglite)

root <- "C:/Users/Gilles Colling/Documents/dev/INLAocc"
df <- read.csv(file.path(root, "benchmark", "param_recovery.csv"),
               stringsAsFactors = FALSE)

method_cols  <- c(INLAocc = "#0072B2", spOccupancy = "#D55E00", Stan = "#CC79A7")
method_shape <- c(INLAocc = 16, spOccupancy = 15, Stan = 17)
df$method <- factor(df$method, levels = c("INLAocc", "spOccupancy", "Stan"))

# Facet by parameter type
df$param_type <- ifelse(grepl("^occ", df$param), "Occupancy", "Detection")
df$param_type <- factor(df$param_type, levels = c("Occupancy", "Detection"))

# Axis range
rng <- range(c(df$true_val, df$estimate), na.rm = TRUE)
pad <- diff(rng) * 0.08
lims <- c(rng[1] - pad, rng[2] + pad)

# Per-method correlation annotation
anno <- do.call(rbind, lapply(split(df, list(df$method, df$param_type)), function(sub) {
  if (nrow(sub) < 3) return(NULL)
  data.frame(
    method     = sub$method[1],
    param_type = sub$param_type[1],
    r          = cor(sub$true_val, sub$estimate),
    stringsAsFactors = FALSE
  )
}))
anno$label <- sprintf("r = %.2f", anno$r)
anno$method <- factor(anno$method, levels = levels(df$method))
anno$param_type <- factor(anno$param_type, levels = levels(df$param_type))
# Position in top-left corner
anno$x <- lims[1] + diff(lims) * 0.05
anno$y <- lims[2] - diff(lims) * 0.05

p <- ggplot(df, aes(x = true_val, y = estimate, colour = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.5,
              linetype = "dashed") +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_text(data = anno, inherit.aes = FALSE,
            aes(x = x, y = y, label = label, colour = method),
            size = 3.3, fontface = "bold", hjust = 0, show.legend = FALSE) +
  facet_grid(method ~ param_type) +
  scale_colour_manual(values = method_cols) +
  scale_shape_manual(values = method_shape) +
  coord_fixed(xlim = lims, ylim = lims) +
  labs(x = "True coefficient", y = "Estimated coefficient",
       colour = NULL, shape = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major  = element_line(colour = "grey92", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(face = "bold", size = 11),
    legend.position   = "none",
    axis.title        = element_text(size = 11),
    plot.margin       = margin(10, 14, 6, 6)
  )

dir.create(file.path(root, "man", "figures"), showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(root, "man", "figures", "parameters.svg"),
       p, width = 6.5, height = 7.5, device = svglite)
ggsave(file.path(root, "man", "figures", "parameters.png"),
       p, width = 6.5, height = 7.5, dpi = 200)

cat("Saved:\n")
cat("  man/figures/parameters.svg / .png\n")
