# =============================================================================
# 02_real_data_analysis.R
#
# Application of Classical PCA vs MacroPCA to a real dataset.
#
# Dataset used: DPOSS (Digitized Palomar Sky Survey) -- included in cellWise
#   - 6215 observations, 10 photometric variables
#   - Known to contain casewise and cellwise outliers
#   - Used as the main illustration in the MacroPCA paper (Section 5)
#
# Tasks:
#   1. Explore and preprocess the data
#   2. Apply Classical PCA
#   3. Apply MacroPCA
#   4. Compare: score plots, residual plots, outlier maps, diagnostic summaries
# =============================================================================

library(cellWise)
library(ggplot2)
library(dplyr)
library(gridExtra)

# =============================================================================
# 1. Load and inspect the DPOSS dataset
# =============================================================================
data("dposs", package = "cellWise")
X_raw <- dposs

cat("=== DPOSS Dataset ===\n")
cat("Dimensions:", nrow(X_raw), "x", ncol(X_raw), "\n")
cat("Variables: ", paste(colnames(X_raw), collapse = ", "), "\n\n")
cat("Missing values per variable:\n")
print(colSums(is.na(X_raw)))
cat("\nBasic summary:\n")
print(summary(X_raw))

# =============================================================================
# 2. Pre-process: standardise columns (MacroPCA expects roughly standardised data)
# =============================================================================
X_scaled <- scale(X_raw)   # mean 0, sd 1 for each variable

# =============================================================================
# 3. Classical PCA
# =============================================================================
k <- 3   # number of components (as in the paper)

pca_classic <- prcomp(X_scaled, center = FALSE, scale. = FALSE)

# Proportion of variance explained
pve <- pca_classic$sdev^2 / sum(pca_classic$sdev^2)
cat("\n=== Classical PCA: Proportion of Variance Explained ===\n")
print(round(cumsum(pve)[1:k], 4))

# Score plot (PC1 vs PC2)
scores_classic <- as.data.frame(pca_classic$x[, 1:2])
colnames(scores_classic) <- c("PC1", "PC2")

p_classic_scores <- ggplot(scores_classic, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.3, size = 0.8, colour = "#E07B54") +
  labs(title = "Classical PCA – Score plot (PC1 vs PC2)",
       x = "PC 1", y = "PC 2") +
  theme_bw(base_size = 12)

print(p_classic_scores)
ggsave("pca_classic_scores.pdf", p_classic_scores, width = 6, height = 5)

# Loadings heatmap
loadings_df <- as.data.frame(pca_classic$rotation[, 1:k])
loadings_df$variable <- rownames(loadings_df)
loadings_long <- tidyr::pivot_longer(loadings_df, -variable,
                                     names_to = "PC", values_to = "loading")

p_loadings <- ggplot(loadings_long, aes(x = PC, y = variable, fill = loading)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "#D7191C", mid = "white", high = "#2C7BB6",
                       midpoint = 0, name = "Loading") +
  labs(title = "Classical PCA – Loadings heatmap",
       x = "Component", y = "Variable") +
  theme_bw(base_size = 12)

ggsave("pca_classic_loadings.pdf", p_loadings, width = 5, height = 5)

# =============================================================================
# 4. MacroPCA
# =============================================================================
cat("\n=== Running MacroPCA (k =", k, ") ===\n")
macro_fit <- MacroPCA(X_scaled, k = k, DDCpars = list(silent = FALSE))

cat("\nMacroPCA summary:\n")
# Number of outlying cells and rows detected
cat("Outlying cells flagged (DDC step):",
    sum(macro_fit$indcells != 0, na.rm = TRUE), "\n")
cat("Casewise outliers flagged:         ",
    length(macro_fit$casewiseOutliers), "\n")

# Proportion of variance explained by MacroPCA components
pve_macro <- macro_fit$eigenvalues / sum(macro_fit$eigenvalues)
cat("Cumulative variance explained:     ",
    round(cumsum(pve_macro)[1:k], 4), "\n")

# --- Score plot ---
scores_macro <- as.data.frame(macro_fit$scores[, 1:2])
colnames(scores_macro) <- c("PC1", "PC2")
scores_macro$outlier <- FALSE
scores_macro$outlier[macro_fit$casewiseOutliers] <- TRUE

p_macro_scores <- ggplot(scores_macro, aes(x = PC1, y = PC2, colour = outlier)) +
  geom_point(alpha = 0.4, size = 0.9) +
  scale_colour_manual(values = c("FALSE" = "#4C8BB5", "TRUE" = "#E07B54"),
                      labels = c("Regular", "Casewise outlier")) +
  labs(title = "MacroPCA – Score plot (PC1 vs PC2)",
       x = "PC 1", y = "PC 2", colour = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p_macro_scores)
ggsave("macropca_scores.pdf", p_macro_scores, width = 6, height = 5)

# --- Outlier map: score distance vs orthogonal distance ---
# Score distance (SD): distance in the PC subspace
# Orthogonal distance (OD): distance from the PC subspace

SD <- macro_fit$SD   # score distances
OD <- macro_fit$OD   # orthogonal distances

cutoff_SD <- macro_fit$cutoffSD
cutoff_OD <- macro_fit$cutoffOD

outlier_type <- case_when(
  SD > cutoff_SD & OD > cutoff_OD ~ "Both",
  SD > cutoff_SD                  ~ "Score outlier",
  OD > cutoff_OD                  ~ "Orthogonal outlier",
  TRUE                            ~ "Regular"
)

outlier_df <- data.frame(SD, OD, type = outlier_type)

p_outlier_map <- ggplot(outlier_df, aes(x = SD, y = OD, colour = type)) +
  geom_point(alpha = 0.5, size = 0.9) +
  geom_vline(xintercept = cutoff_SD, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = cutoff_OD, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(
    "Regular"            = "#AAAAAA",
    "Score outlier"      = "#F0A500",
    "Orthogonal outlier" = "#4C8BB5",
    "Both"               = "#E07B54"
  )) +
  labs(title    = "MacroPCA – Outlier map",
       subtitle = "Dashed lines: 97.5% cut-offs",
       x        = "Score distance (SD)",
       y        = "Orthogonal distance (OD)",
       colour   = "Observation type") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p_outlier_map)
ggsave("macropca_outlier_map.pdf", p_outlier_map, width = 6, height = 5)
message("Saved: macropca_outlier_map.pdf")

# --- Cellwise outlier heatmap (subset of rows for readability) ---
# Show the top 50 most outlying rows
top50 <- order(OD, decreasing = TRUE)[1:50]
cell_flags <- macro_fit$indcells[top50, ]

# Convert to data frame for plotting
cell_df <- as.data.frame(cell_flags)
cell_df$obs <- seq_len(nrow(cell_df))
cell_long  <- tidyr::pivot_longer(cell_df, -obs,
                                  names_to = "variable", values_to = "flag")

p_cell_heat <- ggplot(cell_long, aes(x = variable, y = factor(obs), fill = factor(flag))) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = c("0" = "#F5F5F5", "1" = "#E07B54", "-1" = "#4C8BB5"),
    labels = c("0" = "Normal", "1" = "High outlier", "-1" = "Low outlier"),
    name   = "Cell status"
  ) +
  labs(title    = "Cellwise outlier map (top 50 rows by OD)",
       subtitle = "MacroPCA / DDC step",
       x = "Variable", y = "Observation") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom")

print(p_cell_heat)
ggsave("macropca_cellwise_heatmap.pdf", p_cell_heat, width = 7, height = 6)
message("Saved: macropca_cellwise_heatmap.pdf")

# =============================================================================
# 5. Side-by-side comparison of loadings
# =============================================================================
loadings_macro <- as.data.frame(macro_fit$loadings[, 1:k])
colnames(loadings_macro) <- paste0("PC", 1:k)
loadings_macro$variable <- colnames(X_raw)
loadings_macro$method   <- "MacroPCA"

loadings_classic2 <- as.data.frame(pca_classic$rotation[, 1:k])
colnames(loadings_classic2) <- paste0("PC", 1:k)
loadings_classic2$variable <- colnames(X_raw)
loadings_classic2$method   <- "Classical PCA"

loadings_combined <- rbind(loadings_macro, loadings_classic2)
loadings_long2    <- tidyr::pivot_longer(loadings_combined,
                                         cols      = starts_with("PC"),
                                         names_to  = "PC",
                                         values_to = "loading")

p_load_compare <- ggplot(
  loadings_long2[loadings_long2$PC %in% c("PC1", "PC2"), ],
  aes(x = variable, y = loading, fill = method)
) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("Classical PCA" = "#E07B54", "MacroPCA" = "#4C8BB5")) +
  facet_wrap(~PC, ncol = 1) +
  labs(title = "Loadings comparison: Classical PCA vs MacroPCA",
       x = "Variable", y = "Loading", fill = "Method") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(p_load_compare)
ggsave("loadings_comparison.pdf", p_load_compare, width = 7, height = 6)
message("Saved: loadings_comparison.pdf")

cat("\n=== Real data analysis complete ===\n")
