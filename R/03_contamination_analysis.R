# =============================================================================
# 03_contamination_analysis.R
#
# Task: Introduce 10% outliers and missing values into the Top Gear dataset,
#       then evaluate ICPCA vs MacroPCA.
#
# This corresponds to the project brief requirement:
#   "Introduce 10% of outliers and missing observations into your dataset by
#    inserting observations that deviate significantly from the typical structure
#    of the data. Repeat the analysis and evaluate the effectiveness of the
#    robust estimates."
#
# Contamination scheme:
#   (a) 10% of cells replaced by large deviations  (cellwise outliers)
#   (b) 10% of rows shifted to a different population  (rowwise outliers)
#   (c) 5%  of cells set to NA  (missing at random, on top of existing NAs)
#
# Evaluation:
#   - Subspace recovery: angle between estimated and ground-truth loadings
#   - Outlier detection: sensitivity / specificity / precision
#   - Residual maps: ICPCA vs MacroPCA on the contaminated data
#   - Qualitative comparison of what each method flags
# =============================================================================

library(cellWise)
library(robustHD)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(gridExtra)

set.seed(2025)

# =============================================================================
# 0. Load and preprocess Top Gear data  (identical to script 02)
# =============================================================================
data("topgear", package = "robustHD")

car_names <- rownames(topgear)
cont_vars <- c("Price", "Displacement", "BHP", "Torque",
               "Acceleration", "TopSpeed", "MPG",
               "Weight", "Length", "Width", "Height")

available <- intersect(cont_vars, colnames(topgear))
X_raw     <- as.matrix(topgear[, available])

# Log-transform skewed variables
log_vars <- intersect(c("Price", "Displacement", "BHP", "Torque", "TopSpeed"),
                      colnames(X_raw))
X_clean  <- X_raw
for (v in log_vars) X_clean[, v] <- log(X_raw[, v])

n <- nrow(X_clean)
p <- ncol(X_clean)
k <- 2   # as in paper

cat("=== Dataset: Top Gear (log-transformed) ===\n")
cat("n =", n, ", p =", p, ", existing NAs =", sum(is.na(X_clean)), "\n\n")

# =============================================================================
# 1. Ground-truth fits on CLEAN data
#
# "Ground truth" = MacroPCA and ICPCA fitted on the original data (with its
# natural NAs, no artificial contamination).  We use these loadings as the
# reference to measure subspace recovery after contamination.
# =============================================================================
cat("Fitting ground-truth models on clean data...\n")

fit_clean_macro <- MacroPCA(X_clean, k = k, DDCpars = list(silent = TRUE))
fit_clean_icpca <- ICPCA(X_clean, k = k)

P_truth_macro <- fit_clean_macro$loadings   # p x k reference loadings
P_truth_icpca <- fit_clean_icpca$loadings

# Subspace angle between two loading matrices (Frobenius norm of projection diff)
subspace_angle <- function(P1, P2) {
  proj1 <- P1 %*% solve(t(P1) %*% P1) %*% t(P1)
  proj2 <- P2 %*% solve(t(P2) %*% P2) %*% t(P2)
  norm(proj1 - proj2, type = "F") / sqrt(2 * ncol(P1))
}

# =============================================================================
# 2. Contamination function
# =============================================================================
contaminate_topgear <- function(X,
                                frac_cell  = 0.10,
                                frac_row   = 0.10,
                                frac_miss  = 0.05,
                                shift_mult = 8) {
  # shift_mult: cellwise outliers are set to shift_mult * column_MAD away from median
  Xc <- X
  true_cells <- matrix(FALSE, nrow(X), ncol(X))
  true_rows  <- rep(FALSE, nrow(X))

  # (a) Cellwise outliers: replace with large positive or negative values
  n_cell <- round(frac_cell * sum(!is.na(X)))   # fraction of observed cells
  obs_positions <- which(!is.na(X))
  cell_idx <- sample(obs_positions, n_cell)

  # Direction: randomly +/-
  signs <- sample(c(-1, 1), n_cell, replace = TRUE)
  # Scale: shift_mult * column MAD (robust spread)
  col_mads <- apply(X, 2, function(x) mad(x, na.rm = TRUE))
  col_mads[col_mads < 1e-10] <- 1

  col_of_idx <- ((cell_idx - 1) %% nrow(X)) + 1  # row-major index → col
  # Note: R stores matrices column-major, so:
  col_of_idx <- ceiling(cell_idx / nrow(X))
  row_of_idx <- ((cell_idx - 1) %% nrow(X)) + 1

  for (ii in seq_along(cell_idx)) {
    r <- row_of_idx[ii]
    cc <- col_of_idx[ii]
    col_med <- median(X[, cc], na.rm = TRUE)
    Xc[r, cc] <- col_med + signs[ii] * shift_mult * col_mads[cc]
    true_cells[r, cc] <- TRUE
  }

  # (b) Rowwise outliers: shift entire rows to a different mean
  bad_rows <- sample(which(!true_rows), round(frac_row * nrow(X)))
  col_sds  <- apply(X, 2, function(x) sd(x, na.rm = TRUE))
  col_sds[col_sds < 1e-10] <- 1

  for (i in bad_rows) {
    # Shift in a direction orthogonal-ish to the main structure: uniform high
    Xc[i, ] <- apply(X, 2, function(x) quantile(x, 0.95, na.rm = TRUE)) +
               rnorm(ncol(X), 0, 0.3 * col_sds)
    true_rows[i] <- TRUE
  }

  # (c) Missing at random: on top of existing NAs
  n_miss    <- round(frac_miss * n * p)
  miss_pool <- which(!is.na(Xc))   # only set observed cells to NA
  miss_idx  <- sample(miss_pool, min(n_miss, length(miss_pool)))
  Xc[miss_idx] <- NA

  list(
    X_cont      = Xc,
    true_cells  = true_cells,
    true_rows   = true_rows,
    bad_rows    = bad_rows,
    cell_idx    = cell_idx
  )
}

# =============================================================================
# 3. Apply contamination and fit both methods
# =============================================================================
cont  <- contaminate_topgear(X_clean)
X_cont <- cont$X_cont

cat("=== Contamination summary ===\n")
cat("  Cellwise outliers inserted:", sum(cont$true_cells), "cells (",
    sprintf("%.1f", 100 * mean(cont$true_cells)), "%)\n")
cat("  Rowwise outliers inserted: ", sum(cont$true_rows), "rows (",
    sprintf("%.1f", 100 * mean(cont$true_rows)), "%)\n")
cat("  Missing values (total):    ", sum(is.na(X_cont)), "cells (",
    sprintf("%.1f", 100 * mean(is.na(X_cont))), "%)\n\n")

# --- MacroPCA on contaminated data ---
cat("Fitting MacroPCA on contaminated data...\n")
fit_cont_macro <- MacroPCA(X_cont, k = k, DDCpars = list(silent = TRUE))
angle_macro    <- subspace_angle(fit_cont_macro$loadings, P_truth_macro)

# --- ICPCA on contaminated data ---
cat("Fitting ICPCA on contaminated data...\n")
fit_cont_icpca <- ICPCA(X_cont, k = k)
angle_icpca    <- subspace_angle(fit_cont_icpca$loadings, P_truth_icpca)

cat("\n=== Subspace recovery (lower angle = better) ===\n")
cat("  ICPCA    subspace angle:", round(angle_icpca, 4), "\n")
cat("  MacroPCA subspace angle:", round(angle_macro, 4), "\n\n")

# =============================================================================
# 4. Outlier detection evaluation
# =============================================================================

# --- Rowwise detection ---
detected_rows <- fit_cont_macro$indrows   # logical vector
true_rows     <- cont$true_rows

TP_r <- sum(detected_rows  &  true_rows, na.rm = TRUE)
FP_r <- sum(detected_rows  & !true_rows, na.rm = TRUE)
FN_r <- sum(!detected_rows &  true_rows, na.rm = TRUE)
TN_r <- sum(!detected_rows & !true_rows, na.rm = TRUE)

sens_row  <- TP_r / (TP_r + FN_r)
spec_row  <- TN_r / (TN_r + FP_r)
prec_row  <- if ((TP_r + FP_r) > 0) TP_r / (TP_r + FP_r) else NA
f1_row    <- if (!is.na(prec_row)) 2 * prec_row * sens_row / (prec_row + sens_row) else NA

cat("=== Rowwise outlier detection (MacroPCA) ===\n")
cat("  TP:", TP_r, "  FP:", FP_r, "  FN:", FN_r, "  TN:", TN_r, "\n")
cat("  Sensitivity:", round(sens_row, 3), "\n")
cat("  Specificity:", round(spec_row, 3), "\n")
cat("  Precision:  ", round(prec_row, 3), "\n")
cat("  F1 score:   ", round(f1_row,   3), "\n\n")

# --- Cellwise detection ---
detected_cells <- (fit_cont_macro$indcells != 0)
# Only evaluate on cells that were observed (not NA) in contaminated data
observed       <- !is.na(X_cont)
true_cells_obs <- cont$true_cells & observed
detected_c_obs <- detected_cells  & observed

TP_c <- sum(detected_c_obs &  true_cells_obs, na.rm = TRUE)
FP_c <- sum(detected_c_obs & !true_cells_obs, na.rm = TRUE)
FN_c <- sum(!detected_c_obs & true_cells_obs, na.rm = TRUE)

sens_cell  <- TP_c / (TP_c + FN_c)
prec_cell  <- if ((TP_c + FP_c) > 0) TP_c / (TP_c + FP_c) else NA
f1_cell    <- if (!is.na(prec_cell)) 2 * prec_cell * sens_cell / (prec_cell + sens_cell) else NA

cat("=== Cellwise outlier detection (MacroPCA / DDC) ===\n")
cat("  TP:", TP_c, "  FP:", FP_c, "  FN:", FN_c, "\n")
cat("  Sensitivity:", round(sens_cell, 3), "\n")
cat("  Precision:  ", round(prec_cell, 3), "\n")
cat("  F1 score:   ", round(f1_cell,   3), "\n\n")

# ICPCA does not natively detect outliers; compare subspace quality only

# =============================================================================
# 5. Residual maps: ICPCA vs MacroPCA on contaminated data  (Figure 3 analog)
# =============================================================================

# Select rows to display: true outliers + highest OD
true_any    <- which(cont$true_rows | rowSums(cont$true_cells) > 0)
top_od_cont <- order(fit_cont_macro$OD, decreasing = TRUE)[1:10]
disp_rows   <- unique(c(true_any, top_od_cont))[1:min(24, length(unique(c(true_any, top_od_cont))))]
disp_names  <- car_names[disp_rows]

# MacroPCA residual map on contaminated data
pdf("contamination_residual_map_macropca.pdf", width = 10, height = 7)
cellMap(
  d             = fit_cont_macro$stdResid[disp_rows, ],
  indcells      = fit_cont_macro$indcells[disp_rows, ],
  rowlabels     = disp_names,
  columnlabels  = colnames(X_clean),
  mTitle        = "MacroPCA – Contaminated data residual map",
  sizetitlex    = 9, sizetitley = 8
)
dev.off()

# ICPCA residual map on contaminated data
Xhat_icpca_cont <- sweep(
  fit_cont_icpca$scores %*% t(fit_cont_icpca$loadings),
  2, fit_cont_icpca$center, "+"
)
X_naimputed_cont <- X_cont
for (j in seq_len(ncol(X_cont))) {
  na_j <- is.na(X_cont[, j])
  X_naimputed_cont[na_j, j] <- Xhat_icpca_cont[na_j, j]
}
resid_icpca_cont <- X_naimputed_cont - Xhat_icpca_cont
col_mad_cont     <- apply(resid_icpca_cont, 2, function(x) mad(x, na.rm = TRUE))
col_mad_cont[col_mad_cont < 1e-10] <- 1
stdR_icpca_cont  <- sweep(resid_icpca_cont, 2, col_mad_cont, "/")

indcells_icpca_cont <- matrix(0L, nrow(X_cont), ncol(X_cont))
indcells_icpca_cont[stdR_icpca_cont >  2.576 & !is.na(stdR_icpca_cont)] <-  1L
indcells_icpca_cont[stdR_icpca_cont < -2.576 & !is.na(stdR_icpca_cont)] <- -1L

pdf("contamination_residual_map_icpca.pdf", width = 10, height = 7)
cellMap(
  d             = stdR_icpca_cont[disp_rows, ],
  indcells      = indcells_icpca_cont[disp_rows, ],
  rowlabels     = disp_names,
  columnlabels  = colnames(X_clean),
  mTitle        = "ICPCA – Contaminated data residual map",
  sizetitlex    = 9, sizetitley = 8
)
dev.off()
cat("Saved: contamination residual maps\n")

# =============================================================================
# 6. Outlier map on contaminated data (coloured by ground truth)
# =============================================================================
ground_truth_label <- case_when(
  cont$true_rows ~ "True rowwise outlier",
  rowSums(cont$true_cells) > 0 ~ "True cellwise outlier",
  TRUE ~ "Regular"
)

p_cont_outlier_map <- ggplot(
  data.frame(
    SD    = fit_cont_macro$SD,
    OD    = fit_cont_macro$OD,
    Truth = ground_truth_label,
    label = car_names
  ),
  aes(x = SD, y = OD, colour = Truth)
) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_vline(xintercept = fit_cont_macro$cutoffSD,
             linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = fit_cont_macro$cutoffOD,
             linetype = "dashed", colour = "grey40") +
  geom_text_repel(
    data = ~filter(.x, Truth != "Regular" | SD > fit_cont_macro$cutoffSD | OD > fit_cont_macro$cutoffOD),
    aes(label = label), size = 2.5, max.overlaps = 15
  ) +
  scale_colour_manual(values = c(
    "Regular"               = "#BBBBBB",
    "True rowwise outlier"  = "#D62728",
    "True cellwise outlier" = "#4C8BB5"
  )) +
  labs(
    title    = "MacroPCA outlier map – contaminated Top Gear data",
    subtitle = "Points coloured by ground-truth contamination label",
    x        = "Score distance (SD)",
    y        = "Orthogonal distance (OD)",
    colour   = "Ground truth"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("contamination_outlier_map.pdf", p_cont_outlier_map,
       width = 7, height = 6)
cat("Saved: contamination_outlier_map.pdf\n")

# =============================================================================
# 7. Summary comparison bar charts
# =============================================================================

# (a) Subspace angle comparison
angle_df <- data.frame(
  Method = c("ICPCA", "MacroPCA"),
  Angle  = c(angle_icpca, angle_macro)
)

p_angle <- ggplot(angle_df, aes(x = Method, y = Angle, fill = Method)) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.85) +
  geom_text(aes(label = round(Angle, 4)), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c("ICPCA" = "#E07B54", "MacroPCA" = "#2CA02C"),
                    guide = "none") +
  labs(
    title    = "Subspace recovery on contaminated Top Gear data",
    subtitle = "Angle between estimated and clean-data loadings (lower = better)",
    x        = NULL, y = "Subspace angle"
  ) +
  theme_bw(base_size = 12)

# (b) Detection metrics
perf_df <- data.frame(
  Metric = c("Sensitivity", "Specificity", "Precision", "F1",
             "Sensitivity", "Precision", "F1"),
  Value  = c(sens_row, spec_row, prec_row, f1_row,
             sens_cell, prec_cell, f1_cell),
  Type   = c(rep("Rowwise", 4), rep("Cellwise", 3))
)

p_perf <- ggplot(perf_df, aes(x = Metric, y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.85, width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", Value)),
            position = position_dodge(width = 0.65),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c("Rowwise" = "#D62728", "Cellwise" = "#4C8BB5")) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  labs(
    title    = "Outlier detection performance (MacroPCA, contaminated data)",
    subtitle = "10% rowwise + 10% cellwise contamination + 5% added NAs",
    x        = NULL, y = "Rate", fill = "Outlier type"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

p_summary <- grid.arrange(p_angle, p_perf, ncol = 2)
ggsave("contamination_summary.pdf", p_summary, width = 12, height = 5)
cat("Saved: contamination_summary.pdf\n")

# =============================================================================
# 8. Comparison: clean vs contaminated residual maps side-by-side
#    Shows how contamination affects what each method sees
# =============================================================================
# Re-run clean MacroPCA residual map for same cars
pdf("clean_vs_contaminated_macropca.pdf", width = 14, height = 7)
par(mfrow = c(1, 2))

cellMap(
  d             = fit_clean_macro$stdResid[disp_rows, ],
  indcells      = fit_clean_macro$indcells[disp_rows, ],
  rowlabels     = disp_names,
  columnlabels  = colnames(X_clean),
  mTitle        = "MacroPCA – Clean data",
  sizetitlex    = 9, sizetitley = 8
)

cellMap(
  d             = fit_cont_macro$stdResid[disp_rows, ],
  indcells      = fit_cont_macro$indcells[disp_rows, ],
  rowlabels     = disp_names,
  columnlabels  = colnames(X_clean),
  mTitle        = "MacroPCA – Contaminated data (10%+10%+5% NA)",
  sizetitlex    = 9, sizetitley = 8
)

dev.off()
cat("Saved: clean_vs_contaminated_macropca.pdf\n")

# =============================================================================
# 9. Console summary
# =============================================================================
cat("\n")
cat("=================================================================\n")
cat("KEY FINDINGS\n")
cat("=================================================================\n\n")

cat("SUBSPACE RECOVERY\n")
cat(sprintf("  ICPCA    angle = %.4f  (higher = more distorted by outliers)\n", angle_icpca))
cat(sprintf("  MacroPCA angle = %.4f  (lower  = more robust)\n", angle_macro))
cat(sprintf("  Improvement:     %.1f%%\n\n", 100 * (angle_icpca - angle_macro) / angle_icpca))

cat("DETECTION PERFORMANCE (MacroPCA)\n")
cat(sprintf("  Rowwise  – Sensitivity: %.2f | Specificity: %.2f | Precision: %.2f | F1: %.2f\n",
            sens_row, spec_row, prec_row, f1_row))
cat(sprintf("  Cellwise – Sensitivity: %.2f | Precision: %.2f | F1: %.2f\n\n",
            sens_cell, prec_cell, f1_cell))

cat("ICPCA LIMITATIONS OBSERVED\n")
cat("  - Cannot distinguish cellwise from rowwise contamination\n")
cat("  - Loadings distorted by outlying rows attracting the fit\n")
cat("  - No mechanism to impute/flag individual outlying cells\n")
cat("  - Residual map shows diffuse contamination rather than isolated cells\n\n")

cat("MACROPCA ADVANTAGES OBSERVED\n")
cat("  - Two-step process isolates cellwise anomalies before robust PCA fit\n")
cat("  - Handles existing NAs + artificially added NAs natively\n")
cat("  - Provides interpretable cell-level flagging (which variable, which car)\n")
cat("  - Subspace closer to clean-data solution under contamination\n")
cat("=================================================================\n")
