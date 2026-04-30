# =============================================================================
# 03_contamination_analysis.R
#
# Task 4: Artificially introduce 10% outliers AND missing values into the
# real dataset (DPOSS), then compare Classical PCA vs MacroPCA.
#
# Contamination scheme (as required by project brief):
#   (a) 10% of cells replaced by values that deviate significantly from
#       the data's typical structure (cellwise outliers)
#   (b) 10% of observations shifted as whole rows (casewise outliers)
#   (c) An additional 5% of cells set to NA (missing at random)
#
# We evaluate: subspace angle vs the "ground truth" (fit on clean data),
# outlier detection rate (sensitivity/specificity), and qualitative comparison.
# =============================================================================

library(cellWise)
library(ggplot2)
library(dplyr)
library(gridExtra)

set.seed(2025)

# =============================================================================
# 0. Load and preprocess (same as script 02)
# =============================================================================
data("dposs", package = "cellWise")
X_clean <- scale(dposs)    # standardised clean data
n <- nrow(X_clean)
p <- ncol(X_clean)
k <- 3                     # number of PC components

cat("=== Dataset: DPOSS ===\n")
cat("n =", n, ", p =", p, "\n\n")

# Fit MacroPCA and Classical PCA on clean data → use as "ground truth" loadings
pca_clean   <- prcomp(X_clean, center = FALSE, scale. = FALSE)
macro_clean <- MacroPCA(X_clean, k = k, DDCpars = list(silent = TRUE))

V_truth_class <- pca_clean$rotation[, 1:k]
V_truth_macro <- macro_clean$loadings[, 1:k]

# Helper: subspace angle
subspace_angle <- function(V1, V2) {
  P1 <- V1 %*% solve(t(V1) %*% V1) %*% t(V1)
  P2 <- V2 %*% solve(t(V2) %*% V2) %*% t(V2)
  norm(P1 - P2, type = "F") / sqrt(2 * k)
}

# =============================================================================
# 1. Define contamination function
# =============================================================================
contaminate <- function(X, frac_cell = 0.10, frac_row = 0.10,
                        frac_miss = 0.05, shift_cell = 8, shift_row = 6) {
  Xc <- X
  true_cell_outliers <- matrix(FALSE, nrow(X), ncol(X))
  true_row_outliers  <- rep(FALSE, nrow(X))

  # (a) Cellwise outliers: random cells shifted by a large amount
  n_cell_cont <- round(frac_cell * nrow(X) * ncol(X))
  cell_idx    <- sample(length(X), n_cell_cont)
  Xc[cell_idx] <- Xc[cell_idx] + rnorm(n_cell_cont, mean = shift_cell, sd = 0.5)
  true_cell_outliers[cell_idx] <- TRUE

  # (b) Casewise outliers: whole rows shifted
  bad_rows <- sample(nrow(X), round(frac_row * nrow(X)))
  Xc[bad_rows, ] <- Xc[bad_rows, ] +
    matrix(rnorm(length(bad_rows) * ncol(X), mean = shift_row, sd = 0.5),
           nrow = length(bad_rows))
  true_row_outliers[bad_rows] <- TRUE

  # (c) Missing values (MAR)
  n_miss <- round(frac_miss * nrow(X) * ncol(X))
  miss_idx  <- sample(length(Xc), n_miss)
  Xc[miss_idx] <- NA

  list(
    X_cont             = Xc,
    true_cell_outliers = true_cell_outliers,
    true_row_outliers  = true_row_outliers,
    bad_rows           = bad_rows,
    cell_idx           = cell_idx,
    miss_idx           = miss_idx
  )
}

# =============================================================================
# 2. Apply contamination and run both methods
# =============================================================================
cont <- contaminate(X_clean)
X_cont <- cont$X_cont

cat("Contamination summary:\n")
cat("  Cellwise outliers inserted:  ", sum(cont$true_cell_outliers), "cells (",
    round(100 * mean(cont$true_cell_outliers), 1), "%)\n")
cat("  Casewise outliers inserted:  ", length(cont$bad_rows), "rows (",
    round(100 * length(cont$bad_rows) / n, 1), "%)\n")
cat("  Missing values inserted:     ", sum(is.na(X_cont)), "cells (",
    round(100 * mean(is.na(X_cont)), 1), "%)\n\n")

# --- Classical PCA on contaminated data (requires complete cases only) ---
complete_rows   <- complete.cases(X_cont)
X_cont_complete <- X_cont[complete_rows, ]
pca_cont <- prcomp(X_cont_complete, center = TRUE, scale. = FALSE)

# Compare loadings (on complete cases only, hence approximate)
V_class_cont <- pca_cont$rotation[, 1:k]
angle_class  <- subspace_angle(V_class_cont, V_truth_class)

cat("Classical PCA:\n")
cat("  Complete cases used:    ", sum(complete_rows), "/", n, "\n")
cat("  Subspace angle vs truth:", round(angle_class, 4), "\n\n")

# --- MacroPCA on contaminated data (handles missing values natively) ---
cat("Running MacroPCA on contaminated data...\n")
macro_cont <- MacroPCA(X_cont, k = k, DDCpars = list(silent = TRUE))

V_macro_cont <- macro_cont$loadings[, 1:k]
angle_macro  <- subspace_angle(V_macro_cont, V_truth_macro)

cat("MacroPCA:\n")
cat("  All", n, "observations used (missing values imputed internally)\n")
cat("  Subspace angle vs truth:", round(angle_macro, 4), "\n\n")

# =============================================================================
# 3. Evaluate outlier detection performance
# =============================================================================

# Casewise outlier detection
detected_rows      <- rep(FALSE, n)
detected_rows[macro_cont$casewiseOutliers] <- TRUE
true_rows_all      <- cont$true_row_outliers

TP_row <- sum(detected_rows  &  true_rows_all)
FP_row <- sum(detected_rows  & !true_rows_all)
FN_row <- sum(!detected_rows &  true_rows_all)
TN_row <- sum(!detected_rows & !true_rows_all)

sensitivity_row <- TP_row / (TP_row + FN_row)
specificity_row <- TN_row / (TN_row + FP_row)
precision_row   <- TP_row / (TP_row + FP_row)

cat("=== Casewise Outlier Detection (MacroPCA) ===\n")
cat("  Sensitivity (recall):", round(sensitivity_row, 3), "\n")
cat("  Specificity:         ", round(specificity_row, 3), "\n")
cat("  Precision:           ", round(precision_row,   3), "\n\n")

# Cellwise outlier detection
detected_cells <- (macro_cont$indcells != 0)
# Exclude cells that are NA in the contaminated data (not cells we injected)
known_miss <- is.na(X_cont)
true_cells <- cont$true_cell_outliers & !known_miss
detected_c <- detected_cells & !known_miss

TP_cell <- sum(detected_c & true_cells)
FP_cell <- sum(detected_c & !true_cells)
FN_cell <- sum(!detected_c & true_cells)

sensitivity_cell <- TP_cell / (TP_cell + FN_cell)
precision_cell   <- TP_cell / (TP_cell + FP_cell)

cat("=== Cellwise Outlier Detection (MacroPCA / DDC step) ===\n")
cat("  Sensitivity (recall):", round(sensitivity_cell, 3), "\n")
cat("  Precision:           ", round(precision_cell,   3), "\n\n")

# =============================================================================
# 4. Plots
# =============================================================================

# (a) Outlier map: clean vs contaminated
SD_cont <- macro_cont$SD
OD_cont <- macro_cont$OD

cut_SD <- macro_cont$cutoffSD
cut_OD <- macro_cont$cutoffOD

flag <- case_when(
  true_rows_all                         ~ "True casewise outlier",
  cont$true_cell_outliers[, 1] |        # simplified: flag row if any cell contaminated
    rowSums(cont$true_cell_outliers) > 0 ~ "True cellwise outlier",
  TRUE                                  ~ "Regular"
)

outlier_cont_df <- data.frame(SD = SD_cont, OD = OD_cont, flag = flag)

p_cont_map <- ggplot(outlier_cont_df, aes(x = SD, y = OD, colour = flag)) +
  geom_point(alpha = 0.5, size = 0.9) +
  geom_vline(xintercept = cut_SD, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = cut_OD, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(
    "Regular"               = "#AAAAAA",
    "True casewise outlier" = "#E07B54",
    "True cellwise outlier" = "#4C8BB5"
  )) +
  labs(title    = "MacroPCA outlier map – contaminated data",
       subtitle = "Points coloured by ground truth label",
       x = "Score distance (SD)", y = "Orthogonal distance (OD)",
       colour = "Ground truth") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p_cont_map)
ggsave("contamination_outlier_map.pdf", p_cont_map, width = 6, height = 5)

# (b) Summary bar chart: subspace angle comparison
angle_df <- data.frame(
  method    = c("Classical PCA\n(complete cases)", "MacroPCA\n(all data)"),
  angle     = c(angle_class, angle_macro)
)

p_angle <- ggplot(angle_df, aes(x = method, y = angle, fill = method)) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.85) +
  scale_fill_manual(values = c("#E07B54", "#4C8BB5"), guide = "none") +
  labs(title    = "Subspace recovery on contaminated data",
       subtitle = "Lower angle = closer to clean-data solution",
       x = NULL, y = "Subspace angle") +
  theme_bw(base_size = 12)

print(p_angle)
ggsave("contamination_subspace_angle.pdf", p_angle, width = 5, height = 4)

# (c) Detection performance bar chart
perf_df <- data.frame(
  Metric    = rep(c("Sensitivity", "Specificity/\nPrecision"), 2),
  Value     = c(sensitivity_row, specificity_row,
                sensitivity_cell, precision_cell),
  Type      = rep(c("Casewise", "Cellwise"), each = 2)
)

p_perf <- ggplot(perf_df, aes(x = Metric, y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = c("Casewise" = "#E07B54", "Cellwise" = "#4C8BB5")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Outlier detection performance (MacroPCA)",
       x = NULL, y = "Rate", fill = "Outlier type") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p_perf)
ggsave("contamination_detection_performance.pdf", p_perf, width = 5, height = 4)
message("Saved: contamination_detection_performance.pdf")

cat("\n=== Contamination analysis complete ===\n")

# =============================================================================
# 5. Discussion summary (printed to console for reference)
# =============================================================================
cat("
=== KEY DISCUSSION POINTS ===

1. SUBSPACE RECOVERY
   Classical PCA (complete cases only): angle =", round(angle_class, 4), "
   MacroPCA (all data, handles NAs):    angle =", round(angle_macro, 4), "
   MacroPCA is closer to the ground truth because it:
     - Does not discard missing observations
     - Downweights outlying cells and rows

2. OUTLIER DETECTION
   Casewise: sensitivity =", round(sensitivity_row, 3), ", specificity =", round(specificity_row, 3), "
   Cellwise: sensitivity =", round(sensitivity_cell, 3), ", precision   =", round(precision_cell, 3), "

3. CLASSICAL PCA WEAKNESSES
   - Must discard rows with any NA, losing", n - sum(complete_rows), "observations
   - Its principal components are tilted toward outliers
   - Cannot distinguish cellwise from casewise contamination

4. MACROPCA STRENGTHS
   - Natively handles missing values (no listwise deletion)
   - Two-step approach: first cleans cells (DDC), then applies robust PCA
   - Provides score distances + orthogonal distances for outlier classification

5. MACROPCA LIMITATIONS
   - Computationally heavier than classical PCA
   - DDC step relies on a propagation algorithm that may miss outliers
     in very high-dimensional data (p >> n)
   - Tuning parameters (cut-off quantiles) affect sensitivity/specificity trade-off
")
