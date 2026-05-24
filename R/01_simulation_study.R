# =============================================================================
# 01_simulation_study.R
#
# Partial replication of the simulation study from:
#   Hubert et al. (2019) MacroPCA, Technometrics 61, 459-473.
#
# We replicate:
#   - Figure 7: 20% NAs + 20% cellwise outliers  (MSE vs gamma)
#   - Figure 9: 20% NAs + 10% cellwise + 10% rowwise outliers (MSE vs gamma)
#
# Methods compared: ICPCA, MROBPCA, MacroPCA  (as in the paper)
# Metric: MSE against baseline PCA on clean data  (paper Section 5)
# Data generating process: A09 covariance, n=100, d=200, k=6  (paper Section 5)
#
# NOTE: With d=200 this takes ~20-40 min. Reduce n_sim or d to prototype faster.
# =============================================================================

library(cellWise)   # MacroPCA, ICPCA, MROBPCA, DDC
library(MASS)       # mvrnorm (backup)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(2024)

# =============================================================================
# 1. Data generating process  (Section 5 of the paper)
# =============================================================================

n      <- 100   # observations
d      <- 200   # variables
k      <- 6     # true number of components
n_sim  <- 30    # Monte Carlo replications (paper uses 100; reduce for speed)

cat("Building A09 covariance matrix (d =", d, ")...\n")

# A09 structured correlation: rho_{ij} = (-0.9)^|i-j|
idx    <- matrix(1:d, d, d)
R_A09  <- (-0.9)^abs(idx - t(idx))

# Target eigenvalues: 6 large + 194 small (paper Section 5)
lambda_large  <- c(30, 25, 20, 15, 10, 5)
lambda_small  <- seq(0.098, 0.0015, length.out = d - k)
lambda_target <- c(lambda_large, lambda_small)

# Build Sigma by replacing eigenvalues of R_A09
eig_R     <- eigen(R_A09, symmetric = TRUE)
Sigma     <- eig_R$vectors %*% diag(lambda_target) %*% t(eig_R$vectors)
Sigma     <- (Sigma + t(Sigma)) / 2          # ensure symmetry

# Variables for contamination
sigma_j   <- sqrt(diag(Sigma))               # column SDs, used for cellwise shift
v_kp1     <- eig_R$vectors[, k + 1]          # (k+1)-th eigenvector, for rowwise shift

# Fast clean-data generator using Cholesky
chol_Sig  <- chol(Sigma)                     # upper-triangular
generate_clean <- function(n) {
  matrix(rnorm(n * d), n, d) %*% chol_Sig   # n x d
}

cat("Done. Sigma built.\n")

# =============================================================================
# 2. MSE helper
#
# Baseline: classical PCA on the CLEAN rows of the uncontaminated data X0.
# For each method applied to contaminated data, compute predictions for those
# same clean rows and measure MSE against the baseline predictions.
# Paper eq: MSE = (1/cd) * sum_{i in C} sum_j (xhat_ij - xhat^C_ij)^2
# =============================================================================
compute_predictions <- function(center, loadings, X_data) {
  # X_data: n x d (may contain NAs; missing entries get predicted from subspace)
  # Returns n x d matrix of predicted values
  X_c  <- sweep(X_data, 2, center, "-")      # centre
  # For rows with NAs, use only observed entries to compute scores
  scores <- matrix(NA, nrow(X_data), ncol(loadings))
  for (i in seq_len(nrow(X_data))) {
    obs <- !is.na(X_c[i, ])
    if (sum(obs) >= ncol(loadings)) {
      # Least-squares projection onto observed dimensions
      P_obs    <- loadings[obs, , drop = FALSE]
      scores[i, ] <- solve(t(P_obs) %*% P_obs) %*% t(P_obs) %*% X_c[i, obs]
    } else {
      scores[i, ] <- 0
    }
  }
  Xhat <- sweep(scores %*% t(loadings), 2, center, "+")
  Xhat
}

compute_mse <- function(Xhat_method, Xhat_base, C_rows) {
  # MSE over clean rows C and all d columns
  diff  <- Xhat_method[C_rows, ] - Xhat_base[C_rows, ]
  mean(diff^2, na.rm = TRUE)
}

# =============================================================================
# 3. Contamination functions
# =============================================================================
add_nas <- function(X, frac = 0.20) {
  Xc     <- X
  n_miss <- round(frac * length(X))
  idx    <- sample(length(X), n_miss)
  Xc[idx] <- NA
  Xc
}

add_cellwise <- function(X, frac = 0.20, gamma) {
  Xc     <- X
  n_cont <- round(frac * length(X))
  idx    <- sample(length(X), n_cont)
  # shift by gamma * sigma_j  (paper: replace x_ij with gamma * sigma_j)
  col_idx <- ((idx - 1) %% d) + 1
  Xc[idx] <- gamma * sigma_j[col_idx]
  Xc
}

add_rowwise <- function(X, frac = 0.20, gamma) {
  Xc      <- X
  bad     <- sample(nrow(X), round(frac * nrow(X)))
  # shift from N(gamma * v_{k+1}, Sigma)  (paper Section 5)
  shift   <- gamma * v_kp1                   # d-vector
  for (i in bad) {
    noise      <- matrix(rnorm(d), 1, d) %*% chol_Sig
    Xc[i, ]    <- shift + noise
  }
  list(X = Xc, bad_rows = bad)
}

# =============================================================================
# 4. Run simulation for a given scenario
# Scenario A: 20% NA + 20% cellwise  (Figure 7)
# Scenario B: 20% NA + 10% cellwise + 10% rowwise  (Figure 9)
# =============================================================================
gamma_vals <- c(0, 1, 2, 3, 5, 7, 10, 15, 20)

run_scenario <- function(scenario_name, frac_cell, frac_row) {
  cat("\n=== Scenario:", scenario_name, "===\n")

  results <- expand.grid(
    gamma  = gamma_vals,
    method = c("ICPCA", "MROBPCA", "MacroPCA"),
    mse    = NA_real_,
    stringsAsFactors = FALSE
  )

  for (gi in seq_along(gamma_vals)) {
    gamma <- gamma_vals[gi]
    cat("  gamma =", gamma, "\n")

    mse_icpca    <- numeric(n_sim)
    mse_mrobpca  <- numeric(n_sim)
    mse_macro    <- numeric(n_sim)

    for (s in seq_len(n_sim)) {
      # --- Generate clean data ---
      X0 <- generate_clean(n)

      # --- Baseline: classical PCA on full clean data ---
      pca_base    <- prcomp(X0, center = TRUE, scale. = FALSE)
      center_base <- colMeans(X0)
      scores_base <- X0 %*% pca_base$rotation[, 1:k, drop = FALSE]
      Xhat_base   <- sweep(
        scores_base %*% t(pca_base$rotation[, 1:k, drop = FALSE]),
        2, center_base, "+"
      )

      # --- Contaminate ---
      X_cont  <- X0
      bad_rows <- integer(0)

      if (frac_cell > 0) {
        X_cont <- add_cellwise(X_cont, frac = frac_cell, gamma = gamma)
      }
      if (frac_row > 0) {
        rw       <- add_rowwise(X_cont, frac = frac_row, gamma = gamma)
        X_cont   <- rw$X
        bad_rows <- rw$bad_rows
      }
      # Add NAs last (so we know which cells are outliers vs missing)
      X_cont <- add_nas(X_cont, frac = 0.20)

      # Clean rows for MSE evaluation
      C_rows <- setdiff(seq_len(n), bad_rows)

      # --- ICPCA ---
      fit_icpca <- tryCatch(
        ICPCA(X_cont, k = k),
        error = function(e) NULL
      )
      if (!is.null(fit_icpca)) {
        Xhat_i  <- compute_predictions(fit_icpca$center,
                                       fit_icpca$loadings, X0)
        mse_icpca[s] <- compute_mse(Xhat_i, Xhat_base, C_rows)
      } else {
        mse_icpca[s] <- NA
      }

      # --- MROBPCA ---
      fit_mrob <- tryCatch(
        MROBPCA(X_cont, k = k),
        error = function(e) NULL
      )
      if (!is.null(fit_mrob)) {
        Xhat_m  <- compute_predictions(fit_mrob$center,
                                       fit_mrob$loadings, X0)
        mse_mrobpca[s] <- compute_mse(Xhat_m, Xhat_base, C_rows)
      } else {
        mse_mrobpca[s] <- NA
      }

      # --- MacroPCA ---
      fit_macro <- tryCatch(
        MacroPCA(X_cont, k = k, DDCpars = list(silent = TRUE)),
        error = function(e) NULL
      )
      if (!is.null(fit_macro)) {
        Xhat_p  <- compute_predictions(fit_macro$center,
                                       fit_macro$loadings, X0)
        mse_macro[s] <- compute_mse(Xhat_p, Xhat_base, C_rows)
      } else {
        mse_macro[s] <- NA
      }
    } # end sim loop

    results$mse[results$gamma == gamma & results$method == "ICPCA"]    <- mean(mse_icpca,   na.rm = TRUE)
    results$mse[results$gamma == gamma & results$method == "MROBPCA"]  <- mean(mse_mrobpca, na.rm = TRUE)
    results$mse[results$gamma == gamma & results$method == "MacroPCA"] <- mean(mse_macro,   na.rm = TRUE)
  } # end gamma loop

  results
}

# Run both scenarios
res_fig7 <- run_scenario("Fig7: 20% NA + 20% cellwise",
                         frac_cell = 0.20, frac_row = 0.00)

res_fig9 <- run_scenario("Fig9: 20% NA + 10% cell + 10% row",
                         frac_cell = 0.10, frac_row = 0.10)

# =============================================================================
# 5. Plots  (line plots of avg MSE vs gamma, one curve per method)
# =============================================================================
method_colors <- c(
  "ICPCA"    = "#E07B54",
  "MROBPCA"  = "#4C8BB5",
  "MacroPCA" = "#2CA02C"
)
method_lines <- c(
  "ICPCA"    = "dashed",
  "MROBPCA"  = "dotted",
  "MacroPCA" = "solid"
)

make_mse_plot <- function(results, title_str, y_lim = NULL) {
  p <- ggplot(results, aes(x = gamma, y = mse,
                           colour = method, linetype = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    scale_colour_manual(values = method_colors) +
    scale_linetype_manual(values = method_lines) +
    labs(
      title    = title_str,
      subtitle = paste0("n = ", n, ", d = ", d, ", k = ", k,
                        ", ", n_sim, " replications, A09 covariance"),
      x        = expression(gamma ~ "(contamination distance)"),
      y        = "Average MSE",
      colour   = "Method",
      linetype = "Method"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")

  if (!is.null(y_lim)) p <- p + coord_cartesian(ylim = y_lim)
  p
}

p_fig7 <- make_mse_plot(
  res_fig7,
  "Figure 7 replication: 20% missing + 20% cellwise outliers"
)

p_fig9 <- make_mse_plot(
  res_fig9,
  "Figure 9 replication: 20% missing + 10% cellwise + 10% rowwise outliers"
)

print(p_fig7)
print(p_fig9)

ggsave("sim_figure7_replication.pdf", p_fig7, width = 7, height = 5)
ggsave("sim_figure9_replication.pdf", p_fig9, width = 7, height = 5)

# =============================================================================
# 6. Combined panel (both scenarios side by side)
# =============================================================================
res_fig7$scenario <- "20% NA + 20% cellwise"
res_fig9$scenario <- "20% NA + 10% cell + 10% row"
res_combined      <- rbind(res_fig7, res_fig9)

p_combined <- ggplot(res_combined,
                     aes(x = gamma, y = mse,
                         colour = method, linetype = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = method_colors) +
  scale_linetype_manual(values = method_lines) +
  facet_wrap(~ scenario, scales = "free_y") +
  labs(
    title    = "Simulation study: ICPCA vs MROBPCA vs MacroPCA",
    subtitle = paste0("A09 covariance, n = ", n, ", d = ", d,
                      ", k = ", k, ", ", n_sim, " MC replications"),
    x        = expression(gamma),
    y        = "Average MSE",
    colour   = "Method",
    linetype = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p_combined)
ggsave("sim_combined_panel.pdf", p_combined, width = 11, height = 5)

cat("\n=== Simulation study complete. Saved three PDFs. ===\n")

# Print summary tables
cat("\n--- Figure 7 results (avg MSE) ---\n")
print(pivot_wider(res_fig7[, c("gamma","method","mse")],
                  names_from = method, values_from = mse))

cat("\n--- Figure 9 results (avg MSE) ---\n")
print(pivot_wider(res_fig9[, c("gamma","method","mse")],
                  names_from = method, values_from = mse))
