# =============================================================================
# 01_simulation_study.R
#
# Partial replication of the simulation study from:
#   Hubert et al. (2019) MacroPCA, Technometrics 61, 459-473.
#
# Goal: Compare Classical PCA vs MacroPCA under three contamination scenarios:
#   (A) No contamination (clean data)
#   (B) Cellwise contamination only
#   (C) Casewise (rowwise) contamination only
#   (D) Mixed contamination (cellwise + casewise)
#
# We evaluate: subspace recovery error (angle between true and estimated PC space)
# =============================================================================

library(cellWise)
library(MASS)
library(ggplot2)
library(dplyr)
library(gridExtra)

set.seed(2024)

# =============================================================================
# Simulation parameters (based on paper Section 4)
# =============================================================================
n       <- 100    # number of observations
p       <- 10     # number of variables
q       <- 2      # number of true principal components
n_sim   <- 200    # number of Monte Carlo repetitions

# True loading matrix: first q columns of a random orthonormal basis
set.seed(42)
V_full  <- svd(matrix(rnorm(p * p), p, p))$u
V_true  <- V_full[, 1:q]   # p x q matrix of true loadings

# Variance of scores along each PC
lambda  <- c(10, 5)

# Contamination parameters
epsilon_cell <- 0.10   # fraction of cells contaminated (cellwise)
epsilon_row  <- 0.10   # fraction of rows contaminated (casewise)
shift_cell   <- 10     # how far cellwise outliers are shifted (in SD units)
shift_row    <- 8      # shift for casewise outlier rows

# =============================================================================
# Helper: generate a clean dataset
# =============================================================================
generate_clean <- function(n, p, q, V_true, lambda, sigma_noise = 1) {
  # Scores: n x q matrix
  scores <- matrix(0, n, q)
  for (j in 1:q) scores[, j] <- rnorm(n, 0, sqrt(lambda[j]))
  # Signal: n x p
  signal <- scores %*% t(V_true)
  # Noise
  noise  <- matrix(rnorm(n * p, 0, sigma_noise), n, p)
  X      <- signal + noise
  return(X)
}

# =============================================================================
# Helper: subspace angle (Frobenius-norm-based) between two loading matrices
# Smaller = better subspace recovery
# =============================================================================
subspace_angle <- function(V1, V2) {
  # Projection matrices
  P1 <- V1 %*% solve(t(V1) %*% V1) %*% t(V1)
  P2 <- V2 %*% solve(t(V2) %*% V2) %*% t(V2)
  norm(P1 - P2, type = "F") / sqrt(2 * q)
}

# =============================================================================
# Helper: apply Classical PCA and MacroPCA, return subspace angle to truth
# =============================================================================
run_one_sim <- function(X, V_true, q) {
  results <- list()

  # --- Classical PCA ---
  pca_fit  <- prcomp(X, center = TRUE, scale. = FALSE)
  V_class  <- pca_fit$rotation[, 1:q, drop = FALSE]
  results$classical <- subspace_angle(V_class, V_true)

  # --- MacroPCA ---
  # MacroPCA(X, k = q) performs the two-step robust PCA
  macro_fit <- tryCatch(
    MacroPCA(X, k = q, DDCpars = list(silent = TRUE)),
    error = function(e) NULL
  )
  if (!is.null(macro_fit)) {
    V_macro <- macro_fit$loadings[, 1:q, drop = FALSE]
    results$macropca <- subspace_angle(V_macro, V_true)
  } else {
    results$macropca <- NA
  }

  return(results)
}

# =============================================================================
# Run simulation across four scenarios
# =============================================================================
scenarios <- c("Clean", "Cellwise", "Casewise", "Mixed")

sim_results <- data.frame(
  scenario  = character(),
  method    = character(),
  angle     = numeric(),
  stringsAsFactors = FALSE
)

for (sc in scenarios) {
  message(paste("Running scenario:", sc, "..."))

  for (s in 1:n_sim) {
    X <- generate_clean(n, p, q, V_true, lambda)

    if (sc == "Cellwise" || sc == "Mixed") {
      # Contaminate a fraction of cells uniformly at random
      n_cont <- round(epsilon_cell * n * p)
      idx    <- sample(n * p, n_cont)
      X[idx] <- X[idx] + rnorm(n_cont, mean = shift_cell, sd = 1)
    }

    if (sc == "Casewise" || sc == "Mixed") {
      # Contaminate a fraction of entire rows
      bad_rows <- sample(n, round(epsilon_row * n))
      X[bad_rows, ] <- X[bad_rows, ] +
        matrix(rnorm(length(bad_rows) * p, mean = shift_row, sd = 1),
               nrow = length(bad_rows))
    }

    res <- run_one_sim(X, V_true, q)

    sim_results <- rbind(sim_results,
      data.frame(scenario = sc, method = "Classical PCA", angle = res$classical),
      data.frame(scenario = sc, method = "MacroPCA",      angle = res$macropca)
    )
  }
}

# =============================================================================
# Summarise and plot
# =============================================================================
sim_results$scenario <- factor(sim_results$scenario, levels = scenarios)

# Summary table
summary_tbl <- sim_results %>%
  group_by(scenario, method) %>%
  summarise(
    mean_angle  = round(mean(angle, na.rm = TRUE), 4),
    sd_angle    = round(sd(angle,   na.rm = TRUE), 4),
    .groups = "drop"
  )

cat("\n=== Subspace Recovery: Mean Angle (lower is better) ===\n")
print(as.data.frame(summary_tbl))

# Boxplot
p_box <- ggplot(sim_results, aes(x = scenario, y = angle, fill = method)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  scale_fill_manual(values = c("Classical PCA" = "#E07B54", "MacroPCA" = "#4C8BB5")) +
  labs(
    title    = "Subspace Recovery: Classical PCA vs MacroPCA",
    subtitle = paste0("n=", n, ", p=", p, ", q=", q, ", ",
                      n_sim, " Monte Carlo replications"),
    x        = "Contamination scenario",
    y        = "Subspace angle (lower = better)",
    fill     = "Method"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")

print(p_box)
ggsave("simulation_subspace_angles.pdf", p_box, width = 8, height = 5)
message("Saved: simulation_subspace_angles.pdf")

# =============================================================================
# Additional: visualise one simulated dataset with cellwise outliers
# =============================================================================
set.seed(7)
X_demo <- generate_clean(50, p, q, V_true, lambda)
n_cont <- round(epsilon_cell * 50 * p)
idx    <- sample(50 * p, n_cont)
X_demo[idx] <- X_demo[idx] + rnorm(n_cont, mean = shift_cell, sd = 1)

# Run DDC (step 1 of MacroPCA) to flag outlying cells
ddc_out <- DDC(X_demo, DDCpars = list(silent = TRUE))

cat("\n=== DDC Flagged Cells Summary ===\n")
cat("Number of outlying cells detected:", sum(ddc_out$indcells), "\n")
cat("Number of outlying rows detected: ", length(ddc_out$outlierIndices), "\n")
