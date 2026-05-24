# =============================================================================
# 02_real_data_analysis.R
#
# Replication of Section 3 (and partially Section 4) of:
#   Hubert et al. (2019) MacroPCA, Technometrics 61, 459-473.
#
# Dataset: Top Gear cars (297 cars x 11 continuous variables)
#          from package robustHD (Alfons 2016)
#
# Reproduces:
#   Figure 3 – Residual maps: ICPCA (left) vs MacroPCA (right)
#   Figure 4 – Outlier maps:  ICPCA (left) vs MacroPCA (right)
#   Figure 5 – Online prediction: including vs excluding 24 selected cars
#
# Comparison: ICPCA (classical, handles NAs) vs MacroPCA (robust)
# Components: k = 2  (as in paper)
# =============================================================================

library(cellWise)    # MacroPCA, ICPCA, cellMap
library(robustHD)    # topgear dataset
library(ggplot2)
library(ggrepel)     # non-overlapping labels on outlier map
library(dplyr)
library(gridExtra)

# =============================================================================
# 1. Load Top Gear data
# =============================================================================
data("topgear", package = "robustHD")

cat("=== Top Gear Dataset ===\n")
cat("Dimensions:", nrow(topgear), "rows x", ncol(topgear), "cols\n")
cat("Variables:", paste(colnames(topgear), collapse = ", "), "\n\n")

# Extract car names (row names) and numeric columns
car_names <- rownames(topgear)

# Keep only the 11 continuous variables used in the paper
cont_vars <- c("Price", "Displacement", "BHP", "Torque",
               "Acceleration", "TopSpeed", "MPG",
               "Weight", "Length", "Width", "Height")

# Verify all variables exist (column names may differ slightly by package version)
available <- intersect(cont_vars, colnames(topgear))
missing_v <- setdiff(cont_vars, colnames(topgear))
if (length(missing_v) > 0) {
  cat("NOTE: Variables not found (check column names):", missing_v, "\n")
  cat("Available columns:", paste(colnames(topgear), collapse = ", "), "\n")
}

X_raw <- as.matrix(topgear[, available])
cat("Using", ncol(X_raw), "variables for", nrow(X_raw), "cars.\n")
cat("Missing cells in raw data:", sum(is.na(X_raw)),
    sprintf("(%.1f%%)\n\n", 100 * mean(is.na(X_raw))))

# =============================================================================
# 2. Log-transform skewed variables  (paper Section 3)
# Five variables are right-skewed and log-transformed:
#   Price, Displacement, BHP, Torque, TopSpeed
# =============================================================================
log_vars <- intersect(c("Price", "Displacement", "BHP", "Torque", "TopSpeed"),
                      colnames(X_raw))

X <- X_raw
for (v in log_vars) {
  X[, v] <- log(X_raw[, v])
}

cat("Log-transformed:", paste(log_vars, collapse = ", "), "\n\n")

# =============================================================================
# 3. Fit ICPCA  (classical iterative PCA that handles NAs)
# =============================================================================
k <- 2    # number of components (paper uses k=2 for Top Gear)

cat("Fitting ICPCA (k =", k, ")...\n")
fit_icpca <- ICPCA(X, k = k)

cat("ICPCA: Cumulative variance explained:",
    round(cumsum(fit_icpca$eigenvalues / sum(fit_icpca$eigenvalues))[1:k], 3), "\n")

# =============================================================================
# 4. Fit MacroPCA
# =============================================================================
cat("Fitting MacroPCA (k =", k, ")...\n")
fit_macro <- MacroPCA(X, k = k, DDCpars = list(silent = TRUE))

cat("MacroPCA: Cumulative variance explained:",
    round(cumsum(fit_macro$eigenvalues / sum(fit_macro$eigenvalues))[1:k], 3), "\n")
cat("MacroPCA: Flagged cellwise outliers:",
    sum(fit_macro$indcells != 0, na.rm = TRUE), "\n")
cat("MacroPCA: Flagged casewise outliers:",
    sum(fit_macro$indrows, na.rm = TRUE), "\n\n")

# =============================================================================
# 5. Residual maps  (Figure 3)
#
# Color scheme (paper):
#   Yellow   : regular cell (|r_ij| <= sqrt(chi^2_{1,0.99}) = 2.576)
#   White    : missing value (NA)
#   Orange-Red: positive outlier (r_ij > 2.576)
#   Purple-Blue: negative outlier (r_ij < -2.576)
# Circle on right: OD_i (white = regular, black = large OD)
# =============================================================================

# Cars shown in the paper's Figure 3 (24 selected rows including notable ones)
notable_cars <- c(
  "Bugatti Veyron", "Pagani Huayra",
  "BMW i3", "Chevrolet Volt", "Vauxhall Ampera", "Mitsubishi i-MiEV",
  "Renault Twizy", "Citroen DS5",
  "Land Rover Defender", "Mercedes-Benz G",
  "Ssangyong Rodius"
)

# Find indices; fall back gracefully if some names differ
note_idx  <- which(car_names %in% notable_cars)
# Top OD cars (most outlying by orthogonal distance)
top_od    <- order(fit_macro$OD, decreasing = TRUE)[1:10]
# Union, capped at 24 rows
sel_rows  <- unique(c(note_idx, top_od))[1:min(24, n_distinct(c(note_idx, top_od)))]
sel_names <- car_names[sel_rows]

cat("Selected", length(sel_rows), "cars for residual map.\n")

# --- MacroPCA residual map ---
# stdResid is the standardized residual matrix R_{n,d} (paper eq. after Step 6)
# indcells marks outlying cells

cellMap(
  d             = fit_macro$stdResid[sel_rows, ],
  indcells      = fit_macro$indcells[sel_rows, ],
  rowlabels     = sel_names,
  columnlabels  = colnames(X),
  mTitle        = "MacroPCA – Residual map (selected cars)",
  sizetitlex    = 9,
  sizetitley    = 8
)
# Save to PDF
pdf("residual_map_macropca.pdf", width = 10, height = 7)
cellMap(
  d             = fit_macro$stdResid[sel_rows, ],
  indcells      = fit_macro$indcells[sel_rows, ],
  rowlabels     = sel_names,
  columnlabels  = colnames(X),
  mTitle        = "MacroPCA – Residual map (Figure 3 right)",
  sizetitlex    = 9,
  sizetitley    = 8
)
dev.off()

# --- ICPCA residual map ---
# Construct standardized residuals for ICPCA manually (same color scheme)
# Predictions from ICPCA
Xhat_icpca <- sweep(
  fit_icpca$scores %*% t(fit_icpca$loadings),
  2, fit_icpca$center, "+"
)

# NA-imputed X (replace NAs with ICPCA imputations)
X_naimputed_icpca <- X
for (j in seq_len(ncol(X))) {
  na_j <- is.na(X[, j])
  X_naimputed_icpca[na_j, j] <- Xhat_icpca[na_j, j]
}

resid_icpca  <- X_naimputed_icpca - Xhat_icpca
# Robust column scale (MAD)
col_mad      <- apply(resid_icpca, 2, function(x) mad(x, na.rm = TRUE))
col_mad[col_mad < 1e-10] <- 1
stdR_icpca   <- sweep(resid_icpca, 2, col_mad, "/")

# ICPCA does not flag cells; use threshold ±2.576 as cutoff
indcells_icpca          <- matrix(0L, nrow(X), ncol(X))
indcells_icpca[stdR_icpca >  2.576 & !is.na(stdR_icpca)] <-  1L
indcells_icpca[stdR_icpca < -2.576 & !is.na(stdR_icpca)] <- -1L

pdf("residual_map_icpca.pdf", width = 10, height = 7)
cellMap(
  d             = stdR_icpca[sel_rows, ],
  indcells      = indcells_icpca[sel_rows, ],
  rowlabels     = sel_names,
  columnlabels  = colnames(X),
  mTitle        = "ICPCA – Residual map (Figure 3 left)",
  sizetitlex    = 9,
  sizetitley    = 8
)
dev.off()
cat("Saved: residual_map_macropca.pdf and residual_map_icpca.pdf\n")

# =============================================================================
# 6. Outlier maps: Score Distance (SD) vs Orthogonal Distance (OD)  (Figure 4)
#
# Quadrant classification (paper Section 3):
#   SD <= cSD, OD <= cOD : Regular
#   SD >  cSD, OD <= cOD : Good leverage point
#   SD <= cSD, OD >  cOD : Orthogonal outlier
#   SD >  cSD, OD >  cOD : Bad leverage point
# =============================================================================

make_outlier_map <- function(SD, OD, cSD, cOD, labels, title_str,
                             label_these = NULL) {
  type <- case_when(
    SD > cSD & OD > cOD ~ "Bad leverage point",
    SD > cSD & OD <= cOD ~ "Good leverage point",
    SD <= cSD & OD > cOD ~ "Orthogonal outlier",
    TRUE                 ~ "Regular"
  )

  df <- data.frame(SD, OD, type, label = labels, stringsAsFactors = FALSE)

  # Which cars to label (notable outliers)
  if (is.null(label_these)) {
    label_these <- labels[type != "Regular"]
  }
  df$show_label <- df$label %in% label_these

  ggplot(df, aes(x = SD, y = OD, colour = type)) +
    geom_point(aes(shape = type), size = 1.8, alpha = 0.75) +
    geom_vline(xintercept = cSD, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = cOD, linetype = "dashed", colour = "grey50") +
    geom_text_repel(
      data = subset(df, show_label),
      aes(label = label),
      size = 2.8, max.overlaps = 20, segment.size = 0.3
    ) +
    scale_colour_manual(values = c(
      "Regular"            = "#AAAAAA",
      "Good leverage point" = "#4C8BB5",
      "Orthogonal outlier" = "#E07B54",
      "Bad leverage point" = "#D62728"
    )) +
    scale_shape_manual(values = c(
      "Regular"            = 1,
      "Good leverage point" = 2,
      "Orthogonal outlier" = 16,
      "Bad leverage point" = 17
    )) +
    labs(
      title    = title_str,
      x        = "Score distance (SD)",
      y        = "Orthogonal distance (OD)",
      colour   = NULL, shape = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}

# Cars to label on the outlier map (those mentioned in paper Section 3)
cars_to_label <- c(
  "BMW i3", "Bugatti Veyron", "Pagani Huayra",
  "Vauxhall Ampera", "Chevrolet Volt", "Renault Twizy",
  "Citroen DS5", "Mitsubishi i-MiEV",
  "Land Rover Defender", "Mercedes-Benz G"
)

# MacroPCA outlier map
p_om_macro <- make_outlier_map(
  SD           = fit_macro$SD,
  OD           = fit_macro$OD,
  cSD          = fit_macro$cutoffSD,
  cOD          = fit_macro$cutoffOD,
  labels       = car_names,
  title_str    = "MacroPCA – Outlier map (Figure 4 right)",
  label_these  = cars_to_label
)

# ICPCA outlier map
# Compute SD and OD for ICPCA manually
# OD: ||◦x_i - x̂_i||
OD_icpca  <- sqrt(rowSums((X_naimputed_icpca - Xhat_icpca)^2, na.rm = TRUE))
# SD: robustified Mahalanobis distance in score space (paper eq. 5)
eig_macro <- fit_macro$eigenvalues   # use MacroPCA eigenvalues as reference
scores_icpca_SD <- fit_icpca$scores
SD_icpca <- sqrt(rowSums(
  sweep(scores_icpca_SD^2, 2, fit_icpca$eigenvalues, "/")
))
# Cutoffs using MCD-based approach (approximate with chi-sq quantiles)
cSD_icpca <- sqrt(qchisq(0.99, df = k))
cOD_icpca <- fit_icpca$cutoffOD   # use ICPCA's own cutoff if available
if (is.null(cOD_icpca)) {
  od23 <- OD_icpca^(2/3)
  cOD_icpca <- (median(od23) + mad(od23) * qnorm(0.99))^(3/2)
}

p_om_icpca <- make_outlier_map(
  SD           = SD_icpca,
  OD           = OD_icpca,
  cSD          = cSD_icpca,
  cOD          = cOD_icpca,
  labels       = car_names,
  title_str    = "ICPCA – Outlier map (Figure 4 left)",
  label_these  = cars_to_label
)

# Side-by-side panel
p_outlier_panel <- grid.arrange(p_om_icpca, p_om_macro, ncol = 2)
ggsave("outlier_map_panel.pdf", p_outlier_panel, width = 14, height = 6)
cat("Saved: outlier_map_panel.pdf\n")

# =============================================================================
# 7. Loadings plot  (Figure 11 equivalent for Top Gear)
# =============================================================================
load_df <- data.frame(
  variable = colnames(X),
  MacroPCA_PC1 = fit_macro$loadings[, 1],
  MacroPCA_PC2 = fit_macro$loadings[, 2],
  ICPCA_PC1    = fit_icpca$loadings[, 1],
  ICPCA_PC2    = fit_icpca$loadings[, 2]
)

load_long <- load_df %>%
  pivot_longer(-variable, names_to = "key", values_to = "loading") %>%
  separate(key, into = c("method", "PC"), sep = "_")

p_loadings <- ggplot(
  load_long,
  aes(x = variable, y = loading, fill = method)
) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.85, width = 0.7) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  scale_fill_manual(values = c("MacroPCA" = "#2CA02C", "ICPCA" = "#E07B54")) +
  facet_wrap(~ PC, ncol = 1) +
  labs(
    title = "Loadings: ICPCA vs MacroPCA (Top Gear, k = 2)",
    x     = NULL, y = "Loading", fill = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave("loadings_comparison.pdf", p_loadings, width = 8, height = 7)
cat("Saved: loadings_comparison.pdf\n")

# =============================================================================
# 8. Online prediction demo  (Section 4 / Figure 5)
#
# Exclude the 24 selected cars, fit MacroPCA on the remaining 273 cars,
# then predict the 24 excluded cars using MacroPCApredict + DDCpredict.
# Compare residual maps: in-sample (left) vs out-of-sample (right).
# =============================================================================
cat("\n=== Online prediction (Section 4 / Figure 5) ===\n")

all_rows   <- seq_len(nrow(X))
train_rows <- setdiff(all_rows, sel_rows)
X_train    <- X[train_rows, ]
X_test     <- X[sel_rows, ]

cat("Training on", length(train_rows), "cars, predicting", length(sel_rows), "\n")

# Fit MacroPCA on training set
fit_train <- MacroPCA(X_train, k = k, DDCpars = list(silent = TRUE))

# Predict the 24 excluded cars one-by-one using MacroPCApredict
pred_results <- vector("list", length(sel_rows))
for (i in seq_along(sel_rows)) {
  pred_results[[i]] <- tryCatch(
    MacroPCApredict(
      Xtrain  = X_train,
      Xnew    = X_test[i, , drop = FALSE],
      MacroOut = fit_train
    ),
    error = function(e) NULL
  )
}

# Build standardized residual matrix for the predicted cars
# (each MacroPCApredict result has $stdResid)
stdR_pred <- do.call(rbind, lapply(pred_results, function(r) {
  if (!is.null(r)) r$stdResid else rep(NA, ncol(X))
}))
rownames(stdR_pred) <- sel_names

# Build indcells for predicted cars
indcells_pred <- do.call(rbind, lapply(pred_results, function(r) {
  if (!is.null(r)) r$indcells else rep(0L, ncol(X))
}))

# Side-by-side: in-sample (from fit_macro) vs out-of-sample (from fit_train)
pdf("online_prediction_comparison.pdf", width = 14, height = 7)
par(mfrow = c(1, 2))

cellMap(
  d             = fit_macro$stdResid[sel_rows, ],
  indcells      = fit_macro$indcells[sel_rows, ],
  rowlabels     = sel_names,
  columnlabels  = colnames(X),
  mTitle        = "In-sample (all 297 cars fitted)",
  sizetitlex    = 9, sizetitley = 8
)

cellMap(
  d             = stdR_pred,
  indcells      = indcells_pred,
  rowlabels     = sel_names,
  columnlabels  = colnames(X),
  mTitle        = "Out-of-sample (24 cars predicted)",
  sizetitlex    = 9, sizetitley = 8
)

dev.off()
cat("Saved: online_prediction_comparison.pdf\n")

# =============================================================================
# 9. Summary table of flagged cars
# =============================================================================
macro_type <- case_when(
  fit_macro$SD > fit_macro$cutoffSD & fit_macro$OD > fit_macro$cutoffOD ~
    "Bad leverage point",
  fit_macro$SD > fit_macro$cutoffSD ~
    "Good leverage point",
  fit_macro$OD > fit_macro$cutoffOD ~
    "Orthogonal outlier",
  TRUE ~ "Regular"
)

flagged_df <- data.frame(
  Car          = car_names,
  SD           = round(fit_macro$SD, 2),
  OD           = round(fit_macro$OD, 2),
  Type         = macro_type,
  stringsAsFactors = FALSE
) %>%
  filter(Type != "Regular") %>%
  arrange(desc(OD))

cat("\n=== MacroPCA: Flagged cars (non-regular) ===\n")
print(flagged_df)

cat("\n=== Real data analysis complete ===\n")
