# =============================================================================
# 00_install_packages.R
# Install all required packages for the MacroPCA project
# =============================================================================

required_packages <- c(
  "cellWise",    # MacroPCA, DDC, and the DPOSS dataset
  "ggplot2",     # Plotting
  "dplyr",       # Data manipulation
  "tidyr",       # Data tidying
  "gridExtra",   # Multi-panel plots
  "MASS",        # mvrnorm for simulations
  "robustbase",  # Supplementary robust methods
  "corrplot",    # Correlation plots
  "RColorBrewer" # Color palettes
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing:", pkg))
    install.packages(pkg, dependencies = TRUE)
  } else {
    message(paste("Already installed:", pkg))
  }
}

invisible(lapply(required_packages, install_if_missing))
message("\nAll packages ready.")
