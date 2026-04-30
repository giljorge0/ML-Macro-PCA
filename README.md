# MacroPCA – Project 2 (M4ML 2025-2026)

**Topic 3:** MacroPCA: An all-in-one PCA method allowing for missing values as well as cellwise and rowwise outliers

**Reference paper:**
Hubert, M., Rousseeuw, P.J., Van den Bossche, W. (2019). *MacroPCA: An all-in-one PCA method allowing for missing values as well as cellwise and rowwise outliers.* Technometrics, 61, 459–473.
[Open-access PDF](https://wis.kuleuven.be/stat/robust/papers/publications-2019/hubertetal-macropca-technometrics-open-access.pdf)

**R package:** `cellWise` (CRAN)

---

## Project structure

```
macropca_project/
├── README.md
├── R/
│   ├── 00_install_packages.R          # Package installation
│   ├── 01_simulation_study.R          # Replication of simulation study (Task 2)
│   ├── 02_real_data_analysis.R        # Real data application – clean data (Task 3)
│   └── 03_contamination_analysis.R    # Contamination with outliers & missing values (Task 4)
├── data/
│   └── (datasets downloaded/generated at runtime)
└── report/
    └── report.pdf                     # Final 8-page report (to be added)
```

---

## Tasks covered

| Script | Task |
|--------|------|
| `01_simulation_study.R` | Replication of a simulation comparing classical PCA vs MacroPCA under cellwise and rowwise contamination |
| `02_real_data_analysis.R` | Application to a real dataset: outlier detection, comparison of classical PCA vs MacroPCA, diagnostic plots |
| `03_contamination_analysis.R` | Artificial introduction of 10% outliers + missing values, repeated analysis, comparison with known ground truth |

---

## How to run

```r
# 1. Install dependencies
source("R/00_install_packages.R")

# 2. Run in order
source("R/01_simulation_study.R")
source("R/02_real_data_analysis.R")
source("R/03_contamination_analysis.R")
```

All scripts are self-contained and save plots to the working directory.

---

## Key concepts

- **Casewise outlier:** An entire row (observation) that deviates from the bulk of the data.
- **Cellwise outlier:** An individual cell (a single variable of a single observation) that is anomalous, even if the rest of that row looks normal.
- **MacroPCA:** A two-step robust PCA that first flags and imputes cellwise outliers, then performs a robust PCA robust to casewise outliers, also handling missing values.
- **Classical PCA:** Standard PCA based on the sample covariance matrix – highly sensitive to both types of outliers.

---

## References

- Hubert, M., Rousseeuw, P.J., Van den Bossche, W. (2019). MacroPCA. *Technometrics*, 61, 459–473.
- Raymaekers, J., & Rousseeuw, P.J. (2023). `cellWise`: Analyzing Data with Cellwise Outliers. CRAN.
- Rousseeuw, P.J., Van den Bossche, W. (2018). Detecting deviating data cells. *Technometrics*, 60, 135–145.
- Hubert, M., Raymaekers, J., Rousseeuw, P.J. (2024). Robust discriminant analysis. *WIREs Computational Statistics*.
