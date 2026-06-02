# MacroPCA — Replication and Extension

---

## Overview

This project replicates and extends **MacroPCA** (Hubert, Rousseeuw & Van den Bossche, *Technometrics* 2019), a unified PCA framework that simultaneously handles missing values, cellwise outliers, and rowwise outliers. We compare three methods throughout:

| Method | Handles NAs | Handles cellwise | Handles rowwise | Works when d > n |
|---|---|---|---|---|
| **ICPCA** | ✓ | ✗ | ✗ | ✓ |
| **MROBPCA** | ✓ | ✗ | ✓ | ✗ |
| **MacroPCA** | ✓ | ✓ | ✓ | ✓ |

The project covers three experimental settings:
1. **Simulation** — replication of Figures 7 and 9 from the paper
2. **Top Gear dataset** — replication of Figures 3–5 from the paper, plus a contamination stress-test
3. **World Bank Development Indicators** — original extension to a macroeconomic panel

---

## Repository Structure

```
.
├── macropcatopgear_fixed.ipynb   # Top Gear notebook (simulation + real data + contamination)
├── MPCAfinal_fixed.ipynb         # World Bank notebook (real data + contamination)
├── report.tex                    # LaTeX report source
└── README.md                     # This file
```

Both notebooks are self-contained Kaggle R kernels. Each cell is independent and should be run top to bottom within its notebook.

---

## Notebooks

### `macropcatopgear_fixed.ipynb` — Top Gear

**Cell 1 — Package installation**
Installs `cellWise`, `robustHD`, `ggplot2`, `ggrepel`, `dplyr`, `tidyr`, `gridExtra`, `reshape2`. Safe to skip if packages are already present.

**Cell 2 — Simulation study** *(replicating Figures 7 and 9)*
Generates synthetic data under the A09 covariance structure (n=100, d=200, k=6) and runs 30 Monte Carlo replications across two contamination settings:
- Figure 7: 20% MCAR missing + 20% cellwise outliers, γ ∈ {0,1,2,3,5,7,10,15,20}
- Figure 9: 20% MCAR missing + 10% cellwise + 10% rowwise outliers, same γ grid

MROBPCA is included in the loop but returns NaN throughout — this is expected and documented as an empirical result (MCD requires h > d, which fails here since d=200 > n=100).

Saves:
- `sim_figure7_replication.pdf` — Figure 7 MSE curve
- `sim_figure9_replication.pdf` — Figure 9 MSE curve
- `sim_combined_panel.pdf` — both panels side by side

**Cell 3 — Real data analysis + Online prediction** *(replicating Figures 3, 4, 5)*
Loads the Top Gear dataset from `robustHD`, fits MacroPCA and ICPCA with k=2, then produces all diagnostic plots. Also runs the online prediction experiment from Section 4 of the paper.

Saves:
- `residual_map_macropca.pdf` — MacroPCA residual map, 13 notable cars (Figure 3 right)
- `residual_map_icpca.pdf` — ICPCA residual map, same cars (Figure 3 left)
- `outlier_map_panel.pdf` — ICPCA and MacroPCA outlier maps side by side (Figure 4)
- `loadings_comparison.pdf` — bar chart comparing PC1/PC2 loadings of both methods
- `online_prediction_comparison.pdf` — in-sample vs out-of-sample residual maps (Figure 5); the 24 display cars are held out, MacroPCA is refit on the remaining ~271 cars, and `MacroPCApredict` scores the held-out set

**Cell 4 — Contamination stress-test** *(original extension)*
Injects controlled contamination into Top Gear: 10% cellwise (γ=8), 10% rowwise, and 5% extra MCAR missingness on top of the baseline 3.2%. Evaluates subspace recovery and outlier detection against known ground truth.

Saves:
- `contamination_summary.pdf` — subspace angle bar (ICPCA vs MacroPCA) + detection rate bars
- `contamination_residual_map_macropca.pdf` — MacroPCA residual map on contaminated data
- `contamination_residual_map_icpca.pdf` — ICPCA residual map on contaminated data
- `contamination_outlier_map.pdf` — MacroPCA outlier map coloured by ground-truth label

---

### `MPCAfinal_fixed.ipynb` — World Bank

**Cell 1 — Package installation**
Same packages as the Top Gear notebook.

**Cell 2 — List input files**
Runs `list.files("/kaggle/input", recursive=TRUE)` so you can confirm the dataset path before the main cells execute.

**Cell 3 — Real data analysis + Online prediction**
Loads the World Bank Development Indicators CSV from `/kaggle/input` (auto-detected), aggregates duplicate country entries, drops columns/rows with >50% missingness, log-transforms right-skewed indicators, and fits MacroPCA and ICPCA with k=2. Includes online prediction analogous to the Top Gear analysis.

Saves:
- `residual_map_macropca.pdf` — MacroPCA residual map, top flagged countries
- `residual_map_icpca.pdf` — ICPCA residual map, same selection
- `outlier_map_panel.pdf` — ICPCA and MacroPCA outlier maps side by side
- `loadings_comparison.pdf` — PC1/PC2 loadings comparison
- `online_prediction_comparison.pdf` — in-sample vs out-of-sample residual maps

**Cell 4 — Contamination stress-test**
Same protocol as Top Gear: 10% cellwise + 10% rowwise + 5% extra NAs injected into the World Bank matrix.

Saves:
- `contamination_summary.pdf`
- `contamination_residual_map_macropca.pdf`
- `contamination_residual_map_icpca.pdf`
- `contamination_outlier_map.pdf`

---

## Output Guide (Top Gear PDF pages)

The merged Top Gear output PDF has 11 pages:

| Page | Plot | Experiment |
|---|---|---|
| 1 | Figure 9 replication — MSE vs γ (20% NA + 10% cell + 10% row) | Simulation |
| 2 | Figure 7 replication — MSE vs γ (20% NA + 20% cellwise) | Simulation |
| 3 | Combined simulation panel (both settings side by side) | Simulation |
| 4 | MacroPCA residual map — 13 notable cars (Figure 3 right) | Normal analysis |
| 5 | MacroPCA outlier map on contaminated data, coloured by ground-truth label | Contamination |
| 6 | ICPCA contaminated residual map — diffuse spread across all cars | Contamination |
| 7 | ICPCA residual map — 13 notable cars (Figure 3 left) | Normal analysis |
| 8 | Outlier maps panel: ICPCA (left) and MacroPCA (right) (Figure 4) | Normal analysis |
| 9 | Loadings comparison bar chart — ICPCA vs MacroPCA, PC1 and PC2 | Normal analysis |
| 10 | Subspace angle bar + detection performance bars (sensitivity, precision, F1) | Contamination |
| 11 | MacroPCA contaminated residual map — contamination sharply localised | Contamination |

`online_prediction_comparison.pdf` (Figure 5 replication) is a separate file not included in the merged PDF.

---

## Key Results

### Simulation
| Setting | ICPCA (γ=20) | MROBPCA | MacroPCA (γ=20) |
|---|---|---|---|
| Fig 7: 20% NA + 20% cellwise | MSE = 9.71 | NaN (d > n) | MSE = 0.026 |
| Fig 9: 20% NA + 10% cell + 10% row | MSE = 2.06 | NaN (d > n) | MSE = 0.031 |

MacroPCA's MSE remains below 0.05 across the entire γ grid. MROBPCA fails throughout because MCD requires h > d; with d=200 and n=100 the covariance estimate is singular at every iteration.

### Top Gear (normal analysis)
- MacroPCA retains k=2 components explaining **82.2%** of variance on the clean core; ICPCA reports 89.7% but its loadings are pulled by the hybrid/electric leverage points.
- MacroPCA flags **124 cellwise outliers** and 0 rowwise outliers on the 295-car matrix.
- The notable groups: supercars (Bugatti Veyron, Pagani Huayra) show extreme Price/BHP cells; EVs and hybrids (BMW i3, Chevrolet Volt, Vauxhall Ampera, Mitsubishi i-MiEV, Renault Twizy) show high MPG and zero/missing Displacement; all-terrain vehicles (Land Rover Defender, Mercedes-Benz G-Class) show high Weight/Height and low TopSpeed.

### Top Gear (contamination stress-test)
| Metric | ICPCA | MacroPCA |
|---|---|---|
| Subspace angle | 0.2875 | **0.1504** |
| Cellwise sensitivity | — | 0.95 |
| Cellwise precision | — | 0.46 |
| Cellwise F1 | — | 0.62 |
| Rowwise specificity | — | 1.00 |

MacroPCA reduces the subspace angle by ~48% relative to ICPCA. The moderate precision reflects that the baseline Top Gear data already contains structural cellwise outliers (EVs with zero displacement, supercars with extreme price) which DDC correctly flags but which count as false positives against the synthetic ground truth.

### World Bank
- Matrix: 253 country/aggregate entities × 45 indicators, 13.7% missing.
- MacroPCA flags **1,310 cellwise outliers**, 0 rowwise outliers.
- ICPCA reports 100% variance explained on k=2 because it is dominated by the size leverage of regional aggregates ("World", "IDA & IBRD total", etc.).
- MacroPCA separates size-driven anomalies (raw GDP/population aggregates flagged on economic scale indicators) from governance-driven anomalies (Afghanistan, Zimbabwe, Venezuela flagged on rule-of-law indicators) — a meaningful, domain-interpretable decomposition invisible to ICPCA.

---

## Technical Notes

### Known limitations / design choices

**MROBPCA NaN outputs** are intentional. The simulation is set up with d=200 > n=100, which violates the MCD requirement h > d. Rather than switching to a lower-dimensional simulation (which would differ from the paper), we keep the paper's exact setup and document the failure as a result. The notebook prints a dedicated summary: `Figure 7 setting: 9 of 9 gamma values produced NaN/NA`.

**Online prediction DDC tolerances.** The default DDC implementation in `cellWise` rejects subsets with low column MAD or discrete-looking variables. Both `MacroPCA()` and `MacroPCApredict()` are called with `DDCpars = list(precScale=1e-12, numDiscrete=0, fracNA=1.0, silent=TRUE)` to disable these guards for the held-out sets.

**Car names.** The `robustHD::TopGear` dataset uses integer rownames, not car names. Names are constructed as `paste(Maker, Model)` with a regex fallback, so all labels in the plots are human-readable.

**pretty_cellMap.** `cellWise::cellMap` does not include a continuous colour-scale legend. Both notebooks use a custom `pretty_cellMap()` built on `ggplot2` with `scale_fill_gradient2` (blue–yellow–red, capped at ±4 standardised residual units).

**World Bank CSV path.** The notebook auto-detects any CSV under `/kaggle/input` using `endsWith(tolower(all_files), ".csv")`, preferring files whose path contains "world", "bank", "indicator", "wdi", or "develop". A diagnostic `list.files` block at the top of cells 3 and 4 prints all available files so path issues are immediately visible.

---

## Dependencies

```r
install.packages(c("cellWise", "robustHD", "ggplot2", "ggrepel",
                   "dplyr", "tidyr", "gridExtra", "reshape2", "readr"))
```

All packages are available on CRAN. `cellWise` (≥ 2.1.0) provides `MacroPCA`, `MacroPCApredict`, `ICPCA`, and `DDC`. `robustHD` provides the Top Gear dataset.

---

## Reference

Hubert, M., Rousseeuw, P. J., & Van den Bossche, W. (2019). MacroPCA: An All-in-One PCA Method Allowing for Missing Values as Well as Cellwise and Rowwise Outliers. *Technometrics*, 61(4), 459–473. https://doi.org/10.1080/00401706.2018.1562989
