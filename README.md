# Beyond Keystroke Metrics: A 3D Motion Capture Framework for Smartphone Typing Kinematics

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Language: R](https://img.shields.io/badge/Language-R-276DC3.svg)](https://www.r-project.org/)
[![OSF](https://img.shields.io/badge/OSF-Repository-teal.svg)](https://doi.org/10.17605/OSF.IO/9D4R5)
[![DOI](https://zenodo.org/badge/1216376253.svg)](https://doi.org/10.5281/zenodo.19671897)

This repository contains the R statistical analysis pipeline accompanying the paper:

> **Beyond Keystroke Metrics: A 3D Motion Capture Framework for Smartphone Typing Kinematics**
> *[Authors], [Journal], 2026*

The pipeline extracts and analyzes continuous kinematic descriptors from smartphone micro-oscillations during typing, providing a reproducible framework for characterizing tool-mediated motor behavior. Rather than capturing discrete keystroke events, the framework targets the device itself as the measurement object, deriving dominant frequency, RMS amplitude, and robust peak-to-peak range across both translational and rotational movement domains. All processed datasets, MATLAB extraction code, and supplementary materials are available on the [OSF repository](https://osf.io/9d4r5/overview?view_only=98813d25748944cf87f0b20974073e46).

---

## Repository Structure

```
.
├── BM_BLMM_v2_0.R               # Bayesian MANOVA and Bayesian Linear Mixed Models
├── cfa_script_v2_0.R            # PCA, LDA, and Random Forest analyses
├── figures/
│   ├── BM_BLMM/                 # Output figures from BM_BLMM_v2_0.R
│   └── PCA_LDA_RF/              # Output figures from cfa_script_v2_0.R
├── final_pseudo_draft.csv       # Trial-level kinematic data (from MATLAB pipeline)
├── dfa_aggregated_extended.csv  # Participant-level aggregated feature set
└── README.md
```

> **Note:** The two CSV files are not included in this repository. They are available on the OSF repository linked above.

---

## Scripts

### `BM_BLMM_v2_0.R` — Bayesian MANOVA and Linear Mixed Models

This script implements the primary statistical demonstration of pipeline sensitivity. It takes trial-level kinematic data, applies assumption checks, fits a Bayesian MANOVA across all six kinematic variables, and then fits a series of Bayesian Linear Mixed Models (BLMMs) with explicit prior specification, posterior predictive checks, and LOO model fit indices. The two frequency BLMMs (motion frequency and tilt frequency) are the primary sensitivity demonstration; the amplitude models and interaction terms are treated as secondary and reported in supplementary materials.

**Pipeline:**
1. Data loading, winsorization, z-standardization, and factorization
2. Mardia multivariate normality test and Box's M test
3. Bayesian MANOVA (brms, Student-t family, 8 chains × 4000 iterations)
4. Posterior covariance diagnostics — KL divergence and Kendall correlation
5. Posterior density plots for all six kinematic variables
6. Six Bayesian LMMs (group-only and full interaction models)
7. Prior specification and prior predictive checks
8. Posterior predictive checks and LOO fit indices

**Input:** `final_pseudo_draft.csv`
**Output:** Model summaries, LOO indices, figures saved to `figures/BM_BLMM/`

---

### `cfa_script_v2_0.R` — PCA, LDA, and Random Forest

This script implements the dimensionality reduction and proof-of-concept classification analyses. It takes the aggregated participant-level feature set, screens for collinearity, runs a PCA with Varimax rotation to recover the latent structure of the kinematic feature space, performs an LDA as a descriptive group separation visualization, and fits Random Forest classifiers across three train/test splits as an exploratory supplementary check.

**Pipeline:**
1. Data loading, winsorization, and z-standardization
2. Spearman collinearity screening (cutoff r > 0.95)
3. KMO and Bartlett's test of sphericity
4. Parallel analysis for component retention
5. PCA (unrotated) and Varimax-rotated PCA (4 components, 90.4% variance explained)
6. LDA on first 4 PCA component scores with Wilks' lambda
7. Random Forest classifiers at 50/50, 60/40, and 70/30 splits

**Input:** `dfa_aggregated_extended.csv`
**Output:** Model summaries, confusion matrices, figures saved to `figures/PCA_LDA_RF/`

---

## Requirements

All analyses were run in **R 4.5.1**. The following packages are required:

| Script | Packages |
|--------|----------|
| `BM_BLMM_v2_0.R` | `psych`, `lme4`, `afex`, `DescTools`, `blme`, `performance`, `brms`, `MASS`, `dplyr`, `tidyr`, `ggplot2`, `MVN`, `biotools`, `bayesplot`, `corrplot`, `reshape2`, `patchwork`, `loo`, `rstan` |
| `cfa_script_v2_0.R` | `psych`, `DescTools`, `MASS`, `caret`, `randomForest`, `dplyr`, `tidyr`, `ggplot2`, `MVN`, `biotools`, `corrplot`, `reshape2`, `patchwork` |

Install all packages at once:

```r
packages <- c("psych", "lme4", "afex", "DescTools", "blme", "performance",
              "brms", "MASS", "caret", "randomForest", "dplyr", "tidyr",
              "ggplot2", "MVN", "biotools", "bayesplot", "corrplot",
              "reshape2", "patchwork", "loo", "rstan")

install.packages(packages)
```

> **Note on `brms`:** The Bayesian models require a working Stan installation. See the [rstan installation guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) before running `BM_BLMM_v2_0.R`. The MANOVA is computationally intensive (8 chains × 4000 iterations); expect significant runtime on standard hardware.

---

## Reproducibility Notes

- All processing thresholds (winsorization trim = 0.05, collinearity cutoff r > 0.95, missingness ceiling = 25%, minimum usable window = 2s, filter band 1–15 Hz) are defined once and applied uniformly across all recordings.
- Random seeds are set explicitly throughout (`seed = 003` for BLMMs, `seed = 3` for prior models).
- Every trial inclusion and exclusion is traceable through QC flags in the input CSV, produced by the upstream MATLAB pipeline.
- `figures/` subdirectories are created automatically by the scripts if they do not exist.
- The DN/DI cohort contrast used in the demonstration analyses is an illustrative proof-of-concept label and should not be interpreted as evidence of causal effects of digital exposure. It is confounded with age and generational cohort and is used here purely to demonstrate pipeline sensitivity.

---

## Data Availability

All data, processed datasets, MATLAB extraction code, and supplementary materials are openly available at:

**OSF:** [[https://osf.io/9d4r5/overview?view_only=98813d25748944cf87f0b20974073e46]([https://doi.org/10.17605/OSF.IO/9D4R5]
**Zenodo:** [https://doi.org/10.5281/zenodo.19671898](https://doi.org/10.5281/zenodo.19671898)

---

## Citation

If you use this pipeline in your work, please cite:

```
[Authors] (2026). Beyond Keystroke Metrics: A 3D Motion Capture Framework
for Smartphone Typing Kinematics. [Journal].
https://doi.org/10.5281/zenodo.19671898
```

---

## License

This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.
See the [LICENSE](./LICENSE) file or [https://www.gnu.org/licenses/](https://www.gnu.org/licenses/) for details.
