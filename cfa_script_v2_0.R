# =============================================================================
# MyProject: 3D Motion Capture Analysis in Fine Motor Tasks
# Developing a Methodological Framework with Phone Handling as a Pilot Study.
# Copyright (C) 2026 Author(s)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# Script: PCA, LDA, and Random Forest analyses
# Contains: collinearity screening, KMO/Bartlett checks, parallel analysis,
#           PCA with Varimax rotation, LDA with diagnostics,
#           Random Forest classifiers across three train/test splits.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. PACKAGES
# -----------------------------------------------------------------------------

packages <- c("psych", "DescTools", "MASS", "caret", "randomForest",
              "dplyr", "tidyr", "ggplot2", "MVN", "biotools",
              "corrplot", "reshape2", "patchwork")

lapply(packages, library, character.only = TRUE)

# -----------------------------------------------------------------------------
# 2. OUTPUT DIRECTORY
# -----------------------------------------------------------------------------

dir.create("figures/PCA_LDA_RF", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 3. SHARED PLOT THEME
# -----------------------------------------------------------------------------

theme_paper <- theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 11, color = "grey40"),
    axis.title    = element_text(size = 11),
    axis.text     = element_text(size = 10),
    legend.position = "top",
    legend.title  = element_text(size = 10),
    legend.text   = element_text(size = 9)
  )

# -----------------------------------------------------------------------------
# 4. DATA LOADING AND PREPARATION
# -----------------------------------------------------------------------------

pca_data_full <- read.csv("dfa_aggregated_extended.csv")

# Drop metadata columns (first 6); retain numeric kinematic features only
pca_data <- pca_data_full[, -c(1:6)]

# Winsorize then z-standardize all features
pca_data[] <- lapply(pca_data, winsor, trim = 0.05)
pca_data[] <- lapply(pca_data, scale)

# -----------------------------------------------------------------------------
# 5. COLLINEARITY SCREENING
# -----------------------------------------------------------------------------

cor_result <- corr.test(pca_data, method = "spearman")
cor_mat    <- cor_result$r

# Remove variables with pairwise Spearman r > 0.95
to_remove  <- findCorrelation(cor_mat, cutoff = 0.95, names = TRUE)
pca_reduced <- pca_data[, !(names(pca_data) %in% to_remove)]

# KMO before manual variable selection
KMO(pca_reduced)

# Manual retention of theoretically salient variables (post-KMO review)
vars_to_keep <- c("m_motion_rms", "m_tilt_rms_rad", "max_tilt_freq_Hz",
                  "min_motion_rms", "sd_motion_rms", "sd_motion_ptp",
                  "sd_tilt_freq_Hz", "sd_tilt_rms_rad", "mid_motion_ptp",
                  "mid_tilt_rms_rad", "range_motion_freq_Hz",
                  "range_tilt_freq_Hz", "cv_motion_ptp",
                  "cv_tilt_freq_Hz", "cv_tilt_rms_rad", "cv_tilt_ptp_rad")

pca_reduced <- pca_reduced[, vars_to_keep]

# Final KMO and Bartlett's test on cleaned dataset
KMO(pca_reduced)
psych::cortest.bartlett(pca_reduced)

# -----------------------------------------------------------------------------
# 6. PARALLEL ANALYSIS (COMPONENT RETENTION)
# -----------------------------------------------------------------------------

fa.parallel(pca_reduced, fa = "pc", n.iter = 100, main = "Scree Plot")

# -----------------------------------------------------------------------------
# 7. PCA
# -----------------------------------------------------------------------------

pca1    <- prcomp(pca_reduced, center = TRUE, scale. = TRUE)
pca_rot <- principal(pca_reduced, nfactors = 4, rotate = "varimax")

print(pca_rot$loadings, cutoff = 0.3)

# --- Cumulative variance explained plot ---
var_explained <- pca1$sdev^2 / sum(pca1$sdev^2)
cum_var        <- cumsum(var_explained)

cum_var_df <- data.frame(
  component  = seq_along(cum_var),
  cumulative = cum_var
)

p_cum_var <- ggplot(cum_var_df, aes(x = component, y = cumulative)) +
  geom_point(color = "steelblue", size = 2.5) +
  geom_line(color = "steelblue") +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "firebrick") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    title    = "Cumulative Variance Explained by PCA Components",
    subtitle = "Red dashed line = 80% threshold",
    x        = "Principal Component",
    y        = "Cumulative Proportion of Variance Explained"
  ) +
  theme_paper

ggsave("figures/PCA_LDA_RF/pca_cumulative_variance.png",
       p_cum_var, width = 7, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 8. LDA
# -----------------------------------------------------------------------------

# Attach group label from original full dataset
pca_reduced$group <- factor(pca_data_full$group)

# Assumption checks
lda_mardia <- mvn(pca_reduced[, vars_to_keep], mvn_test = "mardia")
lda_boxM   <- boxM(pca_reduced[, vars_to_keep], pca_reduced$group)

# Fit LDA on first 4 PCA component scores
lda_model <- lda(group ~ pca1$x[, 1:4], data = pca_reduced)

pred <- predict(lda_model)

df_lda <- data.frame(
  Observed     = pca_reduced$group,
  Predicted    = pred$class,
  LD1          = pred$x[, 1],
  Posterior_DI = pred$posterior[, "DI"],
  Posterior_DN = pred$posterior[, "DN"]
) |>
  mutate(
    Confidence = pmax(Posterior_DI, Posterior_DN),
    Correct    = Observed == Predicted
  )

centroids <- df_lda |>
  group_by(Observed) |>
  summarise(LD1 = mean(LD1), .groups = "drop")

boundary <- mean(centroids$LD1)

# --- Plot 1: LD1 group separation ---
p_lda_separation <- ggplot(df_lda, aes(x = LD1, fill = Observed)) +
  geom_density(alpha = 0.4, color = "black") +
  geom_vline(data = centroids, aes(xintercept = LD1, color = Observed),
             linetype = "dashed", linewidth = 1.2) +
  geom_vline(xintercept = boundary, color = "firebrick",
             linetype = "solid", linewidth = 1) +
  labs(
    title    = "LDA Separation along LD1",
    subtitle = "Dashed lines = group centroids | Red line = decision boundary",
    x        = "LD1 discriminant score",
    y        = "Density"
  ) +
  theme_paper

ggsave("figures/PCA_LDA_RF/lda_separation_LD1.png",
       p_lda_separation, width = 7, height = 5, dpi = 300)

# --- Plot 2: Classification confidence by LD1 ---
p_lda_confidence <- ggplot(df_lda, aes(x = LD1, y = Confidence, color = Correct)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey40") +
  geom_vline(xintercept = boundary, linetype = "solid", color = "firebrick") +
  scale_color_manual(
    values = c("TRUE" = "forestgreen", "FALSE" = "firebrick"),
    labels = c("TRUE" = "Correct", "FALSE" = "Incorrect"),
    name   = "Prediction"
  ) +
  labs(
    title    = "Classification Confidence by LD1",
    subtitle = "Points = individual cases | Color = correct vs incorrect prediction",
    x        = "LD1 discriminant score",
    y        = "Model confidence (posterior probability)"
  ) +
  theme_paper

ggsave("figures/PCA_LDA_RF/lda_confidence.png",
       p_lda_confidence, width = 7, height = 5, dpi = 300)

# --- Wilks' lambda and chi-square test ---
eigvals      <- lda_model$svd^2
wilks_lambda <- prod(1 / (1 + eigvals))

n        <- nrow(pca_reduced)
p        <- length(lda_model$scaling[, 1])
n_groups <- length(unique(pca_reduced$group))

chi_sq  <- -((n - 1) - (p + n_groups) / 2) * log(wilks_lambda)
df_chi  <- p * (n_groups - 1)
p_value <- pchisq(chi_sq, df = df_chi, lower.tail = FALSE)

cat("Wilks' lambda:", round(wilks_lambda, 4), "\n")
cat("Chi-square:", round(chi_sq, 3), "| df:", df_chi, "| p =", round(p_value, 4), "\n")

# -----------------------------------------------------------------------------
# 9. RANDOM FOREST
# -----------------------------------------------------------------------------

rf_splits <- list("50_50" = 0.5, "60_40" = 0.6, "70_30" = 0.7)

rf_results <- lapply(names(rf_splits), function(split_name) {
  p <- rf_splits[[split_name]]
  
  train_idx <- createDataPartition(pca_reduced$group, p = p, list = FALSE)
  train     <- pca_reduced[train_idx, ]
  test      <- pca_reduced[-train_idx, ]
  
  rf_model  <- randomForest(group ~ ., data = train,
                             ntree = 500, mtry = 3,
                             proximity = TRUE, importance = TRUE)
  
  preds <- predict(rf_model, newdata = test)
  cm    <- confusionMatrix(preds, test$group)
  
  cat("\n=== Random Forest", split_name, "===\n")
  print(cm)
  
  # Variable importance plot (ggplot)
  imp_df <- as.data.frame(importance(rf_model))
  imp_df$Variable <- rownames(imp_df)
  imp_df <- imp_df[order(imp_df$MeanDecreaseAccuracy, decreasing = TRUE), ]
  imp_df$Variable <- factor(imp_df$Variable, levels = rev(imp_df$Variable))
  
  p_imp <- ggplot(imp_df, aes(x = MeanDecreaseAccuracy, y = Variable)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    geom_col(aes(x = MeanDecreaseGini), fill = "firebrick", alpha = 0.5,
             width = 0.4) +
    labs(
      title    = paste("Variable Importance — RF", split_name),
      subtitle = "Blue = Mean Decrease Accuracy | Red = Mean Decrease Gini",
      x        = "Importance",
      y        = NULL
    ) +
    theme_paper
  
  ggsave(paste0("figures/PCA_LDA_RF/rf_importance_", split_name, ".png"),
         p_imp, width = 7, height = 6, dpi = 300)
  
  # Error over trees plot
  err_df <- data.frame(
    trees = 1:rf_model$ntree,
    OOB   = rf_model$err.rate[, "OOB"],
    DI    = rf_model$err.rate[, "DI"],
    DN    = rf_model$err.rate[, "DN"]
  ) |> tidyr::pivot_longer(-trees, names_to = "class", values_to = "error")
  
  p_err <- ggplot(err_df, aes(x = trees, y = error, color = class)) +
    geom_line(alpha = 0.8) +
    scale_color_manual(values = c("OOB" = "black", "DI" = "steelblue", "DN" = "firebrick"),
                       name = "Class") +
    labs(
      title = paste("Random Forest Error over Trees —", split_name),
      x     = "Number of Trees",
      y     = "Classification Error"
    ) +
    theme_paper
  
  ggsave(paste0("figures/PCA_LDA_RF/rf_error_", split_name, ".png"),
         p_err, width = 7, height = 5, dpi = 300)
  
  list(model = rf_model, confusion = cm, importance = imp_df)
})

names(rf_results) <- names(rf_splits)

# =============================================================================
# License: GNU General Public License v3.0 (GPL-3.0)
# See LICENSE file or https://www.gnu.org/licenses/
# =============================================================================
