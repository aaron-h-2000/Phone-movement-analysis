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
# Script: Bayesian MANOVA and Bayesian Linear Mixed Models (BLMMs)
# Contains: model fitting, posterior checks, covariance diagnostics,
#           posterior density plots, LOO model fit indices.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. PACKAGES
# -----------------------------------------------------------------------------

packages <- c("psych", "lme4", "afex", "DescTools",
              "blme", "performance", "brms", "MASS",
              "dplyr", "tidyr", "ggplot2", "MVN",
              "biotools", "bayesplot", "corrplot",
              "reshape2", "patchwork", "loo", "rstan")

lapply(packages, library, character.only = TRUE)

# -----------------------------------------------------------------------------
# 2. OUTPUT DIRECTORY
# -----------------------------------------------------------------------------

dir.create("figures/BM_BLMM", recursive = TRUE, showWarnings = FALSE)

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

df <- read.csv("final_extracted_data.csv")

# Rows 465-471 are empty placeholder rows appended by the MATLAB pipeline
# export and contain no valid trial data.
df <- df[-c(465:471), ]

# Variables to process
main_vars <- c("motion_freq_Hz", "motion_rms", "motion_ptp",
               "tilt_freq_Hz", "tilt_rms_rad", "tilt_ptp_rad")

# Winsorize then z-standardize all six kinematic variables
df[main_vars] <- lapply(df[main_vars], winsor, trim = 0.05)
df[main_vars] <- lapply(df[main_vars], scale)

# Factorize categorical variables
df$group       <- factor(df$group)
df$gender      <- factor(df$gender)
df$environment <- factor(df$environment)
df$file        <- factor(df$file)

# Subset to numeric DVs only for assumption checks
df_mardia <- df[, main_vars]

# -----------------------------------------------------------------------------
# 5. ASSUMPTION CHECKS
# -----------------------------------------------------------------------------

# Mardia multivariate normality test
mardia_result <- mvn(df_mardia, mvn_test = "mardia")

# Box's M test for homogeneity of covariance matrices across groups
boxM_result <- boxM(df[, main_vars], df$group)

# -----------------------------------------------------------------------------
# 6. BAYESIAN MANOVA
# -----------------------------------------------------------------------------

multi_formula <-
  bf(motion_freq_Hz ~ group * gender * environment + (1 | file)) +
  bf(motion_rms     ~ group * gender * environment + (1 | file)) +
  bf(motion_ptp     ~ group * gender * environment + (1 | file)) +
  bf(tilt_freq_Hz   ~ group * gender * environment + (1 | file)) +
  bf(tilt_rms_rad   ~ group * gender * environment + (1 | file)) +
  bf(tilt_ptp_rad   ~ group * gender * environment + (1 | file)) +
  set_rescor(TRUE)

bigManova <- brm(
  formula      = multi_formula,
  data         = df,
  family       = student(),
  chains       = 8,
  iter         = 4000,
  cores        = 8,
  seed         = 003,
  control      = list(adapt_delta = 0.99, max_treedepth = 15)
)

# -----------------------------------------------------------------------------
# 7. POSTERIOR DRAWS AND COVARIANCE DIAGNOSTICS
# -----------------------------------------------------------------------------

# Draw 10,000 posterior predictive samples
# Structure: draws x observations x DVs
posterior_draws <- posterior_predict(bigManova, ndraws = 10000)

n_draws     <- dim(posterior_draws)[1]
outcome_dim <- dim(posterior_draws)[3]

# Observed covariance matrix
obs_cov <- cov(df_mardia)

# Mean posterior covariance matrix (averaged across all draws)
mean_cov_posterior <- apply(posterior_draws, 1, function(draw_slice) {
  cov(matrix(draw_slice, ncol = outcome_dim))
}, simplify = FALSE) |>
  Reduce("+", x = _) / n_draws

# --- KL divergence between observed and posterior covariance matrices ---
KL_cov <- function(S1, S2) {
  k     <- ncol(S1)
  invS2 <- solve(S2)
  0.5 * (sum(diag(invS2 %*% S1)) - k + log(det(S2) / det(S1)))
}

kl_result  <- KL_cov(obs_cov, mean_cov_posterior)
kendall_r  <- cor(as.vector(obs_cov), as.vector(mean_cov_posterior), method = "kendall")

cat("KL divergence (obs vs posterior cov):", round(kl_result, 3), "\n")
cat("Kendall correlation (obs vs posterior cov):", round(kendall_r, 3), "\n")

# --- Covariance matrix plots ---

plot_cov_matrix <- function(mat, title_text) {
  df_mat <- as.data.frame(as.table(mat))
  ggplot(df_mat, aes(Var1, Var2, fill = Freq)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                         midpoint = 0, name = "Covariance") +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_paper +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_obs_cov  <- plot_cov_matrix(obs_cov,            "Covariance Matrix — Observed Data")
p_post_cov <- plot_cov_matrix(mean_cov_posterior,  "Covariance Matrix — Posterior Means")

ggsave("figures/BM_BLMM/cov_matrix_observed.png",  p_obs_cov,  width = 6, height = 5, dpi = 300)
ggsave("figures/BM_BLMM/cov_matrix_posterior.png", p_post_cov, width = 6, height = 5, dpi = 300)

# -----------------------------------------------------------------------------
# 8. POSTERIOR DENSITY PLOTS (MANOVA)
# -----------------------------------------------------------------------------

# Extract per-draw posterior means for each DV
# Result: matrix of n_draws x 6
post_means <- sapply(1:outcome_dim, function(v) {
  rowMeans(posterior_draws[, , v])
})
colnames(post_means) <- main_vars

# Build a long data frame for ggplot faceting
obs_means <- colMeans(df_mardia)

post_means_long <- as.data.frame(post_means) |>
  tidyr::pivot_longer(cols = everything(), names_to = "variable", values_to = "posterior_mean")

obs_means_df <- data.frame(
  variable  = main_vars,
  obs_mean  = obs_means
)

var_labels <- c(
  motion_freq_Hz = "Motion Frequency (Hz)",
  motion_rms     = "Motion RMS",
  motion_ptp     = "Motion Peak-to-Peak",
  tilt_freq_Hz   = "Tilt Frequency (Hz)",
  tilt_rms_rad   = "Tilt RMS (rad)",
  tilt_ptp_rad   = "Tilt Peak-to-Peak (rad)"
)

post_means_long$variable <- factor(post_means_long$variable,
                                   levels = main_vars,
                                   labels = var_labels)
obs_means_df$variable    <- factor(obs_means_df$variable,
                                   levels = main_vars,
                                   labels = var_labels)

p_posterior_densities <- ggplot(post_means_long, aes(x = posterior_mean)) +
  geom_density(fill = "steelblue", alpha = 0.4, color = "steelblue4") +
  geom_vline(data = obs_means_df, aes(xintercept = obs_mean),
             color = "firebrick", linewidth = 0.8, linetype = "dashed") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(
    title    = "Posterior Mean Densities — Bayesian MANOVA",
    subtitle = "Red dashed line = observed data mean",
    x        = "Posterior mean value",
    y        = "Density"
  ) +
  theme_paper

ggsave("figures/BM_BLMM/posterior_densities_manova.png",
       p_posterior_densities, width = 10, height = 6, dpi = 300)

# -----------------------------------------------------------------------------
# 9. MANOVA MODEL FIT INDICES
# -----------------------------------------------------------------------------

log_lik_matrix  <- log_lik(bigManova, ndraws = 10000)
mean_log_lik    <- mean(rowSums(log_lik_matrix))
manova_waic     <- waic(log_lik_matrix)
manova_loo      <- loo(log_lik_matrix)

cat("Mean log-likelihood:", round(mean_log_lik, 2), "\n")
print(manova_waic)
print(manova_loo)

# -----------------------------------------------------------------------------
# 10. BAYESIAN LINEAR MIXED MODELS (BLMMs)
# -----------------------------------------------------------------------------

# Shared MCMC settings
blmm_args <- list(family = student(), chains = 4, iter = 4000, cores = 8, seed = 003)

# Fit all six BLMMs
fit_blmm <- function(formula) {
  do.call(brm, c(list(formula = formula, data = df), blmm_args))
}

blmm_models <- list(
  motion_freq    = fit_blmm(motion_freq_Hz ~ group + (1|file)),
  tilt_freq      = fit_blmm(tilt_freq_Hz   ~ group + (1|file)),
  tilt_rms       = fit_blmm(tilt_rms_rad   ~ group + (1|file)),
  tilt_rms_inter = fit_blmm(tilt_rms_rad   ~ group * gender * environment + (1|file)),
  tilt_ptp       = fit_blmm(tilt_ptp_rad   ~ group + (1|file)),
  tilt_ptp_inter = fit_blmm(tilt_ptp_rad   ~ group * gender * environment + (1|file))
)
# Summaries
lapply(blmm_models, summary)

# BLMM trace + posterior plots — saved individually
for (model_name in names(blmm_models)) {
  p <- plot(blmm_models[[model_name]], ask = FALSE)
  ggsave(paste0("figures/BM_BLMM/blmm_plot_", model_name, ".png"),
         plot = p[[1]], width = 10, height = 8, dpi = 300)
}

# LOO fit indices for all BLMMs
blmm_loo <- lapply(blmm_models, loo)

for (n in names(blmm_loo)) {
  cat("\n--- LOO:", n, "---\n")
  print(blmm_loo[[n]])
}

# -----------------------------------------------------------------------------
# 11. PRIOR SPECIFICATION AND PRIOR PREDICTIVE CHECKS
# -----------------------------------------------------------------------------

prior_spec <- c(
  set_prior("normal(0, 1)",         class = "b"),
  set_prior("normal(0, 5)",         class = "Intercept"),
  set_prior("student_t(3, 0, 2.5)", class = "sd"),
  set_prior("gamma(2, 0.1)",        class = "nu")
)

prior_formulas <- list(
  motion_freq    = motion_freq_Hz ~ group + (1|file),
  tilt_freq      = tilt_freq_Hz   ~ group + (1|file),
  tilt_rms       = tilt_rms_rad   ~ group + (1|file),
  tilt_rms_inter = tilt_rms_rad   ~ group * gender * environment + (1|file),
  tilt_ptp       = tilt_ptp_rad   ~ group + (1|file),
  tilt_ptp_inter = tilt_ptp_rad   ~ group * gender * environment + (1|file)
)

prior_args <- list(prior = prior_spec, family = student(),
                   sample_prior = "only", chains = 4, iter = 4000,
                   cores = 8, seed = 3)

fit_prior <- function(formula) {
  do.call(brm, c(list(formula = formula, data = df), prior_args))
}

blmm_priors <- lapply(prior_formulas, fit_prior)

# Prior predictive checks — saved as plots
for (model_name in names(blmm_priors)) {
  p <- pp_check(blmm_priors[[model_name]])
  ggsave(paste0("figures/BM_BLMM/prior_ppcheck_", model_name, ".png"),
         plot = p, width = 7, height = 5, dpi = 300)
}

# -----------------------------------------------------------------------------
# 12. POSTERIOR PREDICTIVE CHECKS (BLMMs)
# -----------------------------------------------------------------------------

for (model_name in names(blmm_models)) {
  p <- pp_check(blmm_models[[model_name]])
  ggsave(paste0("figures/BM_BLMM/posterior_ppcheck_", model_name, ".png"),
         plot = p, width = 7, height = 5, dpi = 300)
}

# =============================================================================
# License: GNU General Public License v3.0 (GPL-3.0)
# See LICENSE file or https://www.gnu.org/licenses/
# =============================================================================