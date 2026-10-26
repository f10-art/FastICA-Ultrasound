# ============================================================
# SIMULATED ULTRASOUND SIGNALS: FASTICA
# Complete and corrected analysis pipeline
# ============================================================

# Limpiar el entorno de trabajo antes de empezar
rm(list = ls())

# ============================================================
# 1. PACKAGES
# ============================================================

# Install only if necessary:
# install.packages("fastICA")
# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("signal")

library(fastICA)
library(ggplot2)
library(gridExtra)
library(signal)


# ============================================================
# 2. NORMALIZATION FUNCTIONS
# ============================================================

# Center to mean zero and standard deviation one.
norm_center <- function(x) {
  as.numeric(scale(x, center = TRUE, scale = TRUE))
}

# Rescale to [-1, 1].
norm_pm1 <- function(x) {
  x <- as.numeric(x)
  r <- range(x, na.rm = TRUE)
  
  if (!all(is.finite(r)) || diff(r) == 0) {
    return(rep(0, length(x)))
  }
  
  2 * (x - r[1]) / diff(r) - 1
}

# Align the sign of a component with a reference signal.
# ICA cannot determine the sign of its components.
normalize_sign <- function(x, ref) {
  r <- suppressWarnings(
    cor(x, ref, use = "complete.obs")
  )
  
  if (is.finite(r) && r < 0) {
    return(-x)
  }
  
  x
}

# Calculate common metrics.
calculate_metrics <- function(estimated, reference) {
  
  n <- min(length(estimated), length(reference))
  
  estimated <- as.numeric(estimated[seq_len(n)])
  reference <- as.numeric(reference[seq_len(n)])
  
  valid <- is.finite(estimated) & is.finite(reference)
  
  estimated <- estimated[valid]
  reference <- reference[valid]
  
  if (length(estimated) < 2) {
    return(c(
      Correlation = NA,
      RMSE = NA,
      MAE = NA
    ))
  }
  
  c(
    Correlation = cor(estimated, reference),
    RMSE = sqrt(mean((estimated - reference)^2)),
    MAE = mean(abs(estimated - reference))
  )
}


# ============================================================
# 3. SAMPLING PARAMETERS
# ============================================================

fs <- 50e6
signal_length <- 3500

# Time starts at zero.
t_total <- (seq_len(signal_length) - 1) / fs


# ============================================================
# 4. FILTERING FUNCTIONS
# ============================================================

# Butterworth low-pass filter.
butter_filter <- function(x, fs, cutoff = 2e6, order = 4) {
  Wn <- cutoff / (fs / 2)
  bf <- butter(order, Wn, type = "low")
  as.numeric(filtfilt(bf, x))
}

# Moving average filter.
moving_average_filter <- function(x, n = 5) {
  filt <- stats::filter(
    x,
    rep(1 / n, n),
    sides = 2
  )
  filt <- as.numeric(filt)
  
  bad <- is.na(filt) | is.infinite(filt)
  filt[bad] <- x[bad]
  
  filt
}

# Simple one-dimensional Kalman filter.
kalman_filter <- function(x, Q = 0.01, R = 0.01) {
  n <- length(x)
  x_hat <- numeric(n)
  P <- numeric(n)
  
  x_hat[1] <- x[1]
  P[1] <- 1
  
  for (k in 2:n) {
    x_pred <- x_hat[k - 1]
    P_pred <- P[k - 1] + Q
    
    K <- P_pred / (P_pred + R)
    x_hat[k] <- x_pred + K * (x[k] - x_pred)
    P[k] <- (1 - K) * P_pred
  }
  x_hat
}

# Savitzky-Golay filter.
savitzky_golay_filter <- function(x, p = 3, n = 11) {
  as.numeric(sgolayfilt(x, p = p, n = n))
}

# Replace NA/Inf values with the original samples.
clean_signal <- function(x, original) {
  bad <- is.na(x) | is.infinite(x)
  x[bad] <- original[bad]
  x
}

# Complete filtering pipeline.
full_filter <- function(x, fs) {
  x1 <- butter_filter(x, fs)
  x2 <- moving_average_filter(x1, n = 5)
  x3 <- kalman_filter(x2, Q = 0.01, R = 0.01)
  x4 <- savitzky_golay_filter(x3, p = 3, n = 11)
  
  clean_signal(x4, x)
}


# ============================================================
# 5. FFT PLOTTING FUNCTION
# ============================================================

plot_fft_normalized <- function(
    signal,
    fs,
    title,
    color = "blue",
    max_freq = 4e6) {
  
  signal <- as.numeric(signal)
  n <- length(signal)
  n_half <- floor(n / 2)
  
  fft_values <- abs(fft(signal))[seq_len(n_half)]
  fft_max <- max(fft_values)
  
  if (fft_max > 0) {
    fft_norm <- fft_values / fft_max
  } else {
    fft_norm <- fft_values
  }
  
  freq <- (seq_len(n_half) - 1) * fs / n
  
  df <- data.frame(
    Frequency = freq,
    Magnitude = fft_norm
  )
  
  df <- df[df$Frequency <= max_freq, ]
  
  ggplot(df, aes(x = Frequency, y = Magnitude)) +
    geom_line(color = color, linewidth = 0.5) +
    labs(
      title = title,
      x = "Frequency (Hz)",
      y = "Normalized magnitude"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 9),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7)
    )
}


# ============================================================
# 6. GENERATE THE ORIGINAL SIMULATED SOURCES
# ============================================================

originalUS <- (
  0.5 * sin(2 * pi * 1e6 * t_total) +
    0.3 * sin(2 * pi * 2e6 * t_total)
)

originalRF <- (
  0.7 * sin(2 * pi * 0.5e6 * t_total) +
    0.4 * sin(2 * pi * 1.5e6 * t_total)
)


# ============================================================
# 7. MIX THE ORIGINAL SOURCES AND ADD NOISE
# ============================================================

set.seed(2)

A <- matrix(
  runif(4, min = 0.1, max = 1.0),
  nrow = 2,
  ncol = 2
)

S_original <- cbind(originalUS, originalRF)

mixtures <- S_original %*% t(A)

mixtures <- mixtures +
  matrix(
    rnorm(2 * signal_length, mean = 0, sd = 0.05),
    nrow = signal_length,
    ncol = 2
  )

mixture1 <- as.numeric(mixtures[, 1])
mixture2 <- as.numeric(mixtures[, 2])

cat("\n========================================\n")
cat("ORIGINAL MIXING MATRIX A\n")
cat("========================================\n")
print(A)


# ============================================================
# 8. FILTER AND STANDARDIZE THE MIXTURES
# ============================================================

mixture1_proc <- full_filter(mixture1, fs)
mixture2_proc <- full_filter(mixture2, fs)

mixture1_norm <- norm_center(mixture1_proc)
mixture2_norm <- norm_center(mixture2_proc)

filtered_data <- cbind(
  mixture1_norm,
  mixture2_norm
)

cat("\n========================================\n")
cat("INPUT DIMENSIONS\n")
cat("========================================\n")
cat("Filtered data:", dim(filtered_data), "\n")


# ============================================================
# 9. FASTICA SEPARATION
# ============================================================

set.seed(17)

ica_res <- fastICA(
  filtered_data,
  n.comp = 2,
  row.norm = FALSE,
  fun = "logcosh",
  alpha = 1.5,
  method = "R",
  tol = 0.005
)

S_ica_raw <- as.matrix(ica_res$S)
A_est <- as.matrix(ica_res$A)

# Reconstructed mixtures from FastICA components (X_hat = S %*% t(A))
reconstructed_mixtures <- S_ica_raw %*% t(A_est)
comb1 <- reconstructed_mixtures[, 1]
comb2 <- reconstructed_mixtures[, 2]

cat("\n========================================\n")
cat("FASTICA OUTPUT DIMENSIONS\n")
cat("========================================\n")
cat("Independent components:", dim(S_ica_raw), "\n")
cat("Estimated mixing matrix:", dim(A_est), "\n")


# ============================================================
# 10. MATCH THE ICA COMPONENTS TO THE KNOWN SOURCES
# ============================================================

cor_matrix <- cor(
  S_ica_raw,
  cbind(originalUS, originalRF),
  use = "complete.obs"
)

colnames(cor_matrix) <- c("OriginalUS", "OriginalRF")
rownames(cor_matrix) <- c("ICA_Component1", "ICA_Component2")

cat("\n========================================\n")
cat("CORRELATIONS BEFORE COMPONENT MATCHING\n")
cat("========================================\n")
print(round(cor_matrix, 4))

score_direct <- abs(cor_matrix[1, 1]) + abs(cor_matrix[2, 2])
score_swapped <- abs(cor_matrix[2, 1]) + abs(cor_matrix[1, 2])

if (score_direct >= score_swapped) {
  component_US <- as.numeric(S_ica_raw[, 1])
  component_RF <- as.numeric(S_ica_raw[, 2])
  assignment <- "ICA component 1 -> US; ICA component 2 -> RF"
} else {
  component_US <- as.numeric(S_ica_raw[, 2])
  component_RF <- as.numeric(S_ica_raw[, 1])
  assignment <- "ICA component 2 -> US; ICA component 1 -> RF"
}

S_US <- normalize_sign(component_US, originalUS)
S_RF <- normalize_sign(component_RF, originalRF)

cat("\nSelected assignment:", assignment, "\n")


# ============================================================
# 11. VALIDATION AGAINST THE ORIGINAL SIMULATED SOURCES
# ============================================================

metrics_US <- calculate_metrics(S_US, originalUS)
metrics_RF <- calculate_metrics(S_RF, originalRF)

cat("\n========================================\n")
cat("SOURCE SEPARATION RESULTS\n")
cat("========================================\n")

cat("\nEstimated US component vs original US:\n")
print(round(metrics_US, 6))

cat("\nEstimated RF component vs original RF:\n")
print(round(metrics_RF, 6))

cat("\nCross-correlations after component matching:\n")
cat("US component vs original RF:", round(cor(S_US, originalRF), 4), "\n")
cat("RF component vs original US:", round(cor(S_RF, originalUS), 4), "\n")


# ============================================================
# 12. OPTIONAL LINEAR SCALE-ADJUSTED SOURCE METRICS
# ============================================================

fit_US <- lm(originalUS ~ S_US)
fit_RF <- lm(originalRF ~ S_RF)

S_US_adj <- as.numeric(predict(fit_US))
S_RF_adj <- as.numeric(predict(fit_RF))

metrics_US_adj <- calculate_metrics(S_US_adj, originalUS)
metrics_RF_adj <- calculate_metrics(S_RF_adj, originalRF)

cat("\n========================================\n")
cat("SCALE-ADJUSTED SOURCE METRICS\n")
cat("========================================\n")

cat("\nUS component after linear adjustment:\n")
print(round(metrics_US_adj, 6))

cat("\nRF component after linear adjustment:\n")
print(round(metrics_RF_adj, 6))


# ============================================================
# 14. TEN-WINDOW STABILITY DIAGNOSTIC
# ============================================================

n <- min(length(S_US), length(S_RF), length(originalUS), length(originalRF))

S_US_eval <- S_US[seq_len(n)]
S_RF_eval <- S_RF[seq_len(n)]

US_eval <- originalUS[seq_len(n)]
RF_eval <- originalRF[seq_len(n)]

k <- 10
fold_size <- floor(n / k)

cv_corr_US <- numeric(k)
cv_corr_RF <- numeric(k)

for (i in seq_len(k)) {
  idx_start <- (i - 1) * fold_size + 1
  idx_end <- if (i < k) i * fold_size else n
  idx <- idx_start:idx_end
  
  cv_corr_US[i] <- cor(S_US_eval[idx], US_eval[idx], use = "complete.obs")
  cv_corr_RF[i] <- cor(S_RF_eval[idx], RF_eval[idx], use = "complete.obs")
}

mean_cv_US <- mean(cv_corr_US, na.rm = TRUE)
sd_cv_US <- sd(cv_corr_US, na.rm = TRUE)

mean_cv_RF <- mean(cv_corr_RF, na.rm = TRUE)
sd_cv_RF <- sd(cv_corr_RF, na.rm = TRUE)

cat("\n========================================\n")
cat("LOCALIZED STABILITY: 10 WINDOWS\n")
cat("========================================\n")

cat("\nUS correlations by window:\n")
print(round(cv_corr_US, 4))
cat("US mean correlation =", round(mean_cv_US, 4), "| SD =", round(sd_cv_US, 4), "\n")

cat("\nRF correlations by window:\n")
print(round(cv_corr_RF, 4))
cat("RF mean correlation =", round(mean_cv_RF, 4), "| SD =", round(sd_cv_RF, 4), "\n")


# ============================================================
# 15. RESULTS SUMMARY TABLE
# ============================================================

results_summary <- data.frame(
  Signal = c("US source separation", "RF source separation"),
  Correlation = c(metrics_US["Correlation"], metrics_RF["Correlation"]),
  RMSE_raw = c(metrics_US["RMSE"], metrics_RF["RMSE"]),
  MAE_raw = c(metrics_US["MAE"], metrics_RF["MAE"]),
  RMSE_adjusted = c(metrics_US_adj["RMSE"], metrics_RF_adj["RMSE"]),
  MAE_adjusted = c(metrics_US_adj["MAE"], metrics_RF_adj["MAE"])
)

cat("\n========================================\n")
cat("SUMMARY OF MAIN RESULTS\n")
cat("========================================\n")

results_summary_print <- results_summary
numeric_cols <- vapply(results_summary_print, is.numeric, logical(1))
results_summary_print[numeric_cols] <- lapply(results_summary_print[numeric_cols], round, digits = 6)
print(results_summary_print, row.names = FALSE)


# ============================================================
# 16. PUBLICATION-STYLE PLOT THEME
# ============================================================

theme_pub <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    panel.grid.major = element_line(linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  )


# ============================================================
# 17. DATA FRAME FOR PLOTTING
# ============================================================

S1 <- S_US
S2 <- S_RF

df_signals <- data.frame(
  Time_sec   = t_total,
  Mixture1   = norm_pm1(mixture1_proc),
  Mixture2   = norm_pm1(mixture2_proc),
  Separated1 = norm_pm1(S1),
  Separated2 = norm_pm1(S2),
  Recovered1 = norm_pm1(comb1),
  Recovered2 = norm_pm1(comb2),
  OriginalUS = norm_pm1(originalUS),
  OriginalRF = norm_pm1(originalRF)
)

#=========================================================
# Original sources
#=========================================================

p_originalUS <- ggplot(df_signals, aes(Time_sec, OriginalUS)) +
  geom_line(color = "blue", linewidth = 0.5) +
  labs(title = "(a) Original Source (US)", x = "Time (s)", y = "Amplitude") +
  theme_pub

p_originalRF <- ggplot(df_signals, aes(Time_sec, OriginalRF)) +
  geom_line(color = "blue", linewidth = 0.5) +
  labs(title = "(a) Original Source (RF)", x = "Time (s)", y = "Amplitude") +
  theme_pub

#=========================================================
# Mixtures
#=========================================================

p_mixture1 <- ggplot(df_signals, aes(Time_sec, Mixture1)) +
  geom_line(color = "blue", linewidth = 0.5) +
  labs(title = "(b) Mixture 1", x = "Time (s)", y = "Amplitude") +
  theme_pub

p_mixture2 <- ggplot(df_signals, aes(Time_sec, Mixture2)) +
  geom_line(color = "blue", linewidth = 0.5) +
  labs(title = "(b) Mixture 2", x = "Time (s)", y = "Amplitude") +
  theme_pub

#=========================================================
# FFT
#=========================================================

p_fft_mixture1 <- plot_fft_normalized(mixture1, fs, title = "(c) FFT of Mixture 1", color = "blue")
p_fft_mixture2 <- plot_fft_normalized(mixture2, fs, title = "(c) FFT of Mixture 2", color = "blue")

#=========================================================
# Filtered mixtures
#=========================================================

p_mixture1_filtered <- ggplot(df_signals, aes(Time_sec, Mixture1)) +
  geom_line(color = "blue", linewidth = 0.5) +
  labs(title = "(d) Filtered Mixture 1", x = "Time (s)", y = "Amplitude") +
  theme_pub

p_mixture2_filtered <- ggplot(df_signals, aes(Time_sec, Mixture2)) +
  geom_line(color = "blue", linewidth = 0.5) +
  labs(title = "(d) Filtered Mixture 2", x = "Time (s)", y = "Amplitude") +
  theme_pub

#=========================================================
# FastICA
#=========================================================

p_separated1 <- ggplot(df_signals, aes(Time_sec, Separated1)) +
  geom_line(color = "blue", linewidth = 0.5) +
  labs(title = "(e) Estimated Component 1", x = "Time (s)", y = "Amplitude") +
  theme_pub

p_separated2 <- ggplot(df_signals, aes(Time_sec, Separated2)) +
  geom_line(color = "blue", linewidth = 0.5) +
  labs(title = "(e) Estimated Component 2", x = "Time (s)", y = "Amplitude") +
  theme_pub

#=========================================================
# Reconstruction
#=========================================================

p_recovered1 <- ggplot(df_signals, aes(Time_sec, Recovered1)) +
  geom_line(color = "blue", linewidth = 0.5) +
  labs(title = "(f) Reconstructed Mixture 1", x = "Time (s)", y = "Amplitude") +
  theme_pub

p_recovered2 <- ggplot(df_signals, aes(Time_sec, Recovered2)) +
  geom_line(color = "blue", linewidth = 0.5) +
  labs(title = "(f) Reconstructed Mixture 2", x = "Time (s)", y = "Amplitude") +
  theme_pub

#=========================================================
# Multi-panel figure
#=========================================================

g <- grid.arrange(
  p_originalUS, p_originalRF,
  p_mixture1, p_mixture2,
  p_fft_mixture1, p_fft_mixture2,
  p_mixture1_filtered, p_mixture2_filtered,
  p_separated1, p_separated2,
  p_recovered1, p_recovered2,
  ncol = 2
)


# ============================================================
# 25. SAVE PUBLICATION-QUALITY FIGURES
# ============================================================

ggsave("Figure_Simulated_Results.pdf", g, width = 8.5, height = 12, units = "in", device = cairo_pdf)
ggsave("Figure_Simulated_Results.tiff", g, width = 8.5, height = 12, units = "in", dpi = 600, compression = "lzw", bg = "white")
ggsave("Figure_Simulated_Results.png", g, width = 8.5, height = 12, units = "in", dpi = 600, bg = "white")

cat("\n========================================\n")
cat("SIMULATION COMPLETED\n")
cat("Figures saved:\n")
cat("- Figure_Simulated_Results.pdf\n")
cat("- Figure_Simulated_Results.tiff\n")
cat("- Figure_Simulated_Results.png\n")
cat("========================================\n")


# ============================================================
# 26. ACTIVITY DETECTION METRICS
# ============================================================

evaluar_actividad <- function(estimada, referencia, percentil = 0.75) {
  n <- min(length(estimada), length(referencia))
  
  estimada <- as.numeric(estimada[seq_len(n)])
  referencia <- as.numeric(referencia[seq_len(n)])
  
  valid <- is.finite(estimada) & is.finite(referencia)
  
  estimada <- estimada[valid]
  referencia <- referencia[valid]
  
  abs_estimada <- abs(estimada)
  abs_referencia <- abs(referencia)
  
  umbral_estimada <- as.numeric(quantile(abs_estimada, probs = percentil, na.rm = TRUE))
  umbral_referencia <- as.numeric(quantile(abs_referencia, probs = percentil, na.rm = TRUE))
  
  actividad_estimada <- abs_estimada >= umbral_estimada
  actividad_referencia <- abs_referencia >= umbral_referencia
  
  TP <- sum(actividad_estimada & actividad_referencia)
  TN <- sum(!actividad_estimada & !actividad_referencia)
  FP <- sum(actividad_estimada & !actividad_referencia)
  FN <- sum(!actividad_estimada & actividad_referencia)
  
  precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  sensibilidad <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  especificidad <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  
  data.frame(
    Threshold_estimated = umbral_estimada,
    Threshold_reference = umbral_referencia,
    TP = TP, TN = TN, FP = FP, FN = FN,
    Precision = precision,
    Sensitivity = sensibilidad,
    Specificity = especificidad
  )
}

metrics_activity_US <- evaluar_actividad(estimada = S_US, referencia = originalUS, percentil = 0.75)
metrics_activity_RF <- evaluar_actividad(estimada = S_RF, referencia = originalRF, percentil = 0.75)


# ============================================================
# 27. DISPLAY ACTIVITY DETECTION RESULTS
# ============================================================

cat("\n========================================\n")
cat("ACTIVITY DETECTION METRICS: US\n")
cat("========================================\n")

print(
  transform(
    metrics_activity_US,
    Precision = round(Precision, 4),
    Sensitivity = round(Sensitivity, 4),
    Specificity = round(Specificity, 4)
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("ACTIVITY DETECTION METRICS: RF\n")
cat("========================================\n")

print(
  transform(
    metrics_activity_RF,
    Precision = round(Precision, 4),
    Sensitivity = round(Sensitivity, 4),
    Specificity = round(Specificity, 4)
  ),
  row.names = FALSE
)


# ============================================================
# 28. COMBINED COMPARISON TABLE
# ============================================================

activity_summary <- data.frame(
  Signal = c("US", "RF"),
  Correlation = c(cor(S_US, originalUS), cor(S_RF, originalRF)),
  Precision = c(metrics_activity_US$Precision, metrics_activity_RF$Precision),
  Sensitivity = c(metrics_activity_US$Sensitivity, metrics_activity_RF$Sensitivity),
  Specificity = c(metrics_activity_US$Specificity, metrics_activity_RF$Specificity)
)

cat("\n========================================\n")
cat("SEPARATION AND ACTIVITY DETECTION SUMMARY\n")
cat("========================================\n")

activity_summary[, -1] <- lapply(activity_summary[, -1], round, digits = 4)
print(activity_summary, row.names = FALSE)
