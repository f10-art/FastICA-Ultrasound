# --- Libraries ---  
# This section loads the necessary R packages for the script.  
# Each package provides specific functionalities:  
#   - fastICA: Implements Independent Component Analysis (ICA) using a fast, fixed-point algorithm.  
#   - ggplot2: A powerful and flexible system for creating publication-quality graphics.  
#   - gridExtra: Provides functions to arrange multiple ggplot2 plots into multi-panel composite figures.  
#   - signal: Provides digital signal processing tools (e.g., Butterworth and Savitzky-Golay filters).  

library(fastICA)  
library(ggplot2)  
library(gridExtra)  
library(signal)  

# --- Explicit Normalization and Centering ---
# Centers the signal to zero mean and scales it to unit variance.
# Required by FastICA to ensure numerical stability and proper convergence.
norm_center <- function(x) {
  as.numeric(scale(x, center = TRUE, scale = TRUE))
}

# --- Normalization to Range [-1, 1] ---
# Scales a signal to the interval [-1, 1]. Used primarily for visualization.
norm_pm1 <- function(x) {  
  x <- as.numeric(x)
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) return(rep(0, length(x)))
  (2 * (x - r[1]) / diff(r)) - 1  
}  

# Sampling frequency of the ultrasonic acquisition setup (Hz)
fs <- 50000000  

# --- Butterworth Low-Pass Filter ---
# Applies zero-phase filtering (filtfilt) to prevent phase distortion.
butter_filter <- function(x, fs, cutoff = 2e6, order = 4) {  
  Wn <- cutoff / (fs / 2)  
  bf <- butter(order, Wn, type = "low")  
  filtfilt(bf, x)  
} 

# --- One-Dimensional Kalman Filter ---
# Adaptive state estimation filter to attenuate high-frequency noise.
kalman_filter <- function(x, Q = 0.01, R = 0.01) {  
  n <- length(x)  
  x_hat <- numeric(n)  
  P <- numeric(n)  
  x_hat[1] <- x[1]; P[1] <- 1  
  for (k in 2:n) {  
    x_pred <- x_hat[k-1]  
    P_pred <- P[k-1] + Q  
    K <- P_pred / (P_pred + R)  
    x_hat[k] <- x_pred + K * (x[k] - x_pred)  
    P[k] <- (1 - K) * P_pred  
  }  
  x_hat  
}  

# --- Savitzky–Golay Filter ---
# Polynomial smoothing filter that preserves peak heights and waveform morphology.
savitzky_golay_filter <- function(x, p = 3, n = 11) {  
  sgolayfilt(x, p = p, n = n)  
}  

# --- Moving Average Filter ---
moving_average_filter <- function(x, n = 5) {  
  filt <- stats::filter(x, rep(1 / n, n), sides = 2)  
  filt <- as.numeric(filt)
  filt[is.na(filt)] <- x[is.na(filt)]  
  filt  
}  

# --- Signal Cleaning Helper ---
# Replaces invalid (NA/Inf) numerical values with original raw samples.
clean <- function(x, original) {  
  x[is.na(x) | is.infinite(x)] <- original[is.na(x) | is.infinite(x)]  
  x  
}  

# --- Sequential Preprocessing Filtering Pipeline ---
# Pipeline order: Butterworth -> Moving Average -> Kalman -> Savitzky-Golay
full_filter <- function(x, fs) {  
  x1 <- butter_filter(x, fs)  
  x2 <- moving_average_filter(x1)  
  x3 <- kalman_filter(x2)  
  x4 <- savitzky_golay_filter(x3)  
  clean(x4, x)  
}  

# --- Normalized Fast Fourier Transform Plot ---
plot_fft_normalized <- function(signal, fs, title, color = "blue", max_freq = 4e6) {
  n <- length(signal)
  n_half <- floor(n / 2)
  
  freq <- (seq_len(n_half) - 1) * fs / n
  fft_values <- Mod(fft(signal))[seq_len(n_half)]
  max_fft <- max(fft_values)
  
  fft_norm <- if (max_fft > 0) fft_values / max_fft else rep(0, length(fft_values))
  
  df <- data.frame(Frequency = freq, Magnitude = fft_norm)
  df <- df[df$Frequency <= max_freq, ]
  
  ggplot(df, aes(x = Frequency, y = Magnitude)) +
    geom_line(color = color) +
    labs(title = title, x = "Frequency (Hz)", y = "Normalized Magnitude") +
    theme_minimal() +  
    theme(plot.title = element_text(size = 8, face = "bold"),  
          axis.title = element_text(size = 8),  
          axis.text = element_text(size = 6))  
} 

# --- Data Ingestion ---
mixture1   <- as.numeric(unlist(read.table("mezcla1.txt",   sep = ",", dec = ".", header = FALSE)))
mixture2   <- as.numeric(unlist(read.table("mezcla2.txt",   sep = ",", dec = ".", header = FALSE)))
originalUS <- as.numeric(unlist(read.table("originalUS.txt", sep = ",", dec = ".", header = FALSE)))
originalRF <- as.numeric(unlist(read.table("originalRF.txt", sep = ",", dec = ".", header = FALSE)))

# --- Input Validation Diagnostic ---
cat("\n========================================\n")
cat("INPUT SIGNAL LENGTHS\n")
cat("========================================\n")
cat("Mixture 1   :", length(mixture1), "\n")
cat("Mixture 2   :", length(mixture2), "\n")
cat("Original US :", length(originalUS), "\n")
cat("Original RF :", length(originalRF), "\n")

if (length(mixture1) != length(mixture2)) {
  stop("Mixture 1 and Mixture 2 must have identical lengths.")
}

if (any(!is.finite(mixture1)) || any(!is.finite(mixture2)) ||
    any(!is.finite(originalUS)) || any(!is.finite(originalRF))) {
  stop("Input vectors contain NA, NaN, or Infinite values.")
}

if (length(mixture1) < 11) {
  stop("Signal length is too short for digital filtering.")
}

# --- Signal Filtering and Preprocessing ---  
mixture1_proc <- full_filter(mixture1, fs)  
mixture2_proc <- full_filter(mixture2, fs)  

# --- Standardization and Centering for FastICA Ingestion ---  
mixture1_norm <- norm_center(mixture1_proc)  
mixture2_norm <- norm_center(mixture2_proc)  
filtered_data <- cbind(mixture1_norm, mixture2_norm)  

# --- FastICA Decomposition ---
set.seed(17)  
ica_res <- fastICA(filtered_data, n.comp = 2, row.norm = FALSE, fun = "logcosh", alpha = 1.5, method = "R", tol = 0.005)  

# Extract Raw Independent Components
S1 <- ica_res$S[,1]  
S2 <- ica_res$S[,2]  

# --- Sign Ambiguity Correction against Physical Reference Signals ---
normalize_sign <- function(x, ref) {  
  c_val <- suppressWarnings(cor(x, ref, use = "complete.obs"))
  if (is.finite(c_val) && c_val < 0) return(-x) else return(x)  
}  

# Ensure component polarity matches experimental reference signals
S1 <- normalize_sign(S1, originalUS)  
S2 <- normalize_sign(S2, originalRF)  

# --- Linear Combination / Signal Reconstruction ---
w11 <- 0.7; w12 <- 0.2  
w21 <- -0.9; w22 <- 0.1  
comb1 <- w11 * S1 + w12 * S2  
comb2 <- w21 * S1 + w22 * S2  

# ============================================================
# PERFORMANCE METRICS: ORIGINAL AND LINEARLY ADJUSTED SIGNALS
# Correlation, Root Mean Square Error (RMSE), Mean Absolute Error (MAE)
# ============================================================

# Mixture 1 Reconstructive Evaluation
rmse_original1 <- sqrt(mean((comb1 - mixture1)^2))
mae_original1  <- mean(abs(comb1 - mixture1))
cor_original1  <- cor(comb1, mixture1)

fit1 <- lm(mixture1 ~ comb1)
comb1_adj <- as.numeric(predict(fit1))

rmse_adjusted1 <- sqrt(mean((comb1_adj - mixture1)^2))
mae_adjusted1  <- mean(abs(comb1_adj - mixture1))
cor_adjusted1  <- cor(comb1_adj, mixture1)

# Mixture 2 Reconstructive Evaluation
rmse_original2 <- sqrt(mean((comb2 - mixture2)^2))
mae_original2  <- mean(abs(comb2 - mixture2))
cor_original2  <- cor(comb2, mixture2)

fit2 <- lm(mixture2 ~ comb2)
comb2_adj <- as.numeric(predict(fit2))

rmse_adjusted2 <- sqrt(mean((comb2_adj - mixture2)^2))
mae_adjusted2  <- mean(abs(comb2_adj - mixture2))
cor_adjusted2  <- cor(comb2_adj, mixture2)

# --- Print Reconstruction Performance Diagnostics ---
cat("\n========================================\n")
cat("MIXTURE 1 RECONSTRUCTION METRICS\n")
cat("========================================\n")
cat("Original Correlation :", round(cor_original1, 4), "\n")
cat("Original RMSE        :", round(rmse_original1, 6), "\n")
cat("Original MAE         :", round(mae_original1, 6), "\n")
cat("Adjusted Correlation :", round(cor_adjusted1, 4), "\n")
cat("Adjusted RMSE        :", round(rmse_adjusted1, 6), "\n")
cat("Adjusted MAE         :", round(mae_adjusted1, 6), "\n")

cat("\n========================================\n")
cat("MIXTURE 2 RECONSTRUCTION METRICS\n")
cat("========================================\n")
cat("Original Correlation :", round(cor_original2, 4), "\n")
cat("Original RMSE        :", round(rmse_original2, 6), "\n")
cat("Original MAE         :", round(mae_original2, 6), "\n")
cat("Adjusted Correlation :", round(cor_adjusted2, 4), "\n")
cat("Adjusted RMSE        :", round(rmse_adjusted2, 6), "\n")
cat("Adjusted MAE         :", round(mae_adjusted2, 6), "\n")

# ============================================================
# LOCAL WINDOWED CORRELATION ANALYSIS (k = 10)
# ============================================================
k <- 10

# US Active Signal Window (0.01625 s - 0.01650 s)
n_US <- min(length(S1), length(originalUS))
t_US <- (seq_len(n_US) - 1) / fs
idx_US <- which(t_US >= 0.01625 & t_US <= 0.01650)

S1_US_window <- S1[idx_US]
US_window <- originalUS[idx_US]

# RF Active Signal Window (0.00000 s - 0.01500 s)
n_RF <- min(length(S2), length(originalRF))
t_RF <- (seq_len(n_RF) - 1) / fs
idx_RF <- which(t_RF >= 0 & t_RF <= 0.01500)

S2_RF_window <- S2[idx_RF]
RF_window <- originalRF[idx_RF]

# Contiguous Non-Overlapping Windowed Pearson Correlation Calculator
cor_por_ventanas <- function(x, y, k = 10) {
  n <- min(length(x), length(y))
  if (n < k * 3) {
    stop("Insufficient samples to divide into k windows.")
  }
  
  x <- as.numeric(x[seq_len(n)])
  y <- as.numeric(y[seq_len(n)])
  
  indices <- split(seq_len(n), cut(seq_len(n), breaks = k, labels = FALSE))
  
  vapply(indices, function(idx) {
    valid <- is.finite(x[idx]) & is.finite(y[idx])
    xx <- x[idx][valid]
    yy <- y[idx][valid]
    
    if (length(xx) < 3 || sd(xx) == 0 || sd(yy) == 0) {
      return(NA_real_)
    }
    cor(xx, yy, method = "pearson")
  }, numeric(1))
}

cv_corr_S1_US <- cor_por_ventanas(S1_US_window, US_window, k)
cv_corr_S2_RF <- cor_por_ventanas(S2_RF_window, RF_window, k)

cat("\n========================================\n")
cat("WINDOWED CORRELATION DIAGNOSTICS (ACTIVE INTERVALS)\n")
cat("========================================\n")
cat("\nS1 vs US Burst Interval (0.01625–0.01650 s):\n")
cat("Correlations by window:", round(cv_corr_S1_US, 4), "\n")
cat("Mean Correlation      :", round(mean(cv_corr_S1_US, na.rm = TRUE), 4), "\n")
cat("Standard Deviation    :", round(sd(cv_corr_S1_US, na.rm = TRUE), 4), "\n")

cat("\nS2 vs RF Signal Interval (0–0.01500 s):\n")
cat("Correlations by window:", round(cv_corr_S2_RF, 4), "\n")
cat("Mean Correlation      :", round(mean(cv_corr_S2_RF, na.rm = TRUE), 4), "\n")
cat("Standard Deviation    :", round(sd(cv_corr_S2_RF, na.rm = TRUE), 4), "\n")

# Dataframe for the k = 10 Local Windowed Correlation Plot
df_window_corr <- data.frame(
  Window = factor(rep(1:k, 2)),
  Correlation = c(cv_corr_S1_US, cv_corr_S2_RF),
  Component = rep(c("S1 (Ultrasonic)", "S2 (Radiofrequency)"), each = k)
)

p_window_corr <- ggplot(
  df_window_corr, 
  aes(x = Window, y = Correlation, group = Component, color = Component, shape = Component)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_shape_manual(values = c(16, 17)) +  # 16 = círculo, 17 = triángulo
  scale_y_continuous(limits = c(-1, 1)) +
  labs(
    title = "Temporal Component Stability",
    x = "Window Index (k = 1 to 10)",
    y = "Pearson Correlation Coefficient",
    color = NULL,
    shape = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

# Render individual Windowed Correlation Diagnostic Plot
print(p_window_corr)

# Exportar la figura de correlaciones por ventanas en PDF vectorial de alta resolución
ggsave(
  filename = "temporal_stability.pdf",
  plot = p_window_corr,
  width = 6,
  height = 4.5,
  units = "in",
  dpi = 600
)

# ==============================================================================
# LOCAL ROLLING CORRELATION DIAGNOSTICS (PEAK & MEAN)
# ==============================================================================
t_total <- seq(0, 0.020, by = 1/fs)

t_inicio <- 0.01625
t_final  <- 0.01650
idx_pegar <- which(t_total >= t_inicio & t_total <= t_final)

if (length(idx_pegar) == 0) {
  stop("The designated time window is unavailable.")
}

Separada1_pegada <- rep(0, length(t_total))
t_S1_comprimida <- seq(t_inicio, t_final, length.out = length(S1))

Separada1_pegada[idx_pegar] <- approx(
  x = t_S1_comprimida, y = as.numeric(S1), xout = t_total[idx_pegar], rule = 2
)$y

calc_rolling_cor <- function(x, y, window_size = 50) {
  n <- min(length(x), length(y))
  x <- as.numeric(x[seq_len(n)])
  y <- as.numeric(y[seq_len(n)])
  
  roll_cor <- rep(NA_real_, n)
  half_w <- floor(window_size / 2)
  
  for (i in seq_len(n)) {
    i1 <- max(1, i - half_w)
    i2 <- min(n, i + half_w)
    xx <- x[i1:i2]
    yy <- y[i1:i2]
    valid <- is.finite(xx) & is.finite(yy)
    xx <- xx[valid]; yy <- yy[valid]
    
    if (length(xx) >= 3 && sd(xx) > 1e-12 && sd(yy) > 1e-12) {
      roll_cor[i] <- cor(xx, yy)
    }
  }
  roll_cor
}

S1_local <- as.numeric(S1[idx_US])
US_local <- as.numeric(originalUS[idx_US])
valid_US <- is.finite(S1_local) & is.finite(US_local)
S1_local <- S1_local[valid_US]; US_local <- US_local[valid_US]

rolling_S1_US <- calc_rolling_cor(S1_local, US_local, window_size = 50)
mean_cor_S1_local <- cor(S1_local, US_local)
peak_cor_S1_local <- max(abs(rolling_S1_US), na.rm = TRUE)

S2_local <- as.numeric(S2[idx_RF])
RF_local <- as.numeric(originalRF[idx_RF])
valid_RF <- is.finite(S2_local) & is.finite(RF_local)
S2_local <- S2_local[valid_RF]; RF_local <- RF_local[valid_RF]

rolling_S2_RF <- calc_rolling_cor(S2_local, RF_local, window_size = 100)
mean_cor_S2_local <- cor(S2_local, RF_local)
peak_cor_S2_local <- max(abs(rolling_S2_RF), na.rm = TRUE)

cat("\n========================================\n")
cat("LOCAL ROLLING CORRELATIONS IN VISIBLE WINDOWS\n")
cat("========================================\n")
cat("Ultrasonic Burst Window (0.01625–0.01650 s):\n")
cat("  - Pearson correlation over full window: ", round(mean_cor_S1_local, 4), "\n")
cat("  - Maximum absolute rolling correlation : ", round(peak_cor_S1_local, 4), "\n\n")
cat("Radiofrequency Window (0–0.01500 s):\n")
cat("  - Pearson correlation over full window: ", round(mean_cor_S2_local, 4), "\n")
cat("  - Maximum absolute rolling correlation : ", round(peak_cor_S2_local, 4), "\n")

# ============================================================
# MORPHOLOGY AND AMPLITUDE CALIBRATION DIAGNOSTICS
# ============================================================
normalizar_forma <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0 || any(!is.finite(x))) stop("Vector is empty or contains non-finite values.")
  rango <- max(x) - min(x)
  if (rango == 0) return(rep(0, length(x)))
  2 * (x - min(x)) / rango - 1
}

# 1. Ultrasonic Window Calibration
if (length(idx_US) < 3) stop("Insufficient samples in US window.")
S1_US <- as.numeric(S1[idx_US])
US_ref <- as.numeric(originalUS[idx_US])
tiempo_US <- t_US[idx_US]

df_US_forma <- data.frame(
  Time = rep(tiempo_US, 2),
  Amplitude = c(normalizar_forma(S1_US), normalizar_forma(US_ref)),
  Signal = rep(c("S1 FastICA Component", "US Reference"), each = length(idx_US))
)

p_US_forma <- ggplot(df_US_forma, aes(x = Time, y = Amplitude, color = Signal)) +
  geom_line(linewidth = 0.5) +
  labs(title = "Morphology Alignment: S1 vs US Reference", x = "Time (s)", y = "Normalized Amplitude") +
  theme_minimal()

fit_US <- lm(US_ref ~ S1_US)
S1_US_ajustada <- as.numeric(predict(fit_US))

df_US_ajustada <- data.frame(
  Time = rep(tiempo_US, 2),
  Amplitude = c(S1_US_ajustada, US_ref),
  Signal = rep(c("S1 Linearly Calibrated", "US Reference"), each = length(idx_US))
)

p_US_ajustada <- ggplot(df_US_ajustada, aes(x = Time, y = Amplitude, color = Signal)) +
  geom_line(linewidth = 0.5) +
  labs(title = "Calibrated Scale: S1 vs US Reference", x = "Time (s)", y = "Physical Scale Amplitude") +
  theme_minimal()

cat("\n========================================\n")
cat("S1 vs US COMPARISON DIAGNOSTICS\n")
cat("========================================\n")
cat("Unadjusted Correlation :", round(cor(S1_US, US_ref), 4), "\n")
cat("Post-Adjustment RMSE   :", round(sqrt(mean((S1_US_ajustada - US_ref)^2)), 6), "\n")
cat("Post-Adjustment MAE    :", round(mean(abs(S1_US_ajustada - US_ref)), 6), "\n")
cat("Linear Slope           :", round(coef(fit_US)[2], 6), "\n")
cat("Linear Intercept       :", round(coef(fit_US)[1], 6), "\n")

# 2. Radiofrequency Window Calibration
if (length(idx_RF) < 3) stop("Insufficient samples in RF window.")
S2_RF <- as.numeric(S2[idx_RF])
RF_ref <- as.numeric(originalRF[idx_RF])
tiempo_RF <- t_RF[idx_RF]

df_RF_forma <- data.frame(
  Time = rep(tiempo_RF, 2),
  Amplitude = c(normalizar_forma(S2_RF), normalizar_forma(RF_ref)),
  Signal = rep(c("S2 FastICA Component", "RF Reference"), each = length(idx_RF))
)

p_RF_forma <- ggplot(df_RF_forma, aes(x = Time, y = Amplitude, color = Signal)) +
  geom_line(linewidth = 0.5) +
  labs(title = "Morphology Alignment: S2 vs RF Reference", x = "Time (s)", y = "Normalized Amplitude") +
  theme_minimal()

fit_RF <- lm(RF_ref ~ S2_RF)
S2_RF_ajustada <- as.numeric(predict(fit_RF))

df_RF_ajustada <- data.frame(
  Time = rep(tiempo_RF, 2),
  Amplitude = c(S2_RF_ajustada, RF_ref),
  Signal = rep(c("S2 Linearly Calibrated", "RF Reference"), each = length(idx_RF))
)

p_RF_ajustada <- ggplot(df_RF_ajustada, aes(x = Time, y = Amplitude, color = Signal)) +
  geom_line(linewidth = 0.5) +
  labs(title = "Calibrated Scale: S2 vs RF Reference", x = "Time (s)", y = "Physical Scale Amplitude") +
  theme_minimal()

# Generar grid de calibración y exportar a PDF
g_morphology <- arrangeGrob(p_US_forma, p_RF_forma, p_US_ajustada, p_RF_ajustada, ncol = 2)
grid.arrange(g_morphology)

ggsave(
  filename = "morphology_calibration.pdf",
  plot = g_morphology,
  width = 8.5,
  height = 6,
  units = "in",
  dpi = 600
)

# ============================================================
# DATA PREPARATION FOR COMPOSITE MULTI-PANEL FIGURE
# ============================================================
normalizar_y_completar <- function(x, n) {
  x <- as.numeric(x)
  m <- min(length(x), n)
  x <- x[seq_len(m)]
  resultado <- rep(NA_real_, n)
  valid <- is.finite(x)
  
  if (any(valid)) {
    xv <- x[valid]
    rango <- max(xv) - min(xv)
    if (rango == 0) {
      resultado[which(valid)] <- 0
    } else {
      resultado[which(valid)] <- 2 * (xv - min(xv)) / rango - 1
    }
  }
  resultado
}

df_signals <- data.frame(
  Time_sec = t_total,
  Mixture1_original = normalizar_y_completar(mixture1, length(t_total)),
  Mixture2_original = normalizar_y_completar(mixture2, length(t_total)),
  Mixture1_filtered = normalizar_y_completar(mixture1_proc, length(t_total)),
  Mixture2_filtered = normalizar_y_completar(mixture2_proc, length(t_total)),
  Separated1 = normalizar_y_completar(S1, length(t_total)),
  Separated2 = normalizar_y_completar(S2, length(t_total)),
  Reconstructed1 = normalizar_y_completar(comb1, length(t_total)),
  Reconstructed2 = normalizar_y_completar(comb2, length(t_total)),
  OriginalUS = normalizar_y_completar(originalUS, length(t_total)),
  OriginalRF = normalizar_y_completar(originalRF, length(t_total)),
  Separated1_inserted = Separada1_pegada
)

p_mixture1 <- ggplot(df_signals, aes(x = Time_sec, y = Mixture1_original)) +
  geom_line(color = "blue") +
  labs(title = "(a) Original Mixture 1", x = "Time (s)", y = "Amplitude") +
  theme_minimal() + theme(plot.title = element_text(size = 8, face = "bold"), axis.title = element_text(size = 8), axis.text = element_text(size = 6))

p_mixture2 <- ggplot(df_signals, aes(x = Time_sec, y = Mixture2_original)) +
  geom_line(color = "blue") +
  labs(title = "(a) Original Mixture 2", x = "Time (s)", y = "Amplitude") +
  theme_minimal() + theme(plot.title = element_text(size = 8, face = "bold"), axis.title = element_text(size = 8), axis.text = element_text(size = 6))

p_mixture1_filtered <- ggplot(df_signals, aes(x = Time_sec, y = Mixture1_filtered)) +
  geom_line(color = "blue") +
  labs(title = "(c) Filtered Mixture 1", x = "Time (s)", y = "Amplitude") +
  theme_minimal() + theme(plot.title = element_text(size = 8, face = "bold"), axis.title = element_text(size = 8), axis.text = element_text(size = 6))

p_mixture2_filtered <- ggplot(df_signals, aes(x = Time_sec, y = Mixture2_filtered)) +
  geom_line(color = "blue") +
  labs(title = "(c) Filtered Mixture 2", x = "Time (s)", y = "Amplitude") +
  theme_minimal() + theme(plot.title = element_text(size = 8, face = "bold"), axis.title = element_text(size = 8), axis.text = element_text(size = 6))

p_separated1 <- ggplot(data.frame(Time_sec = t_total, Separated1_inserted = Separada1_pegada), aes(x = Time_sec, y = Separated1_inserted)) +
  geom_line(color = "blue", linewidth = 0.4) +
  coord_cartesian(xlim = c(0, 0.020)) +
  labs(title = "(d) FastICA Component S1 Inserted in US Window", x = "Time (s)", y = "Amplitude") +
  theme_minimal() + theme(plot.title = element_text(size = 8, face = "bold"), axis.title = element_text(size = 8), axis.text = element_text(size = 6))

p_separated2 <- ggplot(df_signals, aes(x = Time_sec, y = Separated2)) + 
  geom_line(color = "blue") + 
  labs(title = "(d) FastICA Component S2", x = "Time (s)", y = "Amplitude") + 
  theme_minimal() + coord_cartesian(xlim = c(0.00, 0.015), ylim = c(-0.3, 0.3)) +
  theme(plot.title = element_text(size = 8, face = "bold"), axis.title = element_text(size = 8), axis.text = element_text(size = 6))

p_reconstructed1 <- ggplot(df_signals, aes(x = Time_sec, y = Reconstructed1)) + 
  geom_line(color = "blue") + labs(title = "(e) Reconstructed Signal 1", x = "Time (s)", y = "Amplitude") + 
  theme_minimal() + theme(plot.title = element_text(size = 8, face = "bold"), axis.title = element_text(size = 8), axis.text = element_text(size = 6))

p_reconstructed2 <- ggplot(df_signals, aes(x = Time_sec, y = Reconstructed2)) + 
  geom_line(color = "blue") + labs(title = "(e) Reconstructed Signal 2", x = "Time (s)", y = "Amplitude") + 
  theme_minimal() + theme(plot.title = element_text(size = 8, face = "bold"), axis.title = element_text(size = 8), axis.text = element_text(size = 6))

p_originalUS <- ggplot(df_signals, aes(x = Time_sec, y = OriginalUS)) + 
  geom_line(color = "blue") + labs(title = "(f) Experimental Reference Signal (US)", x = "Time (s)", y = "Amplitude") + 
  theme_minimal() + theme(plot.title = element_text(size = 8, face = "bold"), axis.title = element_text(size = 8), axis.text = element_text(size = 6))

p_originalRF <- ggplot(df_signals, aes(x = Time_sec, y = OriginalRF)) + 
  geom_line(color = "blue") + labs(title = "(f) Experimental Reference Signal (RF)", x = "Time (s)", y = "Amplitude") + 
  theme_minimal() + theme(plot.title = element_text(size = 8, face = "bold"), axis.title = element_text(size = 8), axis.text = element_text(size = 6))

p_fft_mixture1 <- plot_fft_normalized(mixture1, fs, title = " (b) FFT Original Mixture 1", color = "blue")  
p_fft_mixture2 <- plot_fft_normalized(mixture2, fs, title = " (b) FFT Original Mixture 2", color = "blue")  

# Construir panel compuesto de 12 sub-gráficas
g_full_pipeline <- arrangeGrob(
  p_mixture1, p_mixture2, p_fft_mixture1, p_fft_mixture2, 
  p_mixture1_filtered, p_mixture2_filtered,
  p_separated1, p_separated2, 
  p_reconstructed1, p_reconstructed2, 
  p_originalUS, p_originalRF, ncol = 2
)

grid.arrange(g_full_pipeline)

# Exportar la figura multi-panel de todas las señales en PDF vectorial
ggsave(
  filename = "full_signal_pipeline.pdf",
  plot = g_full_pipeline,
  width = 8.5,
  height = 11,
  units = "in",
  dpi = 600
)

# ============================================================
# ACTIVATION CONCORDANCE EVALUATION (PRECISION, SENSITIVITY, SPECIFICITY)
# ============================================================
evaluar_actividad <- function(componente, referencia, percentil = 0.75) {
  n <- min(length(componente), length(referencia))
  componente <- as.numeric(componente[seq_len(n)])
  referencia <- as.numeric(referencia[seq_len(n)])
  
  valid <- is.finite(componente) & is.finite(referencia)
  componente <- componente[valid]
  referencia <- referencia[valid]
  n <- length(componente)
  
  if (n < 3) stop("Insufficient valid samples for activity assessment.")
  
  umbral_componente <- quantile(abs(componente), probs = percentil)
  umbral_referencia <- quantile(abs(referencia), probs = percentil)
  
  actividad_componente <- abs(componente) >= umbral_componente
  actividad_referencia <- abs(referencia) >= umbral_referencia
  
  TP <- sum(actividad_componente & actividad_referencia, na.rm = TRUE)
  TN <- sum(!actividad_componente & !actividad_referencia, na.rm = TRUE)
  FP <- sum(actividad_componente & !actividad_referencia, na.rm = TRUE)
  FN <- sum(!actividad_componente & actividad_referencia, na.rm = TRUE)
  
  precision    <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  sensibilidad <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  especificidad <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  
  data.frame(
    Precision = precision,
    Sensitivity = sensibilidad,
    Specificity = especificidad,
    TP = TP, TN = TN, FP = FP, FN = FN,
    N = n,
    Threshold_Component = umbral_componente,
    Threshold_Reference = umbral_referencia
  )
}

# --- US Component Evaluation (0.01625 s - 0.01650 s) ---
n_US_metricas <- min(length(S1), length(originalUS))
tiempo_US_metricas <- (seq_len(n_US_metricas) - 1) / fs
idx_US_metricas <- which(tiempo_US_metricas >= 0.01625 & tiempo_US_metricas <= 0.01650)

if (length(idx_US_metricas) == 0) stop("Ultrasonic active window contains no samples.")

metricas_US_original <- evaluar_actividad(
  S1[idx_US_metricas], originalUS[idx_US_metricas], percentil = 0.75
)

# --- RF Component Evaluation (0.00000 s - 0.01500 s) ---
n_RF_metricas <- min(length(S2), length(originalRF))
tiempo_RF_metricas <- (seq_len(n_RF_metricas) - 1) / fs
idx_RF_metricas <- which(tiempo_RF_metricas >= 0 & tiempo_RF_metricas <= 0.01500)

if (length(idx_RF_metricas) == 0) stop("Radiofrequency active window contains no samples.")

metricas_RF_original <- evaluar_actividad(
  S2[idx_RF_metricas], originalRF[idx_RF_metricas], percentil = 0.75
)

# --- Print Classification Concordance Results ---
cat("\n========================================\n")
cat("ACTIVATION CONCORDANCE METRICS\n")
cat("========================================\n")

cat("\nUltrasonic Component (S1 vs US Reference):\n")
cat("  - Precision   :", round(metricas_US_original$Precision, 4), "\n")
cat("  - Sensitivity :", round(metricas_US_original$Sensitivity, 4), "\n")
cat("  - Specificity :", round(metricas_US_original$Specificity, 4), "\n")
cat("  - Confusion Matrix (TP, TN, FP, FN):", 
    metricas_US_original$TP, ",", 
    metricas_US_original$TN, ",", 
    metricas_US_original$FP, ",", 
    metricas_US_original$FN, "\n")

cat("\nRadiofrequency Component (S2 vs RF Reference):\n")
cat("  - Precision   :", round(metricas_RF_original$Precision, 4), "\n")
cat("  - Sensitivity :", round(metricas_RF_original$Sensitivity, 4), "\n")
cat("  - Specificity :", round(metricas_RF_original$Specificity, 4), "\n")
cat("  - Confusion Matrix (TP, TN, FP, FN):", 
    metricas_RF_original$TP, ",", 
    metricas_RF_original$TN, ",", 
    metricas_RF_original$FP, ",", 
    metricas_RF_original$FN, "\n")


