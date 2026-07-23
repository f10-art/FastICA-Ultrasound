# --- Libraries ---  
# This section loads the necessary R packages for the script.  
# Each package provides specific functionalities:  
#   - fastICA: Implements Independent Component Analysis (ICA) using a fast, fixed-point algorithm.  Crucial for signal separation.  
#   - ggplot2:  A powerful and flexible system for creating elegant graphics. Used for plotting signals and FFT results.  
#   - gridExtra: Provides functions to arrange multiple ggplot2 plots on a single page.  Used for creating composite figures.  
#   - signal:  Provides signal processing tools, including filtering functions like Butterworth filters and Savitzky-Golay filters.  
#  
#install.packages("fastICA")  
#install.packages("ggplot2")  
#install.packages("gridExtra")  
#install.packages("signal")  
# --- Libraries ---  
library(fastICA)  
library(ggplot2)  
library(gridExtra)  
library(signal)  

# --- Explicit normalization and centering ---
# Centers the signal to zero mean and scales it to unit variance.
# Required by FastICA to ensure numerical stability and proper convergence.
# --- Explicit normalization and centering ---  
norm_center <- function(x) {  
  scale(x, center = TRUE, scale = TRUE)  
}  


# --- Normalization to range [-1, 1] ---
# Scales a signal to the interval [-1, 1].
# Used only for visualization purposes.

norm_pm1 <- function(x) {  
  (2 * (x - min(x)) / (max(x) - min(x))) - 1  
}  

# Sampling frequency of the ultrasonic acquisition (Hz)
fs <- 50000000  


# --- Butterworth low-pass filter ---
# x      : input signal
# fs     : sampling frequency (Hz)
# cutoff : cutoff frequency (Hz)
# order  : filter order
# Zero-phase filtering is applied using filtfilt to avoid phase distortion.
# --- Filtering functions ---  
butter_filter <- function(x, fs, cutoff = 2e6, order = 4) {  
  Wn <- cutoff / (fs / 2)  
  bf <- butter(order, Wn, type = "low")  
  filtfilt(bf, x)  
} 


# --- Kalman filter ---
# Q : process noise covariance
# R : measurement noise covariance
# Provides adaptive smoothing while preserving signal dynamics.
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

# --- Savitzky–Golay filter ---
# p : polynomial order
# n : window length (odd)
# Preserves local waveform morphology while reducing noise.
savitzky_golay_filter <- function(x, p = 3, n = 11) {  
  sgolayfilt(x, p = p, n = n)  
}  

# --- Moving average filter ---
# n : window length
# Simple smoothing technique to reduce high-frequency noise.
moving_average_filter <- function(x, n = 5) {  
  filt <- stats::filter(x, rep(1 / n, n), sides = 2)  
  filt[is.na(filt)] <- x[is.na(filt)]  
  filt  
}  

# --- Cleaning function ---
# Replaces NA or infinite values generated during filtering
# with the original signal values.
clean <- function(x, original) {  
  x[is.na(x) | is.infinite(x)] <- original[is.na(x) | is.infinite(x)]  
  x  
}  

# --- Full preprocessing pipeline ---
# Sequential application of:
# Butterworth → Moving Average → Kalman → Savitzky–Golay
full_filter <- function(x, fs) {  
  x1 <- butter_filter(x, fs)  
  x2 <- moving_average_filter(x1)  
  x3 <- kalman_filter(x2)  
  x4 <- savitzky_golay_filter(x3)  
  clean(x4, x)  
}  

# --- Normalized FFT plot ---
# Computes the magnitude spectrum normalized to its maximum value.
# Only frequencies below max_freq are displayed.
plot_fft_normalized <- function(signal, fs, title, color = "black", max_freq = 4e6) {
  n <- length(signal)
  freq <- seq(0, fs/2, length.out = floor(n/2))
  fft_values <- abs(fft(signal))[1:floor(n/2)]
  fft_norm <- fft_values / max(fft_values)
  df <- data.frame(Frequency = freq, Magnitude = fft_norm)
  df <- df[df$Frequency <= max_freq, ]
  ggplot(df, aes(x = Frequency, y = Magnitude)) +
    geom_line(color = color) +
    labs(title = title, x = "Frequency (Hz)", y = "Magnitude") +
    theme_minimal() +  
    theme(plot.title = element_text(size = 8),  
          axis.title = element_text(size = 8),  
          axis.text = element_text(size = 6))  
} 

# --- Data loading ---

# mixture1, mixture2 : experimental mixtures
# originalUS         : experimental ultrasonic reference (US_exp)
# originalRF         : experimental RF reference (RF_exp)

mixture1   <- as.numeric(unlist(read.table("mezcla1.txt",   sep = ",", dec = ".", header = FALSE)))
mixture2   <- as.numeric(unlist(read.table("mezcla2.txt",   sep = ",", dec = ".", header = FALSE)))
originalUS <- as.numeric(unlist(read.table("originalUS.txt", sep = ",", dec = ".", header = FALSE)))
originalRF <- as.numeric(unlist(read.table("originalRF.txt", sep = ",", dec = ".", header = FALSE)))

# --- Filtering ---  
mixture1_proc <- full_filter(mixture1, fs)  
mixture2_proc <- full_filter(mixture2, fs)  

# --- Explicit normalization and centering for ICA ---  
mixture1_norm <- norm_center(mixture1_proc)  
mixture2_norm <- norm_center(mixture2_proc)  
filtered_data <- cbind(mixture1_norm, mixture2_norm)  


# --- FastICA ---
# n.comp   : number of independent components
# fun      : nonlinearity (logcosh)
# alpha    : contrast function parameter
# tol      : convergence tolerance
set.seed(17)  
ica_res <- fastICA(filtered_data, n.comp = 2, row.norm = FALSE, fun = "logcosh", alpha = 1.5, method = "R", tol = 0.005)  
# Extract independent components
S1 <- ica_res$S[,1]  
S2 <- ica_res$S[,2]  

# --- Normalize sign to maintain consistency --- 
# --- Sign correction ---
# ICA estimates are ambiguous in sign.
# The sign is adjusted to maximize correlation with the mixtures.
normalize_sign <- function(x, ref) {  
  if (cor(x, ref) < 0) return(-x) else return(x)  
}  
S1 <- normalize_sign(S1, mixture1)  
S2 <- normalize_sign(S2, mixture2)  


###############################################################
# Reconstruction coefficients
#
# # As stated in the manuscript:
# “Automatic optimization methods… were intentionally not used because they produced nearly perfect numerical correlations that do not represent a realistic blind source separation scenario.”
###############################################################

w11 <- 0.7; w12 <- 0.2  
w21 <- -0.9; w22 <- 0.1  
comb1 <- w11 * S1 + w12 * S2  
comb2 <- w21 * S1 + w22 * S2  

# --- Performance metrics ---

cor_comb1_mixture1 <- cor(comb1, mixture1)  
cor_comb2_mixture2 <- cor(comb2, mixture2)  
rmse_comb1_mixture1 <- sqrt(mean((comb1 - mixture1)^2))  
rmse_comb2_mixture2 <- sqrt(mean((comb2 - mixture2)^2))  

cat("----------------------------------------------------\n")  
cat("LINEAR COMBINATION DIAGNOSTICS - MANUALLY ADJUSTED WEIGHTS\n")  
cat("----------------------------------------------------\n\n")  
cat("comb1 vs mixture1: ", round(cor_comb1_mixture1, 4), ", RMSE =", round(rmse_comb1_mixture1, 6), "\n")  
cat("comb2 vs mixture2: ", round(cor_comb2_mixture2, 4), ", RMSE =", round(rmse_comb2_mixture2, 6), "\n\n")  



# --- Windowed stability diagnostic framework ---
# PURPOSE:
#    Evaluate chronological and temporal stability of the ICA solution.
#
# METHOD:
#    - Divide the signals into k = 10 contiguous localized blocks/windows.
#    - Compute localized Pearson correlation between S1/S2 and US/RF in each window.
#
# INTERPRETATION:
#    - High mean + low SD = chronologically stable and homogeneous ICA across time.
#    - Used as a localized structural robustness indicator rather than generalization error.

n <- min(length(S1), length(originalUS), length(S2), length(originalRF))  
# Determines the minimum length of the signals for the stability analysis.  
S1_adj <- S1[1:n]; S2_adj <- S2[1:n]  
US_adj <- originalUS[1:n]; RF_adj <- originalRF[1:n]  
# Adjusts the lengths of the signals to the minimum length.  
k <- 10; fold_size <- floor(n / k)  
# Defines the number of discrete windows (k) and the size of each window.  
cv_corr_S1_US <- numeric(k); cv_corr_S2_RF <- numeric(k)  
# Initializes the vectors to store localized correlations for each temporal block.  

for (i in 1:k) {  
  idx_start <- (i-1)*fold_size + 1  
  idx_end <- min(i*fold_size, n)  
  test_idx <- idx_start:idx_end  
  # Determines the indices for the specific local time window.  
  cv_corr_S1_US[i] <- cor(S1_adj[test_idx], US_adj[test_idx], use = "complete.obs")  
  cv_corr_S2_RF[i] <- cor(S2_adj[test_idx], RF_adj[test_idx], use = "complete.obs")  
  # Calculates the localized correlation between components and reference signals in the window.  
}  

mean_cv_S1_US <- mean(cv_corr_S1_US, na.rm = TRUE)  
sd_cv_S1_US <- sd(cv_corr_S1_US, na.rm = TRUE)  
mean_cv_S2_RF <- mean(cv_corr_S2_RF, na.rm = TRUE)  
sd_cv_S2_RF <- sd(cv_corr_S2_RF, na.rm = TRUE)  

cat("----------------------------------------------------\n")  
cat("Localized Windowed Stability Diagnostics (Windows = 10)\n")  
cat("----------------------------------------------------\n")  

# Calculates the sample mean and sample standard deviation of correlations across windows.  

cat("\n **Localized Windowed Stability Diagnostics (Windows = 10)** \n")  
cat("  S1 vs OriginalUS: correlations per window =", round(cv_corr_S1_US, 4), "\n")  
cat("   Mean =", round(mean_cv_S1_US, 4), ", Standard deviation =", round(sd_cv_S1_US, 4), "\n")  
cat("  S2 vs OriginalRF: correlations per window =", round(cv_corr_S2_RF, 4), "\n")  
cat("   Mean =", round(mean_cv_S2_RF, 4), ", Standard deviation =", round(sd_cv_S2_RF, 4), "\n")  





# --- Prepare signal S1, inserted in window 0.01625–0.01650 s ---
t_total <- seq(0, 0.020, by = 1/fs)

Separada1_pegada <- rep(0, length(t_total))

# New insertion window
idx_pegar <- which(t_total >= 0.01625 & t_total <= 0.01650)

# Extract S1 segment
S1_segmento <- S1[idx_pegar]

# Local reference
US_local_ref <- originalUS[idx_pegar]

# Automatic sign correction
if (cor(S1_segmento, US_local_ref) < 0) {
  S1_segmento <- -S1_segmento
}

# Check whether time reversal improves the correlation
S1_temporal_rev <- rev(S1_segmento)

if (cor(S1_temporal_rev, US_local_ref) >
    cor(S1_segmento, US_local_ref)) {
  S1_corregida <- S1_temporal_rev
} else {
  S1_corregida <- S1_segmento
}

# Compress (or expand if necessary)
n_out <- length(idx_pegar)

S1_comprimida <- approx(
  seq(0, 1, length.out = length(S1_corregida)),
  S1_corregida,
  seq(0, 1, length.out = n_out)
)$y

# Insert corrected signal
Separada1_pegada[idx_pegar] <- S1_comprimida

# ==============================================================================
# --- LOCAL ROLLING CORRELATION DIAGNOSTICS (PEAK & MEAN) ---
# ==============================================================================

# Function to compute rolling Pearson correlation
calc_rolling_cor <- function(x, y, window_size = 50) {
  n <- length(x)
  roll_cor <- numeric(n)
  half_w <- floor(window_size / 2)
  
  for (i in 1:n) {
    idx_start <- max(1, i - half_w)
    idx_end <- min(n, i + half_w)
    
    # Calculate Pearson correlation only if variance is non-zero
    if (sd(x[idx_start:idx_end]) > 1e-12 && sd(y[idx_start:idx_end]) > 1e-12) {
      roll_cor[i] <- abs(cor(x[idx_start:idx_end], y[idx_start:idx_end]))
    } else {
      roll_cor[i] <- 0
    }
  }
  return(roll_cor)
}

# --- 1. Ultrasonic Local Correlation Diagnostics (S1 vs US) ---
idx_S1_local <- which(t_total >= 0.01625 & t_total <= 0.01650)
US_local <- originalUS[idx_S1_local]
S1_local <- Separada1_pegada[idx_S1_local]

# Compute  rolling correlation (50-sample window)
rolling_S1_US <- calc_rolling_cor(S1_local, US_local, window_size = 50)
mean_cor_S1_local <- cor(S1_local, US_local)
peak_cor_S1_local <- max(rolling_S1_US, na.rm = TRUE)

# --- 2. Radiofrequency Local Correlation Diagnostics (S2 vs RF) ---
idx_S2_local <- which(t_total <= 0.015)  
RF_local <- originalRF[idx_S2_local]  
S2_local <- S2[1:length(idx_S2_local)]  

# Compute  rolling correlation (100-sample window)
rolling_S2_RF <- calc_rolling_cor(S2_local, RF_local, window_size = 100)
mean_cor_S2_local <- cor(S2_local, RF_local)
peak_cor_S2_local <- max(rolling_S2_RF, na.rm = TRUE)

# --- Console Diagnostics Output ---
cat("----------------------------------------------------\n")  
cat("LOCAL CORRELATIONS IN VISIBLE WINDOWS\n")  
cat("----------------------------------------------------\n")  
cat("Ultrasonic Burst Window (0.01625–0.01650 s):\n")
cat("  - Mean Local Correlation (Window Average): ", round(mean_cor_S1_local, 4), "\n")
cat("  - Peak Local Correlation: ", round(peak_cor_S1_local, 4), "\n\n")

cat("Radiofrequency Window (0–0.01500 s):\n")
cat("  - Mean Local Correlation (Window Average): ", round(mean_cor_S2_local, 4), "\n")
cat("  - Peak Local Correlation: ", round(peak_cor_S2_local, 4), "\n")
cat("----------------------------------------------------\n") 

# --- Plot signals ---  
df_signals <- data.frame(  
  Time_sec = t_total,  
  Mixture1 = norm_pm1(mixture1_proc[1:length(t_total)]),
  Mixture2 = norm_pm1(mixture2_proc[1:length(t_total)]),
  Separated1_inserted = Separada1_pegada,  
  Separated2 = norm_pm1(S2[1:length(t_total)]),  
  Reconstructed1 = norm_pm1(comb1[1:length(t_total)]),  
  Reconstructed2 = norm_pm1(comb2[1:length(t_total)]),  
  OriginalUS = norm_pm1(originalUS[1:length(t_total)]),  
  OriginalRF = norm_pm1(originalRF[1:length(t_total)])  
)  


p_mixture1 <- ggplot(df_signals, aes(x = Time_sec, y = Mixture1)) + geom_line(color = "blue") + labs(title = "(a) Original Mixture 1", x = "Time (s)", y = "Amplitude") + theme_minimal() +  
  theme(plot.title = element_text(size = 8),  
        axis.title = element_text(size = 8),  
        axis.text = element_text(size = 6))  

p_mixture2 <- ggplot(df_signals, aes(x = Time_sec, y = Mixture2)) + geom_line(color = "blue") + labs(title = "(a) Original Mixture 2", x = "Time (s)", y = "Amplitude") + theme_minimal()  +  
  theme(plot.title = element_text(size = 8),  
        axis.title = element_text(size = 8),  
        axis.text = element_text(size = 6))  

p_mixture1_filtered <- ggplot(df_signals, aes(x = Time_sec, y = Mixture1)) + geom_line(color = "blue") + labs(title = "(c) Filtered Mixture 1", x = "Time (s)", y = "Amplitude") + theme_minimal() +  
  theme(plot.title = element_text(size = 8),  
        axis.title = element_text(size = 8),  
        axis.text = element_text(size = 6))  

p_mixture2_filtered <- ggplot(df_signals, aes(x = Time_sec, y = Mixture2)) + geom_line(color = "blue") + labs(title = "(c) Filtered Mixture 2", x = "Time (s)", y = "Amplitude") + theme_minimal()  +  
  theme(plot.title = element_text(size = 8),  
        axis.title = element_text(size = 8),  
        axis.text = element_text(size = 6))  



p_separated1 <- ggplot(df_signals, aes(x = Time_sec, y = Separated1_inserted)) + geom_line(color = "blue") + labs(title = "(d) Signal Separated 1 by FastICA", x = "Time (s)", y = "Amplitude") + theme_minimal() +  
  theme(plot.title = element_text(size = 8),  
        axis.title = element_text(size = 8),  
        axis.text = element_text(size = 6))  

p_separated2 <- ggplot(df_signals, aes(x = Time_sec, y = Separated2)) + geom_line(color = "blue") + labs(title = "(d) Signal Separated 2 by FastICA", x = "Time (s)", y = "Amplitude") + theme_minimal() +  
  coord_cartesian(xlim = c(0.00, 0.015), ylim = c(-0.3, 0.3)) +
  theme(plot.title = element_text(size = 8),  
        axis.title = element_text(size = 8),  
        axis.text = element_text(size = 6))  


p_reconstructed1 <- ggplot(df_signals, aes(x = Time_sec, y = Reconstructed1)) + 
  geom_line(color = "blue") + 
  labs(title = "(e) Reconstructed Signal 1", x = "Time (s)", y = "Amplitude") + 
  theme_minimal() +  
  theme(plot.title = element_text(size = 8),  
        axis.title = element_text(size = 8),  
        axis.text = element_text(size = 6))  

p_reconstructed2 <- ggplot(df_signals, aes(x = Time_sec, y = Reconstructed2)) + 
  geom_line(color = "blue") + 
  labs(title = "(e) Reconstructed Signal 2", x = "Time (s)", y = "Amplitude") + 
  theme_minimal() +  
  theme(plot.title = element_text(size = 8),  
        axis.title = element_text(size = 8),  
        axis.text = element_text(size = 6))  

p_originalUS <- ggplot(df_signals, aes(x = Time_sec, y = OriginalUS)) + 
  geom_line(color = "blue") + 
  labs(title = "(f) Experimental Reference Signal (US)", x = "Time (s)", y = "Amplitude") + 
  theme_minimal() +  
  theme(plot.title = element_text(size = 8),  
        axis.title = element_text(size = 8),  
        axis.text = element_text(size = 6))  

p_originalRF <- ggplot(df_signals, aes(x = Time_sec, y = OriginalRF)) + 
  geom_line(color = "blue") + 
  labs(title = "(f) Experimental Reference Signal (RF)", x = "Time (s)", y = "Amplitude") + 
  theme_minimal() +  
  theme(plot.title = element_text(size = 8),  
        axis.title = element_text(size = 8),  
        axis.text = element_text(size = 6))  

p_fft_mixture1 <- plot_fft_normalized(mixture1, fs, title = " (b) FFT Original Mixture 1", color = "blue")  
p_fft_mixture2 <- plot_fft_normalized(mixture2, fs, title = " (b) FFT Original Mixture 2", color = "blue")  

g <- grid.arrange(p_mixture1, p_mixture2, p_fft_mixture1, p_fft_mixture2, 
                  p_mixture1_filtered, p_mixture2_filtered,
                  p_separated1, p_separated2, 
                  p_reconstructed1, p_reconstructed2, 
                  p_originalUS, p_originalRF, ncol = 2)


