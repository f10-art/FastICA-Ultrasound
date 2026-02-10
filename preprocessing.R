############################################################
# preprocessing.R
# Preprocessing utilities for FastICA ultrasound project
# Includes: normalization, cleaning, and full filtering pipeline
############################################################

# --- Required library ---
library(signal)

############################################################
# Normalization functions
############################################################

# Standardization: zero mean, unit variance (used for ICA)
norm_center <- function(x) {
  scale(x, center = TRUE, scale = TRUE)
}

# Normalization to [-1, 1] (used ONLY for visualization)
norm_pm1 <- function(x) {
  (2 * (x - min(x)) / (max(x) - min(x))) - 1
}

############################################################
# Filtering functions
############################################################

# Butterworth low-pass filter
butter_filter <- function(x, fs, cutoff = 2e6, order = 4) {
  Wn <- cutoff / (fs / 2)
  bf <- butter(order, Wn, type = "low")
  filtfilt(bf, x)
}

# Moving average filter
moving_average_filter <- function(x, n = 5) {
  filt <- stats::filter(x, rep(1 / n, n), sides = 2)
  filt[is.na(filt)] <- x[is.na(filt)]
  filt
}

# 1D Kalman filter
kalman_filter <- function(x, Q = 0.01, R = 0.01) {
  n <- length(x)
  x_hat <- numeric(n)
  P <- numeric(n)
  x_hat[1] <- x[1]; P[1] <- 1
  
  for (k in 2:n) {
    x_pred <- x_hat[k - 1]
    P_pred <- P[k - 1] + Q
    K <- P_pred / (P_pred + R)
    x_hat[k] <- x_pred + K * (x[k] - x_pred)
    P[k] <- (1 - K) * P_pred
  }
  
  x_hat
}

# Savitzky–Golay filter
savitzky_golay_filter <- function(x, p = 3, n = 11) {
  sgolayfilt(x, p = p, n = n)
}

############################################################
# Cleaning function
############################################################

# Replace NA/Inf values with original signal values
clean_signal <- function(x, original) {
  x[is.na(x) | is.infinite(x)] <- original[is.na(x) | is.infinite(x)]
  x
}

############################################################
# Full preprocessing pipeline
############################################################

# Butterworth → Moving Average → Kalman → Savitzky–Golay → Cleaning
full_filter <- function(x, fs) {
  x1 <- butter_filter(x, fs)
  x2 <- moving_average_filter(x1)
  x3 <- kalman_filter(x2)
  x4 <- savitzky_golay_filter(x3)
  clean_signal(x4, x)
}
