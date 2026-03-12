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
plot_fft_normalized <- function(signal, fs, title, color = "black", max_freq = 5e6) {
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
mixture1 <- as.numeric(unlist(read.table("C:/Users/hp/Desktop/FastICA-Ultrasound/mezcla1.txt", sep = ",", dec = ".", header = FALSE)))  
mixture2 <- as.numeric(unlist(read.table("C:/Users/hp/Desktop/FastICA-Ultrasound/mezcla2.txt", sep = ",", dec = ".", header = FALSE)))  
originalUS <- as.numeric(unlist(read.table("C:/Users/hp/Desktop/FastICA-Ultrasound/originalUS.txt", sep=",", dec=".", header=FALSE)))  
originalRF <- as.numeric(unlist(read.table("C:/Users/hp/Desktop/FastICA-Ultrasound/originalRF.txt", sep=",", dec=".", header=FALSE)))  

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


# --- Manual linear reconstruction ---
# Weights are adjusted manually to approximate the original mixtures.
w11 <- 0.7; w12 <- 0.2  
w21 <- -0.9; w22 <- 0.1  
comb1 <- w11 * S1 + w12 * S2  
comb2 <- w21 * S1 + w22 * S2  

# --- Performance metrics ---
# Global correlation and RMSE are computed
# between reconstructed signals and mixtures.
cor_comb1_mixture1 <- cor(comb1, mixture1)  
cor_comb2_mixture2 <- cor(comb2, mixture2)  
rmse_comb1_mixture1 <- sqrt(mean((comb1 - mixture1)^2))  
rmse_comb2_mixture2 <- sqrt(mean((comb2 - mixture2)^2))  

cat("----------------------------------------------------\n")  
cat("RECONSTRUCTION RESULTS - MANUALLY ADJUSTED WEIGHTS\n")  
cat("----------------------------------------------------\n\n")  
cat("comb1 vs mixture1: ", round(cor_comb1_mixture1, 4), ", RMSE =", round(rmse_comb1_mixture1, 6), "\n")  
cat("comb2 vs mixture2: ", round(cor_comb2_mixture2, 4), ", RMSE =", round(rmse_comb2_mixture2, 6), "\n\n")  

# --- Validation with original signals ---  
cor_S1_US <- cor(S1, originalUS)  
cor_S1_RF <- cor(S1, originalRF)  
cor_S2_US <- cor(S2, originalUS)  
cor_S2_RF <- cor(S2, originalRF)  

cat("Validation with original sources:\n")  
cat("S1 vs OriginalUS: ", round(cor_S1_US, 4), "\n")  
cat("S1 vs OriginalRF: ", round(cor_S1_RF, 4), "\n")  
cat("S2 vs OriginalUS: ", round(cor_S2_US, 4), "\n")  
cat("S2 vs OriginalRF: ", round(cor_S2_RF, 4), "\n")  

# --- k-fold cross-validation ---
# Evaluates robustness of the separation using k = 10 folds.
n <- min(length(S1), length(originalUS), length(S2), length(originalRF))  
# Determines the minimum length of the signals for cross-validation.  
S1_adj <- S1[1:n]; S2_adj <- S2[1:n]  
US_adj <- originalUS[1:n]; RF_adj <- originalRF[1:n]  
# Adjusts the lengths of the signals to the minimum length.  
k <- 10; fold_size <- floor(n / k)  
# Defines the number of folds (k) and the size of each fold.  
cv_corr_S1_US <- numeric(k); cv_corr_S2_RF <- numeric(k)  
# Initializes the vectors to store correlations for each fold.  

for (i in 1:k) {  
  idx_start <- (i-1)*fold_size + 1  
  idx_end <- min(i*fold_size, n)  
  test_idx <- idx_start:idx_end  
  # Determines the indices for the test set in each fold.  
  cv_corr_S1_US[i] <- cor(S1_adj[test_idx], US_adj[test_idx], use = "complete.obs")  
  cv_corr_S2_RF[i] <- cor(S2_adj[test_idx], RF_adj[test_idx], use = "complete.obs")  
  # Calculates the correlation between independent components and original signals in the test set.  
}  

mean_cv_S1_US <- mean(cv_corr_S1_US, na.rm = TRUE)  
sd_cv_S1_US <- sd(cv_corr_S1_US, na.rm = TRUE)  
mean_cv_S2_RF <- mean(cv_corr_S2_RF, na.rm = TRUE)  
sd_cv_S2_RF <- sd(cv_corr_S2_RF, na.rm = TRUE)  
# Calculates the mean and standard deviation of correlations for each fold.  

cat("\n **Cross-Validation (k-fold = 10)** \n")  
cat("  S1 vs OriginalUS: correlations per fold =", round(cv_corr_S1_US, 4), "\n")  
cat("   Mean =", round(mean_cv_S1_US, 4), ", Standard deviation =", round(sd_cv_S1_US, 4), "\n")  
cat("  S2 vs OriginalRF: correlations per fold =", round(cv_corr_S2_RF, 4), "\n")  
cat("   Mean =", round(mean_cv_S2_RF, 4), ", Standard deviation =", round(sd_cv_S2_RF, 4), "\n")  
# Prints the results of cross-validation.  

# --- k-fold cross-validation with ICA re-estimation ---
k <- 10
n <- min(length(mixture1), length(mixture2), length(originalUS), length(originalRF))
fold_size <- floor(n / k)

cv_corr_S1_US <- numeric(k)
cv_corr_S2_RF <- numeric(k)

set.seed(17)

for (i in 1:k) {
  
  # Define fold indices
  idx_start <- (i - 1) * fold_size + 1
  idx_end <- ifelse(i == k, n, i * fold_size)
  idx <- idx_start:idx_end
  
  # Extract fold mixtures
  mix1_fold <- mixture1[idx]
  mix2_fold <- mixture2[idx]
  
  # Filtering
  mix1_f <- full_filter(mix1_fold, fs)
  mix2_f <- full_filter(mix2_fold, fs)
  
  # Normalization
  mix1_n <- norm_center(mix1_f)
  mix2_n <- norm_center(mix2_f)
  fold_data <- cbind(mix1_n, mix2_n)
  
  # ICA re-estimation
  ica_fold <- fastICA(fold_data, n.comp = 2, row.norm = FALSE,
                      fun = "logcosh", alpha = 1.5, method = "R", tol = 0.005)
  
  S1_fold <- ica_fold$S[,1]
  S2_fold <- ica_fold$S[,2]
  
  # Sign correction
  S1_fold <- normalize_sign(S1_fold, mix1_fold)
  S2_fold <- normalize_sign(S2_fold, mix2_fold)
  
  # Extract reference signals for this fold
  US_fold <- originalUS[idx]
  RF_fold <- originalRF[idx]
  
  # Compute correlations
  cv_corr_S1_US[i] <- cor(S1_fold, US_fold, use = "complete.obs")
  cv_corr_S2_RF[i] <- cor(S2_fold, RF_fold, use = "complete.obs")
}

# Summary statistics
mean_cv_S1_US <- mean(cv_corr_S1_US, na.rm = TRUE)
sd_cv_S1_US <- sd(cv_corr_S1_US, na.rm = TRUE)
mean_cv_S2_RF <- mean(cv_corr_S2_RF, na.rm = TRUE)
sd_cv_S2_RF <- sd(cv_corr_S2_RF, na.rm = TRUE)

cat("\n--- ICA k-fold cross-validation (re-estimated ICA) ---\n")
cat("S1 vs US: correlations =", round(cv_corr_S1_US, 4), "\n")
cat("Mean =", round(mean_cv_S1_US, 4), " SD =", round(sd_cv_S1_US, 4), "\n\n")
cat("S2 vs RF: correlations =", round(cv_corr_S2_RF, 4), "\n")
cat("Mean =", round(mean_cv_S2_RF, 4), " SD =", round(sd_cv_S2_RF, 4), "\n")


# --- Prepare signal S1, inserted in window 0.016–0.017 s, rotated + compressed ---  
t_total <- seq(0, 0.020, by = 1/fs)  
# Creates a time vector from 0 to 0.020 seconds with resolution of 1/fs.  
Separada1_pegada <- rep(0, length(t_total))  
# Creates a zero vector with the same length as the time vector.  
idx_pegar <- which(t_total >= 0.016 & t_total <= 0.017)  
# Finds the indices of the time vector corresponding to the interval 0.016–0.017 seconds.  

# Extract the section of S1 to insert  
S1_segmento <- S1[idx_pegar]  
# Extracts the segment of signal S1 corresponding to the identified indices.  

# Rotate the section of S1  
S1_rotada <- rev(S1_segmento)  
# Reverses the order of the extracted segment of S1.  

# Compress the rotated section  
n_out <- length(idx_pegar)  # Same length as the original segment to insert  
S1_comprimida <- approx(seq(0, 1, length.out = length(S1_rotada)), S1_rotada, seq(0, 1, length.out = n_out))$y  
# Compresses or expands the rotated segment to match the length of the original segment.  

Separada1_pegada[idx_pegar] <- S1_comprimida  # Insert the compressed and rotated section  
# Replaces zeros in 'Separada1_pegada' with the compressed and rotated segment.  


# --- Local correlations ---  
idx_S1_local <- which(t_total >= 0.016 & t_total <= 0.017)  
US_local <- originalUS[idx_S1_local]  
S1_local <- Separada1_pegada[idx_S1_local]  
# Defines indices and extracts local sections of S1 and the original US signal.  

idx_S2_local <- which(t_total <= 0.015)  
RF_local <- originalRF[idx_S2_local]  
S2_local <- S2[1:length(idx_S2_local)]  
# Defines indices and extracts local sections of S2 and the original RF signal.  

cat("----------------------------------------------------\n")  
cat("LOCAL CORRELATIONS IN VISIBLE WINDOWS\n")  
cat("----------------------------------------------------\n")  
cat("Local correlation S1 inserted, rotated and compressed vs OriginalUS (0.016–0.017 s): ",  
    round(cor(S1_local, US_local), 4), "\n")  
cat("Local correlation S2 vs OriginalRF (0–0.015 s): ",  
    round(cor(S2_local, RF_local), 4), "\n")  
# Prints local correlations between signals.  

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




