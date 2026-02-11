# --- Libraries ---  
library(fastICA)  
library(ggplot2)  
library(gridExtra)  
library(signal)  

# --- Explicit normalization and centering ---  
norm_center <- function(x) {  
  scale(x, center = TRUE, scale = TRUE)  
}  

# --- Normalization to range [-1, 1] ---  
norm_pm1 <- function(x) {  
  (2 * (x - min(x)) / (max(x) - min(x))) - 1  
}  

fs <- 50e6
signal_length <- 3500
t_total <- seq_len(signal_length) / fs

# --- Filtering functions ---  
butter_filter <- function(x, fs, cutoff = 2e6, order = 4) {  
  Wn <- cutoff / (fs / 2)  
  bf <- butter(order, Wn, type = "low")  
  filtfilt(bf, x)  
}  

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

savitzky_golay_filter <- function(x, p = 3, n = 11) {  
  sgolayfilt(x, p = p, n = n)  
}  

moving_average_filter <- function(x, n = 5) {  
  filt <- stats::filter(x, rep(1 / n, n), sides = 2)  
  filt[is.na(filt)] <- x[is.na(filt)]  
  filt  
}  

clean_signal <- function(x, original) {  
  x[is.na(x) | is.infinite(x)] <- original[is.na(x) | is.infinite(x)]  
  x  
}  

full_filter <- function(x, fs) {  
  x1 <- butter_filter(x, fs)  
  x2 <- moving_average_filter(x1)  
  x3 <- kalman_filter(x2)  
  x4 <- savitzky_golay_filter(x3)  
  clean_signal(x4, x)  
}  

plot_fft_normalized <- function(signal, fs, title, color = "blue", max_freq = 4e6) {
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

# --- Simulated original signals ---  
originalUS <- 0.5*sin(2*pi*1e6*t_total) + 0.3*sin(2*pi*2e6*t_total)
originalRF <- 0.7*sin(2*pi*0.5e6*t_total) + 0.4*sin(2*pi*1.5e6*t_total)

# --- Generate random mixing matrix A ---  
set.seed(2)
A <- matrix(runif(4, min = 0.1, max = 1.0), nrow = 2, ncol = 2)

# --- Automatically generated mixtures using A ---  
S <- cbind(originalUS, originalRF)
mixtures <- S %*% t(A)  # x(t) = A * s(t)

# Add some noise
mixtures <- mixtures + matrix(rnorm(2*signal_length, 0, 0.05), ncol=2)

mixture1 <- mixtures[,1]
mixture2 <- mixtures[,2]

# --- Filtering ---  
mixture1_proc <- full_filter(mixture1, fs)  
mixture2_proc <- full_filter(mixture2, fs)  

# --- Explicitly normalize and center for ICA ---  
mixture1_norm <- norm_center(mixture1_proc)  
mixture2_norm <- norm_center(mixture2_proc)  
filtered_data <- cbind(mixture1_norm, mixture2_norm)  

# --- FastICA ---  
set.seed(17) 
ica_res <- fastICA(filtered_data, n.comp = 2, row.norm = FALSE, fun = "logcosh", alpha = 1.5, method = "R", tol = 0.005)  
S1 <- ica_res$S[,1]  
S2 <- ica_res$S[,2]  

# --- Sign normalization for consistency ---  
normalize_sign <- function(x, ref) {  
  if (cor(x, ref) < 0) return(-x) else return(x)  
}  
S1 <- normalize_sign(S1, mixture1)  
S2 <- normalize_sign(S2, mixture2)  

# --- Manual reconstruction with known weights of A ---  
comb1 <- A[1,1]*S1 + A[1,2]*S2
comb2 <- A[2,1]*S1 + A[2,2]*S2

# --- Correlations ---  
cat("Mixing matrix A:\n")
print(A)

cor_S1_mixture1 <- cor(S1, mixture1)
cor_S1_mixture2 <- cor(S1, mixture2)
cor_S2_mixture1 <- cor(S2, mixture1)
cor_S2_mixture2 <- cor(S2, mixture2)

cat("\n----------------------------------------\n")
cat("CORRELATIONS BETWEEN SEPARATED SIGNALS AND MIXTURES\n")
cat("----------------------------------------\n")
cat("S1 vs mixture1: ", round(cor_S1_mixture1, 4), "\n")
cat("S1 vs mixture2: ", round(cor_S1_mixture2, 4), "\n")
cat("S2 vs mixture1: ", round(cor_S2_mixture1, 4), "\n")
cat("S2 vs mixture2: ", round(cor_S2_mixture2, 4), "\n")

# --- Global correlations and RMSE ---  
cor_comb1_mixture1 <- cor(comb1, mixture1)  
cor_comb2_mixture2 <- cor(comb2, mixture2)  
rmse_comb1_mixture1 <- sqrt(mean((comb1 - mixture1)^2))  
rmse_comb2_mixture2 <- sqrt(mean((comb2 - mixture2)^2))  

cat("----------------------------------------------------\n")  
cat("RECONSTRUCTION RESULTS - ADJUSTED MANUAL WEIGHTS\n")  
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
n <- min(length(S1), length(originalUS), length(S2), length(originalRF))  
S1_adj <- S1[1:n]; S2_adj <- S2[1:n]  
US_adj <- originalUS[1:n]; RF_adj <- originalRF[1:n]  
k <- 10; fold_size <- floor(n / k)  
cv_corr_S1_US <- numeric(k); cv_corr_S2_RF <- numeric(k)  

for (i in 1:k) {  
  idx_start <- (i-1)*fold_size + 1  
  idx_end <- min(i*fold_size, n)  
  test_idx <- idx_start:idx_end  
  cv_corr_S1_US[i] <- cor(S1_adj[test_idx], US_adj[test_idx], use = "complete.obs")  
  cv_corr_S2_RF[i] <- cor(S2_adj[test_idx], RF_adj[test_idx], use = "complete.obs")  
}  

mean_cv_S1_US <- mean(cv_corr_S1_US, na.rm = TRUE)  
sd_cv_S1_US <- sd(cv_corr_S1_US, na.rm = TRUE)  
mean_cv_S2_RF <- mean(cv_corr_S2_RF, na.rm = TRUE)  
sd_cv_S2_RF <- sd(cv_corr_S2_RF, na.rm = TRUE)  

cat("\n **Cross-Validation (k-fold = 10)** \n")  
cat("  S1 vs OriginalUS: correlations per fold =", round(cv_corr_S1_US, 4), "\n")  
cat("   Mean =", round(mean_cv_S1_US, 4), ", Standard deviation =", round(sd_cv_S1_US, 4), "\n")  
cat("  S2 vs OriginalRF: correlations per fold =", round(cv_corr_S2_RF, 4), "\n")  
cat("   Mean =", round(mean_cv_S2_RF, 4), ", Standard deviation =", round(sd_cv_S2_RF, 4), "\n")  


# --- Data frame for plots ---
df_signals <- data.frame(
  Time_sec = t_total,
  Mixture1 = norm_pm1(mixture1_proc),
  Mixture2 = norm_pm1(mixture2_proc),
  Separated1 = S1,
  Separated2 = S2,
  Recovered1 = norm_pm1(comb1),
  Recovered2 = norm_pm1(comb2),
  OriginalUS = norm_pm1(originalUS),
  OriginalRF = norm_pm1(originalRF)
)

# --- Plotting signals ---
p_mixture1 <- ggplot(df_signals, aes(x = Time_sec, y = Mixture1)) +
  geom_line(color = "blue") +
  labs(title = "(b) Mixtures", x = "Time (s)", y = "Amplitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))  

p_mixture2 <- ggplot(df_signals, aes(x = Time_sec, y = Mixture2)) +
  geom_line(color = "blue") +
  labs(title = "(b) Mixtures", x = "Time (s)", y = "Amplitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))  

p_fft_mixture1 <- plot_fft_normalized(mixture1, fs, title = "(c) FFT Mixture 1 Original", color = "blue")  
p_fft_mixture2 <- plot_fft_normalized(mixture2, fs, title = "(c) FFT Mixture 2 Original", color = "blue") 

p_mixture1_filtered <- ggplot(df_signals, aes(x = Time_sec, y = Mixture1)) +
  geom_line(color = "blue") +
  labs(title = "(d) Filtered mixtures", x = "Time (s)", y = "Amplitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))  

p_mixture2_filtered <- ggplot(df_signals, aes(x = Time_sec, y = Mixture2)) +
  geom_line(color = "blue") +
  labs(title = "(d) Filtered mixtures", x = "Time (s)", y = "Amplitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))  

p_separated1 <- ggplot(df_signals, aes(x = Time_sec, y = Separated1)) +
  geom_line(color = "blue") +
  labs(title = "(e) Signals separated by FastICA", x = "Time (s)", y = "Amplitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))  

p_separated2 <- ggplot(df_signals, aes(x = Time_sec, y = Separated2)) +
  geom_line(color = "blue") +
  labs(title = "(e) Signals separated by FastICA", x = "Time (s)", y = "Amplitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))  

p_recovered1 <- ggplot(df_signals, aes(x = Time_sec, y = Recovered1)) +
  geom_line(color = "blue") +
  labs(title = "(f) Reconstructed signals", x = "Time (s)", y = "Amplitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))  

p_recovered2 <- ggplot(df_signals, aes(x = Time_sec, y = Recovered2)) +
  geom_line(color = "blue") +
  labs(title = "(f) Reconstructed signals", x = "Time (s)", y = "Amplitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))  

p_originalUS <- ggplot(df_signals, aes(x = Time_sec, y = OriginalUS)) +
  geom_line(color = "blue") +
  labs(title = "(a) Original sources", x = "Time (s)", y = "Amplitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))  

p_originalRF <- ggplot(df_signals, aes(x = Time_sec, y = OriginalRF)) +
  geom_line(color = "blue") +
  labs(title = "(a) Original sources", x = "Time (s)", y = "Amplitude") +
  theme_minimal() +
  theme(plot.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))  

# --- Arrange plots ---
g <- grid.arrange(p_originalUS, p_originalRF, p_mixture1, p_mixture2, p_fft_mixture1, p_fft_mixture2,
                  p_mixture1_filtered, p_mixture2_filtered, p_separated1, p_separated2,
                  p_recovered1, p_recovered2, ncol = 2)  

ggsave("simulated_signals_fastica.png", g, width = 12, height = 15, dpi = 300)



