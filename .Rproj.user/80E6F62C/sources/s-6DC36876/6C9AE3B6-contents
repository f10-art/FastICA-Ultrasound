############################################################
# simulated_analysis.R
# Full simulated FastICA pipeline for ultrasound signals
############################################################

library(fastICA)
library(ggplot2)
library(gridExtra)

# Load modular utilities
source("R/preprocessing.R")
source("R/utils_plots.R")

############################################################
# Main function
############################################################

run_simulated_analysis <- function(fs = 50e6,
                                   signal_len = 3500,
                                   fig_dir = "figures",
                                   results_dir = "results") {
  
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  
  ############################################################
  # 1. Generate simulated signals
  ############################################################
  
  t <- seq_len(signal_len) / fs
  
  originalUS <- 0.5 * sin(2*pi*1e6*t) +
    0.3 * sin(2*pi*2e6*t)
  
  originalRF <- 0.7 * sin(2*pi*0.5e6*t) +
    0.4 * sin(2*pi*1.5e6*t)
  
  S <- cbind(originalUS, originalRF)
  
  # Random mixing matrix
  set.seed(2)
  A <- matrix(runif(4, 0.1, 1.0), nrow = 2)
  
  mixtures <- S %*% t(A)
  mixtures <- mixtures + matrix(rnorm(2*signal_len, 0, 0.05), ncol = 2)
  
  mixture1 <- mixtures[,1]
  mixture2 <- mixtures[,2]
  
  ############################################################
  # 2. Filtering
  ############################################################
  
  mixture1_f <- full_filter(mixture1, fs)
  mixture2_f <- full_filter(mixture2, fs)
  
  ############################################################
  # 3. Normalization for ICA
  ############################################################
  
  mixture1_n <- norm_center(mixture1_f)
  mixture2_n <- norm_center(mixture2_f)
  
  X <- cbind(mixture1_n, mixture2_n)
  
  ############################################################
  # 4. FastICA
  ############################################################
  
  set.seed(17)
  ica_res <- fastICA(X, n.comp = 2,
                     fun = "logcosh",
                     alpha = 1.5,
                     method = "R",
                     tol = 0.005)
  
  S1 <- ica_res$S[,1]
  S2 <- ica_res$S[,2]
  
  # Sign correction
  normalize_sign <- function(x, ref) {
    if (cor(x, ref) < 0) -x else x
  }
  
  S1 <- normalize_sign(S1, mixture1)
  S2 <- normalize_sign(S2, mixture2)
  
  ############################################################
  # 5. Reconstruction using known A
  ############################################################
  
  comb1 <- A[1,1]*S1 + A[1,2]*S2
  comb2 <- A[2,1]*S1 + A[2,2]*S2
  
  ############################################################
  # 6. Metrics
  ############################################################
  
  metrics <- list(
    cor_S1_US = cor(S1, originalUS),
    cor_S2_RF = cor(S2, originalRF),
    cor_comb1_mixture1 = cor(comb1, mixture1),
    cor_comb2_mixture2 = cor(comb2, mixture2),
    rmse_comb1 = sqrt(mean((comb1 - mixture1)^2)),
    rmse_comb2 = sqrt(mean((comb2 - mixture2)^2))
  )
  
  ############################################################
  # 7. k-fold cross-validation (no re-estimation)
  ############################################################
  
  n <- min(length(S1), length(originalUS))
  k <- 10
  fold_size <- floor(n / k)
  
  cv_S1 <- numeric(k)
  cv_S2 <- numeric(k)
  
  for (i in 1:k) {
    idx <- ((i-1)*fold_size + 1) : min(i*fold_size, n)
    cv_S1[i] <- cor(S1[idx], originalUS[idx], use = "complete.obs")
    cv_S2[i] <- cor(S2[idx], originalRF[idx], use = "complete.obs")
  }
  
  metrics$cv_S1_mean <- mean(cv_S1)
  metrics$cv_S1_sd   <- sd(cv_S1)
  metrics$cv_S2_mean <- mean(cv_S2)
  metrics$cv_S2_sd   <- sd(cv_S2)
  
  ############################################################
  # 8. Build data frame for plotting
  ############################################################
  
  df <- data.frame(
    Time = t,
    Mixture1 = norm_pm1(mixture1_f),
    Mixture2 = norm_pm1(mixture2_f),
    S1 = S1,
    S2 = S2,
    Rec1 = norm_pm1(comb1),
    Rec2 = norm_pm1(comb2),
    OriginalUS = norm_pm1(originalUS),
    OriginalRF = norm_pm1(originalRF)
  )
  
  ############################################################
  # 9. Generate plots
  ############################################################
  
  p1 <- plot_signal(df, "Time", "OriginalUS", "(a) Original US")
  p2 <- plot_signal(df, "Time", "OriginalRF", "(a) Original RF")
  p3 <- plot_signal(df, "Time", "Mixture1", "(b) Mixture 1")
  p4 <- plot_signal(df, "Time", "Mixture2", "(b) Mixture 2")
  p5 <- plot_fft_normalized(mixture1, fs, "(c) FFT Mixture 1")
  p6 <- plot_fft_normalized(mixture2, fs, "(c) FFT Mixture 2")
  p7 <- plot_signal(df, "Time", "S1", "(d) Separated S1")
  p8 <- plot_signal(df, "Time", "S2", "(d) Separated S2")
  p9 <- plot_signal(df, "Time", "Rec1", "(e) Reconstructed 1")
  p10 <- plot_signal(df, "Time", "Rec2", "(e) Reconstructed 2")
  
  g <- grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, ncol = 2)
  
  ggsave(file.path(fig_dir, "simulated_signals_fastica.png"),
         g, width = 12, height = 15, dpi = 300)
  
  ############################################################
  # 10. Return results
  ############################################################
  
  return(list(
    A = A,
    metrics = metrics,
    cv_S1 = cv_S1,
    cv_S2 = cv_S2
  ))
}
