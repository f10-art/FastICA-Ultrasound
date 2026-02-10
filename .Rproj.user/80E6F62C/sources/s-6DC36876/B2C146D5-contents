############################################################
# experimental_analysis.R
# Full experimental FastICA pipeline for ultrasound signals
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

run_experimental_analysis <- function(fs = 50e6,
                                      data_dir = "data",
                                      fig_dir = "figures",
                                      results_dir = "results") {
  
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  
  ############################################################
  # 1. Load experimental data
  ############################################################
  
  mixture1 <- as.numeric(unlist(read.table(file.path(data_dir, "mezcla1.txt"),
                                           sep = ",", dec = ".", header = FALSE)))
  mixture2 <- as.numeric(unlist(read.table(file.path(data_dir, "mezcla2.txt"),
                                           sep = ",", dec = ".", header = FALSE)))
  originalUS <- as.numeric(unlist(read.table(file.path(data_dir, "originalUS.txt"),
                                             sep = ",", dec = ".", header = FALSE)))
  originalRF <- as.numeric(unlist(read.table(file.path(data_dir, "originalRF.txt"),
                                             sep = ",", dec = ".", header = FALSE)))
  
  n <- min(length(mixture1), length(mixture2),
           length(originalUS), length(originalRF))
  
  mixture1 <- mixture1[1:n]
  mixture2 <- mixture2[1:n]
  originalUS <- originalUS[1:n]
  originalRF <- originalRF[1:n]
  
  t <- seq(0, (n-1)/fs, by = 1/fs)
  
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
  # 5. Empirical reconstruction (manual weights)
  ############################################################
  
  w11 <- 0.7; w12 <- 0.2
  w21 <- -0.9; w22 <- 0.1
  
  comb1 <- w11*S1 + w12*S2
  comb2 <- w21*S1 + w22*S2
  
  ############################################################
  # 6. Metrics
  ############################################################
  
  metrics <- list(
    cor_comb1_mixture1 = cor(comb1, mixture1),
    cor_comb2_mixture2 = cor(comb2, mixture2),
    rmse_comb1 = sqrt(mean((comb1 - mixture1)^2)),
    rmse_comb2 = sqrt(mean((comb2 - mixture2)^2)),
    cor_S1_US = cor(S1, originalUS),
    cor_S2_RF = cor(S2, originalRF)
  )
  
  ############################################################
  # 7. k-fold cross-validation (no re-estimation)
  ############################################################
  
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
  # 8. k-fold cross-validation (with ICA re-estimation)
  ############################################################
  
  cv2_S1 <- numeric(k)
  cv2_S2 <- numeric(k)
  
  for (i in 1:k) {
    idx <- ((i-1)*fold_size + 1) : min(i*fold_size, n)
    
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
    ica_fold <- fastICA(fold_data, n.comp = 2,
                        fun = "logcosh",
                        alpha = 1.5,
                        method = "R",
                        tol = 0.005)
    
    S1_f <- normalize_sign(ica_fold$S[,1], mix1_fold)
    S2_f <- normalize_sign(ica_fold$S[,2], mix2_fold)
    
    US_fold <- originalUS[idx]
    RF_fold <- originalRF[idx]
    
    cv2_S1[i] <- cor(S1_f, US_fold, use = "complete.obs")
    cv2_S2[i] <- cor(S2_f, RF_fold, use = "complete.obs")
  }
  
  metrics$cv2_S1_mean <- mean(cv2_S1)
  metrics$cv2_S1_sd   <- sd(cv2_S1)
  metrics$cv2_S2_mean <- mean(cv2_S2)
  metrics$cv2_S2_sd   <- sd(cv2_S2)
  
  ############################################################
  # 9. Local analysis (as in your script)
  ############################################################
  
  t_full <- seq(0, 0.020, by = 1/fs)
  t_full <- t_full[1:n]
  
  # Insert rotated/compressed S1 segment
  idx_local <- which(t_full >= 0.016 & t_full <= 0.017)
  
  S1_seg <- S1[idx_local]
  S1_rot <- rev(S1_seg)
  
  S1_comp <- approx(seq(0,1,length.out=length(S1_rot)),
                    S1_rot,
                    seq(0,1,length.out=length(idx_local)))$y
  
  S1_inserted <- rep(0, length(t_full))
  S1_inserted[idx_local] <- S1_comp
  
  # Local correlations
  metrics$local_S1_US <- cor(S1_inserted[idx_local], originalUS[idx_local])
  
  idx_S2 <- which(t_full <= 0.015)
  metrics$local_S2_RF <- cor(S2[idx_S2], originalRF[idx_S2])
  
  ############################################################
  # 10. Build data frame for plotting
  ############################################################
  
  df <- data.frame(
    Time = t_full,
    Mixture1 = norm_pm1(mixture1_f),
    Mixture2 = norm_pm1(mixture2_f),
    S1_inserted = S1_inserted,
    S2 = norm_pm1(S2),
    Rec1 = norm_pm1(comb1),
    Rec2 = norm_pm1(comb2),
    OriginalUS = norm_pm1(originalUS),
    OriginalRF = norm_pm1(originalRF)
  )
  
  ############################################################
  # 11. Generate plots
  ############################################################
  
  p1 <- plot_signal(df, "Time", "Mixture1", "(a) Original Mixture 1")
  p2 <- plot_signal(df, "Time", "Mixture2", "(a) Original Mixture 2")
  p3 <- plot_fft_normalized(mixture1, fs, "(b) FFT Mixture 1")
  p4 <- plot_fft_normalized(mixture2, fs, "(b) FFT Mixture 2")
  p5 <- plot_signal(df, "Time", "Mixture1", "(c) Filtered Mixture 1")
  p6 <- plot_signal(df, "Time", "Mixture2", "(c) Filtered Mixture 2")
  p7 <- plot_signal(df, "Time", "S1_inserted", "(d) Separated S1")
  p8 <- plot_signal(df, "Time", "S2", "(d) Separated S2")
  p9 <- plot_signal(df, "Time", "Rec1", "(e) Reconstructed 1")
  p10 <- plot_signal(df, "Time", "Rec2", "(e) Reconstructed 2")
  p11 <- plot_signal(df, "Time", "OriginalUS", "(f) Reference US")
  p12 <- plot_signal(df, "Time", "OriginalRF", "(f) Reference RF")
  
  g <- grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, ncol = 2)
  
  ggsave(file.path(fig_dir, "experimental_signals_fastica_real.png"),
         g, width = 12, height = 15, dpi = 300)
  
  ############################################################
  # 12. Return results
  ############################################################
  
  return(list(
    metrics = metrics,
    cv_no_reestimation = list(S1 = cv_S1, S2 = cv_S2),
    cv_reestimation = list(S1 = cv2_S1, S2 = cv2_S2)
  ))
}
