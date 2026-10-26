# ==============================================================================
# FastICA Ultrasound Analysis Application (Shiny Interface)
# File: app.R
# ==============================================================================

library(shiny)
library(fastICA)
library(ggplot2)
library(gridExtra)
library(signal)

# --- Helper Functions & Preprocessing Pipeline ---

norm_center <- function(x) {
  as.numeric(scale(x, center = TRUE, scale = TRUE))
}

norm_pm1 <- function(x) {  
  x <- as.numeric(x)
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) return(rep(0, length(x)))
  (2 * (x - r[1]) / diff(r)) - 1  
}  

butter_filter <- function(x, fs, cutoff = 2e6, order = 4) {  
  Wn <- cutoff / (fs / 2) 
  if (Wn >= 1) Wn <- 0.99
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
  if (n %% 2 == 0) n <- n + 1
  sgolayfilt(x, p = p, n = n)  
}  

moving_average_filter <- function(x, n = 5) {  
  filt <- stats::filter(x, rep(1 / n, n), sides = 2)  
  filt <- as.numeric(filt)
  filt[is.na(filt)] <- x[is.na(filt)]  
  filt  
}  

clean <- function(x, original) {  
  x[is.na(x) | is.infinite(x)] <- original[is.na(x) | is.infinite(x)]  
  x  
}  

full_filter <- function(x, fs, cutoff = 2e6, ma_n = 5, kalman_q = 0.01, kalman_r = 0.01, sg_p = 3, sg_n = 11) {  
  x1 <- butter_filter(x, fs, cutoff = cutoff)  
  x2 <- moving_average_filter(x1, n = ma_n)  
  x3 <- kalman_filter(x2, Q = kalman_q, R = kalman_r)  
  x4 <- savitzky_golay_filter(x3, p = sg_p, n = sg_n)  
  clean(x4, x)  
}  

normalize_sign <- function(x, ref) {  
  c_val <- suppressWarnings(cor(x, ref, use = "complete.obs"))
  if (is.finite(c_val) && c_val < 0) return(-x) else return(x)  
}  

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
    theme(plot.title = element_text(size = 9, face = "bold"),  
          axis.title = element_text(size = 8),  
          axis.text = element_text(size = 7))  
} 

cor_por_ventanas <- function(x, y, k = 10) {
  n <- min(length(x), length(y))
  if (n < k * 3) return(rep(NA_real_, k))
  
  x <- as.numeric(x[seq_len(n)])
  y <- as.numeric(y[seq_len(n)])
  indices <- split(seq_len(n), cut(seq_len(n), breaks = k, labels = FALSE))
  
  vapply(indices, function(idx) {
    valid <- is.finite(x[idx]) & is.finite(y[idx])
    xx <- x[idx][valid]
    yy <- y[idx][valid]
    if (length(xx) < 3 || sd(xx) == 0 || sd(yy) == 0) return(NA_real_)
    cor(xx, yy, method = "pearson")
  }, numeric(1))
}

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

evaluar_actividad <- function(componente, referencia, percentil = 0.75) {
  n <- min(length(componente), length(referencia))
  componente <- as.numeric(componente[seq_len(n)])
  referencia <- as.numeric(referencia[seq_len(n)])
  valid <- is.finite(componente) & is.finite(referencia)
  componente <- componente[valid]
  referencia <- referencia[valid]
  n <- length(componente)
  
  if (n < 3) return(NULL)
  
  umbral_comp <- quantile(abs(componente), probs = percentil)
  umbral_ref  <- quantile(abs(referencia), probs = percentil)
  
  act_comp <- abs(componente) >= umbral_comp
  act_ref  <- abs(referencia) >= umbral_ref
  
  TP <- sum(act_comp & act_ref, na.rm = TRUE)
  TN <- sum(!act_comp & !act_ref, na.rm = TRUE)
  FP <- sum(act_comp & !act_ref, na.rm = TRUE)
  FN <- sum(!act_comp & act_ref, na.rm = TRUE)
  
  precision   <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  sensibilidad <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  especificidad <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  
  data.frame(
    Precision = precision,
    Sensitivity = sensibilidad,
    Specificity = especificidad,
    TP = TP, TN = TN, FP = FP, FN = FN, N = n,
    Thresh_Comp = umbral_comp, Thresh_Ref = umbral_ref
  )
}

normalizar_y_completar <- function(x, n) {
  x <- as.numeric(x)
  m <- min(length(x), n)
  x <- x[seq_len(m)]
  resultado <- rep(NA_real_, n)
  valid <- is.finite(x)
  if (any(valid)) {
    xv <- x[valid]
    rango <- max(xv) - min(xv)
    if (rango == 0) resultado[which(valid)] <- 0
    else resultado[which(valid)] <- 2 * (xv - min(xv)) / rango - 1
  }
  resultado
}

# ==============================================================================
# SHINY USER INTERFACE (UI)
# ==============================================================================

ui <- fluidPage(
  titlePanel("Multi-channel Waveform & Spectral Decomposition (FastICA Platform)"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Preprocessing & Filter Controls"),
      hr(),
      numericInput("fs", "Sampling Frequency (Hz):", value = 50000000, step = 1e6),
      numericInput("cutoff_freq", "Butterworth Cutoff (Hz):", value = 2000000, step = 100000),
      numericInput("ma_window", "Moving Average Window (n):", value = 5, min = 1, max = 51, step = 2),
      numericInput("kalman_q", "Kalman Noise Q:", value = 0.01, min = 0.0001, max = 1, step = 0.005),
      numericInput("kalman_r", "Kalman Noise R:", value = 0.01, min = 0.0001, max = 1, step = 0.005),
      numericInput("sg_p", "Savitzky-Golay Order (p):", value = 3, min = 1, max = 5),
      numericInput("sg_n", "Savitzky-Golay Frame (n):", value = 11, min = 5, max = 51, step = 2),
      numericInput("seed", "Random Seed (FastICA):", value = 17, step = 1),
      hr(),
      sliderInput("window_range", "Active US Window Range (s):", 
                  min = 0, max = 0.020, value = c(0.01625, 0.01650), step = 0.00005),
      hr(),
      h4("Linear Reconstruction Weights"),
      numericInput("w11", "w11:", value = 0.7, step = 0.1),
      numericInput("w12", "w12:", value = 0.2, step = 0.1),
      numericInput("w21", "w21:", value = -0.9, step = 0.1),
      numericInput("w22", "w22:", value = 0.1, step = 0.1),
      hr(),
      actionButton("run", "Run Analysis", class = "btn-primary", style = "width:100%;")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Full Pipeline (12 Panels)", plotOutput("plot_full_pipeline", height = "850px")),
        tabPanel("Temporal Stability (k=10)", 
                 plotOutput("plot_stability", height = "400px"),
                 hr(),
                 verbatimTextOutput("text_stability_metrics")),
        tabPanel("Local Rolling Correlations",
                 h4("Local US Window Rolling Correlation"),
                 plotOutput("p_local_us", height = "300px"),
                 hr(),
                 h4("Local RF Window Rolling Correlation"),
                 plotOutput("p_local_rf", height = "300px")),
        tabPanel("Sensitivity Analysis",
                 h4("Butterworth Cutoff Sensitivity"),
                 plotOutput("p_sens_butter", height = "250px"),
                 hr(),
                 h4("Moving Average Window Sensitivity"),
                 plotOutput("p_sens_ma", height = "250px"),
                 hr(),
                 h4("Kalman Grid Search (RMSE)"),
                 plotOutput("p_sens_kalman", height = "250px")),
        tabPanel("Diagnostics & Metrics", 
                 h4("Reconstruction Metrics"),
                 verbatimTextOutput("text_reconstruction"),
                 hr(),
                 h4("Local Rolling Correlations Summary"),
                 verbatimTextOutput("text_rolling"),
                 hr(),
                 h4("Activation Concordance Metrics (75th Percentile)"),
                 verbatimTextOutput("text_concordance"))
      )
    )
  )
)

# ==============================================================================
# SHINY SERVER LOGIC
# ==============================================================================

server <- function(input, output, session) {
  
  pipeline_data <- eventReactive(input$run, ignoreNULL = FALSE, {
    
    if (!file.exists("mezcla1.txt") || !file.exists("mezcla2.txt") ||
        !file.exists("originalUS.txt") || !file.exists("originalRF.txt")) {
      
      fs <- input$fs
      N <- 1000
      t_total <- seq(0, 0.020, length.out = N)
      originalUS <- ifelse(t_total >= 0.01625 & t_total <= 0.01650, sin(2 * pi * 1e6 * t_total), 0)
      originalRF <- ifelse(t_total <= 0.015, cos(2 * pi * 500e3 * t_total), 0)
      mixture1   <- 0.7 * originalUS + 0.3 * originalRF + rnorm(N, sd = 0.05)
      mixture2   <- 0.2 * originalUS + 0.8 * originalRF + rnorm(N, sd = 0.05)
    } else {
      mixture1   <- as.numeric(unlist(read.table("mezcla1.txt", sep = ",", dec = ".", header = FALSE)))
      mixture2   <- as.numeric(unlist(read.table("mezcla2.txt", sep = ",", dec = ".", header = FALSE)))
      originalUS <- as.numeric(unlist(read.table("originalUS.txt", sep = ",", dec = ".", header = FALSE)))
      originalRF <- as.numeric(unlist(read.table("originalRF.txt", sep = ",", dec = ".", header = FALSE)))
      fs <- input$fs
      t_total <- seq(0, (length(mixture1) - 1) / fs, length.out = length(mixture1))
    }
    
    # Preprocessing
    mixture1_proc <- full_filter(
      mixture1, fs, 
      cutoff   = input$cutoff_freq, 
      ma_n     = input$ma_window, 
      kalman_q = input$kalman_q, 
      kalman_r = input$kalman_r, 
      sg_p     = input$sg_p, 
      sg_n     = input$sg_n
    )
    
    mixture2_proc <- full_filter(
      mixture2, fs, 
      cutoff   = input$cutoff_freq, 
      ma_n     = input$ma_window, 
      kalman_q = input$kalman_q, 
      kalman_r = input$kalman_r, 
      sg_p     = input$sg_p, 
      sg_n     = input$sg_n
    )
    
    mixture1_norm <- norm_center(mixture1_proc)
    mixture2_norm <- norm_center(mixture2_proc)
    filtered_data <- cbind(mixture1_norm, mixture2_norm)
    
    # FastICA
    set.seed(input$seed)
    ica_res <- fastICA(filtered_data, n.comp = 2, row.norm = FALSE, fun = "logcosh", alpha = 1.5, method = "R", tol = 0.005)
    
    S1 <- ica_res$S[,1]
    S2 <- ica_res$S[,2]
    
    S1 <- normalize_sign(S1, originalUS)
    S2 <- normalize_sign(S2, originalRF)
    
    # Linear Reconstruction
    comb1 <- input$w11 * S1 + input$w12 * S2
    comb2 <- input$w21 * S1 + input$w22 * S2
    
    fit1 <- lm(mixture1 ~ comb1)
    comb1_adj <- as.numeric(predict(fit1))
    
    fit2 <- lm(mixture2 ~ comb2)
    comb2_adj <- as.numeric(predict(fit2))
    
    # Windows & Indices
    k <- 10
    
    n_US <- min(length(S1), length(originalUS))
    t_US <- (seq_len(n_US) - 1) / fs
    idx_US <- which(t_US >= input$window_range[1] & t_US <= input$window_range[2])
    if (length(idx_US) < 3) idx_US <- seq_len(n_US)
    
    n_RF <- min(length(S2), length(originalRF))
    t_RF <- (seq_len(n_RF) - 1) / fs
    idx_RF <- which(t_RF >= 0 & t_RF <= 0.01500)
    if (length(idx_RF) < 3) idx_RF <- seq_len(n_RF)
    
    cv_corr_S1_US <- cor_por_ventanas(S1[idx_US], originalUS[idx_US], k)
    cv_corr_S2_RF <- cor_por_ventanas(S2[idx_RF], originalRF[idx_RF], k)
    
    # Inserción de S1 comprimida únicamente dentro del rango t_inicio - t_final
    Separada1_pegada <- rep(0, length(t_total))
    t_inicio <- input$window_range[1]
    t_final  <- input$window_range[2]
    idx_pegar <- which(t_total >= t_inicio & t_total <= t_final)
    if (length(idx_pegar) > 0) {
      t_S1_comp <- seq(t_inicio, t_final, length.out = length(S1))
      Separada1_pegada[idx_pegar] <- approx(x = t_S1_comp, y = as.numeric(S1), xout = t_total[idx_pegar], rule = 2)$y
    }
    
    metricas_US <- evaluar_actividad(S1[idx_US], originalUS[idx_US], percentil = 0.75)
    metricas_RF <- evaluar_actividad(S2[idx_RF], originalRF[idx_RF], percentil = 0.75)
    
    list(
      mixture1 = mixture1, mixture2 = mixture2,
      originalUS = originalUS, originalRF = originalRF,
      mixture1_proc = mixture1_proc, mixture2_proc = mixture2_proc,
      S1 = S1, S2 = S2, comb1 = comb1, comb2 = comb2,
      comb1_adj = comb1_adj, comb2_adj = comb2_adj,
      fit1 = fit1, fit2 = fit2,
      cv_corr_S1_US = cv_corr_S1_US, cv_corr_S2_RF = cv_corr_S2_RF,
      t_total = t_total, Separada1_pegada = Separada1_pegada,
      idx_US = idx_US, idx_RF = idx_RF, t_US = t_US, t_RF = t_RF,
      metricas_US = metricas_US, metricas_RF = metricas_RF,
      fs = fs
    )
  })
  
  processed_results <- reactive({
    pipeline_data()
  })
  
  output$p_local_us <- renderPlot({
    res <- processed_results()
    t_start <- max(0, input$window_range[1] - 0.00125)
    t_end <- min(max(res$t_total), input$window_range[2] + 0.00125)
    idx_focus <- which(res$t_total >= t_start & res$t_total <= t_end)
    
    if (length(idx_focus) > 20) {
      s1_seg <- res$Separada1_pegada[idx_focus]
      us_seg <- res$originalUS[idx_focus]
      time_seg <- res$t_total[idx_focus]
      win <- 20
      n_pts <- length(idx_focus)
      roll_c <- vapply(1:(n_pts - win), function(i) {
        c_v <- suppressWarnings(cor(s1_seg[i:(i+win)], us_seg[i:(i+win)]))
        if(is.na(c_v)) 0 else abs(c_v)
      }, numeric(1))
      df <- data.frame(Time = time_seg[1:(n_pts - win)], Correlation = roll_c)
    } else {
      df <- data.frame(Time = 0, Correlation = 0)
    }
    ggplot(df, aes(x = Time, y = Correlation)) +
      geom_line(color = "#2f855a", linewidth = 0.8) +
      geom_vline(xintercept = c(input$window_range[1], input$window_range[2]), linetype = "dashed", color = "#c53030") +
      labs(x = "Time (s)", y = "|r|") + theme_minimal()
  })
  
  output$p_local_rf <- renderPlot({
    res <- processed_results()
    idx_rf <- which(res$t_total >= 0.000 & res$t_total <= 0.015)
    if (length(idx_rf) > 50) {
      s2_seg <- norm_center(res$S2[idx_rf])
      rf_seg <- norm_center(res$originalRF[idx_rf])
      time_rf <- res$t_total[idx_rf]
      win2 <- 50
      n_pts2 <- length(idx_rf)
      roll_c2 <- vapply(1:(n_pts2 - win2), function(i) {
        c_v <- suppressWarnings(cor(s2_seg[i:(i+win2)], rf_seg[i:(i+win2)]))
        if(is.na(c_v)) 0 else abs(c_v)
      }, numeric(1))
      df <- data.frame(Time = time_rf[1:(n_pts2 - win2)], Correlation = roll_c2)
    } else {
      df <- data.frame(Time = 0, Correlation = 0)
    }
    ggplot(df, aes(x = Time, y = Correlation)) +
      geom_line(color = "#2b6cb0", linewidth = 0.8) +
      labs(x = "Time (s)", y = "|r|") + theme_minimal()
  })
  
  output$p_sens_butter <- renderPlot({
    res <- processed_results()
    cutoffs <- seq(5e5, 1e7, length.out = 5)
    idx_w <- which(res$t_total >= input$window_range[1] & res$t_total <= input$window_range[2])
    if (length(idx_w) < 2) idx_w <- 1:length(res$mixture1)
    cor_b <- vapply(cutoffs, function(co) {
      m1_f <- butter_filter(res$mixture1, res$fs, cutoff = co)
      c_v <- suppressWarnings(cor(m1_f[idx_w], res$originalUS[idx_w]))
      if(is.na(c_v)) 0 else abs(c_v)
    }, numeric(1))
    df <- data.frame(Cutoff_MHz = cutoffs / 1e6, Cor_S1_US = cor_b)
    ggplot(df, aes(x = Cutoff_MHz, y = Cor_S1_US)) +
      geom_line(color = "#c53030", linewidth = 0.8) + geom_point(color = "#742a2a", size = 2) +
      labs(x = "Cutoff (MHz)", y = "|r|") + theme_minimal()
  })
  
  output$p_sens_ma <- renderPlot({
    res <- processed_results()
    ma_wins <- c(3, 5, 9, 13, 17, 21)
    idx_w <- which(res$t_total >= input$window_range[1] & res$t_total <= input$window_range[2])
    if (length(idx_w) < 5) idx_w <- 1:min(1000, length(res$mixture1))
    cor_ma <- vapply(ma_wins, function(w) {
      m1_f <- moving_average_filter(res$mixture1, n = w)
      c_v <- suppressWarnings(cor(as.numeric(m1_f[idx_w]), as.numeric(res$originalUS[idx_w])))
      if (is.na(c_v)) 0 else abs(c_v)
    }, numeric(1))
    df <- data.frame(Window = ma_wins, Cor_S1_US = cor_ma)
    ggplot(df, aes(x = Window, y = Cor_S1_US)) +
      geom_line(color = "#6b46c1", linewidth = 0.8) + geom_point(color = "#44337a", size = 2) +
      labs(x = "Window (n)", y = "|r|") + theme_minimal()
  })
  
  output$p_sens_kalman <- renderPlot({
    res <- processed_results()
    kalman_grid <- expand.grid(Q = c(0.001, 0.01, 0.05), R = c(0.001, 0.01, 0.05))
    
    sub_sample <- res$mixture1[seq(1, min(2000, length(res$mixture1)))]
    kalman_grid$RMSE <- vapply(1:nrow(kalman_grid), function(k) {
      m1_f <- kalman_filter(sub_sample, Q = kalman_grid$Q[k], R = kalman_grid$R[k])
      sqrt(mean((m1_f - sub_sample)^2))
    }, numeric(1))
    
    ggplot(kalman_grid, aes(x = factor(Q), y = factor(R), fill = RMSE)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.2e", RMSE)), color = "white", fontface = "bold", size = 3) +
      scale_fill_gradient(low = "#3182ce", high = "#e53e3e") +
      labs(x = "Q", y = "R") + theme_minimal()
  })
  
  # 1. Full Pipeline Composite Plot (12 Panels)
  output$plot_full_pipeline <- renderPlot({
    d <- pipeline_data()
    
    df_signals <- data.frame(
      Time_sec = d$t_total,
      Mixture1_original = normalizar_y_completar(d$mixture1, length(d$t_total)),
      Mixture2_original = normalizar_y_completar(d$mixture2, length(d$t_total)),
      Mixture1_filtered = normalizar_y_completar(d$mixture1_proc, length(d$t_total)),
      Mixture2_filtered = normalizar_y_completar(d$mixture2_proc, length(d$t_total)),
      Separated2 = normalizar_y_completar(d$S2, length(d$t_total)),
      Reconstructed1 = normalizar_y_completar(d$comb1, length(d$t_total)),
      Reconstructed2 = normalizar_y_completar(d$comb2, length(d$t_total)),
      OriginalUS = normalizar_y_completar(d$originalUS, length(d$t_total)),
      OriginalRF = normalizar_y_completar(d$originalRF, length(d$t_total))
    )
    
    p1 <- ggplot(df_signals, aes(x = Time_sec, y = Mixture1_original)) + geom_line(color = "blue") + labs(title = "(a) Original Mixture 1", x = "Time (s)", y = "Amp") + theme_minimal()
    p2 <- ggplot(df_signals, aes(x = Time_sec, y = Mixture2_original)) + geom_line(color = "blue") + labs(title = "(a) Original Mixture 2", x = "Time (s)", y = "Amp") + theme_minimal()
    p3 <- plot_fft_normalized(d$mixture1, d$fs, title = "(b) FFT Original Mixture 1")
    p4 <- plot_fft_normalized(d$mixture2, d$fs, title = "(b) FFT Original Mixture 2")
    p5 <- ggplot(df_signals, aes(x = Time_sec, y = Mixture1_filtered)) + geom_line(color = "blue") + labs(title = "(c) Filtered Mixture 1", x = "Time (s)", y = "Amp") + theme_minimal()
    p6 <- ggplot(df_signals, aes(x = Time_sec, y = Mixture2_filtered)) + geom_line(color = "blue") + labs(title = "(c) Filtered Mixture 2", x = "Time (s)", y = "Amp") + theme_minimal()
    p7 <- ggplot(data.frame(Time_sec = d$t_total, S1 = d$Separada1_pegada), aes(x = Time_sec, y = S1)) + geom_line(color = "blue", linewidth = 0.4) + coord_cartesian(xlim = c(0, 0.020)) + labs(title = "(d) FastICA S1 Inserted in US Window", x = "Time (s)", y = "Amp") + theme_minimal()
    p8 <- ggplot(df_signals, aes(x = Time_sec, y = Separated2)) + geom_line(color = "blue") + coord_cartesian(xlim = c(0.00, 0.015), ylim = c(-0.3, 0.3)) + labs(title = "(d) FastICA S2", x = "Time (s)", y = "Amp") + theme_minimal()
    p9 <- ggplot(df_signals, aes(x = Time_sec, y = Reconstructed1)) + geom_line(color = "blue") + labs(title = "(e) Reconstructed Signal 1", x = "Time (s)", y = "Amp") + theme_minimal()
    p10 <- ggplot(df_signals, aes(x = Time_sec, y = Reconstructed2)) + geom_line(color = "blue") + labs(title = "(e) Reconstructed Signal 2", x = "Time (s)", y = "Amp") + theme_minimal()
    p11 <- ggplot(df_signals, aes(x = Time_sec, y = OriginalUS)) + geom_line(color = "blue") + labs(title = "(f) Experimental Ref (US)", x = "Time (s)", y = "Amp") + theme_minimal()
    p12 <- ggplot(df_signals, aes(x = Time_sec, y = OriginalRF)) + geom_line(color = "blue") + labs(title = "(f) Experimental Ref (RF)", x = "Time (s)", y = "Amp") + theme_minimal()
    
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, ncol = 2)
  })
  
  # 2. Windowed Stability Plot
  output$plot_stability <- renderPlot({
    d <- pipeline_data()
    k <- 10
    df_window_corr <- data.frame(
      Window = factor(rep(1:k, 2)),
      Correlation = c(d$cv_corr_S1_US, d$cv_corr_S2_RF),
      Component = rep(c("S1 (Ultrasonic)", "S2 (Radiofrequency)"), each = k)
    )
    
    ggplot(df_window_corr, aes(x = Window, y = Correlation, group = Component, color = Component, shape = Component)) +
      geom_line(linewidth = 0.8) + geom_point(size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      scale_shape_manual(values = c(16, 17)) + scale_y_continuous(limits = c(-1, 1)) +
      labs(title = "Temporal Component Stability (k = 10)", x = "Window Index", y = "Pearson Correlation") +
      theme_minimal() + theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5), legend.position = "bottom")
  })
  
  output$text_stability_metrics <- renderText({
    d <- pipeline_data()
    paste0(
      "S1 vs US Burst Interval (", input$window_range[1], "-", input$window_range[2], " s):\n",
      "  - Mean Correlation : ", round(mean(d$cv_corr_S1_US, na.rm = TRUE), 4), "\n",
      "  - Std Deviation    : ", round(sd(d$cv_corr_S1_US, na.rm = TRUE), 4), "\n\n",
      "S2 vs RF Signal Interval (0-0.01500 s):\n",
      "  - Mean Correlation : ", round(mean(d$cv_corr_S2_RF, na.rm = TRUE), 4), "\n",
      "  - Std Deviation    : ", round(sd(d$cv_corr_S2_RF, na.rm = TRUE), 4)
    )
  })
  
  # 3. Text Diagnostic Outputs
  output$text_reconstruction <- renderText({
    d <- pipeline_data()
    cor1 <- cor(d$comb1, d$mixture1); rmse_adj1 <- sqrt(mean((d$comb1_adj - d$mixture1)^2)); mae_adj1 <- mean(abs(d$comb1_adj - d$mixture1))
    cor2 <- cor(d$comb2, d$mixture2); rmse_adj2 <- sqrt(mean((d$comb2_adj - d$mixture2)^2)); mae_adj2 <- mean(abs(d$comb2_adj - d$mixture2))
    
    paste0(
      "MIXTURE 1 RECONSTRUCTION:\n",
      "  - Pearson Correlation : ", round(cor1, 4), "\n",
      "  - Adjusted RMSE         : ", round(rmse_adj1, 6), "\n",
      "  - Adjusted MAE          : ", round(mae_adj1, 6), "\n\n",
      "MIXTURE 2 RECONSTRUCTION:\n",
      "  - Pearson Correlation : ", round(cor2, 4), "\n",
      "  - Adjusted RMSE         : ", round(rmse_adj2, 6), "\n",
      "  - Adjusted MAE          : ", round(mae_adj2, 6)
    )
  })
  
  output$text_rolling <- renderText({
    d <- pipeline_data()
    S1_l <- as.numeric(d$S1[d$idx_US]); US_l <- as.numeric(d$originalUS[d$idx_US])
    r_S1 <- calc_rolling_cor(S1_l, US_l, 50)
    
    S2_l <- as.numeric(d$S2[d$idx_RF]); RF_l <- as.numeric(d$originalRF[d$idx_RF])
    r_S2 <- calc_rolling_cor(S2_l, RF_l, 100)
    
    paste0(
      "Ultrasonic Burst Window (", input$window_range[1], "-", input$window_range[2], " s):\n",
      "  - Mean Pearson Correlation      : ", round(cor(S1_l, US_l), 4), "\n",
      "  - Max Absolute Rolling Correlation: ", round(max(abs(r_S1), na.rm = TRUE), 4), "\n\n",
      "Radiofrequency Window (0-0.01500 s):\n",
      "  - Mean Pearson Correlation      : ", round(cor(S2_l, RF_l), 4), "\n",
      "  - Max Absolute Rolling Correlation: ", round(max(abs(r_S2), na.rm = TRUE), 4)
    )
  })
  
  output$text_concordance <- renderText({
    d <- pipeline_data()
    m_US <- d$metricas_US
    m_RF <- d$metricas_RF
    
    paste0(
      "Ultrasonic Component (S1 vs US Ref):\n",
      "  - Precision   : ", round(m_US$Precision, 4), "\n",
      "  - Sensitivity : ", round(m_US$Sensitivity, 4), "\n",
      "  - Specificity : ", round(m_US$Specificity, 4), "\n",
      "  - Confusion Matrix (TP, TN, FP, FN): ", m_US$TP, ", ", m_US$TN, ", ", m_US$FP, ", ", m_US$FN, "\n\n",
      "Radiofrequency Component (S2 vs RF Ref):\n",
      "  - Precision   : ", round(m_RF$Precision, 4), "\n",
      "  - Sensitivity : ", round(m_RF$Sensitivity, 4), "\n",
      "  - Specificity : ", round(m_RF$Specificity, 4), "\n",
      "  - Confusion Matrix (TP, TN, FP, FN): ", m_RF$TP, ", ", m_RF$TN, ", ", m_RF$FP, ", ", m_RF$FN
    )
  })
}

# Run Shiny Application
shinyApp(ui = ui, server = server)
