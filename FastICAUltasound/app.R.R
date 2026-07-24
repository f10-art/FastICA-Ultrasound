# ==============================================================================
# SCIENTIFIC APPLICATION: BLIND SOURCE SEPARATION ENGINE (FastICA & Ultrasound)
# ==============================================================================

library(shiny)
library(shinythemes)
library(fastICA)
library(ggplot2)
library(gridExtra)
library(signal)

# --- CREACIÓN AUTOMÁTICA DE CARPETA DE SALIDA ---
if (!dir.exists("Figures")) {
  dir.create("Figures")
}

# --- MAXIMUM UPLOAD SIZE CONFIGURATION ---
options(shiny.maxRequestSize = 50 * 1024^2) 

# --- Explicit normalization and centering ---
norm_center <- function(x) {
  as.numeric(scale(x, center = TRUE, scale = TRUE))
}

# --- Normalization to range [-1, 1] ---
norm_pm1 <- function(x) {
  as.numeric((2 * (x - min(x)) / (max(x) - min(x))) - 1)
}

# --- Filtering functions ---
butter_filter <- function(x, fs, cutoff = 2e6, order = 4) {
  Wn <- cutoff / (fs / 2)
  bf <- butter(order, Wn, type = "low")
  as.numeric(filtfilt(bf, x))
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
  as.numeric(x_hat)
}

savitzky_golay_filter <- function(x, p = 3, n = 11) {
  as.numeric(sgolayfilt(x, p = p, n = n))
}

moving_average_filter <- function(x, n = 5) {
  filt <- stats::filter(x, rep(1 / n, n), sides = 2)
  res <- as.numeric(filt)
  res[is.na(res)] <- x[is.na(res)]
  return(as.numeric(res))
}

clean <- function(x, original) {
  x[is.na(x) | is.infinite(x)] <- original[is.na(x) | is.infinite(x)]
  as.numeric(x)
}

# Complete Filter Pipeline
full_filter <- function(x, fs, cutoff, ma_n, kalman_q, kalman_r, sg_p, sg_n) {
  x1 <- butter_filter(x, fs, cutoff = cutoff)
  x2 <- moving_average_filter(x1, n = ma_n)
  x3 <- kalman_filter(x2, Q = kalman_q, R = kalman_r)
  x4 <- savitzky_golay_filter(x3, p = sg_p, n = sg_n)
  clean(x4, x)
}

# --- Normalized FFT plot ---
plot_fft_normalized <- function(signal, fs, title, color = "blue", max_freq = 4e6) {
  n <- length(signal)
  freq <- seq(0, fs/2, length.out = floor(n/2))
  fft_values <- abs(fft(signal))[1:floor(n/2)]
  fft_norm <- as.numeric(fft_values / max(fft_values))
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

# Safe Sign correction
normalize_sign <- function(x, ref) {
  c_val <- suppressWarnings(cor(as.numeric(x), as.numeric(ref), use = "complete.obs"))
  if (is.na(c_val) || is.nan(c_val)) return(as.numeric(x))
  if (c_val < 0) return(as.numeric(-x)) else return(as.numeric(x))
}


# ==============================================================================
# USER INTERFACE (UI)
# ==============================================================================
ui <- fluidPage(
  theme = shinytheme("cosmo"),
  titlePanel("Blind Source Separation Protocol via FastICA & Signal Preprocessing"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      tags$h4("Data Acquisition Framework"),
      fileInput("file_mix1", "Mixture 1 (mezcla1.txt)", accept = ".txt"),
      fileInput("file_mix2", "Mixture 2 (mezcla2.txt)", accept = ".txt"),
      fileInput("file_refUS", "Reference US (originalUS.txt)", accept = ".txt"),
      fileInput("file_refRF", "Reference RF (originalRF.txt)", accept = ".txt"),
      
      hr(),
      tags$h4("Fixed Parameters"),
      numericInput("fs", "Sampling Frequency (Hz):", value = 50000000, min = 1),
      
      hr(),
      tags$h4("Advanced Scientific Filters"),
      sliderInput("cutoff_freq", "Butterworth Cutoff (Hz):", min = 5e5, max = 1e7, value = 2e6, step = 1e5),
      sliderInput("ma_window", "Moving Average Window (n):", min = 3, max = 21, value = 5, step = 2),
      numericInput("kalman_q", "Kalman Q (Process Noise):", value = 0.01, step = 0.001),
      numericInput("kalman_r", "Kalman R (Measurement Noise):", value = 0.01, step = 0.001),
      sliderInput("sg_p", "Savitzky-Golay Poly Order (p):", min = 1, max = 5, value = 3, step = 1),
      sliderInput("sg_n", "Savitzky-Golay Frame Length (n):", min = 5, max = 31, value = 11, step = 2),
      
      hr(),
      tags$h4("Dynamic Windowing Segment"),
      sliderInput("window_range", "Window Range for S1 (s):",
                  min = 0.000, max = 0.020,
                  value = c(0.01625, 0.01650),
                  step = 0.00025),
      
      hr(),
      tags$h4("Manual Weights Framework"),
      numericInput("w11", "Weight w11:", value = 0.7, step = 0.1),
      numericInput("w12", "Weight w12:", value = 0.2, step = 0.1),
      numericInput("w21", "Weight w21:", value = -0.9, step = 0.1),
      numericInput("w22", "Weight w22:", value = 0.1, step = 0.1)
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs",
        tabPanel("Analysis Dashboard", plotOutput("composite_plot", height = "850px")),
        
        tabPanel("Sensitivity & Correlation Diagnostics", 
                 br(),
                 fluidRow(
                   column(6, plotOutput("p_local_us", height = "300px")),
                   column(6, plotOutput("p_local_rf", height = "300px"))
                 ),
                 br(),
                 fluidRow(
                   column(6, plotOutput("p_sens_butter", height = "300px")),
                   column(6, plotOutput("p_sens_ma", height = "300px"))
                 ),
                 br(),
                 fluidRow(
                   column(6, plotOutput("p_sens_kalman", height = "300px")),
                   column(6, plotOutput("p_stability_folds", height = "300px"))
                 ),
                 br(),
                 fluidRow(
                   column(12, plotOutput("p_pipeline_impact", height = "320px"))
                 )
        ),
        
        tabPanel("Statistical Diagnostics & Metrics",
                 tags$h4("Summary Metrics Table"),
                 tableOutput("summary_table"),
                 hr(),
                 tags$h4("Global Linear Combination & Validation Matrix"),
                 verbatimTextOutput("metrics_output"),
                 hr(),
                 tags$h4("Localized Windows & Stability Diagnostics"),
                 verbatimTextOutput("stability_output"))
      )
    )
  )
)


# ==============================================================================
# SERVER LOGIC
# ==============================================================================
server <- function(input, output, session) {
  
  read_signal <- function(file_input) {
    req(file_input)
    as.numeric(unlist(read.table(file_input$datapath, sep = ",", dec = ".", header = FALSE)))
  }
  
  processed_results <- reactive({
    validate(
      need(input$file_mix1, "Please upload Mixture 1."),
      need(input$file_mix2, "Please upload Mixture 2."),
      need(input$file_refUS, "Please upload Reference US."),
      need(input$file_refRF, "Please upload Reference RF.")
    )
    
    fs <- input$fs
    
    mixture1 <- read_signal(input$file_mix1)
    mixture2 <- read_signal(input$file_mix2)
    originalUS <- read_signal(input$file_refUS)
    originalRF <- read_signal(input$file_refRF)
    
    mixture1_proc <- full_filter(mixture1, fs, cutoff = input$cutoff_freq,
                                 ma_n = input$ma_window,
                                 kalman_q = input$kalman_q,
                                 kalman_r = input$kalman_r,
                                 sg_p = input$sg_p,
                                 sg_n = input$sg_n)
    mixture2_proc <- full_filter(mixture2, fs, cutoff = input$cutoff_freq,
                                 ma_n = input$ma_window,
                                 kalman_q = input$kalman_q,
                                 kalman_r = input$kalman_r,
                                 sg_p = input$sg_p,
                                 sg_n = input$sg_n)
    
    mixture1_norm <- norm_center(mixture1_proc)
    mixture2_norm <- norm_center(mixture2_proc)
    filtered_data <- cbind(mixture1_norm, mixture2_norm)
    
    set.seed(17)
    ica_res <- fastICA(filtered_data, n.comp = 2, row.norm = FALSE,
                       fun = "logcosh", alpha = 1.5, method = "R", tol = 0.005)
    
    S1 <- as.numeric(ica_res$S[,1])
    S2 <- as.numeric(ica_res$S[,2])
    
    S1 <- normalize_sign(S1, mixture1)
    S2 <- normalize_sign(S2, mixture2)
    
    # --- Original linear combinations (manual weights) ---
    comb1 <- input$w11 * S1 + input$w12 * S2
    comb2 <- input$w21 * S1 + input$w22 * S2
    
    # --- Adjusted combinations via linear regression (lm) ---
    fit1 <- lm(mixture1 ~ comb1)
    comb1_adj <- as.numeric(predict(fit1))
    
    fit2 <- lm(mixture2 ~ comb2)
    comb2_adj <- as.numeric(predict(fit2))
    
    # --- Time vector and S1 insertion window ---
    t_total <- seq(0, (length(mixture1)-1)/fs, by = 1/fs)
    Separada1_pegada <- rep(0, length(t_total))
    
    idx_pegar <- which(t_total >= input$window_range[1] &
                         t_total <= input$window_range[2])
    idx_pegar <- idx_pegar[idx_pegar <= length(S1)]
    
    if (length(idx_pegar) > 1) {
      S1_segmento <- S1[idx_pegar]
      US_local_ref <- originalUS[idx_pegar]
      
      if (cor(S1_segmento, US_local_ref) < 0) {
        S1_segmento <- -S1_segmento
      }
      
      S1_temporal_rev <- rev(S1_segmento)
      if (cor(S1_temporal_rev, US_local_ref) > cor(S1_segmento, US_local_ref)) {
        S1_corregida <- S1_temporal_rev
      } else {
        S1_corregida <- S1_segmento
      }
      
      n_out <- length(idx_pegar)
      S1_comprimida <- approx(seq(0, 1, length.out = length(S1_corregida)),
                              S1_corregida,
                              seq(0, 1, length.out = n_out))$y
      Separada1_pegada[idx_pegar] <- as.numeric(S1_comprimida)
    }
    
    list(
      mixture1 = mixture1, mixture2 = mixture2,
      mixture1_proc = mixture1_proc, mixture2_proc = mixture2_proc,
      originalUS = originalUS, originalRF = originalRF,
      S1 = S1, S2 = S2,
      comb1 = comb1, comb2 = comb2,
      comb1_adj = comb1_adj, comb2_adj = comb2_adj,
      Separada1_pegada = Separada1_pegada,
      t_total = t_total, fs = fs
    )
  })
  
  
  # --- Render Main Dashboard ---
  output$composite_plot <- renderPlot({
    res <- processed_results()
    t_len <- length(res$t_total)
    
    df_signals <- data.frame(
      Time_sec = res$t_total,
      Mixture1_Orig = norm_pm1(res$mixture1[1:t_len]),
      Mixture2_Orig = norm_pm1(res$mixture2[1:t_len]),
      Mixture1_Filt = norm_pm1(res$mixture1_proc[1:t_len]),
      Mixture2_Filt = norm_pm1(res$mixture2_proc[1:t_len]),
      Separated1_inserted = as.numeric(res$Separada1_pegada),
      Separated2 = norm_pm1(res$S2[1:t_len]),
      Reconstructed1 = norm_pm1(res$comb1[1:t_len]),
      Reconstructed2 = norm_pm1(res$comb2[1:t_len]),
      OriginalUS = norm_pm1(res$originalUS[1:t_len]),
      OriginalRF = norm_pm1(res$originalRF[1:t_len])
    )
    
    p_mixture1 <- ggplot(df_signals, aes(x = Time_sec, y = Mixture1_Orig)) +
      geom_line(color = "blue") +
      labs(title = "(a) Original Mixture 1", x = "Time (s)", y = "Amplitude") +
      theme_minimal() +
      theme(plot.title = element_text(size = 8),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 6))
    
    p_mixture2 <- ggplot(df_signals, aes(x = Time_sec, y = Mixture2_Orig)) +
      geom_line(color = "blue") +
      labs(title = "(a) Original Mixture 2", x = "Time (s)", y = "Amplitude") +
      theme_minimal() +
      theme(plot.title = element_text(size = 8),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 6))
    
    p_mixture1_filtered <- ggplot(df_signals, aes(x = Time_sec, y = Mixture1_Filt)) +
      geom_line(color = "blue") +
      labs(title = "(c) Filtered Mixture 1", x = "Time (s)", y = "Amplitude") +
      theme_minimal() +
      theme(plot.title = element_text(size = 8),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 6))
    
    p_mixture2_filtered <- ggplot(df_signals, aes(x = Time_sec, y = Mixture2_Filt)) +
      geom_line(color = "blue") +
      labs(title = "(c) Filtered Mixture 2", x = "Time (s)", y = "Amplitude") +
      theme_minimal() +
      theme(plot.title = element_text(size = 8),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 6))
    
    p_separated1 <- ggplot(df_signals, aes(x = Time_sec, y = Separated1_inserted)) +
      geom_line(color = "blue") +
      labs(title = "(d) Signal Separated 1 by FastICA", x = "Time (s)", y = "Amplitude") +
      theme_minimal() +
      theme(plot.title = element_text(size = 8),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 6))
    
    p_separated2 <- ggplot(df_signals, aes(x = Time_sec, y = Separated2)) +
      geom_line(color = "blue") +
      labs(title = "(d) Signal Separated 2 by FastICA", x = "Time (s)", y = "Amplitude") +
      theme_minimal() +
      coord_cartesian(xlim = c(0.00, 0.015), ylim = c(-0.3, 0.3)) +
      theme(plot.title = element_text(size = 8),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 6))
    
    p_reconstructed1 <- ggplot(df_signals, aes(x = Time_sec, y = Reconstructed1)) +
      geom_line(color = "blue") +
      labs(title = "(e) Linearly Combined Signal 1", x = "Time (s)", y = "Amplitude") +
      theme_minimal() +
      theme(plot.title = element_text(size = 8),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 6))
    
    p_reconstructed2 <- ggplot(df_signals, aes(x = Time_sec, y = Reconstructed2)) +
      geom_line(color = "blue") +
      labs(title = "(e) Linearly Combined Signal 2", x = "Time (s)", y = "Amplitude") +
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
    
    p_fft_mixture1 <- plot_fft_normalized(res$mixture1, res$fs,
                                          title = " (b) FFT Original Mixture 1",
                                          color = "blue")
    p_fft_mixture2 <- plot_fft_normalized(res$mixture2, res$fs,
                                          title = " (b) FFT Original Mixture 2",
                                          color = "blue")
    
    grid.arrange(p_mixture1, p_mixture2,
                 p_fft_mixture1, p_fft_mixture2,
                 p_mixture1_filtered, p_mixture2_filtered,
                 p_separated1, p_separated2,
                 p_reconstructed1, p_reconstructed2,
                 p_originalUS, p_originalRF,
                 ncol = 2)
  })
  
  # --- Individual Render Plots with Automatic Exporting ---
  
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
      roll_c <- sapply(1:(n_pts - win), function(i) {
        c_v <- suppressWarnings(cor(s1_seg[i:(i+win)], us_seg[i:(i+win)]))
        if(is.na(c_v)) 0 else abs(c_v)
      })
      df <- data.frame(Time = time_seg[1:(n_pts - win)], Correlation = roll_c)
    } else {
      df <- data.frame(Time = 0, Correlation = 0)
    }
    
    p_local_us <- ggplot(df, aes(x = Time, y = Correlation)) +
      geom_line(color = "darkgreen", linewidth = 0.8) +
      geom_vline(xintercept = c(input$window_range[1], input$window_range[2]),
                 linetype = "dashed", color = "darkred") +
      labs(title = "Rolling Pearson Correlation: S1 vs OriginalUS",
           x = "Time (s)", y = "Absolute Pearson |r|") +
      theme_minimal()
    
    ggsave(filename = "Figures/Local_US.png",
           plot = p_local_us, width = 7, height = 4, dpi = 600)
    p_local_us
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
      roll_c2 <- sapply(1:(n_pts2 - win2), function(i) {
        c_v <- suppressWarnings(cor(s2_seg[i:(i+win2)], rf_seg[i:(i+win2)]))
        if(is.na(c_v)) 0 else abs(c_v)
      })
      df <- data.frame(Time = time_rf[1:(n_pts2 - win2)], Correlation = roll_c2)
    } else {
      df <- data.frame(Time = 0, Correlation = 0)
    }
    
    p_local_rf <- ggplot(df, aes(x = Time, y = Correlation)) +
      geom_line(color = "blue", linewidth = 0.8) +
      labs(title = "Rolling Pearson Correlation: S2 vs OriginalRF (0.00 - 0.015 s)",
           x = "Time (s)", y = "Absolute Pearson |r|") +
      theme_minimal()
    
    ggsave(filename = "Figures/Local_RF.png",
           plot = p_local_rf, width = 7, height = 4, dpi = 600)
    p_local_rf
  })
  
  output$p_sens_butter <- renderPlot({
    res <- processed_results()
    cutoffs <- seq(5e5, 1e7, length.out = 5)
    idx_w <- which(res$t_total >= input$window_range[1] &
                     res$t_total <= input$window_range[2])
    cor_b <- sapply(cutoffs, function(co) {
      m1_f <- butter_filter(res$mixture1, res$fs, cutoff = co)
      c_v <- suppressWarnings(cor(m1_f[idx_w], res$originalUS[idx_w]))
      if(is.na(c_v)) 0 else abs(c_v)
    })
    df <- data.frame(Cutoff_MHz = cutoffs / 1e6, Cor_S1_US = cor_b)
    
    p_sens_butter <- ggplot(df, aes(x = Cutoff_MHz, y = Cor_S1_US)) +
      geom_line(color = "red", linewidth = 0.9) +
      geom_point(color = "darkred", size = 2) +
      scale_y_continuous(labels = function(x) sprintf("%.5f", x)) +
      labs(title = "Butterworth Cutoff Frequency vs Local Signal Quality",
           x = "Cutoff Frequency (MHz)",
           y = "Local Cor |r| (M1 vs US)") +
      theme_minimal()
    
    ggsave(filename = "Figures/Butterworth_Sensitivity.png",
           plot = p_sens_butter, width = 7, height = 4, dpi = 600)
    p_sens_butter
  })
  
  output$p_sens_ma <- renderPlot({
    res <- processed_results()
    ma_wins <- c(3, 5, 9, 13, 17, 21)
    
    idx_w <- which(res$t_total >= input$window_range[1] &
                     res$t_total <= input$window_range[2])
    
    if (length(idx_w) < 5) {
      idx_w <- 1:min(1000, length(res$mixture1))
    }
    
    cor_ma <- sapply(ma_wins, function(w) {
      m1_f <- moving_average_filter(res$mixture1, n = w)
      v_filt <- as.numeric(m1_f[idx_w])
      v_ref  <- as.numeric(res$originalUS[idx_w])
      
      c_v <- suppressWarnings(cor(v_filt, v_ref, use = "pairwise.complete.obs"))
      if (is.na(c_v) || is.nan(c_v)) return(0) else return(abs(c_v))
    })
    
    df <- data.frame(Window = ma_wins, Cor_S1_US = as.numeric(cor_ma))
    
    p_sens_ma <- ggplot(df, aes(x = Window, y = Cor_S1_US)) +
      geom_line(color = "purple", linewidth = 0.9) +
      geom_point(color = "darkblue", size = 2) +
      scale_y_continuous(labels = function(x) sprintf("%.5f", x)) +
      labs(title = "Window Size vs Local Signal Correlation",
           x = "Moving Average Window (n)",
           y = "Local Cor |r| (M1 vs US)") +
      theme_minimal()
    
    ggsave(filename = "Figures/MovingAverage_Sensitivity.png",
           plot = p_sens_ma, width = 7, height = 4, dpi = 600)
    p_sens_ma
  })
  
  output$p_sens_kalman <- renderPlot({
    res <- processed_results()
    grid_q <- c(0.001, 0.01, 0.05)
    grid_r <- c(0.001, 0.01, 0.05)
    kalman_grid <- expand.grid(Q = grid_q, R = grid_r)
    kalman_grid$RMSE <- sapply(1:nrow(kalman_grid), function(k) {
      m1_f <- kalman_filter(res$mixture1[1:min(5000, length(res$mixture1))],
                            Q = kalman_grid$Q[k], R = kalman_grid$R[k])
      sqrt(mean((m1_f - res$mixture1[1:min(5000, length(res$mixture1))])^2))
    })
    
    p_sens_kalman <- ggplot(kalman_grid, aes(x = factor(Q), y = factor(R), fill = RMSE)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.2e", RMSE)),
                color = "white", fontface = "bold", size = 3) +
      scale_fill_gradient(low = "#4575b4", high = "#d73027") +
      labs(title = "Kalman Noise Filtering Impact (Q vs R)",
           x = "Process Noise (Q)", y = "Measurement Noise (R)",
           fill = "RMSE") +
      theme_minimal()
    
    ggsave(filename = "Figures/Kalman_Heatmap.png",
           plot = p_sens_kalman, width = 7, height = 4, dpi = 600)
    p_sens_kalman
  })
  
  output$p_stability_folds <- renderPlot({
    res <- processed_results()
    n_pts <- min(length(res$S1), length(res$originalUS),
                 length(res$S2), length(res$originalRF))
    k_folds <- 10
    fold_sz <- floor(n_pts / k_folds)
    df <- data.frame(Fold = 1:k_folds, S1_US = 0, S2_RF = 0)
    for (i in 1:k_folds) {
      idx_s <- (i-1)*fold_sz + 1
      idx_e <- min(i*fold_sz, n_pts)
      v1 <- suppressWarnings(cor(res$S1[idx_s:idx_e], res$originalUS[idx_s:idx_e]))
      v2 <- suppressWarnings(cor(res$S2[idx_s:idx_e], res$originalRF[idx_s:idx_e]))
      df$S1_US[i] <- ifelse(is.na(v1), 0, v1)
      df$S2_RF[i] <- ifelse(is.na(v2), 0, v2)
    }
    
    p_stability_folds <- ggplot(df) +
      geom_line(aes(x = Fold, y = S1_US, color = "S1 vs OriginalUS"),
                linewidth = 0.8) +
      geom_point(aes(x = Fold, y = S1_US, color = "S1 vs OriginalUS")) +
      geom_line(aes(x = Fold, y = S2_RF, color = "S2 vs OriginalRF"),
                linewidth = 0.8) +
      geom_point(aes(x = Fold, y = S2_RF, color = "S2 vs OriginalRF")) +
      scale_x_continuous(breaks = 1:10) +
      scale_color_manual(values = c("S1 vs OriginalUS" = "darkgreen",
                                    "S2 vs OriginalRF" = "blue")) +
      labs(title = "Localized Temporal Correlation",
           x = "Fold Window Number",
           y = "Localized Pearson Correlation",
           color = "Signal Pair") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    ggsave(filename = "Figures/Stability.png",
           plot = p_stability_folds, width = 7, height = 4, dpi = 600)
    p_stability_folds
  })
  
  output$p_pipeline_impact <- renderPlot({
    res <- processed_results()
    idx_w <- which(res$t_total >= input$window_range[1] &
                     res$t_total <= input$window_range[2])
    v1 <- suppressWarnings(cor(res$mixture1[idx_w], res$originalUS[idx_w]))
    m1_s2 <- butter_filter(res$mixture1, res$fs, cutoff = input$cutoff_freq)
    v2 <- suppressWarnings(cor(m1_s2[idx_w], res$originalUS[idx_w]))
    m1_s3 <- moving_average_filter(m1_s2, n = input$ma_window)
    v3 <- suppressWarnings(cor(m1_s3[idx_w], res$originalUS[idx_w]))
    v4 <- suppressWarnings(cor(res$mixture1_proc[idx_w], res$originalUS[idx_w]))
    
    df <- data.frame(
      Stage = factor(c("1. Without Filtering",
                       "2. Butterworth",
                       "3. Butter + MA",
                       "4. Complete Pipeline"),
                     levels = c("1. Without Filtering",
                                "2. Butterworth",
                                "3. Butter + MA",
                                "4. Complete Pipeline")),
      Correlation = abs(c(ifelse(is.na(v1),0,v1),
                          ifelse(is.na(v2),0,v2),
                          ifelse(is.na(v3),0,v3),
                          ifelse(is.na(v4),0,v4)))
    )
    
    p_pipeline_impact <- ggplot(df, aes(x = Stage, y = Correlation, fill = Stage)) +
      geom_bar(stat = "identity", width = 0.55, show.legend = FALSE) +
      geom_text(aes(label = sprintf("%.4f", Correlation)),
                vjust = -0.4, size = 3.5, fontface = "bold") +
      scale_fill_brewer(palette = "Blues") +
      labs(title = "Effect of Successive Preprocessing Stages on Local Signal Quality",
           x = "Preprocessing Pipeline Stage",
           y = "Local Correlation |r| (Filtered vs US)") +
      theme_minimal()
    
    ggsave(filename = "Figures/Pipeline_Impact.png",
           plot = p_pipeline_impact, width = 8, height = 5, dpi = 600)
    p_pipeline_impact
  })
  
  output$summary_table <- renderTable({
    res <- processed_results()
    
    # --- ORIGINAL METRICS ---
    cor_original1 <- cor(res$comb1, res$mixture1)
    cor_original2 <- cor(res$comb2, res$mixture2)
    rmse_original1 <- sqrt(mean((res$mixture1 - res$comb1)^2))
    rmse_original2 <- sqrt(mean((res$mixture2 - res$comb2)^2))
    
    # --- ADJUSTED METRICS (lm-based) ---
    cor_adjusted1 <- cor(res$comb1_adj, res$mixture1)
    cor_adjusted2 <- cor(res$comb2_adj, res$mixture2)
    rmse_adjusted1 <- sqrt(mean((res$mixture1 - res$comb1_adj)^2))
    rmse_adjusted2 <- sqrt(mean((res$mixture2 - res$comb2_adj)^2))
    
    # --- LOCAL WINDOWS ---
    idx_S1_local <- which(res$t_total >= input$window_range[1] &
                            res$t_total <= input$window_range[2])
    US_local <- res$originalUS[idx_S1_local]
    S1_local <- res$Separada1_pegada[idx_S1_local]
    
    idx_S2_local <- which(res$t_total <= 0.015)
    RF_local <- res$originalRF[idx_S2_local]
    S2_local <- res$S2[1:length(idx_S2_local)]
    
    mean_loc_S1 <- cor(S1_local, US_local)
    mean_loc_S2 <- cor(S2_local, RF_local)
    
    data.frame(
      Metric = c(
        "Corr comb1 (original) vs mixture1",
        "Corr comb1 (adjusted) vs mixture1",
        "RMSE comb1 (original)",
        "RMSE comb1 (adjusted)",
        
        "Corr comb2 (original) vs mixture2",
        "Corr comb2 (adjusted) vs mixture2",
        "RMSE comb2 (original)",
        "RMSE comb2 (adjusted)",
        
        "Mean local corr S1 (US Window)",
        "Mean local corr S2 (RF Window)"
      ),
      
      Value = c(
        sprintf("%.4f", cor_original1),
        sprintf("%.4f", cor_adjusted1),
        sprintf("%.6f", rmse_original1),
        sprintf("%.6f", rmse_adjusted1),
        
        sprintf("%.4f", cor_original2),
        sprintf("%.4f", cor_adjusted2),
        sprintf("%.6f", rmse_original2),
        sprintf("%.6f", rmse_adjusted2),
        
        sprintf("%.4f", mean_loc_S1),
        sprintf("%.4f", mean_loc_S2)
      )
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  
  output$metrics_output <- renderText({
    res <- processed_results()
    
    # Original correlations
    cor_original1 <- cor(res$comb1, res$mixture1)
    cor_original2 <- cor(res$comb2, res$mixture2)
    
    # Adjusted correlations (lm-based)
    cor_adjusted1 <- cor(res$comb1_adj, res$mixture1)
    cor_adjusted2 <- cor(res$comb2_adj, res$mixture2)
    
    # RMSE original vs adjusted (script-style)
    rmse_original1 <- sqrt(mean((res$mixture1 - res$comb1)^2))
    rmse_original2 <- sqrt(mean((res$mixture2 - res$comb2)^2))
    
    rmse_adjusted1 <- sqrt(mean((res$mixture1 - res$comb1_adj)^2))
    rmse_adjusted2 <- sqrt(mean((res$mixture2 - res$comb2_adj)^2))
    
    # Local windows
    idx_S1_local <- which(res$t_total >= input$window_range[1] &
                            res$t_total <= input$window_range[2])
    US_local <- res$originalUS[idx_S1_local]
    S1_local <- res$Separada1_pegada[idx_S1_local]
    
    idx_S2_local <- which(res$t_total <= 0.015)
    RF_local <- res$originalRF[idx_S2_local]
    S2_local <- res$S2[1:length(idx_S2_local)]
    
    paste0(
      "----------------------------------------------------\n",
      "LINEAR COMBINATION DIAGNOSTICS - MANUALLY ADJUSTED WEIGHTS\n",
      "----------------------------------------------------\n\n",
      
      "===== MIXTURE 1 =====\n",
      "Correlation original : ", round(cor_original1, 4), "\n",
      "Correlation adjusted : ", round(cor_adjusted1, 4), "\n",
      "RMSE original        : ", round(rmse_original1, 6), "\n",
      "RMSE adjusted        : ", round(rmse_adjusted1, 6), "\n\n",
      
      "===== MIXTURE 2 =====\n",
      "Correlation original : ", round(cor_original2, 4), "\n",
      "Correlation adjusted : ", round(cor_adjusted2, 4), "\n",
      "RMSE original        : ", round(rmse_original2, 6), "\n",
      "RMSE adjusted        : ", round(rmse_adjusted2, 6), "\n\n",
      
      "----------------------------------------------------\n",
      "LOCAL CORRELATIONS IN VISIBLE WINDOWS\n",
      "----------------------------------------------------\n",
      "Local correlation S1 (inserted, aligned) vs OriginalUS (",
      input$window_range[1], "-", input$window_range[2], " s): ",
      round(cor(S1_local, US_local), 4), "\n",
      "Local correlation S2 vs OriginalRF (0–0.015 s): ",
      round(cor(S2_local, RF_local), 4), "\n"
    )
  })
  
  
  # --- Stability Output ---
  output$stability_output <- renderText({
    res <- processed_results()
    
    n <- min(length(res$S1), length(res$originalUS),
             length(res$S2), length(res$originalRF))
    S1_adj <- as.numeric(res$S1[1:n])
    S2_adj <- as.numeric(res$S2[1:n])
    US_adj <- as.numeric(res$originalUS[1:n])
    RF_adj <- as.numeric(res$originalRF[1:n])
    
    k <- 10
    fold_size <- floor(n / k)
    cv_corr_S1_US <- numeric(k)
    cv_corr_S2_RF <- numeric(k)
    
    for (i in 1:k) {
      idx_start <- (i-1)*fold_size + 1
      idx_end <- min(i*fold_size, n)
      test_idx <- idx_start:idx_end
      cv_corr_S1_US[i] <- as.numeric(cor(S1_adj[test_idx],
                                         US_adj[test_idx],
                                         use = "complete.obs"))
      cv_corr_S2_RF[i] <- as.numeric(cor(S2_adj[test_idx],
                                         RF_adj[test_idx],
                                         use = "complete.obs"))
    }
    
    mean_cv_S1_US <- mean(cv_corr_S1_US, na.rm = TRUE)
    sd_cv_S1_US <- sd(cv_corr_S1_US, na.rm = TRUE)
    mean_cv_S2_RF <- mean(cv_corr_S2_RF, na.rm = TRUE)
    sd_cv_S2_RF <- sd(cv_corr_S2_RF, na.rm = TRUE)
    
    paste0(
      "\n **Localized Windowed Stability Diagnostics (Windows = 10)** \n",
      "  S1 vs OriginalUS: correlations per window = ",
      paste(round(cv_corr_S1_US, 4), collapse = "  "), "\n",
      "   Mean = ", round(mean_cv_S1_US, 4),
      ", Standard deviation = ", round(sd_cv_S1_US, 4), "\n",
      "  S2 vs OriginalRF: correlations per window = ",
      paste(round(cv_corr_S2_RF, 4), collapse = "  "), "\n",
      "   Mean = ", round(mean_cv_S2_RF, 4),
      ", Standard deviation = ", round(sd_cv_S2_RF, 4), "\n"
    )
  })
}

shinyApp(ui = ui, server = server)
