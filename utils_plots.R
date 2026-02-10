############################################################
# utils_plots.R
# Plotting utilities for FastICA ultrasound project
############################################################

library(ggplot2)

############################################################
# Normalized FFT plot
############################################################
# Computes the magnitude spectrum normalized to its maximum.
# Only frequencies below max_freq are displayed.
#
# Arguments:
#   signal   : numeric vector
#   fs       : sampling frequency (Hz)
#   title    : plot title
#   color    : line color
#   max_freq : maximum frequency to display (Hz)
############################################################

plot_fft_normalized <- function(signal, fs, title,
                                color = "blue",
                                max_freq = 4e6) {
  
  n <- length(signal)
  freq <- seq(0, fs/2, length.out = floor(n/2))
  fft_values <- abs(fft(signal))[1:floor(n/2)]
  fft_norm <- fft_values / max(fft_values)
  
  df <- data.frame(Frequency = freq, Magnitude = fft_norm)
  df <- df[df$Frequency <= max_freq, ]
  
  ggplot(df, aes(x = Frequency, y = Magnitude)) +
    geom_line(color = color) +
    labs(title = title,
         x = "Frequency (Hz)",
         y = "Magnitude") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text  = element_text(size = 6)
    )
}

############################################################
# Generic time‑domain plot
############################################################
# Used for plotting mixtures, separated signals, reconstructed
# signals, and reference signals.
############################################################

plot_signal <- function(df, time_col, value_col, title,
                        color = "blue") {
  
  ggplot(df, aes_string(x = time_col, y = value_col)) +
    geom_line(color = color) +
    labs(title = title,
         x = "Time (s)",
         y = "Amplitude") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text  = element_text(size = 6)
    )
}
