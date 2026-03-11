# --- Libraries ---

library(ggplot2)
library(gridExtra)

# --- Fold-wise correlation data ---
folds <- 1:10

# Variant 1: without ICA re-estimation
S1_no_re <- c(-1, -1, -1, -1, -0.9889, -0.0748, -0.0735, -0.1247, 0.3226, 0.0438)
S2_no_re <- c(1, 1, 1, 1, 0.9921, 0.9671, 0.9578, 0.9046, 0.3047, 0.3727)

# Variant 2: with ICA re-estimation
S1_re <- c(1, 1, 0.0046, 0.0158, 1, -0.1509, -0.0661, 0.052, 0.1029, 0.0291)
S2_re <- c(0.001, 0.0019, 1, 0.9998, 0.0066, -0.4076, 0.9886, -0.8918, -0.3574, 0.3206)

# --- DataFrames for ggplot ---
df_no_re <- data.frame(Fold = folds,
                       S1 = S1_no_re,
                       S2 = S2_no_re,
                       Variant = "(i) Without ICA re-estimation")

df_re <- data.frame(Fold = folds,
                    S1 = S1_re,
                    S2 = S2_re,
                    Variant = "(ii) With ICA re-estimation")

# --- Combine datasets ---
df_long <- rbind(
  data.frame(Fold = df_no_re$Fold, Correlation = df_no_re$S1, Component = "ŝ₁", Variant = df_no_re$Variant),
  data.frame(Fold = df_no_re$Fold, Correlation = df_no_re$S2, Component = "ŝ₂", Variant = df_no_re$Variant),
  data.frame(Fold = df_re$Fold, Correlation = df_re$S1, Component = "ŝ₁", Variant = df_re$Variant),
  data.frame(Fold = df_re$Fold, Correlation = df_re$S2, Component = "ŝ₂", Variant = df_re$Variant)
)

# --- Plot ---
p <- ggplot(df_long, aes(x = Fold, y = Correlation, color = Component, shape = Component)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~Variant, ncol = 1) +
  scale_color_manual(values = c("ŝ₁" = "#0057B8",
                                "ŝ₂" = "#D00000")) +
  scale_shape_manual(values = c("ŝ₁" = 16, "ŝ₂" = 17)) +
  labs(title = "",
       x = "Fold",
       y = "Pearson Correlation") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9, colour = "gray30"),  # <-- CORRECCIÓN AQUÍ
    axis.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# --- Display ---
print(p)

# --- Save in HIGH RESOLUTION ---
ggsave("k-fold.pdf", p,
       width = 3.5,   
       height = 5.5,
       units = "in")




