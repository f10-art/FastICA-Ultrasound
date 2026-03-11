# --- Librerías ---
library(ggplot2)
library(gridExtra)

# --- Datos de correlación por fold ---
folds <- 1:10

# Variante 1: sin re-estimación ICA
S1_no_re <- c(-1, -1, -1, -1, -0.9889, -0.0748, -0.0735, -0.1247, 0.3226, 0.0438)
S2_no_re <- c(1, 1, 1, 1, 0.9921, 0.9671, 0.9578, 0.9046, 0.3047, 0.3727)

# Variante 2: con re-estimación ICA
S1_re <- c(1, 1, 0.0046, 0.0158, 1, -0.1509, -0.0661, 0.052, 0.1029, 0.0291)
S2_re <- c(0.001, 0.0019, 1, 0.9998, 0.0066, -0.4076, 0.9886, -0.8918, -0.3574, 0.3206)

# --- DataFrames para ggplot ---
df_no_re <- data.frame(Fold = folds,
                       S1 = S1_no_re,
                       S2 = S2_no_re,
                       Variant = "(i) Sin re-estimación ICA")

df_re <- data.frame(Fold = folds,
                    S1 = S1_re,
                    S2 = S2_re,
                    Variant = "(ii) Con re-estimación ICA")

# --- Combinar ---
df_long <- rbind(
  data.frame(Fold = df_no_re$Fold, Correlation = df_no_re$S1, Component = "ŝ₁", Variant = df_no_re$Variant),
  data.frame(Fold = df_no_re$Fold, Correlation = df_no_re$S2, Component = "ŝ₂", Variant = df_no_re$Variant),
  data.frame(Fold = df_re$Fold, Correlation = df_re$S1, Component = "ŝ₁", Variant = df_re$Variant),
  data.frame(Fold = df_re$Fold, Correlation = df_re$S2, Component = "ŝ₂", Variant = df_re$Variant)
)

# --- Gráfico ---
p <- ggplot(df_long, aes(x = Fold, y = Correlation, color = Component, shape = Component)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~Variant, ncol = 1) +
  scale_color_manual(values = c("ŝ₁" = "blue", "ŝ₂" = "orange")) +
  scale_shape_manual(values = c("ŝ₁" = 16, "ŝ₂" = 17)) +
  labs(title = "Correlaciones por fold en validación k-fold",
       x = "Fold",
       y = "Correlación de Pearson") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.text = element_text(size = 8))

# --- Mostrar ---
print(p)
