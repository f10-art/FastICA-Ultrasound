############################################################
# run_all.R
# Master script to reproduce all analyses and figures
# FastICA Ultrasound Project — Sanae Fauzi
############################################################

# --- Load modular components ---
source("R/preprocessing.R")
source("R/utils_plots.R")
source("R/simulated_analysis.R")
source("R/experimental_analysis.R")

# --- Create output directories ---
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

cat("============================================\n")
cat(" FASTICA ULTRASOUND PROJECT — FULL PIPELINE\n")
cat("============================================\n\n")

############################################################
# 1. Run simulated analysis
############################################################

cat("Running simulated analysis...\n")

sim_results <- run_simulated_analysis(
  fs = 50e6,
  signal_len = 3500,
  fig_dir = "figures",
  results_dir = "results"
)

cat("\n--- Simulated Analysis Metrics ---\n")
print(sim_results$metrics)

############################################################
# 2. Run experimental analysis
############################################################

cat("\nRunning experimental analysis...\n")

exp_results <- run_experimental_analysis(
  fs = 50e6,
  data_dir = "data",
  fig_dir = "figures",
  results_dir = "results"
)

cat("\n--- Experimental Analysis Metrics ---\n")
print(exp_results$metrics)

############################################################
# 3. Save results to disk
############################################################

save(sim_results, exp_results,
     file = "results/fastica_results.RData")

cat("\nAll analyses completed successfully.\n")
cat("Figures saved in: figures/\n")
cat("Results saved in: results/fastica_results.RData\n")
cat("============================================\n")
