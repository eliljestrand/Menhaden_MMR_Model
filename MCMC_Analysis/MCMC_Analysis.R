# ==============================================================================
# R Script: ADMB .mcmc Post-Processing & Convergence Diagnostics (No Headers)
# ==============================================================================

library(coda)
library(bayesplot)
library(ggplot2)

# --- 1. Load Data and Metadata ---
# model_name <- "model23" 
# mcmc_file  <- paste0(model_name, ".mcmc")
# names_file <- "names.csv"

# model_name <- "model23Final"
# mcmc_file  <- paste0(model_name, ".mcmc")
# names_file <- "names.csv"

model_name <- "Coston_ME85params"
mcmc_file  <- paste0(model_name, ".mcmc")
names_file <- "names_ault_coston.csv"


# Create a directory with the same name as the model to store outputs
if (!dir.exists(model_name)) {
  dir.create(model_name, recursive = TRUE)
  cat("Created output directory:", model_name, "\n")
}

if (!file.exists(mcmc_file)) {
  stop("The .mcmc file cannot be found. Ensure it is in your working directory.")
}
if (!file.exists(names_file)) {
  stop("The names.csv file cannot be found. Ensure it is in your working directory.")
}

# Read the raw names matrix/vector from the CSV file
raw_names <- read.csv(names_file, header = FALSE, stringsAsFactors = FALSE)
# Flatten the CSV into a clean character vector and trim spaces
col_headers <- trimws(as.character(unlist(raw_names)))

# Load the .mcmc data
raw_mcmc <- read.table(mcmc_file, header = FALSE)

# --- 2. Align Headings and Validate Data ---
if (length(col_headers) != ncol(raw_mcmc)) {
  warning(paste0("Column mismatch! Your names.csv has ", length(col_headers), 
                 " labels, but the .mcmc file has ", ncol(raw_mcmc), " columns.\n",
                 "Assigning generic names fallback to avoid crashing."))
  colnames(raw_mcmc) <- paste0("V", 1:ncol(raw_mcmc))
} else {
  colnames(raw_mcmc) <- col_headers
}

cat("Successfully loaded", nrow(raw_mcmc), "steps with", ncol(raw_mcmc), "tracked variables.\n")


# --- 3. Convert to CODA Object & Apply Burn-in ---
burn_in_pct <- 0.50
n_saved     <- nrow(raw_mcmc)
start_idx   <- floor(n_saved * burn_in_pct) + 1

# Convert into a formal 'mcmc' class for convergence diagnostics
mcmc_chain  <- as.mcmc(raw_mcmc[start_idx:n_saved, , drop = FALSE])
posterior_matrix <- as.matrix(mcmc_chain)


# --- 4. Extract Posterior Values (Median & 95% Credibility Intervals) ---
cat("\n=== EXTRACTED POSTERIOR VALUES ===\n")
posterior_summary_stats <- t(apply(posterior_matrix, 2, function(x) {
  c(Median = median(x), 
    CI_2.5 = quantile(x, probs = 0.025, names = FALSE), 
    CI_97.5 = quantile(x, probs = 0.975, names = FALSE))
}))

print(posterior_summary_stats)
csv_out_path <- file.path(model_name, "extracted_posterior_summaries.csv")
write.csv(posterior_summary_stats, csv_out_path, row.names = TRUE)
cat("Posterior summaries saved to:", csv_out_path, "\n")


# --- 5. Convergence Diagnostics Summary ---
log_out_path <- file.path(model_name, "convergence_diagnostics.txt")
sink(log_out_path)

cat("\n=== POSTERIOR SUMMARY STATS (CODA) ===\n")
print(summary(mcmc_chain))

cat("\n=== EFFECTIVE SAMPLE SIZE (ESS) ===\n")
ess_values <- effectiveSize(mcmc_chain)
print(ess_values)

cat("\n=== GEWEKE CONVERGENCE DIAGNOSTIC ===\n")
geweke_test <- geweke.diag(mcmc_chain)
print(geweke_test)

sink()
cat("Convergence diagnostics text summary saved to:", log_out_path, "\n")


# --- 6. Graphical Analysis (One plot per page in PDFs) ---
n_params <- ncol(posterior_matrix)
param_names <- colnames(posterior_matrix)

# A. Trace Plots
trace_out_path <- file.path(model_name, "traceplots.pdf")
print(paste("Saving Traceplots to", trace_out_path, "..."))
pdf(trace_out_path, width = 8, height = 6)

for (i in 1:n_params) {
  sub_matrix <- posterior_matrix[, i, drop = FALSE]
  p_trace <- mcmc_trace(sub_matrix) + 
    ggtitle(paste0("Trace Plot: ", param_names[i])) + 
    theme_minimal()
  print(p_trace) # Automatically creates a new PDF page
}
dev.off()

# B. Autocorrelation Plots
acf_out_path <- file.path(model_name, "autocorrelation.pdf")
print(paste("Saving Autocorrelation plots to", acf_out_path, "..."))
pdf(acf_out_path, width = 8, height = 6)

for (i in 1:n_params) {
  sub_matrix <- posterior_matrix[, i, drop = FALSE]
  p_acf <- mcmc_acf(sub_matrix, lags = 30) + 
    ggtitle(paste0("Autocorrelation: ", param_names[i])) + 
    theme_minimal()
  print(p_acf)
}
dev.off()

# C. Posterior Densities
dens_out_path <- file.path(model_name, "posterior_densities.pdf")
print(paste("Saving Density plots to", dens_out_path, "..."))
pdf(dens_out_path, width = 8, height = 6)

for (i in 1:n_params) {
  sub_matrix <- posterior_matrix[, i, drop = FALSE]
  p_dens <- mcmc_areas(sub_matrix, prob = 0.5, prob_outer = 0.95) + 
    ggtitle(paste0("Posterior Density: ", param_names[i])) + 
    theme_minimal()
  print(p_dens)
}
dev.off()

cat("\nAnalysis complete. All outputs have been saved to the '", model_name, "' directory.\n", sep = "")