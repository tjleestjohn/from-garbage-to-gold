########################################################################
########################################################################

# Garbage to Gold Simulation in R

# This script simulates the core concepts of the "From Garbage to Gold" (G2G) paper.

########################################################################
########################################################################

##### --- Clear Environment --- #####
rm(list = ls())
gc()

##### --- Installation of Required Packages --- #####

# Define List of Required Packages
packages <- c("data.table", "MASS", "ranger", "pROC", "PRROC", "ggplot2", 
              "clusterGeneration", "stringr", "viridis", "infotheo", 
              "praznik", "rpart", "parallel", "stats", "gridExtra", "grid")

# Identify Missing Packages
missing_pkgs <- packages[!(packages %in% installed.packages()[, "Package"])]

# Install Missing Packages (if any)
if(length(missing_pkgs)) {
  install.packages(missing_pkgs)
}

# Load All Packages & Clean
invisible(lapply(packages, library, character.only = TRUE))
rm(missing_pkgs)
rm(packages)

##### --- Set Seed for Reproducibility --- #####
set.seed(1)

########################################################################
# --- Set Configuration Parameters ---
########################################################################

##### --- S(1) Latent Space Parameters --- #####

# Number of Latent States.
K_PRIMARY_LATENT <- 4

# Min and Max Prevalence Rates of Latent States.
S1_PREVALENCE_MIN <- 0.05
S1_PREVALENCE_MAX <- 0.50

# Latent State Correlation Matrix Parameters.
CORR_SCENARIO_PARAMS <- list(
  "Zero" = list(alphad = Inf, type = "zero"),
  "Low" = list(alphad = 15, type = "rcorr"),
  "Medium" = list(alphad = .5, type = "rcorr"),
  "High" = list(alphad = 0.00000000000000001, type = "rcorr")
)

cat("--- Correlation Matrix Parameters Defined ---\n\n")

##### --- General Simulation Parameters --- #####

# Sample size for Latent Complexity and Spectral Analysis simulations (Parts 1a and 1b).
# This can be changed within simulation prior to Breadth vs Depth (Part 2).
N_SAMPLES <- 12000

# Number of iterations S'(2) Spectral Analysis Sim (Part 1b) runs.
N_SPECTRAL_ITERATIONS <- 5

# Number of iterations Breadth vs Depth Sim (Part 3) and Comparative Spectral Analysis Sim (Part 4) run.
# Additional Iterations of Part 3 can be added within the simulation after the first iteration finishes.
N_SUPER_ITERATIONS <- 5

# Test N to Train N Ratio.
TEST_SET_MULTIPLIER <- 0.25

# Global Artifact Container.
last_run_artifacts <- list()

##### --- S(2) and Y GENERATION PARAMETERS --- #####

# Dynamic S(2) Haystack Sizing.
# Number of Proxies Per Possible Latent Configuration.
POOL_REDUNDANCY_FACTOR <- 1000
M_POOL_SIZE <- K_PRIMARY_LATENT * POOL_REDUNDANCY_FACTOR

# S(2) Feature Selection Budget (e.g., 33.33% of the Haystack).
BUDGET_PCT <- 1/3
M_PREDICTORS <- floor(M_POOL_SIZE * BUDGET_PCT)

# "Dynamic Range" of the S(2) and Y.
PARAMETRIC_MEAN_MIN <- 0.25
PARAMETRIC_MEAN_MAX <- 0.75

# Variable Default Mode (Will be overwritten by Menu)
S2_GENERATION_MODE <- "CONSISTENT"

# Causal Consistency - Interaction Order Limit 
# 1 = Main Effects Only (Linear)
# 2 = Main Effects + Pairwise Interactions
# K_PRIMARY_LATENT = Full Complexity (Main Effects and All Possible Interactions)
MAX_INTERACTION_ORDER <- K_PRIMARY_LATENT

# Percent of S(2) Predictors that are Pure Noise (fully unrelated to S(1) states or configurations)
PURE_NOISE_PERCENTAGE <- 1/4

##### --- Observational Error Parameters for S(2) --- #####

# Set the peak TARGET_CONTRAST_RATIO to a numeric value (e.g., 5.0) to auto-calculate fidelity.
# Set to NA to use the manual OBS_FIDELITY_MIN/MAX defaults below.
TARGET_CONTRAST_RATIO <- NA

# Default manual values (Used if Target Ratio is NA)
OBS_FIDELITY_MIN <- .875
OBS_FIDELITY_MAX <- .925

##### --- Random Forest & Strategy Parameters --- #####
RF_NUM_TREES <- 1000
# % of threads using in ranger() models.
# note: 1 less thread (than specified %) is used in ranger() for system stability.
PCT_THREADS_RANGER <- 1

##### --- Unsupervised P-DCAI Feature Selection Parameters --- #####

# Overwritten by user choice
STRATEGY_CHOICE <- "Random Breadth"

# Weakest link reinforcement smoothing for unsupervised P-DCAI (Spectral)
# Smoothing to prevent "noisy" clusters from dominating weakest link selection
# 1 = no probability dampening (weighted sampling exactly proportional to distance)
# 0.5 = Square root probability smoothing
# 0.0 = no weighted sampling (randomly selects which factor to reinforce with equal probability).
PROB_DAMPENING_EXPONENT <- .75

##### --- Plot/Graph Parameters --- #####
USE_LOG_SCALE_X <- FALSE
PLOT_NO_ERROR_LINE <- FALSE
GENERATE_SPARSITY_GRAPH <- TRUE
GENERATE_AUROC_GRAPH <- FALSE
GENERATE_AUPRC_GRAPH <- FALSE
GENERATE_COND_ENTROPY_GRAPH <- FALSE
ANALYSIS_MODE <- "AUROC/AUPRC"

########################################################################
########################################################################

##### --- Shared Spectral Analysis Helpers --- #####

# 1. Generate Variance-Matched Parallel Analysis Noise Floor
# Correlation Mode (Scale = TRUE)
get_noise_floor <- function(n_rows, col_sds, n_iters = 3) {
  n_cols <- length(col_sds)
  sum_evals <- numeric(n_cols)
  
  for(i in 1:n_iters) {
    # 1. Generate Standard Normal (Var = 1)
    noise_mat <- matrix(rnorm(n_rows * n_cols), nrow = n_rows, ncol = n_cols)
    
    # 2. Standardize Noise (Correlation Strategy)
    noise_scaled <- scale(noise_mat, center = TRUE, scale = TRUE)
    
    # 3. Calculate Eigenvalues
    c_mat <- crossprod(noise_scaled) / (n_rows - 1)
    ev <- eigen(c_mat, symmetric = TRUE, only.values = TRUE)$values
    
    len <- min(length(ev), length(sum_evals))
    sum_evals[1:len] <- sum_evals[1:len] + ev[1:len]
  }
  return(sum_evals / n_iters)
}

# 2. Robust Elbow Finder (Mode-Aware)
find_spectral_elbow <- function(eig_vals, noise_floor) {
  n <- length(eig_vals)
  window_size <- if (n > 55) 50 else max(3, floor(n/3))
  
  if(n < 5) return(1)
  
  # Clamp eigenvalues
  eig_vals <- pmax(eig_vals, 1e-9)
  noise_floor <- pmax(noise_floor, 1e-9)
  
  # --- DYNAMIC SENSITIVITY CONFIGURATION ---
  # Set elbow detection sensitivity by Generation Mode (if distinct criteria are necessary).
  current_mode <- if(exists("S2_GENERATION_MODE")) S2_GENERATION_MODE else "CONSISTENT"
  if (current_mode == "CHAOTIC") {
    # CHAOTIC MODE:
    min_sigma <- 0.025
    sigma_multiplier <- 2.5
  } else {
    # CONSISTENT MODE:
    min_sigma <- 0.025
    sigma_multiplier <- 2.5
  }
  
  # A. Find last component strictly above noise
  check_len <- min(length(eig_vals), length(noise_floor))
  # We add a tiny buffer proportional to sensitivity to ensure we don't just graze the noise
  noise_buffer <- if(current_mode == "CHAOTIC") 1e-10 else 0.01 
  
  above_noise_indices <- which(eig_vals[1:check_len] > (noise_floor[1:check_len] + noise_buffer))
  
  if (length(above_noise_indices) == 0) return(1)
  
  max_signal_idx <- max(above_noise_indices)
  start_k <- min(max_signal_idx + 1, n - 2)
  if (start_k < 3) return(1)
  
  # B. Backward Scan (Triple-Tap)
  log_vals <- log(eig_vals)
  
  for(k in seq(start_k, 3, by=-1)) {
    window_end <- min(n, k + window_size)
    tail_idxs <- k:window_end
    if (length(tail_idxs) < 3) next 
    
    tail_log_vals <- log_vals[tail_idxs]
    fit <- lm(tail_log_vals ~ tail_idxs)
    sigma <- summary(fit)$sigma
    if(is.nan(sigma) || sigma < min_sigma) sigma <- min_sigma
    
    cand_idx <- k - 1
    pred_cand <- predict(fit, newdata=data.frame(tail_idxs=cand_idx))
    actual_cand <- log_vals[cand_idx]
    
    # Check 1: Elbow Shape
    is_elbow_shape <- FALSE
    if(actual_cand > (pred_cand + (sigma_multiplier * sigma))) {
      confirmed <- TRUE
      if(cand_idx > 1) {
        pred_1 <- predict(fit, newdata=data.frame(tail_idxs=cand_idx-1))
        if(log_vals[cand_idx-1] <= (pred_1 + (sigma_multiplier * sigma))) confirmed <- FALSE
      }
      if(confirmed && cand_idx > 2) {
        pred_2 <- predict(fit, newdata=data.frame(tail_idxs=cand_idx-2))
        if(log_vals[cand_idx-2] <= (pred_2 + (sigma_multiplier * sigma))) confirmed <- FALSE
      }
      if(confirmed) is_elbow_shape <- TRUE
    }
    
    # Check 2: Above Noise
    if (is_elbow_shape && cand_idx <= length(noise_floor)) {
      if (eig_vals[cand_idx] > noise_floor[cand_idx]) return(cand_idx)
    }
  }
  return(1)
}

########################################################################
########################################################################

##### --- PART 1a: Latent Sparsity Analysis --- #####

generate_s1 <- function(n_samples, corr_matrix, prevalence_rates) {
  mean_vec <- rep(0, ncol(corr_matrix))
  mvn_data <- mvrnorm(n = n_samples, mu = mean_vec, Sigma = corr_matrix)
  
  if(is.null(dim(mvn_data))) {
    mvn_data <- matrix(mvn_data, ncol = 1)
  }
  
  thresholds <- qnorm(1 - prevalence_rates)
  s1_data <- t(t(mvn_data) > thresholds) * 1
  storage.mode(s1_data) <- "integer"
  return(s1_data)
}

calculate_latent_sparsity <- function(s1_matrix) {
  s1_dt <- as.data.table(s1_matrix)
  freq_dt <- s1_dt[, .N, by = names(s1_dt)]
  probs <- freq_dt$N / nrow(s1_dt)
  
  joint_entropy <- -sum(probs * log2(probs))
  k_rlzd <- nrow(freq_dt)
  min_prevalence <- min(probs)
  
  return(list(entropy = joint_entropy, k_rlzd = k_rlzd, min_prevalence = min_prevalence))
}

run_latent_sparsity_analysis <- function(n_iterations = 100) {
  cat("--- Running Part 1a: Latent Complexity Simulation ---\n")
  results_list <- list()
  
  for (name in names(CORR_SCENARIO_PARAMS)) {
    cat(paste("  Processing Scenario:", name, "\n"))
    iteration_results <- replicate(n_iterations, {
      params <- CORR_SCENARIO_PARAMS[[name]]
      corr_matrix <- if (params$type == "zero") diag(K_PRIMARY_LATENT) else rcorrmatrix(K_PRIMARY_LATENT, alphad = params$alphad)
      avg_abs_corr <- mean(abs(corr_matrix[upper.tri(corr_matrix)]))
      s1_prevalences <- runif(K_PRIMARY_LATENT, S1_PREVALENCE_MIN, S1_PREVALENCE_MAX)
      s1 <- generate_s1(N_SAMPLES, corr_matrix, s1_prevalences) 
      metrics <- calculate_latent_sparsity(s1)
      list(entropy = metrics$entropy, k_rlzd = metrics$k_rlzd, avg_corr = avg_abs_corr, min_prevalence = metrics$min_prevalence)
    }, simplify = FALSE)
    
    entropies <- sapply(iteration_results, function(res) res$entropy)
    k_rlzd_values <- sapply(iteration_results, function(res) res$k_rlzd)
    avg_corrs <- sapply(iteration_results, function(res) res$avg_corr)
    k_eff_values <- 2^entropies
    min_prevalences <- sapply(iteration_results, function(res) res$min_prevalence)
    
    results_list[[name]] <- data.table(
      Scenario = name,
      Mean_K_eff = mean(k_eff_values), SD_K_eff = sd(k_eff_values),
      Mean_K_rlzd = mean(k_rlzd_values), SD_K_rlzd = sd(k_rlzd_values),
      Mean_Avg_Corr = mean(avg_corrs), Mean_Min_Prevalence = mean(min_prevalences)
    )
  }
  cat("Latent Complexity Simulation Complete.\n")
  return(rbindlist(results_list))
}

plot_effective_configurations <- function(results_df, n_iterations) {
  results_df[, Scenario_Label := sprintf("%s\n(Avg. Abs. Corr. = %.2f)", Scenario, Mean_Avg_Corr)]
  results_df$Scenario_Label <- factor(results_df$Scenario_Label, levels = results_df$Scenario_Label)
  
  k_eff_data <- results_df[, .(Scenario_Label, Mean = Mean_K_eff, SD = SD_K_eff)]
  k_eff_data[, Metric := "K_eff (Effective)"]  
  k_rlzd_data <- results_df[, .(Scenario_Label, Mean = Mean_K_rlzd, SD = SD_K_rlzd)]
  k_rlzd_data[, Metric := "K_rlzd (Realized)"]  
  
  plot_data <- rbind(k_eff_data, k_rlzd_data)
  plot_data[, Metric := factor(Metric, levels = c("K_rlzd (Realized)", "K_eff (Effective)"))]
  
  max_y_val <- results_df[Scenario == "Zero", Mean_K_rlzd] + results_df[Scenario == "Zero", SD_K_rlzd]
  theoretical_max <- 2^K_PRIMARY_LATENT
  
  p <- ggplot(plot_data, aes(x = Scenario_Label, y = Mean, fill = Metric)) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = pmax(0, Mean - SD), ymax = Mean + SD), width = 0.25, color = "gray20", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = sprintf("%.1f", Mean)), vjust = 1.5, color = "white", size = 3.5, position = position_dodge(width = 0.9)) +
    annotate("text", x = Inf, y = Inf, label = paste("Theoretical Max:", format(theoretical_max, big.mark=",")), vjust = 2, hjust = 1.1, color = "gray30", size = 3.5) +
    labs(title = "Latent Complexity Simulation",
         subtitle = sprintf("Results are based on %d iterations.", n_iterations),
         x = "Correlation Scenario", y = "Number of S(1) Configurations", fill = "Metric",
         caption = sprintf("Parameters: k = %d, N = %s, S(1) Prevalence Range = [%.2f, %.2f]",
                           K_PRIMARY_LATENT, format(N_SAMPLES, big.mark=","), S1_PREVALENCE_MIN, S1_PREVALENCE_MAX)
    ) +
    scale_fill_viridis_d(option="D", begin=0.3, end=0.8) +
    theme_minimal(base_size = 14) +  
    theme(legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 10)) +
    coord_cartesian(ylim = c(0, max_y_val))
  
  ggsave("effective_configurations_analysis_R.png", p, width = 12, height = 7, dpi = 300)
  print(p)
}

s1_to_int <- function(s1) { 
  if(is.null(dim(s1))) s1 <- matrix(s1, ncol=length(s1))
  # Match expand.grid order (First column is 2^0, Second is 2^1, etc.)
  # This ensures the lookup table aligns with the design matrix structure.
  s1 %*% (2^(0:(ncol(s1)-1))) 
}

########################################################################
########################################################################

##### --- Generating Observed Data  --- #####

create_lookup_tables <- function(k, m,
                                 mean_min = PARAMETRIC_MEAN_MIN,
                                 mean_max = PARAMETRIC_MEAN_MAX,
                                 max_order = MAX_INTERACTION_ORDER,
                                 noise_pct = PURE_NOISE_PERCENTAGE,
                                 mode = S2_GENERATION_MODE,
                                 print_diagnostics = FALSE) { 
  
  n_configs <- 2^k
  m_noise <- floor(m * noise_pct)
  m_signal <- m - m_noise
  
  s1_combinations <- expand.grid(replicate(k, 0:1, simplify = FALSE))
  names(s1_combinations) <- paste0("X", 1:k)
  
  eff_order <- min(k, max(1, max_order))
  
  if (eff_order == 1) {
    formula_str <- as.formula("~ .")
  } else {
    formula_str <- as.formula(paste0("~ .^", eff_order))
  }
  
  design_mat <- model.matrix(formula_str, data = s1_combinations)[, -1, drop = FALSE]
  
  # --- HIERARCHICAL SCALING ---
  term_names <- colnames(design_mat)
  interaction_orders <- stringr::str_count(term_names, ":") + 1
  term_sds <- sqrt(1 / factorial(interaction_orders))
  n_terms <- ncol(design_mat)
  
  apply_copula_mapping <- function(scores, min_bound, range_val, allow_flip = TRUE) {
    
    # Safety: If variance is zero (all scores identical), return mid-point
    # We use sd() checks here instead of max-min
    if (sd(scores) < 1e-9) {
      return(rep(min_bound + (range_val / 2), length(scores)))
    }
    
    # 1. Standardize the scores (Z-Score)
    # This centers the data at 0 and scales it so Standard Deviation is 1.
    z_scores <- as.vector(scale(scores))
    
    # 2. Apply the Probability Integral Transform
    # pnorm() maps the Z-scores to the range [0, 1] using the Gaussian Curve.
    # This flattens the Bell Curve into a Uniform Distribution.
    scaled_0_1 <- pnorm(z_scores)
    
    # --------------------------
    
    # 3. Project to Target Range (e.g., PARAMETRIC_MEAN_MIN to PARAMETRIC_MEAN_MAX)
    final_probs <- min_bound + (scaled_0_1 * range_val)
    
    # 4. Mirror across 0.5 if flipped
    if (allow_flip && runif(1) > 0.5) {
      final_probs <- 1 - final_probs
    }
    
    return(final_probs)
  }
  
  # ==============================================================================
  # 1. Y GENERATION
  # ==============================================================================
  Y_BOUND_MIN <- 0
  Y_BOUND_MAX <- 1
  Y_RANGE <- Y_BOUND_MAX - Y_BOUND_MIN
  
  if (mode == "CONSISTENT") {
    coeffs_y <- numeric(n_terms)
    is_main <- interaction_orders == 1
    is_int <- interaction_orders > 1
    n_main <- sum(is_main); n_int <- sum(is_int); int_scaling <- 6
    
    if (n_main > 0) coeffs_y[is_main] <- rnorm(n_main, mean = 0, sd = 1)
    if (n_int > 0)  coeffs_y[is_int]  <- rnorm(n_int, mean = 0, sd = term_sds[is_int] * int_scaling)
    
    raw_scores_y <- as.vector(design_mat %*% coeffs_y)
    
  } else {
    raw_scores_y <- rnorm(n_configs, mean = 0, sd = 1)
  }
  
  y_probs <- apply_copula_mapping(raw_scores_y, Y_BOUND_MIN, Y_RANGE, allow_flip = FALSE)
  y_lookup <- data.table(config_id = 0:(n_configs - 1), prob = as.vector(y_probs))
  
  # ==============================================================================
  # 2. S(2) GENERATION
  # ==============================================================================
  S2_RANGE <- mean_max - mean_min
  
  if (m_signal > 0) {
    if (mode == "CONSISTENT") {
      coeffs_s2 <- matrix(0, nrow = n_terms, ncol = m_signal)
      is_main <- interaction_orders == 1
      is_int <- interaction_orders > 1
      n_main <- sum(is_main); n_int <- sum(is_int); int_scaling <- .1
      
      for(j in 1:m_signal) {
        if (n_main > 0) coeffs_s2[is_main, j] <- rnorm(n_main, mean = 0, sd = 1)
        if (n_int > 0)  coeffs_s2[is_int, j]  <- rnorm(n_int, mean = 0, sd = term_sds[is_int] * int_scaling)
      }
      raw_scores_s2 <- design_mat %*% coeffs_s2
    } else {
      raw_scores_s2 <- matrix(rnorm(n_configs * m_signal, mean = 0, sd = 1), nrow=n_configs)
    }
    
    final_probs_signal <- apply(raw_scores_s2, 2, function(col) {
      apply_copula_mapping(col, mean_min, S2_RANGE, allow_flip = TRUE)
    })
    
  } else {
    final_probs_signal <- matrix(0, nrow = n_configs, ncol = 0)
  }
  
  # ==============================================================================
  # 3. NOISE GENERATION (PURE ROW RANDOM)
  # ==============================================================================
  noise_indices <- c() # Track which columns are noise
  
  if (m_noise > 0) {
    noise_prob_matrix <- matrix(0, nrow = n_configs, ncol = m_noise)
    
    for(j in 1:m_noise) {
      # 1. Independent Random Draws
      # We create a "Bag of Probabilities" for this variable using the same normal distribution.
      raw_noise <- rnorm(n_configs, mean = 0, sd = 1)
      noise_prob_matrix[, j] <- apply_copula_mapping(raw_noise, mean_min, S2_RANGE, allow_flip = TRUE)
    }
    
    final_probs <- cbind(final_probs_signal, noise_prob_matrix)
    
    # Shuffle columns
    total_cols <- ncol(final_probs)
    mix_idx <- sample(total_cols)
    final_probs <- final_probs[, mix_idx, drop=FALSE]
    
    # Identify where the noise columns ended up
    original_noise_positions <- (m_signal + 1):total_cols
    noise_indices <- which(mix_idx %in% original_noise_positions)
    
  } else {
    final_probs <- final_probs_signal
  }
  
  return(list(y_lookup = y_lookup, 
              s2_prob_matrix = final_probs, 
              noise_indices = noise_indices, 
              intended_targets = NULL, 
              used_sd = 0))
}

generate_truth <- function(s1_input, y_lookup, s2_prob_matrix, noise_indices = NULL) {
  
  if (is.matrix(s1_input)) {
    s1_int <- as.vector(s1_to_int(s1_input))
  } else {
    s1_int <- as.vector(s1_input) 
  }
  
  # 1. Generate Y Truth
  y_probs <- y_lookup$prob[s1_int + 1] 
  y_true <- rbinom(length(y_probs), 1, y_probs)
  
  # 2. Generate S2 (Feature) Truth
  n_samples <- length(s1_int)
  m <- ncol(s2_prob_matrix)
  n_configs_in_matrix <- nrow(s2_prob_matrix)
  
  # Create a matrix of lookup indices
  # Default (Signal): Match the S1 Configuration (Row -> S1 State)
  lookup_indices <- matrix(rep(s1_int + 1, m), nrow = n_samples, ncol = m)
  
  # --- PURE NOISE LOGIC ---
  # For noise columns, we completely ignore S(1).
  # We just sample uniformly from the "Bag of Probabilities" (Rows 1 to n_configs).
  if (!is.null(noise_indices) && length(noise_indices) > 0) {
    for (idx in noise_indices) {
      # replace = TRUE ensures independent random assignment.
      lookup_indices[, idx] <- sample(1:n_configs_in_matrix, n_samples, replace = TRUE)
    }
  }
  
  # Efficient Matrix Lookup
  linear_idxs <- cbind(as.vector(lookup_indices), rep(1:m, each = n_samples))
  prob_vector <- s2_prob_matrix[linear_idxs]
  
  s2_prob_matrix_samples <- matrix(prob_vector, nrow = n_samples, ncol = m)
  
  # Simulate outcomes
  s2_true_matrix <- matrix(rbinom(length(s2_prob_matrix_samples), 1, s2_prob_matrix_samples),
                           nrow = n_samples, ncol = m)
  
  s2_true_dt <- as.data.table(s2_true_matrix)
  setnames(s2_true_dt, paste0("V", 1:ncol(s2_true_dt)))
  
  return(list(y_true = y_true, s2_true_dt = s2_true_dt))
}

apply_observational_error <- function(s2_true_dt, min_fid, max_fid) {
  s2_true_matrix <- as.matrix(s2_true_dt)
  n_rows <- nrow(s2_true_matrix); m <- ncol(s2_true_matrix)
  
  # Single regime using provided bounds
  draws <- runif(m * 2, min = min_fid, max = max_fid)
  
  # Reshape draws into a matrix of parameters (Alpha, Beta) for each feature
  params <- matrix(draws, ncol = 2)
  alpha <- params[, 1]
  beta_spec <- params[, 2]
  
  prevalences <- colMeans(s2_true_matrix)
  fn_rate <- 1 - alpha; fp_rate <- 1 - beta_spec
  weighted_error_rates <- (fn_rate * prevalences) + (fp_rate * (1 - prevalences))
  avg_error_rate <- mean(weighted_error_rates)
  
  alpha_mat <- matrix(alpha, n_rows, m, byrow = TRUE)
  beta_spec_mat <- matrix(beta_spec, n_rows, m, byrow = TRUE)
  p_obs_matrix <- s2_true_matrix * alpha_mat + (1 - s2_true_matrix) * (1 - beta_spec_mat)
  s_prime_2_matrix <- matrix(rbinom(length(p_obs_matrix), 1, p_obs_matrix), nrow = n_rows, ncol = m)
  s_prime_2_dt <- as.data.table(s_prime_2_matrix)
  setnames(s_prime_2_dt, names(s2_true_dt))
  
  return(list(s_prime_2_dt = s_prime_2_dt, avg_error_rate = avg_error_rate))
}

# --- UNIFIED SIMULATION GENERATOR ---
# This ensures Part 1b, Part 3, and Part 4 use Identical Physics.
generate_universe <- function(seed, 
                              n_samples, 
                              k_latent = K_PRIMARY_LATENT, 
                              m_pool = M_POOL_SIZE,
                              s1_min = S1_PREVALENCE_MIN,
                              s1_max = S1_PREVALENCE_MAX,
                              fid_min = OBS_FIDELITY_MIN,
                              fid_max = OBS_FIDELITY_MAX) {
  
  # 1. Set Seed for Exact Reproducibility
  set.seed(seed)
  
  # 2. Generate S(1) Latent Structure
  #    (Used for Latent Sparsity Sim)
  s1_prevalences <- runif(k_latent, s1_min, s1_max)
  # Force medium correlation context (alphad = 0.5) as per paper standard
  medium_corr_matrix <- rcorrmatrix(k_latent, alphad = 0.5) 
  s1_data <- generate_s1(n_samples, medium_corr_matrix, s1_prevalences)
  
  # 3. Generate Physics (Lookups)
  lookups <- create_lookup_tables(k_latent, m_pool)
  
  # 4. Generate Ground Truth (Y and S2)
  truth <- generate_truth(s1_data, lookups$y_lookup, lookups$s2_prob_matrix, lookups$noise_indices)
  
  # 5. Apply Observational Error
  #    (Used for Spectral Sim & Superloop)
  s2_true_dt <- truth$s2_true_dt
  err_data <- apply_observational_error(s2_true_dt, fid_min, fid_max)
  
  # Return everything in a neat package
  return(list(
    seed = seed,
    s1_data = s1_data,
    s1_prevalences = s1_prevalences,
    corr_matrix = medium_corr_matrix,
    lookups = lookups,
    y_true = truth$y_true,
    s2_true = s2_true_dt,
    s2_observed = err_data$s_prime_2_dt,
    avg_error_rate = err_data$avg_error_rate
  ))
}

########################################################################
########################################################################

# --- Fast Spectral Decomposition Helper ---
get_spectral_context <- function(mat) {
  # 1. Center and Scale (Correlation Matrix Strategy)
  mat_centered <- scale(mat, center = TRUE, scale = TRUE)
  mat_centered[is.na(mat_centered)] <- 0
  
  # 2. Calculate Correlation Matrix (M x M)
  c_mat <- crossprod(mat_centered) / (nrow(mat) - 1)
  
  # 3. Eigen Decomposition
  eigen_result <- eigen(c_mat, symmetric = TRUE)
  
  vals <- pmax(eigen_result$values, 1e-9)
  return(list(values = vals, vectors = eigen_result$vectors)) 
}

select_features_spectral_aligned <- function(spectral_map, raw_data_matrix, original_colnames, max_m, top_n = PROTOTYPE_TOP_N) {
  
  k_eff <- ncol(spectral_map) 
  n_features <- nrow(spectral_map)
  if (max_m > length(original_colnames)) max_m <- length(original_colnames)
  
  # --- STEP 1: PROMAX ROTATION (Get Loadings) ---
  if (k_eff > 1) {
    rot_obj <- tryCatch(promax(spectral_map), error = function(e) NULL)
    if (!is.null(rot_obj)) {
      loadings <- unclass(rot_obj$loadings)
    } else {
      loadings <- spectral_map
    }
  } else {
    loadings <- spectral_map
  }
  
  # --- STEP 2: IDENTIFY EMPIRICAL PROTOTYPES (Top N Variables per Factor) ---
  
  # Safety: Ensure we don't ask for more variables than exist
  top_n_prototype <- min(n_features, max(1, top_n))
  
  prototype_vectors <- vector("list", k_eff)
  
  # Standardize Raw Data ONCE for fair distance calculation
  raw_data_scaled <- scale(raw_data_matrix)
  raw_data_scaled[is.na(raw_data_scaled)] <- 0
  feature_profiles <- t(raw_data_scaled) 
  
  for(k in 1:k_eff) {
    ideal_loading_vec <- rep(0, k_eff); ideal_loading_vec[k] <- 1.0
    diffs <- (loadings - rep(ideal_loading_vec, each = n_features))^2
    
    # Select only the Elite Core to define the prototype
    top_indices <- order(rowSums(diffs))[1:top_n_prototype]
    prototype_vectors[[k]] <- colMeans(feature_profiles[top_indices, , drop=FALSE])
  }
  
  # --- STEP 3: RE-SORT QUEUES BASED ON RAW DATA DISTANCE ---
  factor_queues <- vector("list", k_eff)
  for(k in 1:k_eff) {
    proto_mat <- matrix(prototype_vectors[[k]], nrow=n_features, ncol=ncol(feature_profiles), byrow=TRUE)
    diffs_raw <- (feature_profiles - proto_mat)^2
    factor_queues[[k]] <- order(rowSums(diffs_raw))
  }
  
  # --- STEP 4: ADAPTIVE WEAKEST LINK SELECTION ---
  selected_indices_global <- c()
  selected_mask <- logical(n_features) 
  factor_sums <- vector("list", k_eff); for(k in 1:k_eff) factor_sums[[k]] <- rep(0, ncol(feature_profiles))
  factor_counts <- rep(0, k_eff)
  current_errors <- rep(Inf, k_eff)
  queue_pointers <- rep(1, k_eff) 
  
  # === INIT ROUND: Pick Top 5 per Factor ===
  initial_picks_per_factor <- 5
  
  for(k in 1:k_eff) {
    # Inner loop to pick multiple variables per factor
    for(p in 1:initial_picks_per_factor) {
      
      # Safety Check: Stop if we hit the total budget (max_m)
      if(length(selected_indices_global) >= max_m) break
      
      idx <- -1
      while(queue_pointers[k] <= n_features) {
        cand <- factor_queues[[k]][queue_pointers[k]]; queue_pointers[k] <- queue_pointers[k] + 1
        if(!selected_mask[cand]) { idx <- cand; break }
      }
      
      if(idx != -1) {
        selected_indices_global <- c(selected_indices_global, idx); selected_mask[idx] <- TRUE
        
        # Update Cluster Centroid Stats
        factor_sums[[k]] <- factor_sums[[k]] + feature_profiles[idx, ]
        factor_counts[k] <- factor_counts[k] + 1
        
        # Update Error immediately so the Roulette loop starts with correct math
        current_errors[k] <- sum(((factor_sums[[k]] / factor_counts[k]) - prototype_vectors[[k]])^2)
      }
    }
  }
  
  # Roulette Loop
  while(length(selected_indices_global) < max_m) {
    valid_indices <- which(current_errors != -Inf)
    if(length(valid_indices) == 0) break 
    
    valid_errs <- current_errors[valid_indices]
    
    if(sum(valid_errs) == 0) {
      probs <- rep(1/length(valid_errs), length(valid_errs))
    } else {
      flattened_errs <- valid_errs^PROB_DAMPENING_EXPONENT
      probs <- flattened_errs / sum(flattened_errs)
    }
    
    target_k <- sample(valid_indices, size = 1, prob = probs)
    
    idx <- -1
    while(queue_pointers[target_k] <= n_features) {
      cand <- factor_queues[[target_k]][queue_pointers[target_k]]; queue_pointers[target_k] <- queue_pointers[target_k] + 1
      if(!selected_mask[cand]) { idx <- cand; break }
    }
    
    if(idx != -1) {
      selected_indices_global <- c(selected_indices_global, idx); selected_mask[idx] <- TRUE
      factor_sums[[target_k]] <- factor_sums[[target_k]] + feature_profiles[idx, ]
      factor_counts[target_k] <- factor_counts[target_k] + 1
      current_errors[target_k] <- sum(((factor_sums[[target_k]] / factor_counts[target_k]) - prototype_vectors[[target_k]])^2)
    } else {
      current_errors[target_k] <- -Inf
    }
  }
  return(original_colnames[sample(selected_indices_global)])
}

# --- Helper: Outcome-Specific P-DCAI ---
select_features_targeted_residual <- function(y_train, s_prime_mat_train, all_predictor_names, max_m, num_threads = 1) {
  
  if (max_m > length(all_predictor_names)) max_m <- length(all_predictor_names)
  
  # 1. Ensure Integer Storage
  # Reduces memory bandwidth when sub-setting for Ranger
  if(!is.matrix(s_prime_mat_train)) s_prime_mat_train <- as.matrix(s_prime_mat_train)
  storage.mode(s_prime_mat_train) <- "integer"
  
  # Pre-calculate Y formats
  if(is.factor(y_train)) {
    y_num <- as.numeric(y_train) - 1
    y_factor <- y_train
  } else {
    y_num <- as.numeric(y_train)
    y_factor <- factor(y_train)
  }
  
  n_samples_train <- nrow(s_prime_mat_train)
  
  # Pre-Standardize Matrix
  # Create Z-Score matrix ONCE so correlation inside loop is simple Matrix mult.
  cat("            - Pre-calculating Z-score matrix for fast correlation...\n")
  
  X_means <- colMeans(s_prime_mat_train)
  X_sds <- apply(s_prime_mat_train, 2, sd)
  
  # Handle constants
  safe_sds <- X_sds
  safe_sds[safe_sds < 1e-9] <- 1
  
  # Create Z matrix: (X - mu) / sigma
  # This creates a numeric (double) matrix, separate from the integer source
  Z_mat <- scale(s_prime_mat_train, center = X_means, scale = safe_sds)
  # Zero out columns that were constant
  if(any(X_sds < 1e-9)) Z_mat[, which(X_sds < 1e-9)] <- 0
  
  # 2. Initialization: Select feature with highest correlation to Y
  # Cor(X,Y) proportional to crossprod(Z, Y_centered)
  y_centered <- y_num - mean(y_num)
  cor_y <- abs(crossprod(Z_mat, y_centered))
  first_best_idx <- which.max(cor_y)
  
  current_selected_indices <- c(first_best_idx)
  all_indices <- 1:ncol(s_prime_mat_train)
  available_indices <- setdiff(all_indices, current_selected_indices)
  
  cat("            - Targeted Residual Expansion: Running optimized selection loop...\n")
  
  while(length(current_selected_indices) < max_m) {
    if(length(available_indices) == 0) break
    
    # Subset Raw Integer Matrix for Ranger (Memory Efficient)
    x_train_curr <- s_prime_mat_train[, current_selected_indices, drop=FALSE]
    
    # Train RF (Multi-threaded)
    rf_fit <- tryCatch({
      ranger(x = x_train_curr
             , y = y_factor
             , probability = TRUE
             , num.trees = RF_NUM_TREES / 10
             , splitrule = "gini"
             , replace = TRUE
             , sample.fraction = 1
             , min.node.size = 5
             , max.depth = 10
             , mtry = max(1, sqrt(ncol(x_train_curr)))
             , num.threads = num_threads
             , verbose = FALSE)
    }, error = function(e) NULL)
    
    if (is.null(rf_fit)) {
      probs <- rep(0.5, length(y_num)) 
    } else {
      pred_obj <- predict(rf_fit, data = x_train_curr, num.threads = num_threads)$predictions
      probs <- if(ncol(pred_obj) >= 2) pred_obj[, 2] else rep(0.5, length(y_num))
    }
    
    # --- Residual Correlation ---
    residuals <- y_num - probs
    
    # Subset Z-Matrix (Pre-calculated doubles)
    Z_cand <- Z_mat[, available_indices, drop=FALSE]
    
    # Pure BLAS Matrix Multiplication (Vectorized)
    # No means, no SDs, no subtractions computed here.
    cor_proxies <- abs(crossprod(Z_cand, residuals))
    
    best_rel_idx <- which.max(cor_proxies)
    best_absolute_idx <- available_indices[best_rel_idx]
    
    current_selected_indices <- c(current_selected_indices, best_absolute_idx)
    available_indices <- available_indices[-best_rel_idx]
    
    if(length(current_selected_indices) %% 1 == 0) {
      cat(sprintf("\r                                   ... selected %d / %d features", length(current_selected_indices), max_m))
    }
  }
  
  cat("\n")
  return(colnames(s_prime_mat_train)[current_selected_indices])
}

# --- Helper: MRMR (DataFrame Conversion wrapper) ---
helper_build_mrmr_list <- function(y_train, s_prime_mat_train, all_predictor_names, max_m, num_threads = 1) {
  # MRMR requires a dataframe, so we convert ONLY here.
  # We select from the matrix first to save memory if possible, 
  # but here we need the whole pool to select from.
  X_df <- as.data.frame(s_prime_mat_train) 
  Y_fac <- as.factor(y_train)
  
  if (length(unique(Y_fac)) < 2) return(sample(all_predictor_names, max_m))
  
  k_val <- min(max_m, ncol(X_df))
  mrmr_out <- tryCatch({
    praznik::MRMR(X_df, Y_fac, k = k_val, threads = num_threads)
  }, error = function(e) NULL)
  
  if (is.null(mrmr_out)) return(sample(all_predictor_names, max_m))
  
  selected_names <- names(X_df)[mrmr_out$selection]
  
  # Fill if short
  if (length(selected_names) < max_m) {
    remaining <- setdiff(all_predictor_names, selected_names)
    needed <- max_m - length(selected_names)
    fill <- if(length(remaining) >= needed) sample(remaining, needed) else sample(remaining, needed, replace=TRUE)
    selected_names <- c(selected_names, fill)
  }
  return(selected_names)
}

# --- Helper: Fast Eigenvalues for Part 1b ---
calc_eigenvals_raw <- function(data_dt) {
  
  # 1. Use Full Data (No sampling)
  data_matrix <- as.matrix(data_dt)
  
  # 2. NO FILTERING - Keep low variance columns to match Superloop logic
  # col_sds <- apply(data_matrix, 2, sd)
  # non_zero_var_cols <- col_sds > 1e-6
  data_matrix_filtered <- data_matrix # Keep Full
  
  # 3. Fast Calculation
  # Correlation Matrix Strategy
  mat_centered <- scale(data_matrix_filtered, center = TRUE, scale = TRUE)
  
  # Handle NAs from zero-variance columns (scale produces NaNs)
  mat_centered[is.na(mat_centered)] <- 0
  
  # Calculate Correlation Matrix efficiently
  c_mat <- crossprod(mat_centered)
  
  # Eigen decomposition
  eigen_result <- eigen(c_mat, symmetric = TRUE, only.values = TRUE)
  
  # Normalize by (N-1)
  eigenvalues <- eigen_result$values / (nrow(data_matrix_filtered) - 1)
  
  # Clamp to epsilon
  eigenvalues <- pmax(eigenvalues, 1e-9)
  
  return(data.table(Component = 1:length(eigenvalues), Eigenvalue = eigenvalues)) 
}

create_m_schedule <- function(max_m) {
  breaks <- c(0, 20, 100, 250, 500, 1000, 2000, Inf)  
  steps <- c(5, 10, 25, 50, 100, 250, 500)        
  full_schedule <- c()
  for (i in 1:(length(breaks)-1)) {
    start <- breaks[i]; end <- min(breaks[i+1], max_m); step <- steps[i]
    current_start <- if (i == 1) step else start + step
    if (current_start <= end) {
      seq_points <- seq(from = current_start, to = end, by = step)
      full_schedule <- c(full_schedule, seq_points)
    }
  }
  
  full_schedule <- unique(c(full_schedule, 20, 100, max_m))
  full_schedule <- full_schedule[full_schedule <= max_m]
  
  return(sort(unique(full_schedule)))
}

handle_parameter_modification <- function() {
  rerun_needed <- FALSE
  repeat {
    cat("\nWhich parameter would you like to modify?\n")
    cat("1. Modify S(1) Latent Dimensionality (K_PRIMARY_LATENT)\n")
    cat("2. Modify S(1) Prevalence Range (S1_PREVALENCE_MIN/MAX)\n")
    cat("3. Done modifying.\n")
    sub_choice <- readline(prompt = "Enter your choice (1-3): ")
    sub_choice <- trimws(sub_choice)
    
    if (sub_choice == "1") {
      new_K <- as.integer(readline(prompt = paste0("Enter new K_PRIMARY_LATENT (current: ", K_PRIMARY_LATENT, "): ")))
      if (!is.na(new_K) && new_K > 0) { 
        K_PRIMARY_LATENT <<- new_K
        
        # 1. Update Pool (Haystack)
        new_pool <- min(20000, (2^new_K) * POOL_REDUNDANCY_FACTOR)
        M_POOL_SIZE <<- new_pool
        
        # 2. Update Budget (Needle Count)
        new_m <- floor(new_pool * BUDGET_PCT)
        # Ensure m >= 20 for code stability
        if(new_m < 20) new_m <- 20
        M_PREDICTORS <<- new_m
        
        cat(sprintf("   >> Auto-Updated: Pool=%d, Budget=%d (%.0f%%)\n", M_POOL_SIZE, M_PREDICTORS, BUDGET_PCT*100))
        rerun_needed <- TRUE 
      }
    } else if (sub_choice == "2") {
      new_min <- as.numeric(readline(prompt = paste0("Enter new S1_PREVALENCE_MIN (current: ", S1_PREVALENCE_MIN, "): ")))
      new_max <- as.numeric(readline(prompt = paste0("Enter new S1_PREVALENCE_MAX (current: ", S1_PREVALENCE_MAX, "): ")))
      if (!is.na(new_min) && !is.na(new_max) && new_min >= 0 && new_max <= 1 && new_min <= new_max) { 
        S1_PREVALENCE_MIN <<- new_min; S1_PREVALENCE_MAX <<- new_max; rerun_needed <- TRUE 
      }
    } else if (sub_choice == "3") { break }
  }  
  return(rerun_needed)
}

# --- Part 1b: Spectral Analysis ---
run_spectral_analysis_and_plot <- function(k_eff_from_s1, n_iterations = 20, calc_eigenvals_raw_func) {
  cat(sprintf("\n--- Running Part 1b: Spectral Analysis (Medium Corr. Context, %d Iterations) ---\n", n_iterations))
  flush.console()
  
  calc_eigenvals_raw <- calc_eigenvals_raw_func
  all_plot_data_list <- list() 
  
  est_k_obs_values <- numeric(n_iterations)
  true_k_eff_values <- numeric(n_iterations)
  err_rate_obs <- numeric(n_iterations)
  
  k_val <- K_PRIMARY_LATENT
  n_spectral <- N_SAMPLES 
  m_for_spectral <- M_POOL_SIZE 
  
  for (i in 1:n_iterations) {
    if (i %% 1 == 0) cat(sprintf("              - Spectral iteration %d of %d (Calculating noise floors)...\n", i, n_iterations))
    
    # --- USE UNIFIED GENERATOR ---
    univ <- generate_universe(seed = i, n_samples = N_SAMPLES)
    
    iter_metrics <- calculate_latent_sparsity(univ$s1_data)
    true_k_eff_values[i] <- 2^iter_metrics$entropy
    
    # Use the observed matrix directly from the unified generator
    s_prime_obs_dt <- univ$s2_observed
    err_rate_obs[i] <- univ$avg_error_rate
    
    eigen_obs_dt <- calc_eigenvals_raw(s_prime_obs_dt)
    
    if(!is.null(eigen_obs_dt)) {
      # Pass FULL matrix (no filter) to match calc_eigenvals_raw logic
      sds_obs <- apply(s_prime_obs_dt, 2, sd)
      floor_obs <- get_noise_floor(N_SAMPLES, sds_obs, n_iters = 3)
      est_k_obs_values[i] <- find_spectral_elbow(eigen_obs_dt$Eigenvalue, floor_obs)
      
      eigen_obs_dt[, DataSource := "1. S'(2) (Observed)"]
      all_plot_data_list[[length(all_plot_data_list) + 1]] <- eigen_obs_dt
    } else {
      est_k_obs_values[i] <- 1
    }
  } 
  
  full_plot_data <- rbindlist(all_plot_data_list)
  if(nrow(full_plot_data) == 0) return()
  
  agg_plot_data <- full_plot_data[, .(Mean_Eigenvalue = mean(Eigenvalue, na.rm=TRUE)), by = .(DataSource, Component)]
  
  avg_true_k_eff <- mean(true_k_eff_values)
  avg_est_k_obs <- round(mean(est_k_obs_values))
  mean_err_obs_pct <- mean(err_rate_obs) * 100
  lbl_obs <- sprintf("1. S'(2) (Observed Error: %.1f%%)", mean_err_obs_pct)
  agg_plot_data[DataSource == "1. S'(2) (Observed)", DataSource := lbl_obs]
  
  max_component_to_plot <- min(m_for_spectral, max(20, floor(3 * avg_true_k_eff)))
  agg_plot_data_subset <- agg_plot_data[Component <= max_component_to_plot]
  
  color_palette_spectral <- setNames("darkgreen", lbl_obs)
  linetype_spectral <- setNames("solid", lbl_obs)
  
  # Positioning for labels (Restored original logic)
  y_pos_base <- max(agg_plot_data_subset[DataSource == lbl_obs & Component == 1]$Mean_Eigenvalue)
  if(is.infinite(y_pos_base) || is.na(y_pos_base)) y_pos_base <- 1 
  y_pos_keff_text <- y_pos_base * 0.9
  y_pos_est_obs <- y_pos_keff_text * 0.8
  
  display_mode <- if(exists("S2_GENERATION_MODE")) S2_GENERATION_MODE else "CONSISTENT"
  
  p <- ggplot(agg_plot_data_subset, aes(x = Component, y = Mean_Eigenvalue, color = DataSource, linetype = DataSource)) +
    geom_line(linewidth = 1.0) +
    
    # K_eff Line and Label (Restored)
    geom_vline(xintercept = avg_true_k_eff, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = avg_true_k_eff * 1.2, y = y_pos_keff_text, 
             label = sprintf("Effective Latent Complexity (K_eff): %.2f", avg_true_k_eff), 
             color = "black", size = 4, hjust = 0) +
    
    # Observed Estimate Line and Label (Restored)
    geom_vline(xintercept = avg_est_k_obs, linetype = "dotted", color = "darkgreen", linewidth = 1) +
    annotate("text", x = avg_est_k_obs * 1.2, y = y_pos_est_obs, 
             label = sprintf("Effective Signal Rank (r): %d", avg_est_k_obs), 
             color = "darkgreen", size = 4, hjust = 0) +
    
    scale_y_continuous(trans = "log10", labels = scales::number_format(accuracy = 0.001)) + 
    scale_linetype_manual(values = linetype_spectral, name = "Data Source") +
    scale_color_manual(values = color_palette_spectral, name = "Data Source") +
    coord_cartesian(xlim = c(1, max_component_to_plot)) +
    labs(title = paste("S'(2) Spectral Analysis -- Effective Latent Complexity vs Effective Signal Rank"), 
         subtitle = sprintf("Generative Mode: %s -- Results based on %d iteration(s).", display_mode, n_iterations), 
         x = sprintf("Principal Component (First %d)", max_component_to_plot), 
         y = "Average Eigenvalue (log scale)", 
         caption = sprintf("Parameters: k = %d, N = %s, m = %d, Observational Fidelity = [%.3f, %.3f]", 
                           k_val, format(N_SAMPLES, big.mark=","), m_for_spectral, OBS_FIDELITY_MIN, OBS_FIDELITY_MAX)) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right", plot.caption = element_text(hjust = 0, size = 10))
  
  ggsave("spectral_analysis_scree_plot_R.png", p, width = 12, height = 8, dpi = 300)
  print(p)
  
  return(list(k_eff_low = avg_est_k_obs, k_eff_high = avg_est_k_obs))
}

# --- 7. Part 3: Dynamic Strategy Spectral Analysis (Multi-Iteration Average) ---
run_strategy_spectral_comparison <- function(fixed_k_target) {
  # We use the global store instead of the last single artifact
  if (!exists("global_artifact_store") || length(global_artifact_store) == 0) {
    if (length(last_run_artifacts) > 0) {
      # Fallback: If user only ran 1 iteration or didn't save all, use the last one
      cat("Note: Using single last run artifact (Global store empty).\n")
      artifacts_to_process <- list(last_run_artifacts)
    } else {
      cat("Error: No simulation artifacts found.\n")
      return()
    }
  } else {
    artifacts_to_process <- global_artifact_store
  }
  
  n_iters <- length(artifacts_to_process)
  cat(sprintf("\n--- Running Part 3: Strategy Spectral Analysis (Averaging %d Iterations) ---\n", n_iters))
  cat(sprintf("    Benchmarking all strategies at Fixed K = %d\n", fixed_k_target))
  
  if (!requireNamespace("gridExtra", quietly = TRUE)) { stop("Package 'gridExtra' needed.") }
  library(gridExtra); library(grid); library(ggplot2); library(data.table)
  
  strategies <- c("Random Breadth", "MRMR", "Outcome-Specific P-DCAI", "Spectral P-DCAI")
  all_iter_results <- list()
  
  # --- LOOP THROUGH EACH SAVED ITERATION ---
  for (i in 1:n_iters) {
    run_obj <- artifacts_to_process[[i]]
    iter_seed <- run_obj$seed
    
    cat(sprintf("\r    Reconstructing Dataset for Iteration %d (Seed %d)...", i, iter_seed))
    
    # --- USE UNIFIED GENERATOR ---
    # Need full N for training + test reconstruction
    # (Though we only strictly need training for spectral, we match the main loop logic)
    n_train_target <- N_SAMPLES
    n_test_target <- ceiling(n_train_target * TEST_SET_MULTIPLIER)
    n_total_needed <- n_train_target + n_test_target
    
    univ <- generate_universe(seed = iter_seed, n_samples = n_total_needed)
    
    # We only need the Training set for spectral analysis (matching the feature selection step)
    train_indices <- sample(1:n_total_needed, n_train_target)
    
    # Get Observed Error Data
    X_train_full <- as.matrix(univ$s2_observed)[train_indices, ]
    
    # 3. CALCULATE SPECTRA FOR THIS ITERATION
    target_k <- max(1, fixed_k_target)
    
    for(strat in strategies) {
      selected <- NULL
      # Retrieve the specific list saved for THIS iteration
      if(strat == "Random Breadth") selected <- run_obj$random_features
      else if(strat == "Outcome-Specific P-DCAI") selected <- run_obj$deliberate_lists[["Observed"]]
      else if(strat == "Spectral P-DCAI") selected <- run_obj$spectral_lists[["Observed"]]
      else if(strat == "MRMR") selected <- run_obj$mrmr_lists[["Observed"]]
      
      if (!is.null(selected)) {
        # Limit to M predictors
        selected <- selected[1:min(length(selected), M_PREDICTORS)]
        
        # Validation
        available_cols <- intersect(selected, colnames(X_train_full))
        if(length(available_cols) < 2) next
        X_sub <- X_train_full[, available_cols, drop=FALSE]
        
        # NO FILTERING (To match request)
        # col_sds <- apply(X_sub, 2, sd)
        # if(sum(col_sds > 1e-9) < 2) next 
        # X_sub <- X_sub[, col_sds > 1e-9, drop=FALSE]
        
        # Eigenvalues
        X_cent <- scale(X_sub, center=TRUE, scale=TRUE)
        X_cent[is.na(X_cent)] <- 0
        ev <- eigen(crossprod(X_cent), symmetric=TRUE, only.values=TRUE)$values
        ev <- ev / (nrow(X_sub) - 1)
        ev <- pmax(ev, 1e-9)
        
        # Calculate SNR for this specific run
        val_sig <- ev[target_k]
        val_noise <- if((target_k + 1) <= length(ev)) ev[target_k + 1] else 1e-9
        spec_snr <- val_sig / val_noise
        
        energy_sig <- sum(ev[1:target_k])
        energy_noise <- sum(ev[(target_k+1):length(ev)])
        overall_snr <- energy_sig / energy_noise
        
        dt_curve <- data.table(
          Iteration = i,
          Strategy = strat,
          Component = 1:length(ev), 
          Value = ev,
          Spectral_SNR = spec_snr,
          Overall_SNR = overall_snr
        )
        all_iter_results[[length(all_iter_results)+1]] <- dt_curve
      }
    }
  }
  cat("\n")
  
  if(length(all_iter_results) == 0) return()
  
  # --- 4. AGGREGATION & PLOTTING ---
  full_data <- rbindlist(all_iter_results)
  
  # Average Scree Curve
  agg_scree <- full_data[, .(Value = mean(Value)), by = .(Strategy, Component)]
  
  # Average Metrics
  agg_stats <- full_data[, .(Spectral_SNR = mean(unique(Spectral_SNR)), 
                             Overall_SNR = mean(unique(Overall_SNR))), by = .(Strategy)]
  
  # Formatting
  internal_order <- c("Random Breadth", "MRMR", "Outcome-Specific P-DCAI", "Spectral P-DCAI")
  display_labels <- c("Random", "MRMR (Supervised Traditional)", 
                      "Targeted Residual Expansion (Supervised P-DCAI)", 
                      "Spectral (Unsupervised P-DCAI)")
  
  agg_scree <- agg_scree[Strategy %in% internal_order]
  agg_scree[, Strategy := factor(Strategy, levels = internal_order, labels = display_labels)]
  agg_stats <- agg_stats[Strategy %in% internal_order]
  agg_stats[, Strategy := factor(Strategy, levels = internal_order, labels = display_labels)]
  
  strat_colors <- c(
    "Random" = "gray50",
    "MRMR (Supervised Traditional)" = "orange",
    "Targeted Residual Expansion (Supervised P-DCAI)" = "blue",
    "Spectral (Unsupervised P-DCAI)" = "purple"
  )
  
  labels_map <- setNames(
    sprintf("%s\n(Avg SNR = %.2f)", agg_stats$Strategy, agg_stats$Spectral_SNR),
    agg_stats$Strategy
  )
  y_anno <- max(agg_scree$Value) * 0.5
  
  p <- ggplot(agg_scree, aes(x = Component, y = Value, color = Strategy)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = fixed_k_target, linetype = "longdash", color = "black", linewidth = 1) +
    annotate("text", x = fixed_k_target + 0.5, y = y_anno, 
             label = sprintf("Benchmark Elbow (k = %d)", fixed_k_target), 
             color="black", vjust = 0, hjust = 1, size = 4) +
    scale_y_continuous(trans="log10") +
    coord_cartesian(xlim = c(1, fixed_k_target * 2)) + 
    scale_color_manual(values = strat_colors, labels = labels_map) +
    labs(title = paste("Spectral Analysis by Feature Selection Strategy"),
         subtitle = sprintf("Generative Mode: %s -- Avg over %d iteration(s). Analyzing signal strength at Global Haystack Elbow (k = %d).", S2_GENERATION_MODE, n_iters, fixed_k_target),
         y = "Average Eigenvalue (Log Scale)", x = "Principal Component Index") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right", legend.direction = "vertical", legend.text = element_text(size=10),
          plot.subtitle = element_text(size=10, color="gray30"))
  
  ggsave("comparative_snr_analysis_R.png", p, width = 12, height = 7, dpi = 300)
  cat("Generated 'comparative_snr_analysis_R.png' (Averaged Fixed K comparison).\n")
  print(p) 
}

run_one_full_simulation_iteration <- function(iteration_seed, m_schedule, k_eff_latent, inner_threads = 1, save_artifacts = FALSE) {
  # Note: Seed is set INSIDE generate_universe
  cat(sprintf("\n--- Super Loop %d ---\n", iteration_seed))
  flush.console()
  
  # --- 1. Generation (Unified) ---
  n_train_target <- N_SAMPLES
  n_test_target <- ceiling(n_train_target * TEST_SET_MULTIPLIER)
  n_total_needed <- n_train_target + n_test_target
  
  # CALL UNIFIED GENERATOR
  univ <- generate_universe(seed = iteration_seed, n_samples = n_total_needed)
  
  # UNPACK
  s1_medium_corr <- univ$s1_data
  truth_mat <- as.matrix(univ$s2_true)
  obs_mat <- as.matrix(univ$s2_observed)
  
  # Calculate Sparsity Metric
  iter_sparsity_metrics <- calculate_latent_sparsity(s1_medium_corr)
  true_k_eff_iter <- 2^(iter_sparsity_metrics$entropy)
  k_eff_idealized <- max(1, ceiling(true_k_eff_iter))
  
  cat(sprintf("        - K_eff for this Iteration = %.2f (used for Idealized Depth)\n", true_k_eff_iter))
  
  # Split Data
  train_indices <- sample(1:n_total_needed, n_train_target)
  test_indices <- setdiff(1:n_total_needed, train_indices)
  
  y_train_factor <- factor(univ$y_true[train_indices])
  y_train <- univ$y_true[train_indices] # Numeric needed for some helpers
  y_test_truth <- univ$y_true[test_indices]  
  y_test_factor <- factor(y_test_truth)
  
  predictor_names <- colnames(truth_mat)
  shuffled_predictor_names <- sample(predictor_names)
  
  # --- 3. Prepare Error Regimes (Matrices) ---
  error_mats <- list()
  phi_metrics <- list()
  scenarios_to_run <- c("Observed") 
  if (PLOT_NO_ERROR_LINE) scenarios_to_run <- c("No Error", scenarios_to_run)
  
  all_corrs <- cor(univ$y_true, truth_mat)
  
  if("No Error" %in% scenarios_to_run) {
    error_mats[["No Error"]] <- truth_mat
    phi_metrics[["No Error"]] <- list(err = 0, avg_phi = mean(abs(all_corrs)))
  }
  
  # --- DYNAMIC SPECTRAL ANALYSIS (Optimized) ---
  k_eff_empirical_list <- list()
  spectral_contexts <- list() 
  
  # Single Observed Regime
  err_level <- "Observed"
  mat_err <- obs_mat # Already generated
  error_mats[[err_level]] <- mat_err
  phi_metrics[[err_level]] <- list(
    err = univ$avg_error_rate, 
    avg_phi = mean(abs(cor(univ$y_true, mat_err)[1, ]))
  )
  
  mat_err_train <- mat_err[train_indices, ] 
  
  # 1. RUN DECOMPOSITION
  # Use Standard get_spectral_context (Unfiltered, as requested)
  decomp <- tryCatch({ get_spectral_context(mat_err_train) }, error=function(e) NULL)
  
  if(!is.null(decomp)) {
    # 2. Find Elbow
    # Important: Since we are NOT filtering, we calculate floor on full matrix too
    current_sds <- apply(mat_err_train, 2, sd)
    current_floor <- get_noise_floor(nrow(mat_err_train), current_sds, n_iters = 3)
    k_emp <- find_spectral_elbow(decomp$values, current_floor)
    
    # 3. SAVE CLEAN CONTEXT
    valid_k <- min(k_emp, ncol(decomp$vectors))
    spectral_contexts[[err_level]] <- list(
      k_emp = k_emp,
      clean_coords = decomp$vectors[, 1:valid_k, drop=FALSE]
    )
  } else {
    k_emp <- k_eff_idealized
    spectral_contexts[[err_level]] <- NULL
  }
  
  k_eff_empirical_list[[err_level]] <- k_emp
  cat(sprintf("        - Effective Signal Rank for %s Error: %d\n", err_level, k_emp))
  
  # --- 2. Calculate BENCHMARKS ---
  cat("        - Calculating Static Benchmarks (Oracle/Depth)...\n")
  H_Y_uncond <- infotheo::entropy(y_test_factor)
  
  s1_int <- as.vector(s1_to_int(s1_medium_corr))
  oracle_scores_all <- univ$lookups$y_lookup[data.table(config_id=s1_int), on="config_id"]$prob
  oracle_scores_test <- oracle_scores_all[test_indices]
  
  oracle_metrics <- list(
    auroc = as.numeric(roc(y_test_truth, oracle_scores_test, quiet=TRUE)$auc),
    auprc = pr.curve(scores.class0=oracle_scores_test[y_test_truth==1], scores.class1=oracle_scores_test[y_test_truth==0], curve=FALSE)$auc.integral,
    cond_ent = condentropy(y_test_factor, infotheo::discretize(oracle_scores_test))
  )
  
  IDEALIZED_DEPTH_COUNT <- min(k_eff_idealized, 2*M_PREDICTORS)
  X_depth_train_df <- as.data.frame(truth_mat[train_indices, ])
  mrmr_depth <- praznik::MRMR(X_depth_train_df, as.factor(y_train), k = IDEALIZED_DEPTH_COUNT, threads = inner_threads)
  selected_depth_names <- names(X_depth_train_df)[mrmr_depth$selection]
  avg_phi_depth <- mean(abs(all_corrs[1, selected_depth_names]))
  
  s2_depth_train <- truth_mat[train_indices, selected_depth_names, drop=FALSE]
  s2_depth_test <- truth_mat[test_indices, selected_depth_names, drop=FALSE]
  
  model_depth <- ranger(x = s2_depth_train
                        , y = y_train_factor
                        , probability = TRUE
                        , num.trees = RF_NUM_TREES
                        , splitrule = "gini"
                        , replace = TRUE
                        , sample.fraction = 1
                        , min.node.size = 1
                        , max.depth = 0
                        , mtry = max(1, sqrt(ncol(s2_depth_train)))
                        , num.threads = inner_threads
                        , verbose = FALSE)
  
  depth_scores <- predict(model_depth, data = s2_depth_test, num.threads = inner_threads)$predictions[, 2]
  
  depth_metrics <- list(
    auroc = as.numeric(roc(y_test_truth, depth_scores, quiet=TRUE)$auc),
    auprc = pr.curve(scores.class0=depth_scores[y_test_truth==1], scores.class1=depth_scores[y_test_truth==0], curve=FALSE)$auc.integral,
    cond_ent = condentropy(y_test_factor, infotheo::discretize(depth_scores))
  )
  rm(X_depth_train_df, s2_depth_train, s2_depth_test, model_depth); gc()
  
  # --- 4. Strategy Feature Selection ---
  deliberate_lists <- list()
  mrmr_lists <- list()
  spectral_lists <- list()
  
  run_g2g <- STRATEGY_CHOICE %in% c("Outcome-Specific P-DCAI", "All")
  run_unsup_spectral <- STRATEGY_CHOICE %in% c("Spectral P-DCAI", "All")
  run_mrmr <- STRATEGY_CHOICE %in% c("MRMR", "All")
  
  for(lvl in scenarios_to_run) {
    if(run_mrmr) {
      cat(sprintf("        - Building Supervised MRMR list for %s Error...\n", lvl))
      mrmr_lists[[lvl]] <- helper_build_mrmr_list(y_train, error_mats[[lvl]][train_indices, ], predictor_names, M_PREDICTORS, num_threads = inner_threads)
    }
    if(run_g2g) {
      cat(sprintf("        - Building Supervised P-DCAI list for %s Error...\n", lvl))
      deliberate_lists[[lvl]] <- select_features_targeted_residual(y_train, error_mats[[lvl]][train_indices, ], predictor_names, M_PREDICTORS, num_threads = inner_threads)
    }
    
    # SPECTRAL
    if(run_unsup_spectral) {
      cat(sprintf("        - Building Unsupervised P-DCAI list Spectral for %s Error...\n", lvl))
      context_obj <- spectral_contexts[[lvl]]
      if(!is.null(context_obj)) {
      dynamic_prototype_top_n <- max(1, round(M_PREDICTORS / context_obj$k_emp))
      cat(sprintf("          > Dynamic Elite Core Size (m/k_emp): %d\n", dynamic_prototype_top_n))
      spectral_lists[[lvl]] <- select_features_spectral_aligned(
      spectral_map = context_obj$clean_coords,
      raw_data_matrix = error_mats[[lvl]][train_indices, ], 
      original_colnames = predictor_names,
      max_m = M_PREDICTORS,
      top_n = dynamic_prototype_top_n
        )
      } else {
        spectral_lists[[lvl]] <- sample(predictor_names, M_PREDICTORS)
      }
    }
  }
  
  # --- SAVE ARTIFACTS ---
  if(save_artifacts) {
    cat("\n        >>> Saving feature lists and seed for Comparative Spectral Analysis...\n")
    last_run_artifacts <<- list(
      seed = iteration_seed,
      deliberate_lists = if(length(deliberate_lists) > 0) deliberate_lists else NULL,
      spectral_lists = if(length(spectral_lists) > 0) spectral_lists else NULL,
      mrmr_lists = if(length(mrmr_lists) > 0) mrmr_lists else NULL,
      random_features = shuffled_predictor_names,
      true_k_eff = true_k_eff_iter,
      k_eff_empirical = k_eff_empirical_list
    )
  }
  
  # --- 5. The "M" Schedule Loop ---
  results_for_this_iter <- list()
  strategies <- c()
  if(STRATEGY_CHOICE %in% c("Random Breadth", "All")) strategies <- c(strategies, "Random Breadth")
  if(STRATEGY_CHOICE %in% c("Outcome-Specific P-DCAI", "All")) strategies <- c(strategies, "Outcome-Specific P-DCAI")
  if(STRATEGY_CHOICE %in% c("Spectral P-DCAI", "All")) strategies <- c(strategies, "Spectral P-DCAI")
  if(STRATEGY_CHOICE %in% c("MRMR", "All")) strategies <- c(strategies, "MRMR")
  
  if(length(strategies) == 0) cat("\nWARNING: No strategies selected! Check STRATEGY_CHOICE.\n")
  
  total_models <- length(m_schedule) * length(scenarios_to_run) * length(strategies)
  model_counter <- 0
  
  for(m_val in m_schedule) {
    current_mtry <- max(1, sqrt(m_val))
    for(error_level in scenarios_to_run) {
      X_full_train <- error_mats[[error_level]][train_indices, , drop=FALSE]
      X_full_test <- error_mats[[error_level]][test_indices, , drop=FALSE]
      
      for(strat in strategies) {
        model_counter <- model_counter + 1
        cat(sprintf("\r        - Model %d/%d [m=%d, %s, %s]          ", 
                    model_counter, total_models, m_val, error_level, strat))
        flush.console() # Prevent freezing
        
        if(strat == "Random Breadth") feats <- shuffled_predictor_names[1:m_val]
        else if(strat == "Outcome-Specific P-DCAI") feats <- deliberate_lists[[error_level]][1:m_val]
        else if(strat == "Spectral P-DCAI") feats <- spectral_lists[[error_level]][1:m_val]
        else if(strat == "MRMR") feats <- mrmr_lists[[error_level]][1:m_val]
        
        X_sub_train <- X_full_train[, feats, drop=FALSE]
        X_sub_test <- X_full_test[, feats, drop=FALSE]
        
        current_feat_corrs <- cor(as.numeric(y_train), X_sub_train)
        strategy_phi <- mean(abs(current_feat_corrs), na.rm = TRUE)
        
        model_breadth <- ranger(x = X_sub_train
                                , y = y_train_factor
                                , probability = TRUE
                                , num.trees = RF_NUM_TREES
                                , splitrule = "gini"
                                , replace = TRUE
                                , sample.fraction = 1
                                , min.node.size = 1
                                , max.depth = 0
                                , mtry = current_mtry
                                , num.threads = inner_threads
                                , verbose = FALSE)
        
        scores <- predict(model_breadth, data = X_sub_test, num.threads = inner_threads)$predictions[, 2]
        
        res_list <- list(
          iteration = iteration_seed,
          m = m_val,
          error_level = error_level,
          strategy_type = strat,
          auroc = as.numeric(roc(y_test_truth, scores, quiet=TRUE)$auc),
          auprc = pr.curve(scores.class0=scores[y_test_truth==1], scores.class1=scores[y_test_truth==0], curve=FALSE)$auc.integral,
          breadth_cond_entropy = condentropy(y_test_factor, infotheo::discretize(scores)),
          oracle_auroc = oracle_metrics$auroc,
          oracle_auprc = oracle_metrics$auprc,
          oracle_cond_entropy = oracle_metrics$cond_ent,
          depth_auroc = depth_metrics$auroc,
          depth_auprc = depth_metrics$auprc,
          depth_cond_entropy = depth_metrics$cond_ent,
          H_Y_unconditional = H_Y_uncond,
          avg_phi_depth = avg_phi_depth,
          avg_error_rate = phi_metrics[[error_level]]$err,
          avg_phi = phi_metrics[[error_level]]$avg_phi,
          avg_phi_strategy = strategy_phi
        )
        results_for_this_iter[[length(results_for_this_iter) + 1]] <- res_list
      }
    }
  }
  cat("\n")
  cat(sprintf("\n--- Super Loop %d Complete ---\n", iteration_seed))
  return(rbindlist(results_for_this_iter))
}

plot_auroc_results <- function(agg_results, current_iter, total_iters, sparsity_results) {
  plot_data <- agg_results[, .(m, error_level, strategy_type, mean_metric = mean_auroc,  
                               mean_avg_error_rate = mean_avg_error_rate, mean_avg_phi = mean_avg_phi)]
  
  # --- CUSTOM LEGEND ORDER & NAMING ---
  # 1. Define Internal Names (Must match simulation strings exactly)
  internal_order <- c("Random Breadth", "MRMR", "Outcome-Specific P-DCAI", "Spectral P-DCAI")
  
  # 2. Define Display Names (Your desired labels)
  display_labels <- c("Random", "MRMR (Supervised Traditional)", 
                      "Targeted Residual Expansion (Supervised P-DCAI)", 
                      "Spectral (Unsupervised P-DCAI)")
  
  # 3. Apply Factor with Levels/Labels
  plot_data <- plot_data[strategy_type %in% internal_order] # Filter ensures safety
  plot_data[, strategy_type := factor(strategy_type, levels = internal_order, labels = display_labels)]
  
  # 4. Define Colors Mapped to NEW Display Names
  strategy_colors <- c(
    "Random" = "gray70",
    "MRMR (Supervised Traditional)" = "orange",
    "Targeted Residual Expansion (Supervised P-DCAI)" = "blue",
    "Spectral (Unsupervised P-DCAI)" = "purple"
  )
  
  # 5. Define Widths (Optional: Make Random thinner if desired)
  strategy_widths <- c(
    "Random" = 0.8,
    "MRMR (Supervised Traditional)" = 0.8,
    "Targeted Residual Expansion (Supervised P-DCAI)" = 0.8,
    "Spectral (Unsupervised P-DCAI)" = 0.8
  )
  
  # --- Rest of standard logic ---
  all_levels <- c("No Error", "Observed") 
  present_levels <- intersect(all_levels, unique(plot_data$error_level))
  plot_data[, error_level := factor(error_level, levels = present_levels)]
  
  legend_labels_text <- sapply(unique(plot_data$error_level), function(l) {
    row <- plot_data[error_level == l][1]
    if (l == "No Error") sprintf("No Error (True S(2), Avg. Phi = %.3f)", row$mean_avg_phi)
    else sprintf("%s Error (Err: %.1f%%)", l, row$mean_avg_error_rate * 100)
  })
  
  oracle_label_text <- paste("Theoretical Ceiling (S(1) Oracle) =", round(agg_results$mean_oracle_auroc[1], 3))
  k_eff_val <- max(1, ceiling(sparsity_results[Scenario == "Medium", Mean_K_eff]))
  depth_label_text <- sprintf("Idealized Depth (%d vars) = %.3f", k_eff_val * 2, agg_results$mean_depth_auroc[1])
  
  # --- Y-AXIS MODIFICATION START ---
  # Range Min is set to 0.050 below the Idealized Depth Benchmark
  # We clamp at 0.5 to prevent invalid AUROC ranges if depth is very low
  depth_val <- agg_results$mean_depth_auroc[1]
  min_y <- max(0.5, depth_val - 0.050)
  max_y <- agg_results$mean_oracle_auroc[1] + 0.015
  
  p <- ggplot(plot_data, aes(x = m, y = mean_metric, color = strategy_type, linetype = error_level, linewidth = strategy_type)) +
    geom_hline(yintercept = agg_results$mean_oracle_auroc[1], linetype = "solid", color = "black") +  
    geom_hline(yintercept = agg_results$mean_depth_auroc[1], linetype = "solid", color = "black") +
    geom_line() +
    
    scale_color_manual(values = strategy_colors, name = "Strategy") +  
    scale_linewidth_manual(values = strategy_widths, name = "Strategy") +
    scale_linetype_manual(values = c("solid", "longdash"), labels = legend_labels_text, name = "Error Level") +
    
    annotate("text", x = min(plot_data$m), y = agg_results$mean_oracle_auroc[1], label = oracle_label_text, vjust = -0.5, hjust = 0) +
    annotate("text", x = max(plot_data$m), y = agg_results$mean_depth_auroc[1], label = depth_label_text, vjust = 1.5, hjust = 1.0, color="black", size=3.5) +
    
    labs(title = "Breadth Strategy Performance vs. Benchmarks (AUROC)", 
         subtitle = sprintf("Generative Mode: %s -- Avg over %d iteration(s). (Higher is better)", S2_GENERATION_MODE, current_iter), 
         y = "Area Under ROC", x = ifelse(USE_LOG_SCALE_X, "Number of Predictors (m, log)", "Number of Predictors (m)")) +
    theme_minimal(base_size = 14) + theme(legend.position = "right", legend.box = "vertical")
  
  if (USE_LOG_SCALE_X) p <- p + scale_x_log10()
  p <- p + coord_cartesian(ylim = c(min_y, max_y))
  
  ggsave("auroc_simulation_R.png", p, width = 14, height = 8, dpi = 300)
  print(p)
}

plot_auprc_results <- function(agg_results, current_iter, total_iters, sparsity_results) {
  plot_data <- agg_results[, .(m, error_level, strategy_type, mean_metric = mean_auprc,  
                               mean_avg_error_rate = mean_avg_error_rate, mean_avg_phi = mean_avg_phi)]
  
  # --- CUSTOM LEGEND ORDER & NAMING ---
  internal_order <- c("Random Breadth", "MRMR", "Outcome-Specific P-DCAI", "Spectral P-DCAI")
  display_labels <- c("Random", "MRMR (Supervised Traditional)", 
                      "Targeted Residual Expansion (Supervised P-DCAI)", 
                      "Spectral (Unsupervised P-DCAI)")
  
  plot_data <- plot_data[strategy_type %in% internal_order]
  plot_data[, strategy_type := factor(strategy_type, levels = internal_order, labels = display_labels)]
  
  strategy_colors <- c(
    "Random" = "gray70",
    "MRMR (Supervised Traditional)" = "orange",
    "Targeted Residual Expansion (Supervised P-DCAI)" = "blue",
    "Spectral (Unsupervised P-DCAI)" = "purple"
  )
  
  strategy_widths <- c("Random"=0.8, "MRMR (Supervised Traditional)"=0.8, "Targeted Residual Expansion (Supervised P-DCAI)"=0.8, "Spectral (Unsupervised P-DCAI)"=0.8)
  
  # --- Standard Logic ---
  all_levels <- c("No Error", "Observed") 
  present_levels <- intersect(all_levels, unique(plot_data$error_level))
  plot_data[, error_level := factor(error_level, levels = present_levels)]
  
  legend_labels_text <- sapply(unique(plot_data$error_level), function(l) {
    row <- plot_data[error_level == l][1]
    if (l == "No Error") sprintf("No Error (True S(2), Avg. Phi = %.3f)", row$mean_avg_phi)
    else sprintf("%s Error (Err: %.1f%%)", l, row$mean_avg_error_rate * 100)
  })
  
  oracle_label_text <- paste("Theoretical Ceiling (S(1) Oracle) =", round(agg_results$mean_oracle_auprc[1], 3))
  k_eff_val <- max(1, ceiling(sparsity_results[Scenario == "Medium", Mean_K_eff]))
  depth_label_text <- sprintf("Idealized Depth (%d vars) = %.3f", k_eff_val * 2, agg_results$mean_depth_auprc[1])
  
  # --- Y-AXIS MODIFICATION START ---
  # Range Min is set to 0.050 below the Idealized Depth Benchmark
  depth_val <- agg_results$mean_depth_auprc[1]
  min_y <- max(0, depth_val - 0.050)
  max_y <- min(1, agg_results$mean_oracle_auprc[1] + 0.015)
  
  p <- ggplot(plot_data, aes(x = m, y = mean_metric, color = strategy_type, linetype = error_level, linewidth = strategy_type)) +
    geom_hline(yintercept = agg_results$mean_oracle_auprc[1], linetype = "solid", color = "black") +  
    geom_hline(yintercept = agg_results$mean_depth_auprc[1], linetype = "solid", color = "black") +
    geom_line() +
    
    scale_color_manual(values = strategy_colors, name = "Strategy") +  
    scale_linewidth_manual(values = strategy_widths, name = "Strategy") +
    scale_linetype_manual(values = c("solid", "longdash"), labels = legend_labels_text, name = "Error Level") +
    
    annotate("text", x = min(plot_data$m), y = agg_results$mean_oracle_auprc[1], label = oracle_label_text, vjust = -0.5, hjust = 0) +
    annotate("text", x = max(plot_data$m), y = agg_results$mean_depth_auprc[1], label = depth_label_text, vjust = 1.5, hjust = 1.0, color="black", size=3.5) +
    
    labs(title = "Breadth Strategy Performance vs. Benchmarks (AUPRC)", 
         subtitle = sprintf("Generative Mode: %s -- Avg over %d iteration(s). (Higher is better)", S2_GENERATION_MODE, current_iter), 
         y = "Area Under PRC", x = ifelse(USE_LOG_SCALE_X, "Number of Predictors (m, log)", "Number of Predictors (m)")) +
    theme_minimal(base_size = 14) + theme(legend.position = "right", legend.box = "vertical")
  
  if (USE_LOG_SCALE_X) p <- p + scale_x_log10()
  p <- p + coord_cartesian(ylim = c(min_y, max_y))
  
  ggsave("auprc_simulation_R.png", p, width = 14, height = 8, dpi = 300)
  print(p)
}

plot_cond_entropy_results <- function(agg_results, current_iter, total_iters, sparsity_results) {
  plot_data <- agg_results[, .(m, error_level, strategy_type, mean_metric = mean_breadth_cond_entropy,  
                               mean_avg_error_rate = mean_avg_error_rate, mean_avg_phi = mean_avg_phi)]
  
  # --- CUSTOM LEGEND ORDER & NAMING ---
  internal_order <- c("Random Breadth", "MRMR", "Outcome-Specific P-DCAI", "Spectral P-DCAI")
  display_labels <- c("Random", "MRMR (Supervised Traditional)", 
                      "Targeted Residual Expansion (Supervised P-DCAI)", 
                      "Spectral (Unsupervised P-DCAI)")
  
  plot_data <- plot_data[strategy_type %in% internal_order]
  plot_data[, strategy_type := factor(strategy_type, levels = internal_order, labels = display_labels)]
  
  strategy_colors <- c(
    "Random" = "gray70",
    "MRMR (Supervised Traditional)" = "orange",
    "Targeted Residual Expansion (Supervised P-DCAI)" = "blue",
    "Spectral (Unsupervised P-DCAI)" = "purple"
  )
  
  strategy_widths <- c("Random"=0.8, "MRMR (Supervised Traditional)"=0.8, "Targeted Residual Expansion (Supervised P-DCAI)"=0.8, "Spectral (Unsupervised P-DCAI)"=0.8)
  
  # --- Standard Logic ---
  all_levels <- c("No Error", "Observed") 
  present_levels <- intersect(all_levels, unique(plot_data$error_level))
  plot_data[, error_level := factor(error_level, levels = present_levels)]
  
  legend_labels_text <- sapply(unique(plot_data$error_level), function(l) {
    row <- plot_data[error_level == l][1]
    if (l == "No Error") sprintf("No Error (True S(2), Avg. Phi = %.3f)", row$mean_avg_phi)
    else sprintf("%s Error (Err: %.1f%%)", l, row$mean_avg_error_rate * 100)
  })
  
  oracle_label_text <- paste("Theoretical Floor (S(1) Oracle) =", round(agg_results$mean_oracle_cond_entropy[1], 3))
  k_eff_val <- max(1, ceiling(sparsity_results[Scenario == "Medium", Mean_K_eff]))
  depth_label_text <- sprintf("Idealized Depth (%d vars) = %.3f", k_eff_val * 2, agg_results$mean_depth_cond_entropy[1])
  uncond_label <- paste("Max = H(Y) =", round(agg_results$mean_H_Y_unconditional[1], 3))
  
  # Scale: Theoretical Minimum (Oracle) to Theoretical Maximum (Unconditional H(Y))  
  min_y <- agg_results$mean_oracle_cond_entropy[1]
  depth_val <- agg_results$mean_depth_cond_entropy[1]
  max_y <- agg_results$mean_H_Y_unconditional[1] + 0.01
  
  p <- ggplot(plot_data, aes(x = m, y = mean_metric, color = strategy_type, linetype = error_level, linewidth = strategy_type)) +
    geom_hline(yintercept = agg_results$mean_H_Y_unconditional[1], linetype = "solid", color = "black") +  
    geom_hline(yintercept = agg_results$mean_oracle_cond_entropy[1], linetype = "solid", color = "black") +  
    geom_hline(yintercept = agg_results$mean_depth_cond_entropy[1], linetype = "solid", color = "black") +
    geom_line() +
    
    scale_color_manual(values = strategy_colors, name = "Strategy") +  
    scale_linewidth_manual(values = strategy_widths, name = "Strategy") +
    scale_linetype_manual(values = c("solid", "longdash"), labels = legend_labels_text, name = "Error Level") +
    
    annotate("text", x = min(plot_data$m), y = agg_results$mean_H_Y_unconditional[1], label = uncond_label, vjust = -0.5, hjust = 0, size=3.5) +
    annotate("text", x = min(plot_data$m), y = agg_results$mean_oracle_cond_entropy[1], label = oracle_label_text, vjust = 1.5, hjust = 0, size=3.5) +
    annotate("text", x = max(plot_data$m), y = agg_results$mean_depth_cond_entropy[1], label = depth_label_text, vjust = -0.5, hjust = 1.0, color="black", size=3.5) +
    
    labs(title = "Breadth Strategy Performance vs. Benchmarks (Conditional Entropy)", 
         subtitle = sprintf("Generative Mode: %s -- Avg over %d iteration(s). (Lower is better)", S2_GENERATION_MODE, current_iter), 
         y = "Conditional Entropy (nats)", x = ifelse(USE_LOG_SCALE_X, "Number of Predictors (m, log)", "Number of Predictors (m)")) +
    theme_minimal(base_size = 14) + theme(legend.position = "right", legend.box = "vertical")
  
  if (USE_LOG_SCALE_X) p <- p + scale_x_log10()
  p <- p + coord_cartesian(ylim = c(min_y, max_y))
  
  ggsave("cond_entropy_simulation_R.png", p, width = 14, height = 8, dpi = 300)
  print(p)
}

calculate_fidelity_from_contrast <- function(target_ratio, p_min, p_max) {
  # 1. Theoretical Limits
  max_possible_ratio <- p_max / p_min
  
  if (target_ratio >= max_possible_ratio) {
    cat(sprintf("WARNING: Target Ratio %.2f exceeds max possible (%.2f).\n", target_ratio, max_possible_ratio))
    cat("         Setting Fidelity to 1.0 (Perfect).\n")
    return(1.0)
  }
  if (target_ratio <= 1.0) {
    cat("WARNING: Target Ratio must be > 1.0. Setting Fidelity to 0.5 (Random).\n")
    return(0.5)
  }
  
  # 2. Solve for Validity (C)
  numerator <- 0.5 * (target_ratio - 1)
  denominator <- (p_max - 0.5) - target_ratio * (p_min - 0.5)
  
  # Safety check
  if (abs(denominator) < 1e-9) return(0.5)
  
  C <- numerator / denominator
  
  # 3. Convert Validity (C) to Fidelity (F)
  fidelity <- (C + 1) / 2
  return(max(0.5, min(1.0, fidelity)))
}

# --- Main Logic Starts Here ---
run_full_simulation <- function(n_super_iterations = N_SUPER_ITERATIONS) {
  
  if(GENERATE_SPARSITY_GRAPH) {
    sparsity_n_iterations <- 100 
    sparsity_results <- run_latent_sparsity_analysis(n_iterations = sparsity_n_iterations)
    plot_effective_configurations(sparsity_results, n_iterations = sparsity_n_iterations)
    k_eff_medium <- sparsity_results[Scenario == "Medium", Mean_K_eff]
    avg_rarest_prevalence_medium <- sparsity_results[Scenario == "Medium", Mean_Min_Prevalence]
    
    repeat {
      cat("\n--- S(1) PARAMETER REVIEW ---\n")
      cat("Current S(1) Parameters:\n")
      cat(sprintf("  * K_PRIMARY_LATENT (k): %d\n", K_PRIMARY_LATENT))
      cat(sprintf("  * S1_PREVALENCE_MIN:   %.3f\n", S1_PREVALENCE_MIN))
      cat(sprintf("  * S1_PREVALENCE_MAX:   %.3f\n", S1_PREVALENCE_MAX))
      cat("\nResulting Metrics (Medium Corr.):\n")
      cat(sprintf("  * K_eff_S1_Latent: %.1f\n", k_eff_medium))
      cat(sprintf("  * Rarest Prev:          %.6f\n", avg_rarest_prevalence_medium))
      
      cat("\nWould you like to modify the S(1) parameters (k, prevalence range)?\n")
      cat("1. Yes (modify parameters and re-run S(1) analysis).\n")
      cat("2. No, continue with these values.\n")
      cat("3. Quit.\n")
      
      main_choice <- readline(prompt = "Enter choice (1-3): ")
      main_choice <- trimws(main_choice)
      if (main_choice == "2") break 
      if (main_choice == "3") return()
      if (main_choice == "1") {
        if (handle_parameter_modification()) {
          cat("\n--- Re-running Latent Sparsity Analysis... ---\n")
          sparsity_results <- run_latent_sparsity_analysis(n_iterations = sparsity_n_iterations)
          plot_effective_configurations(sparsity_results, n_iterations = sparsity_n_iterations)
          k_eff_medium <- sparsity_results[Scenario == "Medium", Mean_K_eff]
          avg_rarest_prevalence_medium <- sparsity_results[Scenario == "Medium", Mean_Min_Prevalence]
        }
      } else {
        cat("Invalid choice.\n")
      }
    }
  } else {
    cat("\n--- Skipping S(1) Graph, calculating minimal parameters... ---\n")
    quick_res <- run_latent_sparsity_analysis(n_iterations = 10)
    k_eff_medium <- quick_res[Scenario == "Medium", Mean_K_eff]
    avg_rarest_prevalence_medium <- quick_res[Scenario == "Medium", Mean_Min_Prevalence]
    sparsity_results <- quick_res 
  }
  
  # --- OBSERVATIONAL ERROR CONFIGURATION (DYNAMIC WITH HETEROGENEITY) ---
  # Check if a global target is already set
  if (!is.na(TARGET_CONTRAST_RATIO)) {
    cat(sprintf("\n--- USING PRE-SET TARGET CONTRAST RATIO: %.2f ---\n", TARGET_CONTRAST_RATIO))
    calc_fid <- calculate_fidelity_from_contrast(TARGET_CONTRAST_RATIO, PARAMETRIC_MEAN_MIN, PARAMETRIC_MEAN_MAX)
    
    # Define a realistic spread (e.g., +/- 5% fidelity)
    # This creates heterogeneity: some variables are better than the target, some worse.
    fidelity_spread <- 0.05 
    
    # Calculate bounds centered on the target
    final_min <- max(0.5, calc_fid - fidelity_spread)
    final_max <- min(1.0, calc_fid + fidelity_spread)
    
    OBS_FIDELITY_MIN <<- final_min
    OBS_FIDELITY_MAX <<- final_max
    
    cat(sprintf(">> Auto-calculated Center Fidelity: %.4f\n", calc_fid))
    cat(sprintf(">> Applied Heterogeneous Range: [%.4f, %.4f]\n", final_min, final_max))
    
  } else {
    # If no global target, run the interactive menu
    repeat {
      cat("\n--- OBSERVATIONAL ERROR CONFIGURATION ---\n")
      cat(sprintf("Generative Bounds: [%.2f, %.2f] (Max Contrast = %.2f)\n", 
                  PARAMETRIC_MEAN_MIN, PARAMETRIC_MEAN_MAX, PARAMETRIC_MEAN_MAX/PARAMETRIC_MEAN_MIN))
      cat(sprintf("Current Fidelity: [%.4f, %.4f]\n", OBS_FIDELITY_MIN, OBS_FIDELITY_MAX))
      
      cat("1. Keep current fidelity settings.\n")
      cat("2. Set Manual Fidelity (e.g. 0.8 to 0.9).\n")
      cat("3. Set Target Signal-to-Floor Ratio (Dynamic Calculation).\n")
      
      err_choice <- readline(prompt = "Enter choice (1-3): ")
      err_choice <- trimws(err_choice)
      
      if (err_choice == "1") break
      if (err_choice == "2") {
        new_min <- as.numeric(readline("Enter Min Fidelity (0.5 - 1.0): "))
        new_max <- as.numeric(readline("Enter Max Fidelity (0.5 - 1.0): "))
        if(!is.na(new_min) && !is.na(new_max) && new_min <= new_max) {
          OBS_FIDELITY_MIN <<- new_min
          OBS_FIDELITY_MAX <<- new_max
        }
      }
      if (err_choice == "3") {
        cat("\nEnter the desired Contrast Ratio (Signal / Floor).\n")
        target_r <- as.numeric(readline("Target Ratio (e.g., 2.0): "))
        
        if (!is.na(target_r)) {
          calc_fid <- calculate_fidelity_from_contrast(target_r, PARAMETRIC_MEAN_MIN, PARAMETRIC_MEAN_MAX)
          
          C <- 2 * calc_fid - 1
          obs_max <- PARAMETRIC_MEAN_MAX * C + 0.5 * (1 - C)
          obs_min <- PARAMETRIC_MEAN_MIN * C + 0.5 * (1 - C)
          
          cat(sprintf("\n>> Calculated Required Fidelity: %.4f\n", calc_fid))
          cat(sprintf(">> Resulting Observed Range: [%.3f, %.3f] (Ratio: %.2f)\n", obs_min, obs_max, obs_max/obs_min))
          
          OBS_FIDELITY_MIN <<- calc_fid
          OBS_FIDELITY_MAX <<- calc_fid
        }
      }
    }
  }
  
  repeat {
    cat("\n--- GENERATIVE PHYSICS SELECTION ---\n")
    cat("Select how S(2) predictors relate to S(1) latent states:\n")
    cat("1. Causally Consistent (Smooth Manifold): Features depend on S(1) interactions via fixed coefficients.\n")
    cat("2. Chaotic Generative (Lookup Table): Features are random independent draws per configuration.\n")
    cat("3. Exit script.\n")
    
    phys_choice <- readline(prompt = "Enter choice (1-3): ")
    phys_choice <- trimws(phys_choice)
    
    if (phys_choice == "1") { 
      S2_GENERATION_MODE <<- "CONSISTENT"
      cat(">> Selected: CAUSALLY CONSISTENT mode.\n")
      break 
    } else if (phys_choice == "2") { 
      S2_GENERATION_MODE <<- "CHAOTIC"
      cat(">> Selected: CHAOTIC GENERATIVE mode.\n")
      break 
    } else if (phys_choice == "3") { 
      return() 
    } else { 
      cat("Invalid choice.\n") 
    }
  }
  
  cat("\nStarting spectral analysis...\n")
  flush.console()
  cat(sprintf("\n--- Running Part 1b: Spectral Analysis (%d Iterations) ---\n", N_SPECTRAL_ITERATIONS))
  flush.console()
  cat("This will generate the average scree plot and calculate elbow estimates from S'(2) (Effective Signal Rank)...\n")
  flush.console()
  
  # Pass updated function with NO filtering
  dynamic_keffs <- run_spectral_analysis_and_plot(k_eff_medium, n_iterations = N_SPECTRAL_ITERATIONS, calc_eigenvals_raw_func = calc_eigenvals_raw) 
  k_eff_obs_est <- dynamic_keffs$k_eff_low
  
  repeat {
    cat("\n--- N, M & Split REVIEW ---\n")
    cat(sprintf("\nCurrent N (Training) = %s, M_POOL = %s, M_EVAL = %s\n", format(N_SAMPLES, big.mark=","), format(M_POOL_SIZE, big.mark=","), format(M_PREDICTORS, big.mark=",")))
    
    train_pct <- round(100 * (1 / (1 + TEST_SET_MULTIPLIER)), 1)
    test_pct <- round(100 - train_pct, 1)
    cat(sprintf("Split: Train = %.1f%%, Test = %.1f%% (Multiplier = %.2fx)\n", train_pct, test_pct, TEST_SET_MULTIPLIER))
    
    base_target_n <- 1537
    rec_n_rare_low <- ceiling((base_target_n / avg_rarest_prevalence_medium))
    rec_n_keff_s2 <- ceiling((base_target_n * k_eff_medium))
    rec_n_keff_obs <- ceiling((base_target_n * k_eff_obs_est))
    
    cat("\nThe paper suggests N (Training) should be sufficient to characterize the S(1) -> Y link.\n")
    
    cat("\n--- S(1) Latent Sparsity Simulation Results ---\n")
    cat(sprintf("  * K_eff_S1_Latent = %.1f, Rarest S(1) Prev = %.6f\n", k_eff_medium, avg_rarest_prevalence_medium))
    cat(sprintf("\n  1. Set N (Training) to Statistically-based (for P(s_rarest)) (~%s)\n", format(rec_n_rare_low, big.mark=",")))
    cat(sprintf("  2. Set N (Training) to Heuristic (for K_eff_S1_Latent) (~%s) -- RECOMMENDED MINIMUN SAMPLE SIZE\n", format(rec_n_keff_s2, big.mark=",")))
    
    cat("\n--- S'(2) Spectral Analysis Simulation Results ---\n")
    cat(sprintf("  * K_eff_S'(2)_Obs_Est (from %d-run avg) = %d\n", N_SPECTRAL_ITERATIONS, k_eff_obs_est))
    cat(sprintf("\n  3. Set N (Training) to Heuristic (for K_eff_S'(2)_Obs_Est) (~%s)\n", format(rec_n_keff_obs, big.mark=",")))
    
    cat("\n--- Manual Adjustments ---\n")
    cat("  4. Set N (Training) to a custom value.\n")
    cat("  5. Modify M_EVAL_BUDGET (m).\n")
    cat("  6. Modify M_POOL_SIZE (Haystack).\n")
    cat("  7. Modify Test Set Multiplier (Split).\n")
    cat("  8. Continue with current settings.\n")
    cat("  9. Exit script.\n")
    
    n_choice <- readline(prompt = "Enter choice (1-9): ")
    n_choice <- trimws(n_choice)
    if (n_choice == "1") { N_SAMPLES <<- rec_n_rare_low }
    else if (n_choice == "2") { N_SAMPLES <<- rec_n_keff_s2 }
    else if (n_choice == "3") { N_SAMPLES <<- rec_n_keff_obs }
    else if (n_choice == "4") {  
      new_N_custom <- as.integer(readline("Enter N (Training): "))
      if (!is.na(new_N_custom) && new_N_custom > 0) N_SAMPLES <<- new_N_custom
    }
    else if (n_choice == "5") {  
      new_M_custom <- as.integer(readline("Enter m (Evaluation Budget): "))
      if (!is.na(new_M_custom) && new_M_custom > 0) M_PREDICTORS <<- new_M_custom
    }
    else if (n_choice == "6") {  
      new_Pool_custom <- as.integer(readline("Enter Haystack Size (Pool): "))
      if (!is.na(new_Pool_custom) && new_Pool_custom > 0) M_POOL_SIZE <<- new_Pool_custom
    }
    else if (n_choice == "7") {
      new_mult <- as.numeric(readline("Enter Test Set Multiplier (e.g., 0.25 for 80/20, 9 for 10/90): "))
      if (!is.na(new_mult) && new_mult > 0) TEST_SET_MULTIPLIER <<- new_mult
    }
    else if (n_choice == "8") { break }
    else if (n_choice == "9") { return() }
    else { cat("Invalid choice.\n") }
  }
  
  while (M_PREDICTORS < 20) {
    cat("m must be >= 20 for the Outcome-Specific P-DCAI strategy.\n")
    input_m <- readline(prompt = "Enter new m (or 'q' to quit): ")
    input_m <- trimws(input_m)
    if (tolower(input_m) == "q") return()
    M_PREDICTORS <<- as.integer(input_m)
    if(is.na(M_PREDICTORS)) M_PREDICTORS <- 0
  }
  
  repeat {
    cat("\n--- BREADTH STRATEGY SELECTION ---\n")
    cat("1. Random Breadth\n")
    cat("2. Outcome-Specific P-DCAI (Supervised)\n")
    cat("3. Spectral P-DCAI (Unsupervised Prototype-Anchored)\n")
    cat("4. MRMR (Supervised)\n")
    cat("5. All (Comparison Battle) -- RECOMMENDED\n")
    cat("6. Exit script.\n")
    
    sc <- readline("Choice (1-6): ")
    sc <- trimws(sc)
    if (sc == "1") { STRATEGY_CHOICE <<- "Random Breadth"; break }
    if (sc == "2") { STRATEGY_CHOICE <<- "Outcome-Specific P-DCAI"; break }
    if (sc == "3") { STRATEGY_CHOICE <<- "Spectral P-DCAI"; break }
    if (sc == "4") { STRATEGY_CHOICE <<- "MRMR"; break }
    if (sc == "5") { STRATEGY_CHOICE <<- "All"; break }
    if (sc == "6") { return() }
    cat("Invalid choice. Please enter a number between 1 and 6.\n")
  }
  
  repeat {
    cat("\n--- ANALYSIS MODE ---\n1. Conditional Entropy H(Y | Model Score) (Lower is better) -- RECOMMENDED\n2. AUROC (Higher is better)\n3. AUPRC (Higher is better)\n4. Exit script.\n")
    ac <- readline("Choice (1-4): ")
    ac <- trimws(ac)
    if (ac == "1") { ANALYSIS_MODE <<- "ENTROPY"; GENERATE_COND_ENTROPY_GRAPH <<- TRUE; break }
    if (ac == "2") { ANALYSIS_MODE <<- "AUROC"; GENERATE_AUROC_GRAPH <<- TRUE; break }
    if (ac == "3") { ANALYSIS_MODE <<- "AUPRC"; GENERATE_AUPRC_GRAPH <<- TRUE; break }
    if (ac == "4") { return() }
    cat("Invalid choice.\n")
  }
  
  cat("\n--- Running Part 2 Simulation ---\n")
  if (require("parallel")) total_cores_available <- parallel::detectCores() else total_cores_available <- 1
  n_inner_threads <- max(1, floor(total_cores_available*PCT_THREADS_RANGER) - 1)
  m_schedule <- create_m_schedule(M_PREDICTORS)
  all_results_dt <- data.table()
  current_total_iterations <- 0
  current_batch_size <- n_super_iterations
  
  # Initialize a global storage list
  global_artifact_store <<- list()
  
  repeat {
    if (current_batch_size <= 0) break
    start_iter <- current_total_iterations + 1
    end_iter <- current_total_iterations + current_batch_size
    
    for (i in start_iter:end_iter) {
      
      one_iter_results <- run_one_full_simulation_iteration(
        iteration_seed = i, 
        m_schedule = m_schedule, 
        k_eff_latent = k_eff_medium, 
        inner_threads = n_inner_threads, 
        save_artifacts = TRUE
      )
      
      all_results_dt <- rbindlist(list(all_results_dt, one_iter_results))
      
      # Store the artifact for this iteration
      global_artifact_store[[i]] <<- last_run_artifacts
      
      agg_results <- all_results_dt[, .(
        mean_auroc = mean(auroc), mean_auprc = mean(auprc),
        mean_oracle_auroc = mean(oracle_auroc), mean_depth_auroc = mean(depth_auroc),
        mean_oracle_auprc = mean(oracle_auprc), mean_depth_auprc = mean(depth_auprc),
        mean_breadth_cond_entropy = mean(breadth_cond_entropy),
        mean_oracle_cond_entropy = mean(oracle_cond_entropy),
        mean_depth_cond_entropy = mean(depth_cond_entropy),
        mean_H_Y_unconditional = mean(H_Y_unconditional),
        mean_avg_phi_depth = mean(avg_phi_depth),  
        mean_avg_error_rate = mean(avg_error_rate),
        mean_avg_phi = mean(avg_phi),
        mean_avg_phi_strategy = mean(avg_phi_strategy)
      ), by = .(m, error_level, strategy_type)]
      
      if(GENERATE_AUROC_GRAPH) plot_auroc_results(agg_results, i, end_iter, sparsity_results)
      if(GENERATE_AUPRC_GRAPH) plot_auprc_results(agg_results, i, end_iter, sparsity_results)
      if(GENERATE_COND_ENTROPY_GRAPH) plot_cond_entropy_results(agg_results, i, end_iter, sparsity_results)
    }
    
    current_total_iterations <- end_iter  
    cat(sprintf("\n--- Batch Complete. Total iterations run: %d ---\n", current_total_iterations))
    
    # Recalculate Aggregates for Final Report
    agg_results <- all_results_dt[, .(
      mean_auroc = mean(auroc), mean_auprc = mean(auprc),
      mean_oracle_auroc = mean(oracle_auroc), mean_depth_auroc = mean(depth_auroc),
      mean_oracle_auprc = mean(oracle_auprc), mean_depth_auprc = mean(depth_auprc),
      mean_breadth_cond_entropy = mean(breadth_cond_entropy),
      mean_oracle_cond_entropy = mean(oracle_cond_entropy),
      mean_depth_cond_entropy = mean(depth_cond_entropy),
      mean_H_Y_unconditional = mean(H_Y_unconditional),
      mean_avg_phi_depth = mean(avg_phi_depth),  
      mean_avg_error_rate = mean(avg_error_rate),
      mean_avg_phi = mean(avg_phi),
      mean_avg_phi_strategy = mean(avg_phi_strategy)
    ), by = .(m, error_level, strategy_type)]
    
    menu_loop <- TRUE
    while(menu_loop) {
      cat("\n--- Main Menu ---\n1. Run more simulation iterations.\n2. Next.\n3. Exit script.\n")
      choice <- readline(prompt = "Enter choice (1-3): ")
      choice <- trimws(choice)
      if (choice == "1") {
        add_iters_char <- readline(prompt = "How many *additional* iterations? ")
        add_iters <- as.integer(add_iters_char)
        if (!is.na(add_iters) && add_iters > 0) { current_batch_size <- add_iters; menu_loop <- FALSE }
        else { cat("Invalid input.\n") }
      } else if (choice == "2") { current_batch_size <- 0; menu_loop <- FALSE }
      else if (choice == "3") { return() }
      else { cat("Invalid choice.\n") }
    }  
    if (current_batch_size == 0) break  
  }  
  
  cat("\n--- Final Aggregation for Post-Simulation Analysis ---\n")
  agg_results <- all_results_dt[, .(
    mean_auroc = mean(auroc), mean_auprc = mean(auprc),
    mean_oracle_auroc = mean(oracle_auroc), mean_depth_auroc = mean(depth_auroc),
    mean_oracle_auprc = mean(oracle_auprc), mean_depth_auprc = mean(depth_auprc),
    mean_breadth_cond_entropy = mean(breadth_cond_entropy),
    mean_oracle_cond_entropy = mean(oracle_cond_entropy),
    mean_depth_cond_entropy = mean(depth_cond_entropy),
    mean_H_Y_unconditional = mean(H_Y_unconditional),
    mean_avg_phi_depth = mean(avg_phi_depth),  
    mean_avg_error_rate = mean(avg_error_rate),
    mean_avg_phi = mean(avg_phi),
    mean_avg_phi_strategy = mean(avg_phi_strategy)
  ), by = .(m, error_level, strategy_type)]
  
  repeat {
    cat("\n--- Post-Simulation Analysis Menu ---\n")
    cat(paste("Total iterations completed:", current_total_iterations, "\n"))
    available_options <- list()
    if (!GENERATE_COND_ENTROPY_GRAPH) available_options[[length(available_options) + 1]] <- list(text = "Generate Conditional Entropy graph.", action = "ENTROPY")
    if (!GENERATE_AUROC_GRAPH) available_options[[length(available_options) + 1]] <- list(text = "Generate AUROC graph.", action = "AUROC")
    if (!GENERATE_AUPRC_GRAPH) available_options[[length(available_options) + 1]] <- list(text = "Generate AUPRC graph.", action = "AUPRC")
    
    available_options[[length(available_options) + 1]] <- list(text = "Generate Strategy Spectral Comparison Graph (Runs new analysis).", action = "SPECTRAL_COMP")
    
    option_counter <- 1
    if (length(available_options) > 0) {
      for (i in 1:length(available_options)) {
        cat(paste(option_counter, ". ", available_options[[i]]$text, "\n", sep=""))
        option_counter <- option_counter + 1
      }
    }
    toggle_option_num <- option_counter
    cat(paste(toggle_option_num, ". Toggle X-axis scale and regenerate all active graphs.\n", sep=""))
    exit_option_num <- option_counter + 1
    cat(paste(exit_option_num, ". Exit script.\n", sep=""))
    post_choice <- as.integer(readline(prompt = paste("Enter choice (1-", exit_option_num, "): ", sep="")))
    if (is.na(post_choice)) { cat("Invalid input.\n"); next }
    if (post_choice > 0 && post_choice < option_counter) {
      chosen_action <- available_options[[post_choice]]$action
      if (chosen_action == "AUROC") {
        GENERATE_AUROC_GRAPH <<- TRUE
        plot_auroc_results(agg_results, current_total_iterations, current_total_iterations, sparsity_results)
      } else if (chosen_action == "AUPRC") {
        GENERATE_AUPRC_GRAPH <<- TRUE
        plot_auprc_results(agg_results, current_total_iterations, current_total_iterations, sparsity_results)
      } else if (chosen_action == "ENTROPY") {
        GENERATE_COND_ENTROPY_GRAPH <<- TRUE
        plot_cond_entropy_results(agg_results, current_total_iterations, current_total_iterations, sparsity_results)
      } else if (chosen_action == "SPECTRAL_COMP") {
        k_haystack <- last_run_artifacts$k_eff_empirical[["Observed"]]
        cat(sprintf("\nUsing Haystack Detected Elbow: k = %d\n", k_haystack))
        run_strategy_spectral_comparison(k_haystack)
      }
    } else if (post_choice == toggle_option_num) {
      USE_LOG_SCALE_X <<- !USE_LOG_SCALE_X
      cat("Toggled X-axis scale to:", ifelse(USE_LOG_SCALE_X, "Logarithmic", "Linear"), "\n")
      if (GENERATE_AUROC_GRAPH) plot_auroc_results(agg_results, current_total_iterations, current_total_iterations, sparsity_results)
      if (GENERATE_AUPRC_GRAPH) plot_auprc_results(agg_results, current_total_iterations, current_total_iterations, sparsity_results)
      if (GENERATE_COND_ENTROPY_GRAPH) plot_cond_entropy_results(agg_results, current_total_iterations, current_total_iterations, sparsity_results)
    } else if (post_choice == exit_option_num) {
      cat("Exiting script.\n"); return() 
    } else {
      cat(paste("Invalid choice.\n"))
    }
  }  
  cat(sprintf("\n--- Simulation Complete ---\n"))
}

if (sys.nframe() == 0) {
  run_full_simulation()
}