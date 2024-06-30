# Load necessary libraries
library(rms)
library(ggplot2)
library(parallel)
library(pbapply)
library(dplyr)

# Set the random seed for reproducibility
set.seed(1)

# Define a function to calculate bias for given sample sizes with variance heterogeneity
calculate_bias <- function(n1, n2, var1, var2) {
  y1 <- rnorm(n1, 0, sqrt(var1))
  y2 <- rnorm(n2, 1, sqrt(var2))
  
  group <- c(rep(1, n1), rep(2, n2))
  y <- c(y1, y2)
  
  conc <- (mean(rank(y)[group == 2]) - (n2 + 1) / 2) / n1
  
  f <- orm(y ~ group)
  or <- exp(coef(f)['group'])
  capprox <- or^0.66 / (1 + or^0.66)
  
  bias <- (capprox - conc) / conc * 100
  
  return(bias)
}

# Define a function to run multiple simulations for given sample sizes and variances
run_simulations <- function(total_n, ratio, num_runs, var1, var2) {
  n1 <- round(total_n * ratio / (1 + ratio))
  n2 <- total_n - n1
  replicate(num_runs, calculate_bias(n1, n2, var1, var2))
}

# Define parameters
num_runs <- 100  # Number of simulation runs for each sample size
ratios <- c(1, 5, 10, 20)  # Ratios of n1 to n2
total_sample_sizes <- c(100, 500)  # Sequence of total sample sizes
var2 <- 1  # Fixed variance for group 2
var_ratios <- c(1, 2, 5, 10)  # Variance ratios for group 1 compared to group 2 (var1 = var2 * ratio)

# Set up parallel processing
num_cores <- min(4, detectCores())
cl <- makeCluster(num_cores)

# Export necessary objects and load required packages on each core
clusterExport(cl, c("calculate_bias", "run_simulations", "num_runs", "ratios", "total_sample_sizes", "var2", "var_ratios"))
clusterEvalQ(cl, library(rms))

# Run the simulation in parallel with progress bar
bias_values <- pblapply(seq_along(ratios), function(i) {
  ratio <- ratios[i]
  lapply(total_sample_sizes, function(total_n) {
    lapply(var_ratios, function(vr) run_simulations(total_n, ratio, num_runs, var2 * vr, var2))
  })
}, cl = cl)

# Stop the cluster
stopCluster(cl)

# Process results
bias_data <- data.frame()
for (i in seq_along(ratios)) {
  ratio <- ratios[i]
  for (j in seq_along(total_sample_sizes)) {
    total_n <- total_sample_sizes[j]
    for (k in seq_along(var_ratios)) {
      vr <- var_ratios[k]
      mean_bias <- mean(bias_values[[i]][[j]][[k]])
      temp_data <- data.frame(
        TotalSampleSize = total_n,
        MeanBias = mean_bias,
        Ratio = ratio,
        VarRatio = vr
      )
      temp_data$n1 <- round(temp_data$TotalSampleSize * ratio / (1 + ratio))
      temp_data$n2 <- temp_data$TotalSampleSize - temp_data$n1
      bias_data <- rbind(bias_data, temp_data)
    }
  }
}

# Convert Ratio and VarRatio to factors for better plotting
bias_data$Ratio <- factor(bias_data$Ratio, levels = ratios)
bias_data$VarRatio <- factor(bias_data$VarRatio, levels = var_ratios)


# Plot the results using facet_grid
# Create a custom labeller for the facets
custom_labeller <- function(variable, value) {
  paste("var1:var2 Variance Ratio:", value)
}

# Plot the results using facet_wrap to create a 2x2 grid layout
plot.panel <- ggplot(bias_data, aes(x = TotalSampleSize, y = MeanBias, color = Ratio)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  labs(title = "Bias in C-index Approximation with the Proportional Odds Model",
       subtitle = paste("Based on", num_runs, "simulations per scenario"),
       x = "Total Sample Size",
       y = "Mean Bias (%)",
       color = "n1:n2 Sample Size Ratio") +
  theme_minimal(base_size = 14) +
  facet_wrap(~ VarRatio, ncol = 2, labeller = custom_labeller) +  # Create a 2x2 grid layout with custom labels
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(breaks = c(-6, -5, -4, -3, -2, -1, 0, 1)) +  # Set y-axis ticks explicitly
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),  # Make legend title bold
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

# Display the plot
print(plot.panel)

# Print summary statistics
summary_stats <- bias_data %>%
  group_by(Ratio, VarRatio) %>%
  summarise(
    Mean_Bias = mean(MeanBias),
    SD_Bias = sd(MeanBias),
    Min_Bias = min(MeanBias),
    Max_Bias = max(MeanBias)
  )
print(summary_stats)

# Save the plot
ggsave("c_index_bias_var_heterogeneity_plot.jpg", plot.panel, width = 10, height = 10, dpi = 300)

# Save the data
write.csv(bias_data, "c_index_bias_var_heterogeneity_data.csv", row.names = FALSE)

############ WITH COVERAGE PROBABILITIES ###############
# Load necessary libraries
library(rms)
library(ggplot2)
library(parallel)
library(pbapply)
library(dplyr)

# Set the random seed for reproducibility
set.seed(1)

# Define a function to calculate bias and CI coverage for given sample sizes with variance heterogeneity
calculate_bias_and_coverage <- function(n1, n2, var1, var2) {
  tryCatch({
    # Generate data
    y1 <- rnorm(n1, 0, sqrt(var1))
    y2 <- rnorm(n2, 1, sqrt(var2))
    
    group <- c(rep(1, n1), rep(2, n2))
    y <- c(y1, y2)
    
    # Calculate sample concordance
    conc <- (mean(rank(y)[group == 2]) - (n2 + 1) / 2) / n1
    
    # Normal theory-based CI for concordance
    var_conc <- (conc * (1 - conc) + (n2 - 1) * (0.5 - conc)^2 + (n1 - 1) * (0.5 - (1 - conc))^2) / (n1 * n2)
    se_conc <- sqrt(var_conc)
    ci_lower_conc_normal <- conc - 1.96 * se_conc
    ci_upper_conc_normal <- conc + 1.96 * se_conc
    
    # Bootstrap confidence interval for concordance
    boot_conc <- replicate(1000, {
      boot_sample <- sample(length(y), replace = TRUE)
      boot_y <- y[boot_sample]
      boot_group <- group[boot_sample]
      (mean(rank(boot_y)[boot_group == 2]) - (sum(boot_group == 2) + 1) / 2) / sum(boot_group == 1)
    })
    ci_conc_boot <- quantile(boot_conc, c(0.025, 0.975))
    
    # Calculate C-index approximation
    f <- orm(y ~ group)
    or <- exp(coef(f)['group'])
    capprox <- or^0.66 / (1 + or^0.66)
    
    # Calculate standard error of log odds ratio
    se_log_or <- sqrt(diag(vcov(f)))['group']
    
    # Calculate standard error for C-index approximation using the delta method
    dcapprox_dor <- 0.66 * or^(-0.34) / (1 + or^0.66)^2
    se_capprox <- abs(dcapprox_dor) * or * se_log_or
    
    # Calculate 95% confidence interval for C-index approximation
    ci_lower_approx <- capprox - 1.96 * se_capprox
    ci_upper_approx <- capprox + 1.96 * se_capprox
    
    # Calculate theoretical true concordance
    delta <- 1 / sqrt((var1 + var2) / 2)
    true_conc <- pnorm(delta / sqrt(2))
    
    # Check coverage for each method
    coverage_approx <- (true_conc >= ci_lower_approx) && (true_conc <= ci_upper_approx)
    coverage_conc_normal <- (true_conc >= ci_lower_conc_normal) && (true_conc <= ci_upper_conc_normal)
    coverage_conc_boot <- (true_conc >= ci_conc_boot[1]) && (true_conc <= ci_conc_boot[2])
    
    bias <- (capprox - true_conc) / true_conc * 100
    
    return(list(
      bias = bias, 
      coverage_approx = coverage_approx, 
      coverage_conc_normal = coverage_conc_normal,
      coverage_conc_boot = coverage_conc_boot,
      capprox = capprox, 
      conc = conc, 
      true_conc = true_conc,
      ci_lower_approx = ci_lower_approx, 
      ci_upper_approx = ci_upper_approx,
      ci_lower_conc_normal = ci_lower_conc_normal,
      ci_upper_conc_normal = ci_upper_conc_normal,
      ci_lower_conc_boot = ci_conc_boot[1],
      ci_upper_conc_boot = ci_conc_boot[2]
    ))
  }, error = function(e) {
    message("Error in calculate_bias_and_coverage: ", e$message)
    return(list(error = e$message))
  })
}

# Define a function to run multiple simulations for given sample sizes and variances
run_simulations <- function(total_n, ratio, num_runs, var1, var2) {
  n1 <- round(total_n * ratio / (1 + ratio))
  n2 <- total_n - n1
  results <- replicate(num_runs, calculate_bias_and_coverage(n1, n2, var1, var2), simplify = FALSE)
  
  bias_winP <- mean(sapply(valid_results, function(x) x$bias_winP), na.rm = TRUE)
  # Add debugging output
  cat("Mean bias_winP:", bias_winP, "\n")
  
  # Filter out error results
  valid_results <- results[!sapply(results, function(x) !is.null(x$error))]
  
  if (length(valid_results) == 0) {
    return(list(error = "All simulations failed"))
  }
  
  bias <- mean(sapply(valid_results, function(x) x$bias))
  coverage_approx <- mean(sapply(valid_results, function(x) x$coverage_approx))
  coverage_conc_normal <- mean(sapply(valid_results, function(x) x$coverage_conc_normal))
  coverage_conc_boot <- mean(sapply(valid_results, function(x) x$coverage_conc_boot))
  
  return(list(
    bias = bias, 
    coverage_approx = coverage_approx, 
    coverage_conc_normal = coverage_conc_normal, 
    coverage_conc_boot = coverage_conc_boot,
    num_valid = length(valid_results),
    num_total = num_runs
  ))
}

# Define parameters
num_runs <- 100  # Number of simulation runs for each sample size
ratios <- c(1, 5, 10, 20)  # Ratios of n1 to n2
total_sample_sizes <- c(100, 500)  # Sequence of total sample sizes
var2 <- 1  # Fixed variance for group 2
var_ratios <- c(1, 2, 5, 10)  # Variance ratios for group 1 compared to group 2 (var1 = var2 * ratio)

# Set up parallel processing
num_cores <- min(4, detectCores())
cl <- makeCluster(num_cores)

# Export necessary objects and load required packages on each core
clusterExport(cl, c("calculate_bias_and_coverage", "run_simulations", "num_runs", "ratios", "total_sample_sizes", "var2", "var_ratios"))
clusterEvalQ(cl, library(rms))

# Run the simulation in parallel with progress bar
simulation_results <- pblapply(seq_along(ratios), function(i) {
  ratio <- ratios[i]
  lapply(total_sample_sizes, function(total_n) {
    lapply(var_ratios, function(vr) {
      result <- run_simulations(total_n, ratio, num_runs, var2 * vr, var2)
      if (!is.null(result$error)) {
        message(sprintf("Error in simulation: n=%d, ratio=%f, var_ratio=%f, Error: %s", 
                        total_n, ratio, vr, result$error))
      }
      return(result)
    })
  })
}, cl = cl)

# Stop the cluster
stopCluster(cl)

# Process results
results_data <- data.frame()
for (i in seq_along(ratios)) {
  ratio <- ratios[i]
  for (j in seq_along(total_sample_sizes)) {
    total_n <- total_sample_sizes[j]
    for (k in seq_along(var_ratios)) {
      vr <- var_ratios[k]
      sim_result <- simulation_results[[i]][[j]][[k]]
      if (!is.null(sim_result$error)) {
        message(sprintf("Skipping failed simulation: n=%d, ratio=%f, var_ratio=%f", total_n, ratio, vr))
        next
      }
      temp_data <- data.frame(
        TotalSampleSize = total_n,
        MeanBias = sim_result$bias,
        CoverageProbabilityApprox = sim_result$coverage_approx,
        CoverageProbabilityConcNormal = sim_result$coverage_conc_normal,
        CoverageProbabilityConcBoot = sim_result$coverage_conc_boot,
        Ratio = ratio,
        VarRatio = vr,
        ValidSimulations = sim_result$num_valid,
        TotalSimulations = sim_result$num_total
      )
      temp_data$n1 <- round(temp_data$TotalSampleSize * ratio / (1 + ratio))
      temp_data$n2 <- temp_data$TotalSampleSize - temp_data$n1
      results_data <- rbind(results_data, temp_data)
    }
  }
}

# Convert Ratio and VarRatio to factors for better plotting
results_data$Ratio <- factor(results_data$Ratio, levels = ratios)
results_data$VarRatio <- factor(results_data$VarRatio, levels = var_ratios)

# Custom labeller for the facets
custom_labeller <- function(variable, value) {
  paste("var1:var2 Variance Ratio:", value)
}

# Plot the results for coverage probability (all three methods)
plot_coverage <- ggplot(results_data) +
  geom_line(aes(x = TotalSampleSize, y = CoverageProbabilityApprox, color = Ratio, linetype = "Approx C-index"), size = 1) +
  geom_line(aes(x = TotalSampleSize, y = CoverageProbabilityConcNormal, color = Ratio, linetype = "Normal Theory"), size = 1) +
  geom_line(aes(x = TotalSampleSize, y = CoverageProbabilityConcBoot, color = Ratio, linetype = "Bootstrap"), size = 1) +
  geom_point(aes(x = TotalSampleSize, y = CoverageProbabilityApprox, color = Ratio), size = 1.5) +
  geom_point(aes(x = TotalSampleSize, y = CoverageProbabilityConcNormal, color = Ratio), size = 1.5) +
  geom_point(aes(x = TotalSampleSize, y = CoverageProbabilityConcBoot, color = Ratio), size = 1.5) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  labs(title = "Coverage Probability Comparison",
       subtitle = paste("Based on", num_runs, "simulations per scenario"),
       x = "Total Sample Size",
       y = "Coverage Probability",
       color = "n1:n2 Sample Size Ratio",
       linetype = "Method") +
  theme_minimal(base_size = 14) +
  facet_wrap(~ VarRatio, ncol = 2, labeller = custom_labeller) +
  scale_color_brewer(palette = "Dark2") +
  scale_y_continuous(limits = c(0.5, 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )

# Print summary statistics
summary_stats <- results_data %>%
  group_by( Ratio, VarRatio, TotalSampleSize,) %>%
  summarise(
    Mean_Bias = mean(MeanBias),
    Mean_Coverage_Approx = mean(CoverageProbabilityApprox),
    Mean_Coverage_Conc_Normal = mean(CoverageProbabilityConcNormal),
    Mean_Coverage_Conc_Boot = mean(CoverageProbabilityConcBoot),
    Mean_Valid_Simulations = mean(ValidSimulations),
    .groups = 'drop'
  ) %>%
  arrange(TotalSampleSize, Ratio, VarRatio)

print(summary_stats)

# Save the plot
ggsave("c_index_coverage_comparison_plot.jpg", plot_coverage, width = 12, height = 10, dpi = 300)

# Save the data
write.csv(results_data, "c_index_bias_coverage_comparison_data.csv", row.names = FALSE)


####### Version 3 ##########

# Load necessary libraries
library(rms)
library(ggplot2)
library(parallel)
library(pbapply)
library(dplyr)

# Set the random seed for reproducibility
set.seed(1)

#' Calculate Win Probability (Concordance) and Confidence Interval
#' This function calculates the win probability (concordance) between two groups
#' and provides its confidence interval.
#' # Assuming 'mydata' is your dataset with 'outcome' and 'treatment' columns
#' result <- win_prob(mydata, response = "outcome", group = "treatment")

win_prob <-
  function(data,
           response = NULL,
           group = NULL,
           alpha = 0.05,
           beta = 0.2,
           group.ratio = 1,
           sample.size = FALSE,
           print.tables = FALSE,
           dec = 3) {
    if (is.null(response)) {
      response <- names(data)[1]
    }
    
    if (is.null(group)) {
      group <- names(data)[2]
    }
    
    group_levels <- levels(factor(data[, group]))
    
    if (length(group_levels) != 2) {
      stop("The group has to contain 2, and only 2 levels.")
    }
    
    if (!is.numeric(group.ratio)) {
      stop("Group ratio must be a numeric")
    }
    
    if (!is.logical(sample.size)) {
      stop("Sample size must be a logical")
    }
    
    # Internal helper function for calculating ranks
    freq_rank <- function(data, x = "Freq") {
      lapply(data, function(i) {
        rank <- c()
        n_i <- nrow(i)
        for (j in seq_len(n_i)) {
          if (j < n_i) {
            rank[j] <-
              ((i[j, x] + 1) / 2 + (sum(i[seq_len(n_i)[(j + 1):n_i], x])))
          }
          if (j == n_i) {
            rank[j] <- (i[j, x] + 1) / 2
          }
        }
        cbind(i, rank)
      })
    }
    
    overall <-
      freq_rank(list(data.frame(xtabs(data = data[c(response)]))))[[1]]
    
    tbl <- xtabs(data = data[c(response, group)])
    
    tbl_df <- data.frame(tbl)
    
    prop_df <- data.frame(proportions(tbl, group))
    
    df <- cbind(tbl_df, prop = prop_df[, "Freq"])
    
    list_cum <- split(df, df[, group])
    
    list_cum <- lapply(list_cum, function(i) {
      data.frame(i, overall_rank = overall$rank)
    })
    
    list_cum <- freq_rank(list_cum)
    
    sum_a <- sum(df$Freq[df$rtreat == group_levels[1]])
    sum_b <- sum(df$Freq[df$rtreat == group_levels[2]])
    
    list_cum[[1]]$win_frac <-
      with(list_cum[[1]], overall_rank - rank) / sum_b
    list_cum[[2]]$win_frac <-
      with(list_cum[[2]], overall_rank - rank) / sum_a
    
    winP_a <- sum(with(list_cum[[1]], prop * win_frac))
    winP_b <- sum(with(list_cum[[2]], prop * win_frac))
    
    var_win_frac_a <-
      sum(with(list_cum[[1]], prop * win_frac ^ 2)) - winP_a ^ 2
    var_win_frac_b <-
      sum(with(list_cum[[2]], prop * win_frac ^ 2)) - winP_b ^ 2
    
    var_win_prob <- var_win_frac_a / sum_a + var_win_frac_b / sum_b
    
    se_win_prob <- sqrt(var_win_prob)
    
    ci_up <-
      exp(log(winP_a / (1 - winP_a)) - qnorm(1 - alpha / 2) * se_win_prob /
            (winP_a / (1 - winP_a))) /
      (1 + exp(log(winP_a / (1 - winP_a)) - qnorm(1 - alpha / 2) *
                 se_win_prob / (winP_a / (1 - winP_a))))
    
    ci_lo <-
      exp(log(winP_b / (1 - winP_b)) + qnorm(1 - alpha / 2) * se_win_prob /
            (winP_b / (1 - winP_b))) /
      (1 + exp(log(winP_b / (1 - winP_b)) + qnorm(1 - alpha / 2) *
                 se_win_prob / (winP_b / (1 - winP_b))))
    
    test_stat <-
      abs(log(winP_a / (1 - winP_a))) / (se_win_prob / (winP_a * (1 - winP_a)))
    p_val <- 2 * (1 - pnorm(test_stat))
    
    nnt <- 1 / (winP_a - 0.5)
    
    ss_n <- NA
    
    if (sample.size) {
      ss_n <-
        ceiling((group.ratio + 1) / group.ratio *
                  (qnorm(1 - alpha / 2) + qnorm(1 - beta)) ^ 2 *
                  (var_win_frac_a + group.ratio * var_win_frac_b) /
                  (winP_a * (1 - winP_a) * log(winP_a / (1 - winP_a))) ^ 2
        )
    }
    
    out <- list(
      list_cum = list_cum,
      group_levels = group_levels,
      sum_a = sum_a,
      sum_b = sum_b,
      winP_a = winP_a,
      winP_b = winP_b,
      var_win_frac_a = var_win_frac_a,
      var_win_frac_b = var_win_frac_b,
      var_win_prob = var_win_prob,
      se_win_prob = se_win_prob,
      conf.int = c(ci_lo, ci_up),
      test_stat = test_stat,
      p_val = p_val,
      nnt = nnt,
      ss_n = ss_n,
      param.record = list(
        data = data,
        response = response,
        group = group,
        alpha = alpha,
        beta = beta,
        group.ratio = group.ratio,
        sample.size = sample.size,
        print.tables = print.tables,
        dec = dec
      )
    )
    class(out) <- c("win_Prob", class(out))
    return(out)
  }



# Define simulation parameters
num_runs <- 100  # Number of simulation runs for each scenario
ratios <- c(1, 5, 10, 20)  # Ratios of n1 to n2
total_sample_sizes <- c(100, 500)  # Sequence of total sample sizes
var2 <- 1  # Fixed variance for group 2
var_ratios <- c(1, 2, 5, 10)  # Variance ratios for group 1 compared to group 2 (var1 = var2 * ratio)

# Define the calculate_bias_and_coverage function
calculate_bias_and_coverage <- function(n1, n2, var1, var2) {
  tryCatch({
    # Generate data
    y1 <- rnorm(n1, 0, sqrt(var1))
    y2 <- rnorm(n2, 1, sqrt(var2))
    
    group <- c(rep(1, n1), rep(2, n2))
    y <- c(y1, y2)
    
    # Create a data frame for win_prob function
    sim_data <- data.frame(y = y, group = factor(group))
    
    # Calculate sample concordance
    conc <- (mean(rank(y)[group == 2]) - (n2 + 1) / 2) / n1
    
    # Normal theory-based CI for concordance
    var_conc <- (conc * (1 - conc) + (n2 - 1) * (0.5 - conc)^2 + (n1 - 1) * (0.5 - (1 - conc))^2) / (n1 * n2)
    se_conc <- sqrt(var_conc)
    ci_lower_conc_normal <- conc - 1.96 * se_conc
    ci_upper_conc_normal <- conc + 1.96 * se_conc
    
    # Bootstrap confidence interval for concordance
    boot_conc <- replicate(1000, {
      boot_sample <- sample(length(y), replace = TRUE)
      boot_y <- y[boot_sample]
      boot_group <- group[boot_sample]
      (mean(rank(boot_y)[boot_group == 2]) - (sum(boot_group == 2) + 1) / 2) / sum(boot_group == 1)
    })
    ci_conc_boot <- quantile(boot_conc, c(0.025, 0.975))
    
    # Calculate win probability and its CI
    wp_result <- win_prob(sim_data, response = "y", group = "group")
    winP <- wp_result$winP_a
    ci_lower_winP <- wp_result$conf.int[1]
    ci_upper_winP <- wp_result$conf.int[2]
    
    # Calculate C-index approximation using the proportional odds model
    f <- orm(y ~ group)
    or <- exp(coef(f)['group'])
    capprox <- or^0.66 / (1 + or^0.66)
    
    # Calculate CI for C-index approximation from the proportional odds model
    se_log_or <- sqrt(diag(vcov(f)))['group']
    dcapprox_dor <- 0.66 * or^(-0.34) / (1 + or^0.66)^2
    se_capprox <- abs(dcapprox_dor) * or * se_log_or
    ci_lower_approx <- capprox - 1.96 * se_capprox
    ci_upper_approx <- capprox + 1.96 * se_capprox
    
    # Calculate theoretical true concordance
    delta <- 1 / sqrt((var1 + var2) / 2)
    true_conc <- pnorm(delta / sqrt(2))
    
    # Check coverage for each method
    coverage_approx <- (true_conc >= ci_lower_approx) && (true_conc <= ci_upper_approx)
    coverage_winP <- (true_conc >= ci_lower_winP) && (true_conc <= ci_upper_winP)
    coverage_conc_normal <- (true_conc >= ci_lower_conc_normal) && (true_conc <= ci_upper_conc_normal)
    coverage_conc_boot <- (true_conc >= ci_conc_boot[1]) && (true_conc <= ci_conc_boot[2])
    
    bias_approx <- (capprox - true_conc) / true_conc * 100
    bias_winP <- (winP - true_conc) / true_conc * 100
    bias_conc <- (conc - true_conc) / true_conc * 100
    
    return(list(
      bias_approx = bias_approx,
      bias_winP = bias_winP,
      bias_conc = bias_conc,
      coverage_approx = coverage_approx, 
      coverage_winP = coverage_winP,
      coverage_conc_normal = coverage_conc_normal,
      coverage_conc_boot = coverage_conc_boot,
      capprox = capprox,
      winP = winP,
      conc = conc,
      true_conc = true_conc
    ))
  }, error = function(e) {
    message("Error in calculate_bias_and_coverage: ", e$message)
    return(list(error = e$message))
  })
}

# Define the run_simulations function
run_simulations <- function(total_n, ratio, num_runs, var1, var2) {
  n1 <- round(total_n * ratio / (1 + ratio))
  n2 <- total_n - n1
  results <- replicate(num_runs, calculate_bias_and_coverage(n1, n2, var1, var2), simplify = FALSE)
  
  valid_results <- results[!sapply(results, function(x) !is.null(x$error))]
  
  if (length(valid_results) == 0) {
    return(list(error = "All simulations failed"))
  }
  
  bias_approx <- mean(sapply(valid_results, function(x) x$bias_approx))
  bias_winP <- mean(sapply(valid_results, function(x) x$bias_winP))
  bias_conc <- mean(sapply(valid_results, function(x) x$bias_conc))
  coverage_approx <- mean(sapply(valid_results, function(x) x$coverage_approx))
  coverage_winP <- mean(sapply(valid_results, function(x) x$coverage_winP))
  coverage_conc_normal <- mean(sapply(valid_results, function(x) x$coverage_conc_normal))
  coverage_conc_boot <- mean(sapply(valid_results, function(x) x$coverage_conc_boot))
  
  return(list(
    bias_approx = bias_approx,
    bias_winP = bias_winP,
    bias_conc = bias_conc,
    coverage_approx = coverage_approx,
    coverage_winP = coverage_winP,
    coverage_conc_normal = coverage_conc_normal,
    coverage_conc_boot = coverage_conc_boot,
    num_valid = length(valid_results),
    num_total = num_runs
  ))
}

# Set up parallel processing
num_cores <- min(4, detectCores())
cl <- makeCluster(num_cores)

# Export necessary objects and load required packages on each core
clusterExport(cl, c("calculate_bias_and_coverage", "run_simulations", "win_prob", "num_runs", "ratios", "total_sample_sizes", "var2", "var_ratios"))
clusterEvalQ(cl, {
  library(rms)
  library(dplyr)
})

# Run the simulation in parallel with progress bar
simulation_results <- pblapply(seq_along(ratios), function(i) {
  ratio <- ratios[i]
  lapply(total_sample_sizes, function(total_n) {
    lapply(var_ratios, function(vr) {
      result <- run_simulations(total_n, ratio, num_runs, var2 * vr, var2)
      if (!is.null(result$error)) {
        message(sprintf("Error in simulation: n=%d, ratio=%f, var_ratio=%f, Error: %s", 
                        total_n, ratio, vr, result$error))
      }
      return(result)
    })
  })
}, cl = cl)

# Stop the cluster
stopCluster(cl)

# Process results
results_data <- data.frame()
for (i in seq_along(ratios)) {
  ratio <- ratios[i]
  for (j in seq_along(total_sample_sizes)) {
    total_n <- total_sample_sizes[j]
    for (k in seq_along(var_ratios)) {
      vr <- var_ratios[k]
      sim_result <- simulation_results[[i]][[j]][[k]]
      if (is.null(sim_result) || !is.null(sim_result$error)) {
        message(sprintf("Skipping failed simulation: n=%d, ratio=%f, var_ratio=%f", total_n, ratio, vr))
        next
      }
      
      # Create a list of the data we want to include
      result_list <- list(
        TotalSampleSize = total_n,
        Ratio = ratio,
        VarRatio = vr,
        BiasApprox = sim_result$bias_approx,
        BiasWinP = sim_result$bias_winP,
        BiasConc = sim_result$bias_conc,
        CoverageApprox = sim_result$coverage_approx,
        CoverageWinP = sim_result$coverage_winP,
        CoverageConcNormal = sim_result$coverage_conc_normal,
        CoverageConcBoot = sim_result$coverage_conc_boot,
        ValidSimulations = sim_result$num_valid,
        TotalSimulations = sim_result$num_total
      )
      
      # Remove any NULL elements from the list
      result_list <- result_list[!sapply(result_list, is.null)]
      
      # Convert the list to a data frame
      temp_data <- as.data.frame(result_list)
      
      # Add n1 and n2 calculations
      temp_data$n1 <- round(temp_data$TotalSampleSize * ratio / (1 + ratio))
      temp_data$n2 <- temp_data$TotalSampleSize - temp_data$n1
      
      # Combine with existing results
      if (nrow(results_data) == 0) {
        results_data <- temp_data
      } else {
        results_data <- rbind(results_data, temp_data)
      }
    }
  }
}

# Print the structure of results_data to check
print(str(results_data))

# Convert Ratio and VarRatio to factors for better plotting
results_data$Ratio <- factor(results_data$Ratio, levels = ratios)
results_data$VarRatio <- factor(results_data$VarRatio, levels = var_ratios)

# Calculate summary statistics
summary_stats <- results_data %>%
  group_by(TotalSampleSize, Ratio, VarRatio) %>%
  summarise(
    Mean_Bias_Approx = mean(BiasApprox, na.rm = TRUE),
    Mean_Bias_WinP = mean(BiasWinP, na.rm = TRUE),
    Mean_Bias_Conc = mean(BiasConc, na.rm = TRUE),
    Mean_Coverage_Approx = mean(CoverageApprox, na.rm = TRUE),
    Mean_Coverage_WinP = mean(CoverageWinP, na.rm = TRUE),
    Mean_Coverage_Conc_Normal = mean(CoverageConcNormal, na.rm = TRUE),
    Mean_Coverage_Conc_Boot = mean(CoverageConcBoot, na.rm = TRUE),
    Mean_Valid_Simulations = mean(ValidSimulations, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(TotalSampleSize, Ratio, VarRatio)

# Print summary statistics
print(summary_stats)

# Save results
write.csv(results_data, "c_index_simulation_results.csv", row.names = FALSE)
write.csv(summary_stats, "c_index_summary_statistics.csv", row.names = FALSE)

# Print a message to indicate completion
cat("\nSimulation completed. Results saved to CSV files.\n")





######### REDO ######

win_prob <-
  function(data,
           response = NULL,
           group = NULL,
           alpha = 0.05,
           beta = 0.2,
           group.ratio = 1,
           sample.size = FALSE,
           print.tables = FALSE,
           dec = 3) {
    if (is.null(response)) {
      response <- names(data)[1]
    }
    
    if (is.null(group)) {
      group <- names(data)[2]
    }
    
    group_levels <- levels(factor(data[, group]))
    
    if (length(group_levels) != 2) {
      stop("The group has to contain 2, and only 2 levels.")
    }
    
    if (!is.numeric(group.ratio)) {
      stop("Group ratio must be a numeric")
    }
    
    if (!is.logical(sample.size)) {
      stop("Sample size must be a logical")
    }
    
    # Internal helper function for calculating ranks
    freq_rank <- function(data, x = "Freq") {
      lapply(data, function(i) {
        rank <- c()
        n_i <- nrow(i)
        for (j in seq_len(n_i)) {
          if (j < n_i) {
            rank[j] <-
              ((i[j, x] + 1) / 2 + (sum(i[seq_len(n_i)[(j + 1):n_i], x])))
          }
          if (j == n_i) {
            rank[j] <- (i[j, x] + 1) / 2
          }
        }
        cbind(i, rank)
      })
    }
    
    overall <-
      freq_rank(list(data.frame(xtabs(data = data[c(response)]))))[[1]]
    
    tbl <- xtabs(data = data[c(response, group)])
    
    tbl_df <- data.frame(tbl)
    
    prop_df <- data.frame(proportions(tbl, group))
    
    df <- cbind(tbl_df, prop = prop_df[, "Freq"])
    
    list_cum <- split(df, df[, group])
    
    list_cum <- lapply(list_cum, function(i) {
      data.frame(i, overall_rank = overall$rank)
    })
    
    list_cum <- freq_rank(list_cum)
    
    sum_a <- sum(df$Freq[df$rtreat == group_levels[1]])
    sum_b <- sum(df$Freq[df$rtreat == group_levels[2]])
    
    list_cum[[1]]$win_frac <-
      with(list_cum[[1]], overall_rank - rank) / sum_b
    list_cum[[2]]$win_frac <-
      with(list_cum[[2]], overall_rank - rank) / sum_a
    
    winP_a <- sum(with(list_cum[[1]], prop * win_frac))
    winP_b <- sum(with(list_cum[[2]], prop * win_frac))
    
    var_win_frac_a <-
      sum(with(list_cum[[1]], prop * win_frac ^ 2)) - winP_a ^ 2
    var_win_frac_b <-
      sum(with(list_cum[[2]], prop * win_frac ^ 2)) - winP_b ^ 2
    
    var_win_prob <- var_win_frac_a / sum_a + var_win_frac_b / sum_b
    
    se_win_prob <- sqrt(var_win_prob)
    
    ci_up <-
      exp(log(winP_a / (1 - winP_a)) - qnorm(1 - alpha / 2) * se_win_prob /
            (winP_a / (1 - winP_a))) /
      (1 + exp(log(winP_a / (1 - winP_a)) - qnorm(1 - alpha / 2) *
                 se_win_prob / (winP_a / (1 - winP_a))))
    
    ci_lo <-
      exp(log(winP_b / (1 - winP_b)) + qnorm(1 - alpha / 2) * se_win_prob /
            (winP_b / (1 - winP_b))) /
      (1 + exp(log(winP_b / (1 - winP_b)) + qnorm(1 - alpha / 2) *
                 se_win_prob / (winP_b / (1 - winP_b))))
    
    test_stat <-
      abs(log(winP_a / (1 - winP_a))) / (se_win_prob / (winP_a * (1 - winP_a)))
    p_val <- 2 * (1 - pnorm(test_stat))
    
    nnt <- 1 / (winP_a - 0.5)
    
    ss_n <- NA
    
    if (sample.size) {
      ss_n <-
        ceiling((group.ratio + 1) / group.ratio *
                  (qnorm(1 - alpha / 2) + qnorm(1 - beta)) ^ 2 *
                  (var_win_frac_a + group.ratio * var_win_frac_b) /
                  (winP_a * (1 - winP_a) * log(winP_a / (1 - winP_a))) ^ 2
        )
    }
    
    out <- list(
      list_cum = list_cum,
      group_levels = group_levels,
      sum_a = sum_a,
      sum_b = sum_b,
      winP_a = winP_a,
      winP_b = winP_b,
      var_win_frac_a = var_win_frac_a,
      var_win_frac_b = var_win_frac_b,
      var_win_prob = var_win_prob,
      se_win_prob = se_win_prob,
      conf.int = c(ci_lo, ci_up),
      test_stat = test_stat,
      p_val = p_val,
      nnt = nnt,
      ss_n = ss_n,
      param.record = list(
        data = data,
        response = response,
        group = group,
        alpha = alpha,
        beta = beta,
        group.ratio = group.ratio,
        sample.size = sample.size,
        print.tables = print.tables,
        dec = dec
      )
    )
    class(out) <- c("win_Prob", class(out))
    return(out)
  }




#' @title Prints win_prob results
#' @param x win_prob results.
#' @param ... ignored for now
#' @return Prints win_prob statistics. 
#' @export
print.win_Prob <- function (x, ...) {
  args <- list(...)
  
  cat("\t Zou et al's winP (doi: 10.1161/STROKEAHA.121.037744) \n\n")
  cat(
    sprintf(
      "Probability of a random observation in %s group 
      will have a higher response score than a random
      observation in %s group:\n\n",
      x$group_levels[2],
      x$group_levels[1]
    )
  )
  
  cat("  ")
  
  cat(
    sprintf(
      "\t   winP: %2.3f (%2.3f, %2.3f)      p=%1.4f",
      x$winP_a,
      x$conf.int[1],
      x$conf.int[2],
      x$p_val
    )
  )
  
  cat("\n")
  
  cat("--------------------------------------------\n\n")
  
  cat(sprintf("The numbers needed to treat (NNT) are: %s\n\n",
              ceiling(x$nnt)))
  
  cat("\n")
  
  if (x$param.record$sample.size) {
    cat("--------------------------------------------\n\n")
    
    cat(
      sprintf(
        "\tWith %s/%s ratio = %s and beta = %s
              the sample size needed is: %s\n\n",
        x$group_levels[1],
        x$group_levels[2],
        x$param.record$group.ratio,
        x$param.record$beta,
        x$ss_n
      )
    )
  }
  
  cat("\n")
  
  if (x$param.record$print.tables) {
    cat("--------------------------------------------\n\n")
    
    lc <- lapply(x$list_cum, function(i, t = c("prop", "win_frac")) {
      i[t] <- round(i[t], x$param.record$dec)
      i[, !names(i) == x$param.record$group]
    })
    
    for (i in x$group_levels) {
      tab <- knitr::kable(lc[[i]], row.names = FALSE)
      cat(sprintf("Results for the %s group:\n", i))
      cat(sprintf("\t%s\n",
                  tab))
      cat("\n")
    }
  }
  return(invisible(x))
}
# Load necessary librarie# Load necessary libraries
library(dplyr)

# Set the random seed for reproducibility
set.seed(1)

# Create a simulated dataset
n1 <- 100  # Sample size for group 1
n2 <- 100  # Sample size for group 2
response_group1 <- rnorm(n1, mean = 0, sd = 1)  # Response variable for group 1
response_group2 <- rnorm(n2, mean = 1.5, sd = 1)  # Response variable for group 2

# Combine the data into a single data frame
simulated_data <- data.frame(
  response = c(response_group1, response_group2),
  group = factor(c(rep(1, n1), rep(2, n2)))
)

# Print the first few rows of the simulated dataset
head(simulated_data)

# Test the win_prob function with the simulated dataset
win_prob_result <- win_prob(simulated_data, response = "response", group = "group")

# Print the result
print(win_prob_result)

win_prob_result$conf.int

print.win_Prob(win_prob_result)





# Define the function
calculate_winP <- function(data, group_var, post_var) {
  
  # Calculate overall rank
  data <- data[order(data[[group_var]], -data[[post_var]]), ]
  data$rpost <- rank(-data[[post_var]], ties.method = "average")
  
  # Calculate group-specific rank
  data$grpost <- ave(data$rpost, data[[group_var]], FUN = function(x) rank(-x, ties.method = "average"))
  
  # Calculate win fractions
  freqcnt <- aggregate(id ~ get(group_var), data, length)
  names(freqcnt) <- c(group_var, "count")
  freqcnt[[group_var]] <- 1 - freqcnt[[group_var]]
  freqcnt <- freqcnt[order(freqcnt[[group_var]]), ]
  
  data <- merge(data, freqcnt, by = group_var)
  data$winf <- (data$rpost - data$grpost) / data$count
  
  # Perform regression on win fractions
  model1 <- lm(winf ~ as.factor(data[[group_var]]), data = data)
  
  # Extract coefficients and standard errors
  estMW <- summary(model1)$coefficients["as.factor(data[[group_var]])1", ]
  
  # Calculate Win Probability and Other Metrics
  WinP <- estMW["Estimate"] / 2 + 0.5
  lgt <- log(WinP / (1 - WinP))
  lgtSe <- estMW["Std. Error"] / (WinP * (1 - WinP))
  ln_Odds_l <- lgt - 1.96 * lgtSe
  ln_Odds_u <- lgt + 1.96 * lgtSe
  
  Odds <- exp(lgt)
  Odds_l <- exp(ln_Odds_l)
  Odds_u <- exp(ln_Odds_u)
  
  WinP_l <- 1 / (1 + exp(-ln_Odds_l))  # logistic function
  WinP_u <- 1 / (1 + exp(-ln_Odds_u))  # logistic function
  test <- (lgt - 0) / lgtSe
  pvalue <- 2 * (1 - pnorm(abs(test)))
  SomersD <- 2 * WinP - 1
  D_L <- 2 * WinP_l - 1  # lower limit for Somers D
  D_U <- 2 * WinP_u - 1  # upper limit for Somers D
  
  # Create a data frame to store the results
  WinP_result <- data.frame(
    WinP = WinP,
    WinP_l = WinP_l,
    WinP_u = WinP_u,
    test = test,
    pvalue = pvalue,
    Odds = Odds,
    Odds_l = Odds_l,
    Odds_u = Odds_u,
    SomersD = SomersD,
    D_L = D_L,
    D_U = D_U
  )
  
  # Return the results
  return(WinP_result)
}

# Test the function with the simulated dataset
data <- data.frame(
  id = 1:61,
  group = c(rep(0, 27), rep(1, 34)),
  pre = c(18, 27, 16, 17, 15, 20, 16, 28, 28, 25, 24, 16, 26, 21, 21, 22, 26, 19, 22, 16, 21, 20, 17, 22, 19, 21, 18, 21, 27, 15, 24, 15, 17, 20, 18, 28, 21, 18, 27.46, 19, 20, 16, 21, 23, 23, 24, 25, 22, 20, 20, 25, 18, 26, 20, 17, 22, 22, 23, 17, 22, 26),
  post = c(17, 26, 17, 14, 12, 19, 13, 26, 26, 9, 14, 19, 13, 7, 18, 18, 19, 19, 20, 7, 19, 16, 15, 20, 16, 7, 19, 13, 8, 8, 14, 15, 9, 7, 8, 11, 7, 8, 22, 14, 13, 17, 19, 11, 16, 16, 20, 15, 7, 12.13, 15, 17, 1, 27, 20, 12, 15.38, 11, 15, 7, 24)
)

# Test the function
win_prob_result <- calculate_winP(data, group_var = "group", post_var = "post")

# Print the result
print(win_prob_result)












######## AGAIN #####################################


# Define the function
calculate_winP_pre_post <- function(data, group_var, post_var, pre_var = NULL) {
 
  required_packages <- c("dplyr", "sandwich", "lmtest", "broom")
  
  # Install missing packages
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
  
  
  # Convert data to win fractions
  
  if (!is.null(pre_var)) {
    # Calculate overall rank
    data <- data %>%
      arrange(!!sym(group_var), desc(!!sym(pre_var)), desc(!!sym(post_var))) %>%
      mutate(rpre = rank(-!!sym(pre_var), ties.method = "average"),
             rpost = rank(-!!sym(post_var), ties.method = "average"))
    
    # Calculate group-specific rank
    data <- data %>%
      group_by(!!sym(group_var)) %>%
      mutate(grpre = rank(-!!sym(pre_var), ties.method = "average"),
             grpost = rank(-!!sym(post_var), ties.method = "average")) %>%
      ungroup()
    
    # Calculate win fractions
    freqcnt <- data %>%
      group_by(!!sym(group_var)) %>%
      summarise(count = n()) %>%
      mutate(!!sym(group_var) := 1 - !!sym(group_var)) %>%
      arrange(!!sym(group_var))
    
    data <- data %>%
      left_join(freqcnt, by = group_var) %>%
      mutate(winf = (rpost - grpost) / count,
             winf_pre = (rpre - grpre) / count)
    
    # Perform regression on win fractions
    formula2 <- as.formula(paste("winf ~ winf_pre +", group_var))
    model2 <- lm(formula2, data = data)
  } else {
    # Calculate overall rank
    data <- data %>%
      arrange(!!sym(group_var), desc(!!sym(post_var))) %>%
      mutate(rpost = rank(-!!sym(post_var), ties.method = "average"))
    
    # Calculate group-specific rank
    data <- data %>%
      group_by(!!sym(group_var)) %>%
      mutate(grpost = rank(-!!sym(post_var), ties.method = "average")) %>%
      ungroup()
    
    # Calculate win fractions
    freqcnt <- data %>%
      group_by(!!sym(group_var)) %>%
      summarise(count = n()) %>%
      mutate(!!sym(group_var) := 1 - !!sym(group_var)) %>%
      arrange(!!sym(group_var))
    
    data <- data %>%
      left_join(freqcnt, by = group_var) %>%
      mutate(winf = (rpost - grpost) / count)
    
    # Perform regression on win fractions
    formula1 <- as.formula(paste("winf ~", group_var))
    model2 <- lm(formula1, data = data)
  }
  
  # Extract coefficients and standard errors with HC2 adjustment
  hc_se <- coeftest(model2, vcov = vcovHC(model2, type = "HC2"))
  
  # Get the coefficient for the group variable
  estMW <- tidy(hc_se) %>%
    filter(term == group_var)
  
  # Check if estMW is empty
  if (nrow(estMW) == 0) {
    stop("No coefficients found for the specified group variable. Check your input data.")
  }
  
  # Calculate Win Probability and Other Metrics
  WinP <- estMW$estimate / 2 + 0.5
  lgt <- log(WinP / (1 - WinP))
  lgtSe <- estMW$std.error / (WinP * (1 - WinP))
  ln_Odds_l <- lgt - 1.96 * lgtSe
  ln_Odds_u <- lgt + 1.96 * lgtSe
  
  Odds <- exp(lgt)
  Odds_l <- exp(ln_Odds_l)
  Odds_u <- exp(ln_Odds_u)
  
  WinP_l <- 1 / (1 + exp(-ln_Odds_l))  # logistic function
  WinP_u <- 1 / (1 + exp(-ln_Odds_u))  # logistic function
  test <- (lgt - 0) / lgtSe
  pvalue <- 2 * (1 - pnorm(abs(test)))
  SomersD <- 2 * WinP - 1
  D_L <- 2 * WinP_l - 1  # lower limit for Somers D
  D_U <- 2 * WinP_u - 1  # upper limit for Somers D
  
  # Create a data frame to store the results
  WinP_result <- data.frame(
    WinP = WinP,
    WinP_l = WinP_l,
    WinP_u = WinP_u,
    test = test,
    pvalue = pvalue,
    Odds = Odds,
    Odds_l = Odds_l,
    Odds_u = Odds_u,
    SomersD = SomersD,
    D_L = D_L,
    D_U = D_U
  )
  
  # Print the results
  return(WinP_result)
}

# Example usage with pre_var
data <- data.frame(
  id = 1:61,
  group = c(rep(0, 27), rep(1, 34)),
  pre = c(18, 27, 16, 17, 15, 20, 16, 28, 28, 25, 24, 16, 26, 21, 21, 22, 26, 19, 22, 16, 21, 20, 17, 22, 19, 21, 18, 21, 27, 15, 24, 15, 17, 20, 18, 28, 21, 18, 27.46, 19, 20, 16, 21, 23, 23, 24, 25, 22, 20, 20, 25, 18, 26, 20, 17, 22, 22, 23, 17, 22, 26),
  post = c(17, 26, 17, 14, 12, 19, 13, 26, 26, 9, 14, 19, 13, 7, 18, 18, 19, 19, 20, 7, 19, 16, 15, 20, 16, 7, 19, 13, 8, 8, 14, 15, 9, 7, 8, 11, 7, 8, 22, 14, 13, 17, 19, 11, 16, 16, 20, 15, 7, 12.13, 15, 17, 1, 27, 20, 12, 15.38, 11, 15, 7, 24)
)

result_with_pre <- calculate_winP_pre_post(data, group_var = "group", pre_var = "pre", post_var = "post")
print(result_with_pre)

# Example usage without pre_var
result_without_pre <- calculate_winP_pre_post(data, group_var = "group", post_var = "post")
print(result_without_pre)




