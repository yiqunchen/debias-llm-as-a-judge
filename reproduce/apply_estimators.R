# Apply debiasing estimators to LLM-as-a-Judge results
# Uses 10% labeled sample to correct for judge bias

library(tidyverse)
library(jsonlite)
library(debiasLLMReporting)

set.seed(2026)

# Configuration
LABEL_RATIO <- 0.10  # 10% labeled
ALPHA <- 0.10        # 90% CI

# Color scheme matching simulation plots
METHOD_COLORS <- c(
  "PPI" = "#1B9E77",
  "PPI++" = "#D95F02",
  "Rogan-Gladen" = "#7570B3",
  "Joint MLE" = "#E7298A",
  "Naive" = "#525252",
  "EIF" = "#E6AB02",
  "True" = "black"
)

#' Apply all estimators to a single dataset
#'
#' @param judge_df Data frame with human_winner and judge_pick columns
#' @param label_ratio Fraction to use as labeled calibration set
#' @param alpha Significance level for CIs
#' @return Data frame with estimator results
apply_all_estimators <- function(judge_df, label_ratio = 0.10, alpha = 0.10) {
  # Filter to valid rows (judge made a pick)
  df <- judge_df %>%
    filter(!is.na(judge_pick)) %>%
    mutate(
      # Y = 1 if human picked model_a, 0 otherwise
      Y = as.integer(human_winner == "model_a"),
      # Yhat = 1 if judge picked A, 0 otherwise
      Yhat = as.integer(judge_pick == "A")
    )

  N <- nrow(df)
  m <- max(4L, round(label_ratio * N))  # Calibration size
  n <- N - m  # Test size

  # Random split
  idx_cal <- sample(N, m)
  idx_test <- setdiff(seq_len(N), idx_cal)

  # Calibration (labeled)
  Y_cal <- df$Y[idx_cal]
  Yhat_cal <- df$Yhat[idx_cal]

  # Test (unlabeled)
  Yhat_test <- df$Yhat[idx_test]

  # True win rate for model_a (using all data)
  theta_true <- mean(df$Y)

  # Class counts in calibration
  m0 <- sum(Y_cal == 0)
  m1 <- sum(Y_cal == 1)

  # Confusion matrix estimates
  q0_hat <- mean(Yhat_cal[Y_cal == 0] == 0)  # Specificity
  q1_hat <- mean(Yhat_cal[Y_cal == 1] == 1)  # Sensitivity
  p_hat <- mean(Yhat_test)  # Test positive rate

  results <- list()

  # 1. Naive (no correction)
  naive_theta <- p_hat
  naive_var <- p_hat * (1 - p_hat) / n
  z_alpha <- qnorm(1 - alpha / 2)
  results$naive <- tibble(
    method = "Naive",
    theta_hat = naive_theta,
    ci_lower = pmax(naive_theta - z_alpha * sqrt(naive_var), 0),
    ci_upper = pmin(naive_theta + z_alpha * sqrt(naive_var), 1)
  )

  # 2. PPI
  if (m0 >= 2 && m1 >= 2) {
    ppi_est <- ppi_point_and_ci(Y_L = Y_cal, f_L = Yhat_cal, f_U = Yhat_test,
                                 alpha = alpha)
    results$ppi <- tibble(
      method = "PPI",
      theta_hat = ppi_est$theta,
      ci_lower = ppi_est$ci_lower,
      ci_upper = ppi_est$ci_upper
    )

    # 3. PPI++
    ppi_pp_est <- ppi_pp_point_and_ci(Y_L = Y_cal, f_L = Yhat_cal,
                                      f_U = Yhat_test, alpha = alpha)
    results$ppi_pp <- tibble(
      method = "PPI++",
      theta_hat = ppi_pp_est$theta,
      ci_lower = ppi_pp_est$ci_lower,
      ci_upper = ppi_pp_est$ci_upper
    )

    # 4. Rogan-Gladen
    rg_est <- rg_point_and_ci(
      p_hat = p_hat,
      q0_hat = q0_hat,
      q1_hat = q1_hat,
      n = n,
      m0 = m0,
      m1 = m1,
      alpha = alpha
    )
    results$rg <- tibble(
      method = "Rogan-Gladen",
      theta_hat = rg_est$theta,
      ci_lower = rg_est$ci_lower,
      ci_upper = rg_est$ci_upper
    )

    # 5. EIF
    eif_est <- eif_point_and_ci(
      Y_cal = Y_cal,
      Yhat_cal = Yhat_cal,
      Yhat_test = Yhat_test,
      alpha = alpha
    )
    results$eif <- tibble(
      method = "EIF",
      theta_hat = eif_est$theta,
      ci_lower = eif_est$ci_lower,
      ci_upper = eif_est$ci_upper
    )

    # 6. Joint MLE
    mle_result <- tryCatch({
      mle_fit <- mle_point_and_ci(
        Y_cal = Y_cal,
        Yhat_cal = Yhat_cal,
        Yhat_test = Yhat_test,
        alpha = alpha
      )
      tibble(
        method = "Joint MLE",
        theta_hat = mle_fit$theta,
        ci_lower = mle_fit$ci_lower,
        ci_upper = mle_fit$ci_upper
      )
    }, error = function(e) {
      tibble(
        method = "Joint MLE",
        theta_hat = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_
      )
    })
    results$mle <- mle_result
  }

  bind_rows(results) %>%
    mutate(
      theta_true = theta_true,
      N = N,
      m = m,
      n = n,
      m0 = m0,
      m1 = m1,
      q0_hat = q0_hat,
      q1_hat = q1_hat
    )
}


# Main analysis
main <- function() {
  results_dir <- "data/judge_results"
  output_dir <- "results/llm_judge_analysis"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Get all result files (JSON format)
  files <- list.files(results_dir, pattern = "\\.json$", full.names = TRUE)
  files <- files[!grepl("checkpoint", files)]  # Exclude checkpoints

  if (length(files) == 0) {
    stop("No judge result files found in ", results_dir)
  }

  all_results <- list()

  for (f in files) {
    fname <- basename(f)
    # Parse filename: dataset_judge.json
    parts <- str_match(fname, "(.+)_(gpt-.+)\\.json")
    if (is.na(parts[1])) next

    dataset_name <- parts[2]
    judge_model <- parts[3]

    message("Processing: ", dataset_name, " with ", judge_model)

    df <- fromJSON(f) %>% as_tibble()

    # Apply estimators
    est_results <- apply_all_estimators(df, label_ratio = LABEL_RATIO, alpha = ALPHA)

    est_results <- est_results %>%
      mutate(
        dataset = dataset_name,
        judge = judge_model
      )

    all_results[[fname]] <- est_results
  }

  combined <- bind_rows(all_results)

  # Save results
  write_csv(combined, file.path(output_dir, "estimator_results.csv"))

  # Create plots for each dataset
  datasets <- unique(combined$dataset)

  for (ds in datasets) {
    ds_data <- combined %>% filter(dataset == ds)

    # Get true value
    theta_true <- ds_data$theta_true[1]

    # Plot with large fonts, no subtitle
    p <- ggplot(ds_data, aes(x = method, y = theta_hat, color = method)) +
      geom_point(size = 4) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, linewidth = 1.2) +
      geom_hline(yintercept = theta_true, linetype = "solid", color = "black", linewidth = 1) +
      facet_wrap(~judge, nrow = 1) +
      scale_color_manual(values = METHOD_COLORS) +
      labs(
        title = paste("Win Rate Estimation:", gsub("_", " ", ds)),
        x = "Method",
        y = "Estimated Win Rate (Model A)",
        color = "Method"
      ) +
      theme_bw(base_size = 18) +
      theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        strip.text = element_text(size = 16),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20)
      ) +
      coord_cartesian(ylim = c(0, 1))

    ggsave(file.path(output_dir, paste0(ds, "_ci_comparison.png")),
           p, width = 10, height = 6, dpi = 300)

    message("Saved plot for ", ds)
  }

  # Summary table
  summary_table <- combined %>%
    select(dataset, judge, method, theta_hat, ci_lower, ci_upper, theta_true, q0_hat, q1_hat) %>%
    mutate(
      ci_width = ci_upper - ci_lower,
      covers_true = theta_true >= ci_lower & theta_true <= ci_upper,
      bias = theta_hat - theta_true
    )

  print(summary_table)
  write_csv(summary_table, file.path(output_dir, "summary_table.csv"))

  message("\nResults saved to ", output_dir)
}

#' Run coverage simulation over many random splits
#'
#' @param n_sims Number of simulation repetitions
run_coverage_simulation <- function(n_sims = 100) {
  results_dir <- "data/judge_results"
  output_dir <- "results/llm_judge_analysis"

  files <- list.files(results_dir, pattern = "\\.json$", full.names = TRUE)
  files <- files[!grepl("checkpoint", files)]

  all_sim_results <- list()

  for (f in files) {
    fname <- basename(f)
    parts <- str_match(fname, "(.+)_(gpt-.+)\\.json")
    if (is.na(parts[1])) next

    dataset_name <- parts[2]
    judge_model <- parts[3]

    message("Running ", n_sims, " simulations for: ", dataset_name, " + ", judge_model)

    df <- fromJSON(f) %>% as_tibble() %>%
      filter(!is.na(judge_pick)) %>%
      mutate(
        Y = as.integer(human_winner == "model_a"),
        Yhat = as.integer(judge_pick == "A")
      )

    theta_true <- mean(df$Y)
    N <- nrow(df)
    m <- max(4L, round(LABEL_RATIO * N))
    n <- N - m

    sim_results <- map_dfr(seq_len(n_sims), function(sim_id) {
      # Random split
      idx_cal <- sample(N, m)
      idx_test <- setdiff(seq_len(N), idx_cal)

      Y_cal <- df$Y[idx_cal]
      Yhat_cal <- df$Yhat[idx_cal]
      Yhat_test <- df$Yhat[idx_test]

      m0 <- sum(Y_cal == 0)
      m1 <- sum(Y_cal == 1)

      if (m0 < 2 || m1 < 2) return(NULL)

      q0_hat <- mean(Yhat_cal[Y_cal == 0] == 0)
      q1_hat <- mean(Yhat_cal[Y_cal == 1] == 1)
      p_hat <- mean(Yhat_test)

      results <- list()

      # Naive
      naive_theta <- p_hat
      naive_var <- p_hat * (1 - p_hat) / n
      z_alpha <- qnorm(1 - ALPHA / 2)
      results$naive <- tibble(
        method = "Naive",
        theta_hat = naive_theta,
        ci_lower = pmax(naive_theta - z_alpha * sqrt(naive_var), 0),
        ci_upper = pmin(naive_theta + z_alpha * sqrt(naive_var), 1)
      )

      # PPI
      tryCatch({
        ppi_est <- ppi_point_and_ci(Y_L = Y_cal, f_L = Yhat_cal, f_U = Yhat_test, alpha = ALPHA)
        results$ppi <- tibble(method = "PPI", theta_hat = ppi_est$theta,
                              ci_lower = ppi_est$ci_lower, ci_upper = ppi_est$ci_upper)
      }, error = function(e) NULL)

      # PPI++
      tryCatch({
        ppi_pp_est <- ppi_pp_point_and_ci(Y_L = Y_cal, f_L = Yhat_cal, f_U = Yhat_test, alpha = ALPHA)
        results$ppi_pp <- tibble(method = "PPI++", theta_hat = ppi_pp_est$theta,
                                  ci_lower = ppi_pp_est$ci_lower, ci_upper = ppi_pp_est$ci_upper)
      }, error = function(e) NULL)

      # Rogan-Gladen
      tryCatch({
        rg_est <- rg_point_and_ci(p_hat = p_hat, q0_hat = q0_hat, q1_hat = q1_hat,
                                  n = n, m0 = m0, m1 = m1, alpha = ALPHA)
        results$rg <- tibble(method = "Rogan-Gladen", theta_hat = rg_est$theta,
                             ci_lower = rg_est$ci_lower, ci_upper = rg_est$ci_upper)
      }, error = function(e) NULL)

      # EIF
      tryCatch({
        eif_est <- eif_point_and_ci(Y_cal = Y_cal, Yhat_cal = Yhat_cal, Yhat_test = Yhat_test, alpha = ALPHA)
        results$eif <- tibble(method = "EIF", theta_hat = eif_est$theta,
                              ci_lower = eif_est$ci_lower, ci_upper = eif_est$ci_upper)
      }, error = function(e) NULL)

      # Joint MLE
      tryCatch({
        mle_fit <- mle_point_and_ci(Y_cal = Y_cal, Yhat_cal = Yhat_cal, Yhat_test = Yhat_test, alpha = ALPHA)
        if (!is.na(mle_fit$ci_lower)) {
          results$mle <- tibble(method = "Joint MLE", theta_hat = mle_fit$theta,
                                ci_lower = mle_fit$ci_lower, ci_upper = mle_fit$ci_upper)
        }
      }, error = function(e) NULL)

      bind_rows(results) %>% mutate(sim_id = sim_id)
    })

    # Compute coverage with SE
    coverage_summary <- sim_results %>%
      mutate(
        covers = theta_true >= ci_lower & theta_true <= ci_upper,
        ci_width = ci_upper - ci_lower
      ) %>%
      group_by(method) %>%
      summarise(
        n_sims = n(),
        coverage = mean(covers, na.rm = TRUE),
        mean_ci_width = mean(ci_width, na.rm = TRUE),
        sd_ci_width = sd(ci_width, na.rm = TRUE),
        mean_bias = mean(theta_hat - theta_true, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        # SE for coverage (binomial proportion)
        se_coverage = sqrt(coverage * (1 - coverage) / n_sims),
        # SE for CI width
        se_ci_width = sd_ci_width / sqrt(n_sims),
        dataset = dataset_name,
        judge = judge_model,
        theta_true = theta_true
      )

    all_sim_results[[fname]] <- coverage_summary
  }

  combined_coverage <- bind_rows(all_sim_results)

  # Save CSV
  write_csv(combined_coverage, file.path(output_dir, "coverage_simulation.csv"))

  # Create coverage vs CI width plot (Pareto frontier)
  combined_coverage <- combined_coverage %>%
    mutate(
      dataset_short = gsub("claude-opus-4-20250514_vs_", "", dataset),
      label = paste0(dataset_short, "\n", judge)
    )

  p_pareto <- ggplot(combined_coverage, aes(x = mean_ci_width, y = coverage, color = method)) +
    geom_errorbar(aes(ymin = pmax(coverage - 1.96 * se_coverage, 0),
                      ymax = pmin(coverage + 1.96 * se_coverage, 1)),
                  width = 0, linewidth = 0.5, alpha = 0.7) +
    geom_errorbarh(aes(xmin = pmax(mean_ci_width - 1.96 * se_ci_width, 0),
                       xmax = mean_ci_width + 1.96 * se_ci_width),
                   height = 0, linewidth = 0.5, alpha = 0.7) +
    geom_point(size = 4) +
    geom_hline(yintercept = 1 - ALPHA, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_text(aes(label = paste0(method, " (", round(100*coverage), "%)")),
              hjust = -0.05, vjust = 0.5, size = 2.3, show.legend = FALSE) +
    facet_wrap(~ label, ncol = 3) +
    scale_color_manual(values = METHOD_COLORS) +
    scale_y_continuous(limits = c(-0.05, 1.05), labels = scales::percent,
                       breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(limits = c(0, 1.15)) +
    labs(
      x = "Mean CI Width",
      y = "Empirical Coverage",
      color = "Method"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      strip.text = element_text(size = 9, face = "bold")
    )

  ggsave(file.path(output_dir, "coverage_vs_width.png"), p_pareto, width = 14, height = 10, dpi = 300)

  message("\nCoverage simulation saved to ", output_dir)
  print(combined_coverage %>% select(dataset, judge, method, coverage, mean_ci_width))
}

main()
run_coverage_simulation(n_sims = 1000)
