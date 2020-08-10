summarise_results <- function(fit, conceptual_matrix, incidence_df, pop_df, 
                              specs_list, unique_params, actual_K,
                              scenario = "perfect information") {
  
  age_groups <- pop_df$group
  
  posterior_df            <- as.data.frame(fit)
  samples                 <- rstan::extract(fit)
  param_samples           <- samples$params
  colnames(param_samples) <- unique_params
  
  rho_df <- NULL
  
  if(scenario == "underreporting") {
    rho_df <- as.data.frame(samples$rho)
    colnames(rho_df) <- "rho"
  }
  
  summary_k <- summarise_k(param_samples, conceptual_matrix, age_groups,
                           rho_df = rho_df, actual_K) 
  
  sim_data     <- extract_incidences(posterior_df, age_groups)
  real_data    <- incidence_df %>% mutate(source = "syn data") 
  g_comparison <- g_compare_ts(sim_data, real_data, intervals = TRUE)
  
  metrics <- accuracy_metrics(sim_data, real_data, age_groups)
  
  log_lik <- mean(posterior_df$log_lik)
  
  ts_summary <- list(g_comparison = g_comparison,
                     MSE          = metrics$avg_MSE,
                     MASE         = metrics$avg_MASE,
                     log_lik      = log_lik)
  
  R_0_obj <- posterior_R0(param_samples, conceptual_matrix, pop_df)
  
  output <- list(
    summary_k    = summary_k,
    ts_summary   = ts_summary,
    R_0          = R_0_obj)
  
  if(scenario == "underreporting") {
    
    bounds  <- quantile(rho_df$rho, c(0.025, 0.975))
    rho_lb  <- bounds[[1]] # lower bound
    rho_ub  <- bounds[[2]] # upper bound
    
    rho_list <- list(mean        = mean(rho_df$rho) ,
                     lower_bound = rho_lb,
                     upper_bound = rho_ub)
    
    output$rho_hat <- rho_list
  }
  
  output
}

extract_incidences <- function(posterior_df, age_groups) {
  
  incidence_sims <- imap_dfr(age_groups, function(ag, index, posterior_df) {
    y_var <- str_glue("incidence{index}")
    extract_timeseries_var(y_var, posterior_df) %>% 
      mutate(cohort = ag) %>% 
      dplyr::select(-variable)
  }, posterior_df = posterior_df)
  
  sim_data <- incidence_sims %>% group_by(cohort, time) %>% 
    summarise(lower_bound = quantile(value, c(0.025, 0.975))[[1]],
              upper_bound = quantile(value, c(0.025, 0.975))[[2]],
              y = mean(value),
              .groups = "drop") %>% mutate(source = "sim data")
  
  sim_data
}

summarise_k <- function(param_samples, conceptual_matrix, age_groups,
                        rho_df = NULL, actual_K) {
  
  param_means <- apply(param_samples, 2, function(col) mean(col))
  
  # Matrix K means
  nc    <- length(age_groups)
  k_hat       <- get_mean_k_hat(param_samples, conceptual_matrix, nc, 
                                age_groups)
  
  interval_df <- get_k_intervals(param_samples, conceptual_matrix)
  
  g_k_hat <- draw_WAIFW(k_hat, "", interval_df, precision = 1)
  
  param_samples_df <- as.data.frame(param_samples)
  
  if(!is.null(rho_df)) {
    param_samples_df <- bind_cols(param_samples_df, rho_df)
  }
  g_pairs          <- pairs_posterior(param_samples_df)
  
  MSE_k <- mse(k_hat, actual_K)
  
  list(k_hat   = k_hat,
       MSE_k   = MSE_k,
       g_k_hat = g_k_hat,
       g_pairs = g_pairs)
}

posterior_R0 <- function(param_samples, conceptual_matrix, pop_df,
                         tau_I = 2) {
  
  nc <- nrow(pop_df)
  
  r_noughts <- sapply(1:nrow(param_samples), function(i) {
    row      <- param_samples[i, ]
    
    
    k_matrix <- sapply(conceptual_matrix, 
                    function(kij, row) row[[kij]], row = row) %>%
      matrix(nrow = nc, byrow = TRUE)
    
    R0_from_K(k_matrix, pop_df, tau_I)
  })
  
  credible_interval <- quantile(r_noughts, c(0.025, 0.975))
  
  g_rNougths <- draw_density(r_noughts, specs_list)
  
  list(g = g_rNougths,
       mean        = mean(r_noughts),
       lower_bound = credible_interval[[1]],
       upper_bound = credible_interval[[2]])
}

consolidate_R0 <-function(smr_list) {
  
  data.frame(R0          = smr_list$R_0$mean, 
             lower.bound = smr_list$R_0$lower_bound, 
             upper.bound = smr_list$R_0$upper_bound)
}

extract_K_error <- function(smr_list) smr_list$summary_k$MSE_k

get_mean_k_hat <- function(param_samples, cm, nc, age_groups) {
  
  param_means <- apply(param_samples, 2, function(col) mean(col))
  
  k_hat <- sapply(cm, function(kij, row) row[[kij]], row = param_means) %>%
    matrix(nrow = nc, byrow = T)
  
  colnames(k_hat) <- rownames(k_hat) <- age_groups
  
  k_hat
}

get_k_intervals <- function(param_samples, conceptual_matrix) {
  map_df(conceptual_matrix, function(kij, param_samples) {
    
    vals              <- param_samples[, kij]
    credible_interval <- quantile(vals, c(0.025, 0.975))
    
    data.frame(lower.interval = round(credible_interval[1], 1), 
               upper.interval = round(credible_interval[2], 1))
    
  }, param_samples = param_samples)
}


  


