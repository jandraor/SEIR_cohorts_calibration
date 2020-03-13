summarise_results <- function(fit, conceptual_matrix, incidence_df, pop_sizes, 
                              specs_list, unique_params, cumulative = FALSE,
                              stochastic = FALSE) {
  
  age_groups <- c("00-04", "05-14", "15-44", "45+")
  
  posterior_df <- as.data.frame(fit)
  samples                 <- rstan::extract(fit)
  param_samples           <- samples$params
  colnames(param_samples) <- unique_params
  param_medians           <- apply(param_samples, 2, function(col) median(col))
  
  
  WAIFW_medians <- sapply(
    conceptual_matrix, function(Bij, row) row[[Bij]], row = param_medians) %>%
    matrix(nrow = 4, byrow = T)
  
  normalised_WAIFW           <- WAIFW_medians
  colnames(normalised_WAIFW) <- rownames(normalised_WAIFW) <- age_groups
  
  
  interval_df <- map_df(conceptual_matrix, function(Bij, param_samples) {
   
    vals              <- param_samples[, Bij]
    credible_interval <- rethinking::HPDI(vals, 0.95)
    
    data.frame(lower.interval = round(credible_interval[1], 0), 
               upper.interval = round(credible_interval[2],0 ))
   
  }, param_samples = param_samples)
  
  g_WAIFW <- draw_WAIFW(normalised_WAIFW, "", interval_df)
  
  #=============================================================================
  age_groups <- c("00-04", "05-14", "15-44", "45+")
  translation_df <- data.frame(index_group = 1:4, 
                               age_group = c("00-04", "05-14", "15-44", "45+"))
  
  if(cumulative == FALSE) {
    sim_data <- extract_mean_incidences(posterior_df, stochastic)
  }
  
  if(cumulative == TRUE) {
    init_vals <- get_inits(fit)[[1]]$y0
    
    sim_data <- map2_df(age_groups, c(3, 7, 11, 15), function(ag, i) {
      
      infected_index  <- paste0(",", i,"]")
      recovered_index <- paste0(",", i + 1,"]")
      
      col_means    <- posterior_df %>%
        select(ends_with(infected_index), ends_with(recovered_index)) %>% 
        colMeans()
      
      extracted_df <- str_extract_all(names(col_means), "\\d+") %>% 
        map_df(function(extracted_obj) {
          data.frame(stringsAsFactors = FALSE,
                     time = extracted_obj[[1]], stock_id = extracted_obj[[2]])
        })
      
      tbl_colmeans <- tibble(variable = names(col_means), value = col_means,
                             time = extracted_df$time, 
                             stock_id = extracted_df$stock_id)
      
      incidence_sim <- tbl_colmeans %>% group_by(time) %>% 
        summarise(ever_infected = sum(value)) %>% ungroup() %>% 
        mutate(time = as.numeric(time), age_group = ag)   %>% 
        bind_rows(data.frame(stringsAsFactors = FALSE, 
                             time = 0, ever_infected = init_vals[i], 
                             age_group = ag)) %>% arrange(time) %>% 
        mutate(value = ever_infected - lag(ever_infected, 
                                           default = ever_infected[1])) %>% 
        slice(-1) %>% mutate(value = round(value, 0)) %>% 
        select(-ever_infected)
    }) %>% mutate(source = "sim data")
  }
  
  
  real_data <- incidence_df %>% rename(value = incidence) %>% 
    mutate(source = "syn data") 
  
  if(!"age_group" %in% names(real_data)) {
    real_data <- left_join(real_data, translation_df, by = "index_group")
  }
  
  comparison_data <- bind_rows(sim_data, real_data)
  
  g_comparison <- ggplot(comparison_data, aes(x = time, y = value)) +
    geom_line(aes(group = source, colour = source)) +
    scale_colour_manual(values = c("lightgrey", "blue")) +
    facet_wrap(~ age_group) +
    theme_test()
  
  #=============================================================================
  MSE_per_ag <- map_dbl(age_groups, function(ag, comparison_data) {
    
    ag_data  <- comparison_data %>% filter(age_group == ag)
    syn_data <- ag_data %>% filter(source == "syn data") %>% pull(value)
    sim_data <- ag_data %>% filter(source == "sim data") %>% pull(value)
    MSE(sim_data, syn_data)
    
  }, comparison_data = comparison_data)
  #=============================================================================
 
  r_noughts <- sapply(1:nrow(param_samples), function(i) {
    row <- param_samples[i, ]
    WAIFW <- sapply(conceptual_matrix, 
                    function(Bij, row) row[[Bij]], row = row) %>%
      matrix(nrow = 4, byrow = TRUE)
    next_gen_matrix <- WAIFW / 1e5 * pop_sizes * 2
    eigensystem     <- eigen(next_gen_matrix)
    max(abs(eigensystem$values))
  })
  
  credible_interval <- HPDI(r_noughts, prob = 0.95)
  
  g_rNougths <- draw_density(r_noughts, specs_list)
  
  output <- list(
    g_WAIFW = g_WAIFW,
    WAIFW = normalised_WAIFW,
    g_comparison = g_comparison,
    g_rNougths = g_rNougths,
    mean_rNought = mean(r_noughts),
    lower_bound  = credible_interval[[1]],
    upper_bound  = credible_interval[[2]],
    MSE = sum(MSE_per_ag))
  
  if("p" %in% colnames(posterior_df)) {
    output$p_hat    <- mean(posterior_df$p)
    bounds          <- rethinking::HPDI(posterior_df$p, 0.95)
    output$p_hat_lb <- bounds[[1]] # lower bound
    output$p_hat_ub <- bounds[[2]] # upper bound
  }
  
  output
}

extract_mean_incidences <- function(posterior_df, stochastic = FALSE) {
  age_groups <- c("00-04", "05-14", "15-44", "45+")
  
  incidences      <- posterior_df %>% select(contains("incidence"))
  
  if(stochastic) {
    p_values   <- posterior_df %>% pull(p)
    incidences <- incidences * p_values 
  }
  
  mean_incidences <- colMeans(incidences)
  
  names_mi <- names(mean_incidences)
  
  sim_data <- map_df(age_groups, function(ag) {
    i <- which(ag == age_groups)
    regex <- stringr::regex(paste0("incidence", i, "\\["))
    data.frame(time = 1:50, value = mean_incidences[grepl(regex, names_mi)], 
               age_group = ag, stringsAsFactors = F)
  }) %>% mutate(source = "sim data")
}
