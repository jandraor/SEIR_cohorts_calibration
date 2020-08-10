draw_multiple_densities <- function(densities_df, stats_df,
                                    x_pos, text_size = 5) {
  g <- ggplot(densities_df, aes(x = x, y = y)) +
    geom_line() +
    geom_ribbon(aes(ymin = 0, ymax = y, fill = area)) +
    facet_wrap(~ param, scales = "free") +
    scale_fill_manual(values = c(NA, "lightgrey"))  +
    geom_vline(aes(xintercept = mean_value), stats_df, color = "#1261A0",
               linetype ="dashed", size = 1) +
    geom_text(aes(x = x_pos, y = y_pos_mean, label = mean_label),
              stats_df, colour = "#1261A0", size = text_size) +
    geom_text(aes(x = x_pos, y = y_pos_median, label = median_label),
              stats_df, colour = "#1261A0", size = text_size) +
    geom_text(aes(x = x_pos, y = y_pos_interval, label = interval_label),
               stats_df, colour = "grey", size = text_size) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
}

generate_graph_inputs <- function(stan_fit, x_pos, y_pos, y_pos_median, 
                                  y_pos_interval, pars, rename_pars = NULL) {
  
  posterior_df <- as.data.frame(stan_fit) 
  
  for(rename_list in rename_pars) {
    posterior_df <- posterior_df %>% 
      rename(!!rename_list$new := !!rename_list$old)
  }
  
  if("recoveryTime" %in% pars) {
    posterior_df <- mutate(posterior_df, recoveryTime = 1 / recoveryProportion)
  }
  
  if("latent_period" %in% pars) {
    posterior_df <- mutate(posterior_df, 
                           latent_period = 1 / incubation_proportion)
  }
     
  params_df <- select(posterior_df, pars)
  
  credible_intervals <- apply(params_df, 2, HPDI, prob = 0.95) %>% t() %>% 
    as.data.frame() %>% 
    rename(lower_interval = "|0.95", upper_interval = "0.95|") %>% 
    mutate(param = rownames(.))
  
  means   <- apply(params_df, 2, mean)
  medians <- apply(params_df, 2, median)
  
  stats_df <- data.frame(stringsAsFactors = FALSE,
                         param = names(means), 
                         mean_value = means,
                         median_value = medians,
                         lower_interval = credible_intervals$lower_interval,
                         upper_interval = credible_intervals$upper_interval) %>% 
    mutate(mean_label     = paste0("mean   = ", round(mean_value, 3)),
           median_label   = paste0("median = ", round(median_value, 3)),
           interval_label = paste0("[ ", round(lower_interval, 3), ", ", 
                                   round(upper_interval, 3), " ]"),
           x_pos          = x_pos,
           y_pos_mean     = y_pos,
           y_pos_median   = y_pos_median,
           y_pos_interval = y_pos_interval)
  
  densities <- apply(params_df, 2, density) %>% lapply(function(densityObj){
    data.frame(x = densityObj$x, y = densityObj$y)
  }) %>% bind_rows(.id = "param") %>% 
    inner_join(credible_intervals) %>% 
    mutate(area = x >= lower_interval & x <= upper_interval)
  
  list(densities = densities,
       stats_df  = stats_df)
}

draw_density <- function(data_vector, g_params) {
  credible_interval <- quantile(data_vector, c(0.025, 0.975))
  mean_param        <- mean(data_vector)
  median_param      <- median(data_vector)
  
  hist.y <- density(data_vector) %$% 
    data.frame(x = x, y = y) %>% 
    mutate(area = x >= credible_interval[1] & x <= credible_interval[2])
  
  g1 <- ggplot(hist.y, aes(x = x)) + 
    geom_line(aes(y = y)) +
    geom_text(aes(x = g_params$x_pos, y = g_params$ypos_mean), 
              label = paste0("mean   = ", round(mean_param, 4)),
              colour = "#1261A0", size = g_params$text_size) +
    geom_text(aes(x = g_params$x_pos, y = g_params$ypos_median), 
              label = paste0("median = ", round(median_param, 4)),
              colour = "#009999", size = g_params$text_size) +
    annotate("text", x = g_params$x_pos, y = g_params$ypos_interval, 
             label = paste0("[",
                            round(credible_interval[1], 4),",",
                            round(credible_interval[2], 4), "]"),
             size = g_params$text_size,
             colour = "grey") +
    geom_ribbon(aes(ymin = 0, ymax = y, fill = area)) +
    scale_fill_manual(values = c(NA, "lightgrey")) + 
    geom_vline(aes(xintercept = mean_param),
               color = "#1261A0", linetype ="dashed", size = 1) +
    geom_vline(aes(xintercept = median_param),
               color = "#009999", linetype ="dashed", size = 1) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(color = "#404040", size = 8)) +
    labs(x = g_params$xlabel, title = g_params$title)
}

draw_WAIFW <- function(WAIFW, subtitle, interval_df = NULL,
                       precision = 0) {
  library(reshape2)

  WAIFW_df <- WAIFW %>% t() %>% melt()
  
  if(!is.null(interval_df)) {
    WAIFW_df <- bind_cols(WAIFW_df, interval_df)
  }
  
  g <- ggplot(data = WAIFW_df, aes(x = Var1, 
                              y = ordered(Var2, levels = rev(sort(unique(Var2)))), 
                              fill = value)) + 
    geom_tile() +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    geom_text(aes(label = round(value, precision)), colour = "white", 
              size = 2) +
    theme_minimal() + 
    labs(y ="", x = "",subtitle = subtitle) +
    theme(legend.position = "none",
          plot.subtitle = element_text(color = "#404040", size = 8),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6))
  
  if(!is.null(interval_df)) {
    g <- g + geom_text(
      aes(label = paste("[", lower.interval," ,", upper.interval, "]")),
      nudge_y = -0.2, size = 1.5, colour = "white")
  }
  
  g
}

ts_errors <- function(summaries_optim, inits) {
  metrics_list <- purrr::map(summaries_optim, "metrics")
  
  MASEs    <- metrics_list %>% map_dbl("avg_MASE")
  MSEs     <- metrics_list %>% map_dbl("avg_MSE")
  log_liks <- metrics_list %>% map_dbl("log_lik")
  
  df_params_list <- list(
    list(label = "MASE",
         vals  = MASEs),
    list(label = "MSE",
         vals  = MSEs),
    list(label = "Log lik",
         vals  = log_liks))
  
  metric_df <- map_df(df_params_list, function(df_params, inits) {
    
    if(df_params$label != "Log lik") var_best <- min(df_params$vals)
    
    if(df_params$label == "Log lik") var_best <- max(df_params$vals)
    
    
    
    tibble(init = inits, value = df_params$vals) %>% 
      mutate(is_Best = value == var_best,
             metric = df_params$label)
  }, inits = inits)
  
  MASE_df <- metric_df %>% filter(metric == "MASE") %>% 
    arrange(desc(value))
  
  metric_df <- mutate(metric_df, init = factor(init, levels = MASE_df$init))
  
  ggplot(metric_df, aes(x = init, y = value)) +
    facet_wrap(~ metric, scales = "free", nrow = 1) +
    coord_flip() +
    geom_lollipop(aes(colour = is_Best)) +
    scale_colour_manual(values = c("grey", "steelblue")) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x = "Init id",
         y = "Value",
         title = "Predicted incidence accuracy") 
}

draw_inits_comparison <- function(summaries_optim, actual_R0, inits) {
  
  MSE_K  <- map_dbl(summaries_optim, "MSE_K")
  MSE_R0     <- map_dbl(summaries_optim, function(summary, actual_R0) {
    MSE(actual_R0, summary$R_nought)
  }, actual_R0 = actual_R0) 
  
  df_params_list <- list(
    list(label = "K (MSE)",
         vals  = MSE_K),
    list(label = "R0 (MSE)",
         vals  = MSE_R0))
  
  MSEs_df <- map_df(df_params_list, function(df_params, inits) {
    var_min <- min(df_params$vals)
    
    tibble(init = inits, MSE = df_params$vals) %>% 
      mutate(is_Min = MSE == var_min,
             variable = df_params$label)
  }, inits = inits)
  
  ggplot(MSEs_df, aes(x = as.factor(init), y = MSE)) +
    facet_wrap(~ variable, scales = "free", nrow = 1) +
    coord_flip() +
    geom_lollipop(aes(colour = is_Min)) +
    scale_colour_manual(values = c("grey", "steelblue")) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x = "Init id",
         y = "Error") 
  

}

# Draw distance comparison graph
draw_dcg <- function(df, limits, actual_val) {
  
  g <- ggplot(df, aes(x = x, y = y)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower.bound, ymax = upper.bound), width =.1) +
    geom_hline(yintercept = actual_val, linetype = "dashed") +
    scale_y_continuous(limits = limits) +
    facet_wrap(~ method) +
    theme_test() +
    theme(legend.text  = element_text(size = 3)) +
    labs(x = "Structure", y = "Reporting probability")
}

# Compare time-series to data points

g_compare_ts <- function(sim_data, real_data, intervals = TRUE,
                         scales = "fixed", xlabel = "Days") {
  g <- ggplot(sim_data, aes(x = time, y = y)) +
    geom_line(colour = "steelblue", alpha = 0.9, size = 0.25) +
    geom_point(data = real_data, size = 0.5, colour = "grey30",
               alpha = 0.8) +
    scale_y_continuous(labels = comma) +
    facet_wrap(~ cohort, scales = scales) 
  
  if(isTRUE(intervals)) {
    g <- g + geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound),
                         alpha = 0.5, fill = "steelblue")
  }
  
  g <- g + 
    labs(x = xlabel, y = "Incidence") +
    theme_pubr()
  
  g
}

# ===================Pairs======================================================

dens_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    stat_density2d(aes(fill=..density..), geom = "tile", contour = FALSE) +
    scale_fill_viridis_c()
  p
}

cor_fun <- function(data, mapping, method = "pearson", ndp = 2, sz=5, 
                    stars=TRUE, ...){
  
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  corr <- cor.test(x, y, method=method)
  est <- corr$estimate
  lb.size <- sz* abs(est) 
  
  palette <- gradient_n_pal(c("lightgrey", "black"))
  
  if(stars){
    stars <- c("***", "**", "*", "")[findInterval(corr$p.value, 
                                                  c(0, 0.001, 0.01, 0.05, 1))]
    
    lbl  <- paste0(round(est, ndp), stars)
    cor_colour <-  palette(abs(est))
  }else{
    lbl <- round(est, ndp)
  }
  
  ggplot(data = data, mapping = mapping) + 
    annotate("text", x = mean(x, na.rm = TRUE), y = mean(y, na.rm=TRUE), 
             label=lbl, size = 3, colour = cor_colour,...)+
    theme(panel.grid = element_blank())
}

pairs_posterior <- function(posterior, strip_text = 3) {
  ggpairs(posterior, lower = list(continuous = dens_fn),
          upper = list(continuous = cor_fun)) +
    theme_pubr() +
    theme(axis.text = element_text(size = 4),
          strip.text = element_text(size = strip_text))
}

#===============================================================================

g_time_comparison <- function(t_df) {
  ggplot(t_df, aes(x = matrix, y = time)) +
    geom_lollipop(colour = "steelblue") +
    scale_y_continuous() +
    coord_flip() +
    geom_text(aes(label = round(time, 0)), nudge_y = 10, size = 3) +
    scale_colour_manual(values = c("grey", "steelblue")) +
    facet_grid(scenario ~ method) +
    theme_test() +
    theme(legend.position = "none") +
    labs(x = "Structure", y = "Time [Minutes]",
         title = "Run time")
}




