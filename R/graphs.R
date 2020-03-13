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
  credible_interval <- HPDI(data_vector, prob = 0.95)
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

draw_WAIFW <- function(WAIFW, subtitle, interval_df = NULL) {
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
    geom_text(aes(label = round(value)), colour = "white", size = 2) +
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

draw_inits_comparison <- function(summaries_optim, actual_R0, inits) {
  
  MSEs       <- map_dbl(summaries_optim, "MSE")
  MSE_WAIFW  <- map_dbl(summaries_optim, "MSE_WAIFW")
  MSE_R0     <- map_dbl(summaries_optim, function(summary, actual_R0) {
    MSE(actual_R0, summary$R_nought)
  }, actual_R0 = actual_R0) 
  
  df_params_list <- list(
    list(label = "Timeseries",
         vals  = MSEs),
    list(label = "WAIFW",
         vals  = MSE_WAIFW),
    list(label = "R0",
         vals  = MSE_R0))
  
  MSEs_df <- map_df(df_params_list, function(df_params, inits) {
    var_min <- min(df_params$vals)
    
    tibble(init = inits, MSE = df_params$vals) %>% 
      mutate(is_Min = MSE == var_min,
             variable = df_params$label)
  }, inits = inits)
  
  g_MSEs <- ggplot(MSEs_df, aes(x = as.factor(init), y = MSE)) +
    facet_wrap(~ variable, scales = "free", nrow = 1) +
    coord_flip() +
    geom_lollipop(aes(colour = is_Min)) +
    scale_colour_manual(values = c("grey", "steelblue")) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x = "Init values")
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


