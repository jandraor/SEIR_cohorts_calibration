summarise_optim_fit <- function(
  optim_fit, plot_title = "", incidence_data, conceptual_matrix, pop_df, 
  actual_K, g_options = NULL) {
  
  incidence_data <- incidence_data %>% mutate(cohort = as.character(cohort))
  
  translation_df <- incidence_data %>% dplyr::select(index, cohort) %>% unique()
  
  age_groups <- translation_df$cohort
  
  #=============================================================================
  K_analysis <- run_K_analysis(optim_fit, conceptual_matrix, age_groups, 
                               actual_K, plot_title, pop_df)
  
  R_nought <- K_analysis$R_nought
  #=============================================================================
  
  actual_df <- incidence_data %>% mutate(source = "syn data")
  
  sym_df    <- optim_fit$fit_df %>% left_join(translation_df, by = "index")
  
  comparison_df <- bind_rows(actual_df, sym_df)
    
  
  #=============================================================================
  metrics         <- accuracy_metrics(sym_df, actual_df, age_groups)
  metrics$log_lik <- optim_fit$fit$value * -1 # because of optim minimisation
  #=============================================================================
  rounded_R0 <- round(R_nought, 2)
  
  g_comparison <- ggplot(comparison_df, aes(x = time, y = y)) +
    geom_line(aes(group = source, colour = source)) +
    scale_colour_manual(values = c("steelblue", "grey")) +
    facet_wrap(~cohort) +
    theme_test() +
    labs(title = bquote(list(.(plot_title), R[0]== .(rounded_R0))),
         x = "Day", 
         y = "Incidence [New cases per day]")
  
  if(!is.null(g_options)) {
    
    if(g_options$stripSize) {
      g_comparison <- g_comparison +
        theme(strip.text.x = element_text(size = g_options$stripSize))
    }
    
    if(!is.null(g_options[["legendPosition"]])) {
      g_comparison <- g_comparison +
        theme(legend.position = g_options$legendPosition)
    }
    
    if(!is.null(g_options[["axisTextXSize"]])) {
      g_comparison <- g_comparison +
        theme(axis.text.x = element_text(size = g_options$axisTextXSize))
    }
    
    if(!is.null(g_options[["axisTextYSize"]])) {
      g_comparison <- g_comparison +
        theme(axis.text.y = element_text(size = g_options$axisTextYSize))
    }
    axis.title=element_text(size=14,face="bold")
    
    if(!is.null(g_options[["axisTitleSize"]])) {
      g_comparison <- g_comparison +
       theme(axis.title = element_text(size = g_options$axisTitleSize))
    }
    
    if(!is.null(g_options[["legendSize"]])) {
      legend_title <- g_options$legendSize
      legend_text  <- legend_title * 0.75
      g_comparison <- g_comparison +
        theme(legend.title = element_text(size = legend_title),
              legend.text  = element_text(size = legend_text))
    }
    
    if(!is.null(g_options[["legendBoxMargin"]]) & 
       g_options[["legendBoxMargin"]] == FALSE) {
      
      g_comparison <- g_comparison +
        theme(legend.box.margin = margin(-12,-12,-12,-12))
    }
    
    if(!is.null(g_options[["titleSize"]])) {
      g_comparison <- g_comparison +
        theme(plot.title = element_text(size = g_options$titleSize))
    }
  }
  
  output <- list(
    g_comparison = g_comparison,
    K_hat        = K_analysis$K_hat,
    g_K_hat      = K_analysis$g_K_hat,
    SMAPE_K      = K_analysis$SMAPE_K,
    metrics      = metrics,
    R_nought     = K_analysis$R_nought)
  
  params_names <- names(optim_fit$fit$par)
  
  if("rho" %in% params_names) {
    output$rho_hat <- optim_fit$fit$par[["rho"]]
  }
  
  output
}

run_K_analysis <- function(optim_fit, conceptual_matrix, age_groups, 
                               actual_K, plot_title, pop_df) {
  
  K_hat <- sapply(conceptual_matrix, function(kij, row) row[[kij]],
                            row = optim_fit$fit$par) %>%
    matrix(nrow = 4, byrow = T)
  
  colnames(K_hat) <- rownames(K_hat) <- age_groups
  
  SMAPE_K <- smape(as.vector(K_hat), actual_K)
  
  g_K_hat <- draw_WAIFW(K_hat, plot_title, precision = 2)
  
  R_nought <- R0_from_K(K_hat, pop_df)
  
  list(R_nought        = R_nought,
       K_hat           = K_hat,
       g_K_hat         = g_K_hat,
       SMAPE_K         = SMAPE_K)
}

find_best_fits <- function(summaries, n_inits) {
  
  vals   <- purrr::map(summaries, "metrics") %>% 
    map_dbl("log_lik")
  
  names(vals)   <- 1:n_inits
  sorted_vals          <- sort(vals, decreasing = TRUE)
  top_6                <- sorted_vals[1:6]
  top_id               <- which(vals %in% top_6)
  
  list(top_id       = top_id,
       pos_best_fit = which.max(vals))
}
