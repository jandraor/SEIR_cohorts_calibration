summarise_optim_fit <- function(
  optim_fit, plot_title = "", incidence_data, conceptual_matrix, pop_cohorts, 
  actual_WAIFW, g_options = NULL) {
  
  age_groups <- c("00-04", "05-14", "15-44", "45+")
  
  #=============================================================================
  WAIFW_analysis <- run_WAIFW_analysis(optim_fit, conceptual_matrix, age_groups, 
                               actual_WAIFW, plot_title, pop_cohorts)
  
  R_nought <- WAIFW_analysis$R_nought
  #=============================================================================
  
  translation_df <- tibble(
    index_group = 1:4,
    var_incidence = paste0("incidence", 1:4),
    age_group = age_groups)
  
  actual_df <- incidence_data %>% left_join(translation_df) %>% 
    select(-var_incidence, -index_group) %>% mutate(source = "syn_data")
  
  comparison_df <- bind_rows(actual_df, optim_fit$fit_df)
  
  #=============================================================================
  
  MSE_per_ag <- map_dbl(age_groups, function(ag, comparison_df) {
    
    ag_data  <- comparison_df %>% filter(age_group == ag)
    syn_data <- ag_data %>% filter(source == "syn_data") %>% pull(incidence)
    sim_data <- ag_data %>% filter(source == "sim data") %>% pull(incidence)
    MSE(sim_data, syn_data)
    
  }, comparison_df = comparison_df)
  #=============================================================================
  rounded_R0 <- round(R_nought, 2)
  
  g_comparison <- ggplot(comparison_df, aes(x = time, y = incidence)) +
    geom_line(aes(group = source, colour = source)) +
    scale_colour_manual(values = c("steelblue", "grey")) +
    facet_wrap(~age_group) +
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
    WAIFW = WAIFW_analysis$WAIFW_estimates,
    g_WAIFW = WAIFW_analysis$g_WAIFW,
    MSE_WAIFW = WAIFW_analysis$MSE_WAIFW,
    MSE = sum(MSE_per_ag),
    R_nought = WAIFW_analysis$R_nought)
  
  params_names <- names(optim_fit$fit$par)
  
  if("p" %in% params_names) {
    output$p_hat <- optim_fit$fit$par[["p"]]
  }
  
  output
}

run_WAIFW_analysis <- function(optim_fit, conceptual_matrix, age_groups, 
                               actual_WAIFW, plot_title, pop_cohorts) {
  
  WAIFW_estimates <- sapply(conceptual_matrix, function(Bij, row) row[[Bij]],
                            row = optim_fit$fit$par) %>%
    matrix(nrow = 4, byrow = T)
  
  colnames(WAIFW_estimates) <- rownames(WAIFW_estimates) <- age_groups
  
  MSE_WAIFW <- MSE(as.vector(WAIFW_estimates), actual_WAIFW)
  
  g_WAIFW <- draw_WAIFW(WAIFW_estimates, plot_title)
  
  #  2 is the recovery delay
  next_gen_matrix <- WAIFW_estimates / 1e5 * pop_cohorts * 2 
  eigensystem     <- eigen(next_gen_matrix)
  R_nought        <- max(Re(eigensystem$values))
  
  
  list(R_nought        = R_nought,
       WAIFW_estimates = WAIFW_estimates,
       g_WAIFW         = g_WAIFW,
       MSE_WAIFW       = MSE_WAIFW)
}