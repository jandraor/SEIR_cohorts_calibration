accuracy_metrics <- function(sim_data, real_data, age_groups) {
  comparison_data <- bind_rows(sim_data, real_data)
  
  accuracy_per_ag <- purrr::map(age_groups, function(ag, comparison_data) {
    
    ag_data  <- comparison_data %>% filter(cohort == ag)
    pred_y   <- ag_data %>% filter(source == "syn data") %>% pull(y)
    actual_y <- ag_data %>% filter(source == "sim data") %>% pull(y)
    
    list(cohort = ag,
         MSE = mse(pred_y, actual_y),
         MASE = mase(pred_y, actual_y))
    
  }, comparison_data = comparison_data)
  
  avg_MSE  <- purrr::map_dbl(accuracy_per_ag, "MSE") %>% mean()
  avg_MASE <- purrr::map_dbl(accuracy_per_ag, "MASE") %>% mean()
  
  list(avg_MSE  = avg_MSE,
       avg_MASE = avg_MASE)
}