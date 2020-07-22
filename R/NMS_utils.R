
get_matrix_inits <- function(n_inits, conceptual_matrix) {
  
  unique_elements <- unique(conceptual_matrix)
  
  init_matrix_list <- lapply(seq_len(n_inits), function(i) {
    random_inits <- runif(length(unique_elements), 0, 10) %>% round(1)
    names(random_inits) <- unique_elements
    random_inits
  })
  
  init_matrix_list
}

optimise_SEIR2 <- function(inits, matrix_type, data_list, sim_specs, 
                           var_cache) {
  
  filepath <- str_glue("./object_fits/optim/optim_det_{matrix_type}.rds")
  
  if(var_cache == FALSE) {
    
    output      <- optimise_SEIR(inits, matrix_type, data_list, sim_specs)
    saveRDS(output, filepath)
  } else {
    output    <- readRDS(filepath)
  }
  
  output
}
