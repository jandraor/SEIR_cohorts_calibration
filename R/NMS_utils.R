
get_matrix_inits <- function(n_inits, conceptual_matrix, 
                             reporting_fraction = FALSE) {
  
  unique_elements <- unique(conceptual_matrix)
  
  init_matrix_list <- lapply(seq_len(n_inits), function(i, rep_f) {
    random_inits <- runif(length(unique_elements), 0, 10) %>% round(1)
    names_inits  <- unique_elements
    
    if(rep_f) {
      rep_f_init   <- runif(1) %>% round(2) # Reporting fraction init
      random_inits <- c(random_inits, rep_f_init)
      names_inits <- c(names_inits, "rho")
    }
    
    names(random_inits) <- names_inits
    random_inits
  }, rep_f = reporting_fraction)
  
  init_matrix_list
}

run_optim <- function(init_K, loglik) {
  
  optim_fit <- optim(par     = init_K, 
                     fn      = loglik, 
                     method  = "Nelder-Mead",
                     control = list(maxit = 3000))
}

optimise_SEIR2 <- function(inits, matrix_type, data_list, sim_specs, 
                           var_cache, pop_size, 
                           scenario = "perfect information") {
  
  if(scenario == "perfect information") {
    prefix <- "det"
  }
  
  if(scenario == "underreporting") {
    prefix <- "und"
  }
  
  filepath <- str_glue("./object_fits/optim/optim_{prefix}_{matrix_type}.rds")
  
  if(var_cache == FALSE) {
    
    output      <- optimise_SEIR(inits, matrix_type, data_list, sim_specs, 
                                 pop_size, scenario)
    saveRDS(output, filepath)
  } else {
    output    <- readRDS(filepath)
  }
  
  output
}

generate_inc_func <- function(mdl, sim_specs) {
  
  ds_consts      <- mdl$deSolve_components$consts %>% sort()
  ds_func        <- mdl$deSolve_components$func
  ds_stocks      <- mdl$deSolve_components$stocks
  constant_names <- names(ds_consts)
  
  # Create time vector
  start_time   <- sim_specs$start_time
  stop_time    <- sim_specs$stop_time
  step_size    <- sim_specs$step_size
  integ_method <- sim_specs$integ_method
  
  simtime <- seq(start_time, stop_time, by = step_size)
  
  function(pars) {
    
    for(i in seq_along(pars)) {
      par_name             <- names(pars[i])
      ds_consts[par_name]  <- pars[[i]]
    }
    
    o <- ode(y     = ds_stocks, times  = simtime, func  = ds_func,  
             parms = ds_consts, method = integ_method)
    
    o_df <- data.frame(o) %>% filter(time - trunc(time) == 0)
    
    incidence_df <- o_df %>% dplyr::select(time, starts_with("C")) %>% 
      pivot_longer(-time) %>% 
      group_by(name) %>% 
      mutate(y_hat = round(value - lag(value), 0)) %>% ungroup() %>% 
      filter(time != 0) %>% 
      mutate(variable = str_replace(name, "C", "y_hat")) %>% 
      dplyr::select(-name, -value) %>% 
      pivot_wider(names_from = variable, values_from = y_hat)
  }
}

get_mdl <- function(matrix_type, pop_size, stock_list) {
  
  model_file     <- paste0("./deterministic_models/4_cohorts_SEIR_matrix_",
                           matrix_type, ".stmx")
  
  const_list     <- list(latent_period = 1,
                         recovery_time = 2,
                         population    = pop_size)
  
  read_xmile(model_file, const_list = const_list,
             stock_list = stock_list)
}

best_fit <- function(optim_fit, get_incidence) {
  
  fit_df <- get_incidence(optim_fit$par) %>%
    pivot_longer(-time, names_to = "cohort", values_to = "y") %>% 
    mutate(index = as.numeric(str_replace(cohort, "y_hat", "")),
           source = "sim data") %>% 
    dplyr::select(-cohort)
  
  list(fit = optim_fit,
       fit_df = fit_df)
}
