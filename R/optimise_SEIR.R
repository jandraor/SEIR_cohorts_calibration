
#' Optimise SEIR
#'
#' @param init_K List of init values
#' @param matrix_type A character. It is an upper-case letter indicating the 
#' matrix structure assumed for the simulation model.
#' @param data_list A list in which each element is a vector of incidence data
#' for a specific cohort.
optimise_SEIR <- function(init_K, matrix_type, data_list, sim_specs) {
  
  model_file     <- paste0("./deterministic_models/4_cohorts_SEIR_matrix_",
                           matrix_type, ".stmx")
  
  const_list     <- list(latent_period = 1,
                         recovery_time = 2)
  
  mdl            <- read_xmile(model_file, const_list = const_list,
                               stock_list = sim_specs$stock_list)
  
  get_incidence  <- generate_inc_func(mdl, sim_specs)
  
  poisson.loglik <- generate_ll(get_incidence)
  
  tic()
  
  optim_result <- mclapply(init_K, run_optim, mc.cores = 4,
                           loglik = poisson.loglik)
  
  toc(quiet = TRUE, log = TRUE)
  
  log.lst <- tic.log(format = FALSE)
  
  translation_df <- tibble(
    index_group = 1:4,
    var_incidence = paste0("incidence", 1:4),
    age_group = age_groups <- c("00-04", "05-14", "15-44", "45+"))
  
  result <- lapply(optim_result, best_fit, get_incidence = get_incidence)
  
  list(result    = result,
       time      = log.lst)
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

run_optim <- function(init_K, loglik) {
  
  optim_fit <- optim(par     = init_K, 
                     fn      = loglik, 
                     method  = "Nelder-Mead",
                     control = list(maxit = 3000))
}

# Generate log-likelihood function
generate_ll <- function(get_incidence) {
  
  function (pars) {
    
    if(any(sign(pars) == -1)) return(Inf)
    
    sim_incidences <- get_incidence(pars)
    
    n_cohorts      <- ncol(sim_incidences) - 1
    
    # Log-likelihoods
    ll <- vector(mode = "numeric", length = n_cohorts)
    
    for(i in seq_len(n_cohorts)) {
      ll[[i]] <- -1 * sum(dpois(lambda = sim_incidences[[i + 1]] + 0.0000001, 
                                x = data_list[[i]], 
                                log = TRUE))
    }
    
    sum(ll)
  }
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
