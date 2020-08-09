
#' Optimise SEIR
#'
#' @param init_K List of init values
#' @param matrix_type A character. It is an upper-case letter indicating the 
#' matrix structure assumed for the simulation model.
#' @param data_list A list in which each element is a vector of incidence data
#' for a specific cohort.
#' @param pop_size A number.
optimise_SEIR <- function(init_K, matrix_type, data_list, sim_specs, pop_size,
                          scenario) {
  
  stock_list    <- sim_specs$stock_list
  mdl           <- get_mdl(matrix_type, pop_size, stock_list)  
  get_incidence <- generate_inc_func(mdl, sim_specs)
  
  if(scenario == "perfect information") {
    poisson.loglik <- generate_ll(get_incidence, data_list)
  }
  
  if(scenario == "underreporting") {
    poisson.loglik <- generate_ll_und(get_incidence, data_list)
  }
  
  tic.clearlog()
  tic()
  
  optim_result <- mclapply(init_K, run_optim, mc.cores = 4,
                           loglik = poisson.loglik)
  
  toc(quiet = TRUE, log = TRUE)
  
  log.lst <- tic.log(format = FALSE)
  
  result <- lapply(optim_result, best_fit, get_incidence = get_incidence)
  
  list(result    = result,
       time      = log.lst)
}

# Generate log-likelihood function
generate_ll <- function(get_incidence, data_list) {
  
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

# Generate log-likelihood function
generate_ll_und <- function(get_incidence, data_list) {
  
  function (pars) {
    
    if(any(sign(pars) == -1)) return(Inf)
    
    p <- pars[length(pars)]
    
    if(p > 1) return(Inf)
    
    SEIR_params    <- pars[1:length(pars) - 1]
    sim_incidences <- get_incidence(SEIR_params)
    n_cohorts      <- ncol(sim_incidences) - 1
    
    # Log-likelihoods
    ll <- vector(mode = "numeric", length = n_cohorts)
    
    for(i in seq_len(n_cohorts)) {
      ll[[i]] <- -1 * sum(dpois(lambda = p * (sim_incidences[[i + 1]] + 0.0000001), 
                                x = data_list[[i]], 
                                log = TRUE))
    }
    
    sum(ll)
  }
}