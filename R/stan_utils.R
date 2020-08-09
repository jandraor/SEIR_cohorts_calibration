create_stan_file <- function(stan_text, filename) {
  file_connection <- file(filename)
  writeLines(stan_text, file_connection)
  close(file_connection)
}

# n_constants is the number of unknown parameters (except init values) in the
# ODE Model

fit_stan_model <- function(n_constants, incidence_data,
                           warm_up, iterations, seed = NULL,
                           adapt_delta = 0.8,
                           chains = 3,
                           stan_file = "SIR.stan",
                           n_difeq = 3) {
  
  stan_d = list(n_obs    = length(incidence_data),
                n_params = n_constants, 
                n_difeq  = n_difeq, # number of differential equations
                y        = incidence_data,
                t0       = 0,
                ts       = 1:length(incidence_data))
  
  # Test / debug the model:
  test <- stan(stan_file,
               data = stan_d,
               chains = 1, 
               iter = 10,
               verbose = FALSE,
               refresh = 0)
  
  if(is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1)
  }
  
  # Fit and sample from the posterior
  stan_fit <- stan(fit    = test,
                   data   = stan_d,
                   chains = chains,
                   warmup = warm_up,
                   iter   = iterations,
                   cores  = 3,
                   seed    = seed,
                   control = list(adapt_delta = adapt_delta))
}

generate_stan_model <- function(prior_susceptible = "", 
                                             likelihood = "normal",
                                             n_params) {
  
  if(n_params == 1 & likelihood == "normal") {
    model <- paste("model {", 
                   "  params[1] ~ uniform(0, 500);",
                   "  sigma ~ cauchy( 0 , 2 );",
                   "  y ~ normal(incidence, sigma); ",
                   "}", sep = "\n")
  }
  
  if(n_params == 1 & likelihood == "poisson") {
    model <- paste("model {", 
                   "  params[1] ~ uniform(0, 500);",
                   "  y ~ poisson(incidence); ",
                   "}", sep = "\n")
  }
  
  if(n_params == 2 & likelihood == "poisson") {
    model <- paste("model {", 
                   "  params[1] ~ uniform(0, 500);",
                   "  params[2] ~ uniform(0, 20);",
                   "  y ~ poisson(incidence); ",
                   "}", sep = "\n")
  }
  
  if(n_params == 3) {
    if(likelihood == "normal") {
      
      model <- paste("
model {
  params[1] ~ uniform(0, 500);
  params[2] ~ uniform(0, 20);",
                     prior_susceptible,
                     "sigma ~ cauchy( 0 , 2 );
  y ~ normal(incidence, sigma); 
}
", sep = "\n  ")
      
    }
  
    if(likelihood == "poisson") {
    
    model <- paste("
model {
  params[1] ~ uniform(0, 500);
  params[2] ~ uniform(0, 20);",
                   prior_susceptible,
                   "y ~ poisson(incidence);
}
  ", sep = "\n  ")
    }
  }
  
  model
}

generate_transformed_parameters <- function(model_type, initValues) {
  
  if(model_type == "SIR") {
    
    S_text <- paste0("  y0[1] = ", initValues$Susceptible, ";")
    I_text <- paste0("  y0[2] = ", initValues$Infected, ";")
    R_text <- paste0("  y0[3] = ", initValues$Recovered, ";")
    
    SIR_text <- paste(S_text, I_text, R_text, sep = "\n")
    
    stan_transformed_parameters <- paste(
      "transformed parameters{
         real y_hat[n_obs, n_difeq]; // Output from the ODE solver
         real y0[n_difeq]; // Initial conditions
         real incidence[n_obs];
      ", SIR_text, 
      "  y_hat = integrate_ode_rk45(SIR, y0, t0, ts, params, x_r, x_i);
      
         incidence[1] =  y_hat[1, 2] + y_hat[1, 3] - y0[2] - y0[3];
         
         for (i in 1:n_obs-1) {
           incidence[i + 1] = y_hat[i + 1, 2] + y_hat[i + 1, 3] - y_hat[i, 2] - y_hat[i, 3] + 0.000001; 
         }
      }", sep = "\n")
  }
  
  if(model_type == "SEIR") {
    
    S_text <- paste0("  y0[1] = ", initValues$Susceptible, ";")
    E_text <- paste0("  y0[2] = ", initValues$Exposed, ";")
    I_text <- paste0("  y0[3] = ", initValues$Infected, ";")
    R_text <- paste0("  y0[4] = ", initValues$Recovered, ";") 
    
    SEIR_text <- paste(S_text, E_text, I_text, R_text, sep = "\n")
    
    stan_transformed_parameters <- paste(
      "transformed parameters{
         real y_hat[n_obs, n_difeq]; // Output from the ODE solver
         real y0[n_difeq]; // Initial conditions
         real incidence[n_obs];
      ", SEIR_text, 
      "  y_hat = integrate_ode_rk45(SEIR, y0, t0, ts, params, x_r, x_i);
      
         incidence[1] =  y_hat[1, 3] + y_hat[1, 4] - y0[3] - y0[4];
         
         for (i in 1:n_obs-1) {
           incidence[i + 1] = y_hat[i + 1, 3] + y_hat[i + 1, 4] - y_hat[i, 3] - y_hat[i, 4] + 0.000001; 
         }
      }", sep = "\n")
  }
  
  stan_transformed_parameters
}

generate_model_text <- function(model_lines) {
  text_body <- sapply(model_lines, function(model_line) {
    paste0("  ", model_line, ";")
  }) %>% paste(collapse = "\n")
  model_text <- paste(
    "model {",
    text_body,
    "}", sep = "\n")
}

generate_data_block <- function(distribution) {
  stan_data <- NULL
  
  if(distribution == "poisson") {
    
    stan_data <- paste(
      "data {",
      "  int<lower = 1> n_obs; // Number of weeks sampled",
      "  int<lower = 1> n_params; // Number of model parameters",
      "  int<lower = 1> n_difeq; // Number of differential equations in the system",
      "  int y[n_obs];",
      "  real t0; // Initial time point (zero)",
      "  real ts[n_obs]; // Time points that were sampled",
      "}", sep = "\n")
  }
  
  if(distribution == "normal") {
    stan_data <- paste(
      "data {",
      "  int<lower = 1> n_obs; // Number of weeks sampled",
      "  int<lower = 1> n_params; // Number of model parameters",
      "  int<lower = 1> n_difeq; // Number of differential equations in the system",
      "  real y[n_obs]; ",
      "  real t0; // Initial time point (zero)",
      "  real ts[n_obs]; // Time points that were sampled",
      "}", sep = "\n") 
  }
  
  stan_data
}

get_stock_inits <- function(mdl, unknown = NA) {
  stocks           <- mdl$deSolve_components$stocks
  stan_init_stocks <- vector(mode = "character", length = length(stocks))
  
  for(i in seq_along(stocks)) {
    stock      <- stocks[i]
    stock_name <- names(stock)
    
    stock_val  <- ifelse(stock_name %in% unknown, 
                         paste0(stock_name, "0"), stocks[i])
    
    stan_init_stocks[[i]] <- str_glue("  y0[{i}] = {stock_val};")
  }
  
  stock_init_text <- paste(stan_init_stocks, collapse = "\n")
}


run_stan <- function(filename, var_cache, cache_file, arg_list) {
  
  
  if(var_cache == FALSE) {
    # # Test / debug the model:
    # test <- stan(filename, data = arg_list$data, 
    #              chains = 1, iter = 10,
    #              verbose = TRUE, refresh = 0)
    
    set_cmdstan_path(file.path(Sys.getenv("HOME"), ".cmdstanr", 
                               "cmdstan-2.24.0-rc1"))
    mod <- cmdstan_model(filename)
    #---------------------------------------------------------------------------
    tic.clearlog()
    tic()
    # Fit and sample from the posterior
    # arg_list$fit <- test
    # stan_fit     <- do.call("stan", arg_list)
    
    fit <- mod$sample(data            = arg_list$data,
                      seed            = arg_list$seed,
                      chains          = arg_list$chains,
                      parallel_chains = arg_list$cores,
                      iter_warmup     = arg_list$warmup,
                      iter_sampling   = arg_list$iter,
                      refresh         = 5,
                      save_warmup     = TRUE)
  
    
    
    toc(log = TRUE, quiet = FALSE)
    log.lst <- tic.log(format = FALSE)
    #---------------------------------------------------------------------------
    stan_fit <- rstan::read_stan_csv(fit$output_files())
    output <- list(stan_fit = stan_fit, time = log.lst)
    saveRDS(output, cache_file)
  } else {
    output    <- readRDS(cache_file)
  }
  
  output
}

