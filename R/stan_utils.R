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

generate_stan_function <- function(model_type, params = NULL) {
  stan_function <- NULL
  
  
  if(model_type == "SIR") {
  
  SIR_fun <- "
functions {
  real[] SIR(real t,
            real[] y, // stocks
            real[] params,
            real[] x_r,
            int[] x_i) {
              
      real dydt[3];
      real population;
      real probability;
      real contactsPerInfected;
      real IR; // Infection rate
      real RR; // Recovery rate
      
      population          = y[1] + y[2] + y[3];
      probability         = y[1] / population;
      contactsPerInfected = y[2] * params[1]; // effective contacts
      
      IR                  = contactsPerInfected * probability;
      RR                  = y[2] * params[2]; // recovery proportion
      
      dydt[1] = -IR;
      dydt[2] = IR - RR;
      dydt[3] = RR;
      
      return dydt;
    }
}
"
    if(!is.null(params$recovery_delay)) {
      SIR_fun <- SIR_fun %>% 
        stringr::str_replace("params\\[2\\]", 
                             as.character(params$recovery_delay)) 
    }
  
   
    
  stan_function <- SIR_fun
  }

  if(model_type == "SEIR") {
    SEIR_fun <- paste(
      "functions {",
      "  real[] SEIR(real t,",
      "              real[] y, // stocks",
      "              real[] params,",
      "              real[] x_r,",
      "              int[] x_i) {",
      "  real dydt[4];",
      "  real population;",
      "  real probability;",
      "  real contactsPerInfected;",
      "  real IR; // Infection rate",
      "  real InR; // Incidence rate",
      "  real RR; // Recovery rate",
      "  population = y[1] + y[2] + y[3] + y[4];",
      "  probability = y[1] / population;",
      "  contactsPerInfected = y[3] * params[1]; // effective contacts",
      "  IR                  = contactsPerInfected * probability;",
      "  InR                 = y[2] * params[2]; // latent period",
      "  RR                  = y[3] * params[3]; // recovery delay",
      "  dydt[1] = -IR;",
      "  dydt[2] = IR - InR;",
      "  dydt[3] = InR - RR;",
      "  dydt[4] = RR;",
      "  return dydt;",
      "  }",
      "}",
      sep = "\n")
    
    params_modified <- 0;
    
    if(!is.null(params$effective_contacts)){
      # Do something
      params_modified <- params_modified + 1
    }
    
    if(!is.null(params$latent_period)){
      SEIR_fun <- SEIR_fun %>% 
        stringr::str_replace("params\\[2\\]", 
                             as.character(params$latent_period)) 
    }
    
    if(is.null(params$latent_period) & params_modified > 0){
      SEIR_fun <- SEIR_fun %>% 
        stringr::str_replace("params\\[2\\]", "params\\[", 
                             params_modified, "\\]") 
    }
    
    if(!is.null(params$recovery_delay)){
      SEIR_fun <- SEIR_fun %>% 
        stringr::str_replace("params\\[3\\]", 
                             as.character(params$recovery_delay)) 
    }
    
    if(is.null(params$recovery_delay) & params_modified > 0){
      SEIR_fun <- SEIR_fun %>% 
        stringr::str_replace("params\\[2\\]", "params\\[", 
                             params_modified, "\\]") 
    }
    
    
    stan_function <- SEIR_fun
  }

  stan_function
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

generate_transformed_data_block <- function() {
  stan_transformed_data <- paste(
    "transformed data {",
    "  real x_r[0];",
    "  int x_i[0];",
    "}", sep = "\n")
}

generate_parameters_block <- function(distribution, stocks = NULL) {
  stan_parameters <- NULL
  stock_params    <- NULL
  
  if(!is.null(stocks)) {
    stock_params <- purrr::map_chr(stocks, function(stock) {
      paste0("  real<lower = 0> ", stock,";")
    }) %>% paste(collapse = "\n")
  }
  
  if(distribution == "poisson") {
    stan_parameters <- paste(
      "parameters {",
      "  real<lower = 0> params[n_params]; // Model parameters",
      stock_params,
      "}", sep = "\n")
  }
  
  if(distribution == "normal") {
    stan_parameters <- paste(
      "parameters {",
      "  real<lower = 0> params[n_params]; // Model parameters",
      "  real<lower = 0> sigma;",
      stock_params,
      "}", sep = "\n")
    
  }
  
  stan_parameters
}


