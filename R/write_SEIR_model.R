#' @param filename A string that indicates the file path where the Stan file 
#' will be saved.
write_SEIR_model <- function(matrix_type, filename, stock_list, params_prior,
                             pop_size, scenario = "perfect information",
                             nc = 4, recovery_time = 2,
                             latent_period = 1) {
  
  model_file   <- str_glue("./deterministic_models/{nc}_cohorts_SEIR_matrix_{matrix_type}.stmx")
  ODE_fun_name <- str_glue("SEIR_matrix_{matrix_type}")
  
  mdl            <- read_xmile(model_file, stock_list = stock_list)
  constants      <- sd_constants(mdl)
  constant_names <- constants$name
  # serial interval parameters & population
  fixed_pars     <- c("recovery_time", "latent_period", "population")
  params         <- constant_names[!constant_names %in% fixed_pars] %>% sort()
  override_consts <- list(recovery_time = recovery_time,
                          latent_period = latent_period,
                          population    = pop_size)
  
  stan_fun <- stan_ode_function(model_file, ODE_fun_name,
                                pars = params, const_list = override_consts)
  
  stock_inits <- sd_stocks(mdl)$init_value
  stock_names <- sd_stocks(mdl)$name
  C_indexes   <- which(str_detect(stock_names, "C"))
  
  stan_data  <- stan_data_block(nc)
  
  stan_params <- get_stan_params(scenario)
  
  stan_init_stocks <- get_stock_inits(mdl)
  
  fun_exe_line <- str_glue("  y_hat = ode_rk45({ODE_fun_name}, y0, t0, ts, params);") 
  
  stan_tp <- get_tp(stan_init_stocks, fun_exe_line, scenario, nc, C_indexes)
  
  likelihood <- str_glue("  y{1:nc}     ~ poisson(incidence{1:nc})") %>% 
    as.vector()
    
  
  stan_model_text <- c(params_prior, likelihood)
  stan_model  <- generate_model_text(stan_model_text)
  
  stan_text   <- paste(stan_fun, stan_data, stan_params,
                       stan_tp, stan_model, sep = "\n")
  
  create_stan_file(stan_text, filename)
  
  list(params      = params,
       stock_inits = stock_inits)
}

get_stan_params <- function(scenario) {
  
  if(scenario == "perfect information") {
    stan_params <- paste(
      "parameters {",
      "  real<lower = 0> params[n_params]; // Model parameters",
      "}", sep = "\n")
  }
  
  if(scenario == "underreporting") {
    stan_params <- paste(
      "parameters {",
      "  real<lower = 0> params[n_params]; // Model parameters",
      "  real<lower = 0, upper = 1> rho;",
      "}", sep = "\n")
  }
  
 stan_params
}

get_tp <- function(stan_init_stocks, fun_exe_line, scenario, nc, C_indexes) {
  
  inc_declaration <- str_glue("  real incidence{1:nc}[n_obs];") %>% 
    paste(collapse = "\n")
  
  if(scenario == "perfect information") {
    
    inc_t1 <- str_glue("  incidence{1:nc}[1] =  y_hat[1, {C_indexes}]  - y0[{C_indexes}];") %>% 
      paste(collapse = "\n")
    
    inc_tn <- str_glue("    incidence{1:nc}[i + 1] = y_hat[i + 1, {C_indexes}] - y_hat[i, {C_indexes}] + 0.00001;") %>% 
      paste(collapse = "\n")
    
    stan_tp <- paste(
      "transformed parameters{",
      "  vector[n_difeq] y_hat[n_obs]; // Output from the ODE solver",
      inc_declaration,
      fun_exe_line,
      inc_t1,
      "  for (i in 1:n_obs-1) {",
      inc_tn,
      "   }",
      "}", sep = "\n")
  }
  
  if(scenario == "underreporting") {
    
    inc_t1 <- str_glue("  incidence{1:nc}[1] =  rho * (y_hat[1, {C_indexes}]  - y0[{C_indexes}]);") %>% 
      paste(collapse = "\n")
    
    inc_tn <- str_glue("    incidence{1:nc}[i + 1] = rho * (y_hat[i + 1, {C_indexes}] - y_hat[i, {C_indexes}] + 1e-5);") %>% 
      paste(collapse = "\n")
    
    stan_tp <- paste(
      "transformed parameters{",
      "  vector[n_difeq] y_hat[n_obs]; // Output from the ODE solver",
      inc_declaration,
      fun_exe_line,
      inc_t1,
      "  for (i in 1:n_obs-1) {",
      inc_tn,
      "   }",
      "}", sep = "\n")
  }
  stan_tp
}

stan_data_block <- function(n_datasets) {
  data_lines <- str_glue("  int y{1:n_datasets}[n_obs];") %>% 
    paste(collapse = "\n")
  
  stan_data  <- paste(
    "data {",
    "  int<lower = 1> n_obs; // Number of days sampled",
    "  int<lower = 1> n_params; // Number of model parameters",
    "  int<lower = 1> n_difeq; // Number of differential equations in the system",
    data_lines,
    "  real t0; // Initial time point (zero)",
    "  real ts[n_obs]; // Time points that were sampled",
    "  vector[n_difeq] y0;",
    "}", sep = "\n")
}
