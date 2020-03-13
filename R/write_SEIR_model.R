write_SEIR_model <- function(matrix_type, filename, stocks, params_prior) {
  
  if(matrix_type == "sym") {
    model_file   <- "./deterministic_models/4_cohorts_SEIR_matrix_sym.stmx"
    ODE_fun_name <- "SEIR_matrix_sym"
  }
  
  if(matrix_type == "A") {
    model_file   <- "./deterministic_models/4_cohorts_SEIR_matrix_A.stmx"
    ODE_fun_name <- "SEIR_matrix_A"
  }
  
  if(matrix_type == "B") {
    model_file   <- "./deterministic_models/4_cohorts_SEIR_matrix_B.stmx"
    ODE_fun_name <- "SEIR_matrix_B"
  }
  
  if(matrix_type == "C") {
    model_file     <- "./deterministic_models/4_cohorts_SEIR_matrix_C.stmx"
    ODE_fun_name <- "SEIR_matrix_C"
  }
  
  if(matrix_type == "D") {
    model_file     <- "./deterministic_models/4_cohorts_SEIR_matrix_D.stmx"
    ODE_fun_name <- "SEIR_matrix_D"
  }
  
  mdl            <- read_xmile(model_file)
  constant_names <- map_chr(mdl$description$constants, "name")
  # serial interval parameters
  si_pars        <- c("recovery_time", "latent_period")
  params         <- constant_names[!constant_names%in% si_pars ] %>% sort()
  
  override_consts <- list(
    list(name = "recovery_time", value = 2),
    list(name = "latent_period", value = 1))
  
  stan_fun   <- create_stan_function(model_file, ODE_fun_name, 
                                     pars = params,
                                     override.consts = override_consts)
  
  stan_data  <- paste(
    "data {",
    "  int<lower = 1> n_obs; // Number of weeks sampled",
    "  int<lower = 1> n_params; // Number of model parameters",
    "  int<lower = 1> n_difeq; // Number of differential equations in the system",
    "  int y1[n_obs];",
    "  int y2[n_obs];",
    "  int y3[n_obs];",
    "  int y4[n_obs];",
    "  real t0; // Initial time point (zero)",
    "  real ts[n_obs]; // Time points that were sampled",
    "}", sep = "\n")
  
  stan_td     <- generate_transformed_data_block()
  stan_params <- generate_parameters_block("poisson")
  
  counter <- 0
  
  stan_init_stocks <- sapply(stocks, function(stock) {
    counter <<- counter + 1
    paste0("  y0[", counter, "] = ", stock, ";")
  }) %>% paste(collapse = "\n")
  
  fun_exe_line <- paste0("  y_hat = integrate_ode_rk45(",
                         ODE_fun_name,", y0, t0, ts, params, x_r, x_i);")
  
  stan_tp <- paste(
    "transformed parameters{",
    "  real y_hat[n_obs, n_difeq]; // Output from the ODE solver",
    "  real y0[n_difeq]; // Initial conditions",
    "  real incidence1[n_obs];",
    "  real incidence2[n_obs];",
    "  real incidence3[n_obs];",
    "  real incidence4[n_obs];",
    stan_init_stocks,
    fun_exe_line,
    "  incidence1[1] =  y_hat[1, 3] + y_hat[1, 4] - y0[3] - y0[4];",
    "  incidence2[1] =  y_hat[1, 7] + y_hat[1, 8] - y0[7] - y0[8];",
    "  incidence3[1] =  y_hat[1, 11] + y_hat[1, 12] - y0[11] - y0[12];",
    "  incidence4[1] =  y_hat[1, 15] + y_hat[1, 16] - y0[15] - y0[16];",
    "  for (i in 1:n_obs-1) {",
    "    incidence1[i + 1] = y_hat[i + 1, 3] + y_hat[i + 1, 4] - y_hat[i, 3] - y_hat[i, 4] + 0.00001;",
    "    incidence2[i + 1] = y_hat[i + 1, 7] + y_hat[i + 1, 8] - y_hat[i, 7] - y_hat[i, 8] + 0.00001;",
    "    incidence3[i + 1] = y_hat[i + 1, 11] + y_hat[i + 1, 12] - y_hat[i, 11] - y_hat[i, 12] + 0.00001;",
    "    incidence4[i + 1] = y_hat[i + 1, 15] + y_hat[i + 1, 16] - y_hat[i, 15] - y_hat[i, 16] + 0.00001;",
    "   }",
    "}", sep = "\n")
  
  stan_model_text <- c(
    params_prior,
    "  y1     ~ poisson(incidence1)",
    "  y2     ~ poisson(incidence2)",
    "  y3     ~ poisson(incidence3)",
    "  y4     ~ poisson(incidence4)"
  )
  stan_model  <- generate_model_text(stan_model_text)
  
  stan_text   <- paste(stan_fun, stan_data, stan_td, stan_params,
                       stan_tp, stan_model, sep = "\n")
  
  create_stan_file(stan_text, filename)
  
  list(params = params)
}