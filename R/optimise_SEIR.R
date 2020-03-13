optimise_SEIR <- function(matrix_type, data_list, init.WAIFW = 20) {
  
  age_0_4_data     <- data_list$age_0_4_data
  age_5_14_data    <- data_list$age_5_14_data
  age_15_44_data   <- data_list$age_15_44_data
  age_45_over_data <- data_list$age_45_over_data
  
  model_file     <- paste0("./deterministic_models/4_cohorts_SEIR_matrix_",
                           matrix_type, ".stmx")
  mdl            <- read_xmile(model_file)
  ds_consts      <- mdl$deSolve_components$consts %>% sort()
  ds_func        <- mdl$deSolve_components$func
  ds_stocks      <- output_gsi$stocks
  constant_names <- names(ds_consts)
  si_pars        <- c("recovery_time", "latent_period")
  params         <- constant_names[!constant_names%in% si_pars ] %>% sort()
  
  # Create time vector
  simtime <- seq(0, 50, by = 1 / 128)
  
  generate_inc_func <- function() {
    
    function(pars) {
      ds_consts["latent_period"] <- 1
      ds_consts["recovery_time"] <- 2
      
      for(i in seq(length(params))) {
        param             <- params[[i]]
        ds_consts[param]  <- pars[[i]]
      }
      
      o <- ode(
        y     = ds_stocks, times = simtime,
        func  = ds_func,
        parms = ds_consts, method = "rk4") %>% data.frame() %>% 
        filter(time - trunc(time) == 0) %>% 
        mutate(eI1 = I1 + R1, 
               incidence1 = eI1 - lag(eI1, default = eI1[1]),
               incidence1 = round(incidence1),
               eI2 = I2 + R2, 
               incidence2 = eI2 - lag(eI2, default = eI2[1]),
               incidence2 = round(incidence2),
               eI3 = I3 + R3, 
               incidence3 = eI3 - lag(eI3, default = eI3[1]),
               incidence3 = round(incidence3),
               eI4 = I4 + R4, 
               incidence4 = eI4 - lag(eI4, default = eI4[1]),
               incidence4 = round(incidence4)) %>% 
        select(time, contains("incidence")) %>% 
        filter(time != 0)
    }
  }
  
  get_incidence <- generate_inc_func()
  
  poisson.loglik <-  function (params) {
    
    if(any(sign(params) == -1)) return(Inf)

    sim_incidences <- get_incidence(params)
    
    lik1 <- -1 * sum(dpois(lambda = sim_incidences$incidence1 + 0.0000001, x = age_0_4_data, 
                           log = TRUE))
    
    lik2 <- -1 * sum(dpois(lambda = sim_incidences$incidence2 + 0.0000001, x = age_5_14_data, 
                           log = TRUE))
    
    lik3 <- -1 * sum(dpois(lambda = sim_incidences$incidence3 + 0.0000001, x = age_15_44_data, 
                           log = TRUE))
    
    lik4 <- -1 * sum(dpois(lambda = sim_incidences$incidence4 + 0.0000001 , x = age_45_over_data, 
                           log = TRUE))
    
    sum(lik1, lik2, lik3, lik4)
  }
  
  optim_fit <- optim(par = rep(init.WAIFW, length(params)), 
                     fn = poisson.loglik, 
                     method = "Nelder-Mead")
  
  names(optim_fit$par) <- params
  
  translation_df <- tibble(
    index_group = 1:4,
    var_incidence = paste0("incidence", 1:4),
    age_group = age_groups <- c("00-04", "05-14", "15-44", "45+"))
  
  fit_df <- get_incidence(optim_fit$par) %>%
    pivot_longer(-time, names_to = "var_incidence", values_to = "incidence") %>% 
    mutate(source = "sim data") %>% left_join(translation_df) %>% 
    select(-var_incidence, index_group)
  
  list(fit = optim_fit,
       fit_df = fit_df)
}
  
    
