generate_synthetic_incidence <- function(output_sp) {
  
  model_file <- "./deterministic_models/4_cohorts_SEIR_matrix_full.stmx"
  model_structure <- read_xmile(model_file)
  syn_WAIFW <- output_sp$synthetic_WAIFW * 1e5
  constants <- names(model_structure$deSolve_components$consts)
  
  for(i in 1:nrow(syn_WAIFW)) {
    
    for(j in 1:ncol(syn_WAIFW)) {
      param <- paste0("B", i, j)
      pos   <- which(param == constants)
      model_structure$deSolve_components$consts[[pos]] <- syn_WAIFW[i, j]
    }
  }
  
  pos_lp <- which("latent_period" == constants)
  model_structure$deSolve_components$consts[[pos_lp]] <- 1
  
  pos_rt <- which("recovery_time" == constants)
  model_structure$deSolve_components$consts[[pos_rt]] <- 2
  
  #-------------------------------------------------------------------------------
  # Initial values of stocks
  #-------------------------------------------------------------------------------
  stock_names <- names(model_structure$deSolve_components$stocks)
  
  syn_pop <- output_sp$syn_pop %>% mutate(index = row_number(),
                                          init_I = c(0, 0, 10, 0),
                                          init_S = syn_pop - init_I)
  for(i in 1:nrow(syn_pop)) {
    row      <- syn_pop[i, ]
    S_cohort <- paste0("S", row$index)
    S_value  <- row$init_S
    pos_S    <- which(S_cohort == stock_names)
    model_structure$deSolve_components$stocks[pos_S] <- S_value
    I_cohort <- paste0("I", row$index)
    I_value  <- row$init_I
    pos_I    <- which(I_cohort == stock_names)
    model_structure$deSolve_components$stocks[pos_I] <- I_value
  }
  
  
  # Create the start time, finish time, and time step
  START  <- 0
  FINISH <- 50
  STEP   <- 1 / 128
  
  # Create time vector
  simtime <- seq(START, FINISH, by = STEP )
  
  o <- data.frame(ode(y = model_structure$deSolve_components$stocks,
                       times = simtime,
                       func = model_structure$deSolve_components$func,
                       parms = model_structure$deSolve_components$consts, 
                       method = "rk4"))
  
  generate_incidence_df <- function(index, df) {
    Infected  <- paste0("I", index)
    Recovered <- paste0("R", index)
    
    incidence_data <- df %>% select(time, c(Infected, Recovered)) %>% 
      rename(Infected := !!Infected, Recovered := !!Recovered) %>% 
      filter(time - trunc(time) == 0) %>% 
      mutate(everInfected = Infected + Recovered,
             incidence = everInfected - lag(everInfected, 
                                            default = everInfected[1]),
             incidence = round(incidence)) %>% 
      select(time, incidence) %>% 
      filter(time != 0) %>% 
      mutate(index_group = index)
  }
  
  generate_cml_incidence_df <- function(index, df) {
    Infected  <- paste0("I", index)
    Recovered <- paste0("R", index)
    
    incidence_data <- df %>% select(time, c(Infected, Recovered)) %>% 
      rename(Infected := !!Infected, Recovered := !!Recovered) %>% 
      filter(time - trunc(time) == 0) %>% 
      mutate(cml_incidence = round(Infected + Recovered, 0)) %>% 
      select(time, cml_incidence) %>% 
      filter(time != 0) %>% 
      mutate(index_group = index)
  }
  
  incidence_df  <- map_df(1:4, generate_incidence_df, df = o)
  
  cumulative_df <- map_df(1:4, generate_cml_incidence_df , df = o)
  
  age_groups <- c("00-04", "05-14", "15-44", "45+")
  names(age_groups) <- as.character(1:4)
  
  g_incidences <- ggplot(incidence_df, aes(x = time, y = incidence)) +
    geom_point(size = 0.5, shape = 1, colour = "#008080", fill = "#008080") +
    geom_line(colour ="#008080", alpha = 0.5) +
    facet_wrap(~index_group, ncol = 2, 
               labeller = labeller(index_group = age_groups)) +
    geom_vline(xintercept = 25, linetype = "dashed", alpha = 0.3) +
    theme_test()
  
  list(incidence_df = incidence_df,
       cumulative_df = cumulative_df,
       g_incidences = g_incidences,
       constants    = constants,
       stocks       = model_structure$deSolve_components$stocks)
}

