generate_synthetic_incidence <- function(output_sp, stop_time, 
                                         reporting_fraction, seed) {
  
  pop_df    <- output_sp$population
  n_cohorts <- nrow(pop_df)
  K_matrix  <- output_sp$K_matrix
  
  
  #-------------------------------------------------------------------------------
  # Values of constants
  #-------------------------------------------------------------------------------
  
  pad_width <- trunc(log(n_cohorts, base = 10)) + 1
  
  # number of separate elements
  n_sep <- (1 + n_cohorts) * (n_cohorts / 2)
  
  k_params <- vector(mode = "list", length = n_sep)
  k_names  <- vector(mode = "character", length = n_sep)  
  
  row_index <- 1
  cnt       <- 1 # counter
  
  for(i in 1:n_cohorts) {
    
    for(j in row_index:n_cohorts) {
      k_params[[cnt]] <- K_matrix[[i, j]]
      S               <- str_pad(i, width = pad_width, pad = "0")
      I               <- str_pad(j, width = pad_width, pad = "0")
      k_names[[cnt]]  <- str_glue("k{S}{I}")
      cnt             <- cnt + 1
    }
    
    row_index <- row_index + 1
  }
  
  names(k_params) <- k_names
  
  const_list      <- c(k_params, list(latent_period = 1,
                                      recovery_time = 2,
                                      population    = sum(pop_df$population)))
  
  
  #-----------------------------------------------------------------------------
  # Initial values of stocks
  #-----------------------------------------------------------------------------
  stock_names <- paste0(c("S", "E", "I", "R", "C"), rep(1:n_cohorts, each = 5))
  
  stock_matrix <- matrix(0, nrow = 1, ncol = length(stock_names))
  stock_names  -> colnames(stock_matrix)
  stock_df     <- as.data.frame(stock_matrix)
  
  infected_cohort <- 3
  infected_value  <- 10
  
  for(i in seq_len(n_cohorts)) {
    S_stock <- paste0("S", i)
    S_val   <- pop_df[i, "population"]
    
    if(i == infected_cohort) {
      I_stock <- paste0("I", i)
      C_stock <- paste0("C", i)
      stock_df[, I_stock] <- infected_value
      stock_df[, C_stock] <- infected_value
      S_val   <- S_val - infected_value
    }
    stock_df[, S_stock] <- S_val
  }
  
  stock_list <- as.list(stock_df)
  
  #-----------------------------------------------------------------------------
  model_file      <- str_glue("./deterministic_models/{n_cohorts}_cohorts_SEIR_matrix_sym.stmx")
  model_structure <- read_xmile(model_file, const_list = const_list,
                                stock_list = stock_list)
  

  # Create the start time, finish time, and time step
  START  <- 0
  FINISH <- stop_time
  STEP   <- 1 / 32
  
  # Create time vector
  simtime <- seq(START, FINISH, by = STEP )
  
  o <- data.frame(ode(y = model_structure$deSolve_components$stocks,
                       times = simtime,
                       func = model_structure$deSolve_components$func,
                       parms = model_structure$deSolve_components$consts, 
                       method = "rk4"))
  
  cumulative_df <- o %>% dplyr::select(time, starts_with("C")) %>% 
    filter(time - trunc(time) == 0)
  
  age_groups <- output_sp$population$group
  
  incidence_df <- cumulative_df %>% pivot_longer(-time) %>% 
    group_by(name) %>% 
    mutate(y = round(value - lag(value), 0)) %>% ungroup() 
  
  incidence_df <- incidence_df %>%
    mutate(index = as.numeric(str_replace(name, "C", "")),
           cohort = cut(index, 
                        breaks = length(age_groups), 
                        labels = age_groups)) %>% 
    dplyr::select(-name, -value) %>% filter(!is.na(y)) %>% 
    mutate(Scenario = "Perfect information")
  
  data_prf <- incidence_df %>% dplyr::select(time, y, cohort) %>% 
    pivot_wider(names_from = cohort, values_from = y) %>% 
    dplyr::select(-time) %>% as.list()
  
  # Underreported data
  set.seed(seed)
  data_und <- lapply(data_prf, function(incidence_data) {
    rbinom(length(incidence_data), incidence_data, 0.8)
  })
  
  data_list <- list(prf = data_prf,
                    und = data_und)
  
  und_incidence_df <- map2_df(data_und, as.factor(age_groups), make_df,
                              stop_time = stop_time)
  
  tidy_list <- list(prf = incidence_df,
                    und = und_incidence_df)
  
  g_incidences <- ggplot(incidence_df, aes(x = time, y = y)) +
    geom_col(fill = "#008080") +
    scale_y_continuous(labels = comma) +
    facet_wrap(~cohort) +
    labs(x = "Days", y = "Incidence [Cases / Day]") +
    theme_pubr()
  
  list(tidy_list     = tidy_list,
       data_list     = data_list,     
       cumulative_df = cumulative_df,
       g_incidences  = g_incidences,
       const_list    = const_list,
       stock_list    = stock_list)
}

make_df <- function(inc_data, ag, stop_time) {
  tibble(time = 1:stop_time, y = inc_data, cohort = as.character(ag),
         index = as.numeric(ag), Scenario = "Underreporting")
}

