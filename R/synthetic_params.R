produce_synthetic_params <- function(pop_size, age_limits) {
  
  age_groups <- c("00-04", "05-14", "15-44", "45+")
  
  cm_object <- contact_matrix(polymod, age.limits = seq(0, 70, 5),
                              countries = "Finland",
                              survey.pop = "Finland",
                              symmetric = TRUE)
    
  raw_M_matrix             <- cm_object$matrix %>% t()
 
  raw_age_cohorts <- c("00-04", "05-09", "10-14", "15-19", "20-24", "25-29",
                       "30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                       "60-64", "65-69", "70+") 
  
  dimnames(raw_M_matrix) <- list(raw_age_cohorts, raw_age_cohorts)
  
  g_raw_M_matrix <- draw_WAIFW(raw_M_matrix, "Finland's social contact matrix")

  
  # Contact matrix object
  cm_object <- contact_matrix(polymod, age.limits = age_limits,
                              countries = "Finland",
                              survey.pop = "Finland",
                              symmetric = TRUE)
  
  # Synthetic population
  agg_pop    <- cm_object$demography %>% 
    rename(group = age.group) %>% 
    mutate(group = format_age_group(as.character(group)),
           population = round(pop_size * proportion, 0)) 
  
  age_groups <- agg_pop$group
  
  # Social contact matrix corrected for reciprocity
  corrected_M_matrix           <- cm_object$matrix %>% t()
  dimnames(corrected_M_matrix) <- list(age_groups, age_groups)


  
  # Normalised age-specific contact rate matrix
  symmetric_C_matrix <- corrected_M_matrix / agg_pop$proportion
  
  g_c_M_matrix <- draw_WAIFW(corrected_M_matrix, "Aggregated contacts")
  
  infectious_period <- 2 # days  
  infectivity       <- 0.1
  
  
  # Normalised age-specific effective contact rate
  K_matrix <- symmetric_C_matrix * infectivity
  
  g_K_matrix <- draw_WAIFW(K_matrix, "K matrix", precision = 1)
  
  #=============================================================================
  # Theoretical R0
  #=============================================================================
  next_generation_matrix <- corrected_M_matrix * infectious_period * infectivity
  eigensystem            <- eigen(next_generation_matrix)
  theoretical_R0         <- max(Re(eigensystem$values))
  
  syn_WAIFW <- K_matrix / pop_size
  
  g_syn_WAIFW <- draw_WAIFW(syn_WAIFW, "WAIFW matrix", precision = 6)
  
  list(g_original_contacts = g_raw_M_matrix,
       g_c_M_matrix        = g_c_M_matrix,
       synthetic_WAIFW     = syn_WAIFW,
       g_syn_WAIFW         = g_syn_WAIFW,
       population          = agg_pop,
       theoretical_R0      = theoretical_R0,
       K_matrix            = K_matrix,
       g_K_matrix          = g_K_matrix)
}

format_age_group <- function(ag_vector) {
  pattern <- regex("\\[(\\d+),(\\d+)\\)")
  
  new_ag_vector <- vector(length = length(ag_vector), mode = "character") 
  
  for(i in seq_along(ag_vector)) {
    current_ag <- ag_vector[[i]]
    
    if(str_detect(current_ag, pattern)) {
      output_sm   <- str_match(current_ag, pattern)
      lower_bound <- output_sm[[2]] %>% str_pad(width = 2, pad = "0")
      upper_bound <- (as.numeric(output_sm[[3]]) - 1) %>% 
        str_pad(width = 2, pad = "0")
      current_ag  <- paste(lower_bound, upper_bound, sep = "-")
    }
    
    new_ag_vector[[i]] <- current_ag
  }
  new_ag_vector
}

  
  