produce_synthetic_params <- function() {
  
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
  cm_object <- contact_matrix(polymod, age.limits = c(0, 5, 15, 45),
                              countries = "Finland",
                              survey.pop = "Finland",
                              symmetric = TRUE)
  
  # Social contact matrix corrected for reciprocity
  corrected_M_matrix           <- cm_object$matrix %>% t()
  dimnames(corrected_M_matrix) <- list(age_groups, age_groups)

  # Synthetic population
  age_groups <- c("00-04", "05-14", "15-44", "45+")
  agg_pop    <- cm_object$demography %>% 
    rename(group = age.group) %>% 
    mutate(group = age_groups) 
  
  # Normalised matrix
  symmetric_C_matrix <- corrected_M_matrix / agg_pop$proportion
  
  syn_pop_size         <- 1e4
  synthetic_population <- agg_pop %>% 
    mutate(syn_pop = round(syn_pop_size * proportion)) %>% 
    select(group, syn_pop)
  
  g_c_M_matrix <- draw_WAIFW(corrected_M_matrix, "Aggregated contacts")
  
  infectious_period <- 2 # days  
  infectivity       <- 0.1
  
  
  # Normalised age-specific effective contact rate
  necr_matrix <- symmetric_C_matrix * infectivity
  
  #=============================================================================
  # Theoretical R0
  #=============================================================================
  next_generation_matrix <- corrected_M_matrix * infectious_period * infectivity
  eigensystem            <- eigen(next_generation_matrix)
  theoretical_R0         <- max(Re(eigensystem$values))
  
  syn_WAIFW <- necr_matrix / syn_pop_size 
  
  g_syn_WAIFW <- draw_WAIFW(syn_WAIFW * 1e5, "")
  
  list(g_original_contacts = g_raw_M_matrix,
       g_c_M_matrix = g_c_M_matrix,
       synthetic_WAIFW = syn_WAIFW,
       g_syn_WAIFW = g_syn_WAIFW,
       syn_pop = synthetic_population,
       theoretical_R0 = theoretical_R0)
}


  
  