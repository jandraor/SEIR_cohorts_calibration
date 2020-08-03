get_cm_sym <- function(size) {
  
  cm <- vector(mode = "character", length = size ** 2)
  k  <- 1
  
  for(i in seq_len(size)) {
    
    for(j in seq_len(size)) {
      
      k_value <- str_glue("k{i}{j}")
      
      if(i > j) k_value <- str_glue("k{j}{i}")
      cm[[k]] <- k_value
      k       <- k + 1
    }
  }
  cm
}