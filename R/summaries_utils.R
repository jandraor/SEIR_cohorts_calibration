R0_from_K <- function(K_hat, pop_df) {
  tau_I <- 2 # Infectious period
  
  #Next Generation matrix
  NGM               <- K_hat * pop_df$proportion * tau_I
  eigensystem       <- eigen(NGM)
  max(Re(eigensystem$values))
}