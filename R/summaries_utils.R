
#' Estimate R0 from K matrix
#'
#' @param K_hat A matrix
#' @param pop_df A data frame
#' @param tau_I A number that corresponds to the infectious period.
#'
#' @return A number.
R0_from_K <- function(K_hat, pop_df, tau_I = 2) {
  
  #Next Generation matrix
  NGM               <- K_hat * pop_df$proportion * tau_I
  eigensystem       <- eigen(NGM)
  max(Re(eigensystem$values))
}

calculate_time <- function(t_list) {
  t_obj <- t_list[[1]]
  (t_obj$toc - t_obj$tic) / 60
}