lambda_equation <- function(n_cohorts, susceptible_index, pad_width) {
  equation_list <- vector(mode = "character", length = n_cohorts)
  
  for(i in seq_len(n_cohorts)) {
    S <- str_pad(susceptible_index, width = pad_width, pad = "0")
    I <- str_pad(i, width = pad_width, pad = "0")
    equation_list[[i]] <- str_glue("k{S}{I}*I{i}")
  }
  
  paste(equation_list, collapse = " + ")
}
