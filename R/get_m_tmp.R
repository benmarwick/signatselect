get_m_tmp <-
function(cov_g_inv_diag, cov_g_inv_above, cov_g_inv_below, cov_b_diag) {
  n <- length(cov_b_diag)
  diag_fun <- function(i){1+cov_g_inv_diag[i] * cov_b_diag[i]}
  above_fun <- function(i){cov_g_inv_above[i] * cov_b_diag[i+1]}
  below_fun <- function(i){cov_g_inv_below[i] * cov_b_diag[i]}
  newdiag <- sapply(1:n,diag_fun)
  newabove <- sapply(1:(n-1),above_fun)
  newbelow <- sapply(1:(n-1),below_fun)
  return(list(diag=newdiag,above=newabove,below=newbelow))
}
