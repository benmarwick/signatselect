get_logz <-
function(mu_b, mu_g, cov_b_diag, cov_g, log_det_g, cov_g_inv_diag, cov_g_inv_above, cov_g_inv_below){
  # Get the number of remaining time points
  l <- length(mu_b)
  
  # Generate a temporary covariance matrix from the tridiagonal matrix
  m_tmp <- get_m_tmp(cov_g_inv_diag, cov_g_inv_above, cov_g_inv_below, cov_b_diag)
  
  # Get the determinant of the temporary matrix
  det_m_tmp <- det_tridiag(m_tmp)
  
  # Combine the determinant of the temporary matrix with that of the actual covariance matrix
  log_det <- log(det_m_tmp) + log_det_g
  
  # Invert the temporary matrix
  m_tmp_inv <- invert_tridiag(m_tmp$diag,m_tmp$above,m_tmp$below)
  
  # Get a vector of the difference between the predicted and empirical frequencies
  mu12<- mu_g - mu_b
  
  # Multiply the inverted temporary matrix with the inverted tridiagonal covariance matrix
  m_tmp2_inv <- m_tmp_inv %*% vec2mat_tridiag(cov_g_inv_diag,cov_g_inv_above,cov_g_inv_below)
  
  # Multiply the output of the previous row by the difference between the predicted and empirical frequencies
  y <- m_tmp2_inv %*% mu12
  
  # Multiply the output of the previous row by the difference between the predicted and empirical frequencies
  exparg <- mu12 %*% y
  
  # Calculate 2 * the output
  numerator <- -(c(exparg) + l * log(2 * pi) + c(log_det))
  return(numerator/2)
}
