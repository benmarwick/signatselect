get_log_internal_integral <-
function(lgth, mu_b, cov_b_diag, tvec, s, alpha, x0){
  # Get the time vector minus the first time point
  tvec1 <- tvec[-1]
  
  # Calculate the logistic predicted frequencies
  muvec_g <- sapply(tvec1, function(x){get_g(x0,s,x)})
  
  # Get the predicted mean values for each time point
  mvec <- get_Mvec(lgth, x0, s, tvec)
  
  # Get the predicted standard deviations for each time point
  sigmavec <- get_sigmavec(lgth, x0, s, tvec)

  # Generate a covariance matrix from the population size, standard deviation and  predicted means
  covresults <- get_cov_mat_g(alpha, sigmavec, mvec)
  
  # Extract the components from the function above
  
  # Actual covariance matrix
  cov_mat_g <- covresults[[1]]
  
  # The log determinant of the covariance matrix
  log_det_g <- covresults[[2]]
  
  # The three vectors of the tridiagonal matrix underlying the covariance matrix
  cov_g_inv_diag <- covresults[[3]]
  cov_g_inv_above <- covresults[[4]]
  cov_g_inv_below <- covresults[[5]]
  
  if (alpha == Inf){
    # With infinite population size, just use the multivariate gaussian proabilities
    logz <- multivariate_gaussian_diag_logpdf(mu_b,cov_b_diag,muvec_g)
  } else {
    # Otherwise, use another function to predict the log-likelihood components
    logz <- get_logz(mu_b,muvec_g,cov_b_diag,cov_mat_g,log_det_g,cov_g_inv_diag,cov_g_inv_above,cov_g_inv_below)
  }
  return(logz)
}
