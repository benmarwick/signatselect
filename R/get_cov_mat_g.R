get_cov_mat_g <-
function(alpha, sigmavec, mvec){
  len <- length(sigmavec)
  cov_mat_inv_diag <- rep.int(NaN,len)
  cov_mat_inv_above <- rep.int(NaN,len-1)
  cov_mat_inv_below <- rep.int(NaN,len-1)
  log_det_g <- -len * log(alpha)
  for (i in 1:(len-1)) {
    log_det_g <- log_det_g + 2 * log(sigmavec[i])
    cov_mat_inv_diag[i] <- alpha/sigmavec[i]^2 + alpha * mvec[i+1]^2/sigmavec[i+1]^2
    cov_mat_inv_above[i] <- -alpha*mvec[i+1]/sigmavec[i+1]^2
    cov_mat_inv_below[i] <- -alpha*mvec[i+1]/sigmavec[i+1]^2
  }
  log_det_g <- log_det_g + 2 * log(sigmavec[len])
  cov_mat_inv_diag[len]<-alpha/sigmavec[len]^2
  cov_mat <- invert_tridiag(cov_mat_inv_diag,cov_mat_inv_above,cov_mat_inv_below)
  return(list(cov_mat, log_det_g, cov_mat_inv_diag, cov_mat_inv_above, cov_mat_inv_below))
}
