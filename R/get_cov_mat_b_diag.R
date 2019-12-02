get_cov_mat_b_diag <-
function(bvec,nvec){
  len <- length(bvec)
  return(mapply(function(b,n){return((b+1)*(n-b+1)/(n+2)^2/(n+3))},bvec,nvec))
}
