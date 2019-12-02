get_mean_vec_b <-
function(bvec, nvec){
  f <- function(b, n){(b + 1)/(n + 2)}
  return(mapply(f, bvec, nvec))
}
