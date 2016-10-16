get_sigmavec <-
function(lgth, x0, s, tvec) {
  gvec <- sapply(tvec, function(x){get_g(x0,s,x)})
  sigmafun <- function(i) {
    dt <- tvec[i+1] - tvec[i]
    return(sqrt(get_var(gvec[i], s, dt)))
  }
  return(unlist(sapply(1:(lgth-1), sigmafun)))
}
