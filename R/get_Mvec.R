get_Mvec <-
function(lgth, x0, s, tvec){
  gvec <- sapply(tvec, function(x){get_g(x0,s,x)})
  f <- function(i){
    dt <- tvec[i+1] - tvec[i]
    return(get_M(gvec[i], s, dt))
  }
  return(sapply(1:(lgth),f))
}
