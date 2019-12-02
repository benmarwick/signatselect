det_tridiag <-
function(m){
  n <- length(m$diag)
  det_tridiag_r <- function(d, i) {
    if (i == 0) {
      return(d)
    } else if (i == n) {
      return(det_tridiag_r(m$diag[i],i-1))
    } else {
      x <- m$diag[i] * d - m$above[i] * m$below[i]
      return(det_tridiag_r(x, i-1))
    }
  }
  return(det_tridiag_r(0,n))
}
