#' @importFrom limSolve Solve.tridiag
invert_tridiag <-
function(diag,abovediag,belowdiag){
  n <- length(diag)
  invm <- matrix(NaN, n, n)
  for (i in 1:n){
    barr <- rep(0, n)
    barr[i]<-1
    invm[i,]<-Solve.tridiag(belowdiag,diag,abovediag,barr)
  }
  return(invm)
}
