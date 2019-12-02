get_second_deriv <-
function(myf,x,prec,maxiter){
  dx <- 2 * prec
  fx <- myf(x)
  fxplus <- myf(x+dx)
  fxminus <- myf(x-dx)
  return((fxminus-2 * fx + fxplus) / (dx^2))
}
