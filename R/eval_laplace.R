#' @importFrom stats dnorm optimize t.test
#' 
eval_laplace <-
function(f, a, b, prec, maxiter, verbose) {
  # Use the negative of the function for minimisation
  funcneg <- function(x){return(-f(x))}
  
  # Get a small search interval for optimisation
  minsearch <- find_min_interval(funcneg, a, b, 0.1, 0.001, 2, verbose)
  x1 <- minsearch[[1]]
  x2 <- minsearch[[2]]
  guess <- minsearch[[3]]
  if (verbose) {
    print(paste0('Maximum of f(x) is in the interval[',x1,',',x2,'], guess=',guess,'\n'))
  }
  
  # Use the R optimisation function to find best value within the interval
  optout <- optimize(funcneg, lower=x1, upper=x2, tol=prec)
  xmax <- optout$minimum
  fmax <- optout$objective
  if (verbose) {
    print(paste0("Maximum of the f(x) function found: x0 = ",xmax,"\n"))
  }
  
  # Calculate the second derivative at optimal value
  f2max <- get_second_deriv(f,xmax,prec,maxiter)
  
  # Evaluate the validity of the second derivative and return the Laplace approximation if the second derivative is valid.
  if (abs(f2max) < prec) {
    stop("Laplace method failed because the second derivative appears to be zero at the maximum")
  } else if (f2max > 0) {
    browser()
    stop("Laplace method failed because the second derivative appears to be positive at the maximum (non-sense!)")
  } else {
    return(-fmax + (log(2*pi)-log(-f2max))/2)
  }
}
