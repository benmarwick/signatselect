#' Get the log-likelihood of data while incorporating sampling error
#'
#' \code{get_llh} uses Laplace approximation to compute the log-likelihood of the gaussian approximation assuming the data is drawn from a series of binomial choices.
#'
#' @param prior a function that returns the prior distribution (default: function(x){1})
#' @param lgth number of data points in the time series
#' @param tvec time points for the time series (start with 0)
#' @param nvec total number of samples for each time point
#' @param bvec number of successful samples at each time point
#' @param prec minimum precision for optimisation
#' @param maxiter maximum number of iterations for optimisation
#' @param verbose print intermediate output (TRUE or FALSE)
#' @param s selection coefficient
#' @param alpha population size, bounded [0,Inf)
#' @param f0 initial frequency for the underlying logistic change
#' @return A log-likelihood estimate.
#' @export
get_llh <-
function(prior, lgth, tvec, nvec, bvec, prec, maxiter, verbose, s, alpha, f0) {
  # Takes a prior distribution, number of time points, time vector, success vectore, number of tokens vectore,
  # precision, maximum number of iterations, whether or not to print intermediate outpus,
  # selection coefficient, population size and initial frequency (current initial frequency does not do anything)
  
  # Get the mean predicted values from the gaussian approximation
  muvec_b <- get_mean_vec_b(bvec, nvec)
  
  # Get the diagnal of the covariance matrix
  cov_mat_b_diag = get_cov_mat_b_diag(bvec,nvec)
  
  # Extract the initial predicted mean and s.d.
  mu_b_0 <- muvec_b[1]
  sigma_b_0 <- cov_mat_b_diag[1]
  
  # Generate vectors with the remaining means and diagnal of the covariance matrix
  mu_b <- muvec_b[-c(1)]
  cov_b_diag <- cov_mat_b_diag[-c(1)]
  
  # Get a scaling constant based on the number of tokens
  logsc <- get_log_scale_constant(nvec[-c(1)])
  
  # Generate functions for the Laplace approximation of the integral
  func0 <- function(x0){return((dnorm(x0,mu_b_0,sigma_b_0,log=T)+log(prior(x0))))}
  func1 <- function(x0){return(get_log_internal_integral(lgth, mu_b,cov_b_diag, tvec, s, alpha, x0))}
  func <- function(x0){
    return(func0(x0)+func1(x0))}
  
  # Calculate the Laplace approximation of the integral
  logint <- eval_laplace(func, prec, 1-prec, prec, maxiter, verbose)
  
  # Return the log-likelihood (the Laplace approximation plus the scaling constant)
  return(logint + logsc)
}
