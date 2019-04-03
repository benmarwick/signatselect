#' Get the log-likelihood of data assuming sample frequencies are actual frequencies
#'
#' \code{get_llh_freq} computes the Gaussian approximation assuming that the actual frequencies are the sampled frequencies
#'
#' @param prior a function that returns the prior distribution (default: function(x){1})
#' @param lgth number of data points in the time series
#' @param prec minimum precision for optimisation
#' @param tvec time points for the time series (start with 0)
#' @param nuvec percent of successful samples per time point
#' @param s selection coefficient
#' @param alpha population size, bounded [0,Inf)
#' @param f0 initial frequency for the underlying logistic change
#' @return A log-likelihood estimate.
#' @export
get_llh_freq <-
function(prior, lgth, prec, tvec, nuvec, s, alpha, f0) {
  # Takes a prior distribution, the number of time points, a precision, the time vector, the frequency vector,
  # the selection coefficient, the population size and the initial frequency
  if (alpha == Inf) {
    # If alpha is infinite, use the logistic model exclusively
    
    # Extract the predicted frequencies based on the logistic with specified selection coefficient and initial probability
    gvec <- sapply(tvec, function(x){get_g(f0,s,x)})
    
    # Function for capturing the difference between the predicted frequencies and the observed frequencies 
    gfun <- function(g, nu, c) {
      (abs(g - nu) <= prec) && c
    }
    
    # Return different log-likelihoods depending on how similar the predicted and actual frequencies are
    if (Reduce2(gfun, gvec, nuvec, TRUE) && (prior(f0) > 0)) {
      return(Inf)
    } else {
      return(-Inf)
    }
  } else {
    # For when there is a finite population size
    
    # Generate a vector of functions that predicted values for each time point based on the previous time point
    muvec <- get_muvec(lgth, f0, s, tvec)
    
    # Generate a vector of predicted s.d. for each time point
    sigmavec <- get_sigmavec(lgth, f0, s, tvec)
    
    # Get the log-likelihood by looking up the probability of getting
    # the observed frequencies at each time point given the predicted
    # mean and s.d. and summing the log-probabilities (i.e., 
    # multiplying the probabilities)
    myfun_r <- function(acc, k) {
      if (k == 0) {
        acc + (log(prior(f0)))
      } else {
        x <- nuvec[k+1]
        mu <- muvec[[k]](nuvec[k])
        sigma <- sigmavec[k]/sqrt(alpha)
        g <- dnorm(x,mu,sigma,TRUE)
        return(myfun_r(acc+g,k-1))
      }
    }
    return(myfun_r(0,lgth-1))
  }
}
