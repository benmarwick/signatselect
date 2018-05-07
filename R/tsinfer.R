#' Infer selection/drift parameters for time-series allele-frequency data
#'
#' \code{tsinfer} computes the selection coefficient, population size and initial frequency for the logistic using a Gaussian approximation of the Kimura expression for time-series allele-frequency data.
#'
#' @param tvec Time coordinates for the time-series (start with 0)
#' @param bvec Number of new form per time point
#' @param nvec Total number of samples per time point
#' @param maxiter Maximum number of iterations
#' @param prec Precision for optimisation
#' @param iffreq whether to assume that the sample frequencies are the population frequencies (default=FALSE)
#' @param ifneut whether to compute only the neutral model (default=FALSE)
#' @param iffixedf0 whether to use the initial sample frequency as the initial frequency of the logistic component of the model (default=FALSE)
#' @param verbose whether to print intermediate output (detault=FALSE)
#' @param mins minimum s value to consider (default=-2)
#' @param mins maximum s value to consider (default=2)
#' @param minalpha minimum alpha value to consider (default=10)
#' @param minalpha maximum alpha value to consider (default=1e8)
#' @param minf0 minimum f0 value to consider in log-odds (default=-10)
#' @param maxf0 maximum f0 value to consider in log-odds (default=10)
#' @return A list with the neutral and non-neutral parameter values and associated log-likelihoods
#' @export
tsinfer <-
function(tvec, bvec, nvec,
                    maxiter=200,
                    prec=1.0e-06,
                    iffreq=FALSE,
                    ifneut=FALSE,
		    iffixedf0=FALSE,
                    verbose=FALSE,
		    mins=-2,
		    maxs=2,
		    minalpha=10,
		    maxalpha=1e8,
		    minf0=-10,
		    maxf0=10) {
  my_likelihood <- function(x) {
    # Likelihood function
    
    # Break the initial vector into the 3 parameters (undoing any transformations)
    s <- x[1] 
    alpha <- exp(x[2])
    f0 <- 1/(1+exp(-x[3]))
    
    # Set the log-likelihood infinitely high, in case there is an error in calculating the log-likelihood
    ll <- Inf
    
    if (alpha <= 0) { # Make sure the alpha is valid
      return(-Inf)
    } else if (iffreq) { # Use different log-likelihood functions depending on whether iffreq is TRUE
      try(ll <- -get_llh_freq(prior, l, prec, tvec, nuvec, s, alpha, f0))
    } else {
      try(ll <- -get_llh(prior, l, tvec, nvec, bvec, prec, maxiter, verbose, s, alpha, f0))
    }
    #print(ll)
    return(ll)
  }
  # Initialise output
  output <- list(s.0 = 0, alpha.0=NA, f0.0 = NA, LL.0=NA, s=NA,alpha=NA,f0=NA,LL=NA)
  # Convert bvec and nvec into vector of sample frequencies
  nuvec <- bvec / nvec
  
  # Define a meaningless prior
  prior <- function(x){return(1)}
  
  # Get the number of time points
  l <- length(nuvec)
  
  # Calculate an initial guess for s
  s_guess <- max(c(mins,get_rough_s_guess(l, tvec, nuvec, prec)))
 
  # Calculate an initial guess for alpha
    alpha_guess <- log(max(c(minalpha,get_rough_alpha_guess(l, nvec, bvec, tvec,s=s_guess, minalpha, maxalpha))))

  # Set initial probability equal to initial frequency
  f0_guess <- nuvec[1]
  
  # Convert initial probability into log-odds
  f0_guess <- log(f0_guess/(1-f0_guess))

  # Create bounds for optimisation
  lowervec <- c(max(mins,s_guess-.5),max(log(minalpha),alpha_guess/2),minf0)
  uppervec <- c(min(maxs,s_guess+.5),min(log(maxalpha),alpha_guess*2),maxf0)
  lowervec2 <- c(mins,log(minalpha),minf0)
  uppervec2 <- c(maxs,log(maxalpha),maxf0)

  # Fix the inital frequncy if iffixedf0 == TRUE
  if (iffixedf0) {
    lowervec[3] <- f0_guess
    uppervec[3] <- f0_guess
    lowervec2[3] <- f0_guess
    uppervec2[3] <- f0_guess
  }

  if (!ifneut) {
    res1 <- nloptr(c(s_guess,alpha_guess,f0_guess),my_likelihood,lb=lowervec,ub=uppervec,
                     opts=list(algorithm="NLOPT_GN_MLSL_LDS",local_opts=list(algorithm='NLOPT_LN_SBPLX'),xtol_rel=.1,maxeval=500))
     # Use a local optimisation method to refine results from global search
    mod <- nloptr(res1$sol,my_likelihood,lb=lowervec2,ub=uppervec2,
                     opts=list(algorithm="NLOPT_LN_SBPLX",xtol_rel=prec,maxeval=maxiter))
    output$s <- mod$sol[1]
    output$alpha <- exp(mod$sol[2])
    output$f0 <- 1/(1+exp(-mod$sol[3]))
    output$LL <- mod$obj
  }
  # Neutral model: set s to 0
  s_guess<-0
   
  # Generate guess of alpha using s=0
  alpha_guess <- log(max(c(minalpha,get_rough_alpha_guess(l, nvec, bvec, tvec,s=s_guess,minalpha, maxalpha))))
   
  # Use bounds to force s to be 0
  lowervec[1] <- 0
  uppervec[1] <- 0
  lowervec2[1] <- 0
  uppervec2[1] <- 0

  # f0 and alpha cannot both be estimated for neutral models; so fix f0
  lowervec[3] <- f0_guess 
  uppervec[3] <- f0_guess
  lowervec2[3] <- f0_guess
  uppervec2[3] <- f0_guess
  # Initial optimisation does global search of narrow bounds around initial guesses to try and
  # escape local optima
  res1 <- nloptr(c(s_guess,alpha_guess,f0_guess),my_likelihood,lb=lowervec,ub=uppervec,
                   opts=list(algorithm="NLOPT_GN_MLSL_LDS",local_opts=list(algorithm='NLOPT_LN_SBPLX'),xtol_rel=.1,maxeval=500))
  # Use a local optimisation method to refine results from global search
  mod <- nloptr(res1$sol,my_likelihood,lb=lowervec2,ub=uppervec2,
                   opts=list(algorithm="NLOPT_LN_SBPLX",xtol_rel=prec,maxeval=maxiter))
  output$s.0 <- mod$sol[1]
  output$alpha.0 <- exp(mod$sol[2])
  output$f0.0 <- 1/(1+exp(-mod$sol[3]))
  output$LL.0 <- mod$obj
  return(output)
 }
