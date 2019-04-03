#' Infer selection/drift parameters for time-series allele-frequency data
#'
#' \code{tsinfer} computes the selection coefficient, population size and initial frequency for the logistic using a Gaussian approximation of the Kimura expression for time-series variant-frequency data. 
#' 
#' Essential arguments are tvec (time point labels for time series starting at 0), bvec (number of new variants at each time point), and nvec (total number of samples at each time point).
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
#' @param maxs maximum s value to consider (default=2)
#' @param minalpha minimum alpha value to consider (default=10)
#' @param maxalpha maximum alpha value to consider (default=1e8)
#' @param minf0 minimum f0 value to consider in log-odds (default=-10)
#' @param maxf0 maximum f0 value to consider in log-odds (default=10)
#' @return A list with the neutral and non-neutral parameter values and associated log-likelihoods
#' The output of the execution is the following list:
#' s = selection coefficient for non-neutral model
#' alpha = population size for non-neutral model
#' f0 = initial frequency for the logistic in non-neutral model
#' LL = log-likelihood of non-neutral model
#' s.0 = 0 (selection coefficient for the neutral model)
#' alpha.0 = population size for the neutral model
#' f0.0 = initial frequency for the logistic in neutral model
#' LL.0 = log-likelihood of neutral model
#' 
#' @export
tsinfer <- function (tvec,
                     bvec,
                     nvec,
                     maxiter = 200,
                     prec = 1e-06,
                     iffreq = FALSE,
                     ifneut = FALSE,
                     iffixedf0 = FALSE,
                     verbose = FALSE,
                     mins = -2,
                     maxs = 2,
                     minalpha = 10,
                     maxalpha = 1e+08,
                     minf0 = -10,
                     maxf0 = 10) {
  my_likelihood <- function(x) {
    s <- x[1]
    alpha <- exp(x[2])
    f0 <- 1 / (1 + exp(-x[3]))
    ll <- Inf
    if (alpha <= 0) {
      return(-Inf)
    }
    else if (iffreq) {
      try(ll <- -get_llh_freq(prior, l, prec, tvec, nuvec,
                                        s, alpha, f0))
    }
    else {
      try(ll <- -get_llh(prior,
                                   l,
                                   tvec,
                                   nvec,
                                   bvec,
                                   prec,
                                   maxiter,
                                   verbose,
                                   s,
                                   alpha,
                                   f0))
    }
    return(ll)
  }
  output <- list(
    s.0 = 0,
    alpha.0 = NA,
    f0.0 = NA,
    LL.0 = NA,
    s = NA,
    alpha = NA,
    f0 = NA,
    LL = NA
  )
  nuvec <- bvec / nvec
  prior <- function(x) {
    return(1)
  }
  l <- length(nuvec)
  s_guess <- max(c(mins, get_rough_s_guess(l, tvec, nuvec,
                                                     prec)))
  alpha_guess <- log(max(c(
    minalpha,
    get_rough_alpha_guess(l,
                                    nvec, bvec, tvec, s = s_guess, minalpha, maxalpha)
  )))
  f0_guess <- nuvec[1]
  f0_guess <- log(f0_guess / (1 - f0_guess))
  lowervec <- c(max(mins, s_guess - 0.5), max(log(minalpha),
                                              alpha_guess / 2), minf0)
  uppervec <- c(min(maxs, s_guess + 0.5), min(log(maxalpha),
                                              alpha_guess * 2), maxf0)
  lowervec2 <- c(mins, log(minalpha), minf0)
  uppervec2 <- c(maxs, log(maxalpha), maxf0)
  if (iffixedf0) {
    lowervec[3] <- f0_guess
    uppervec[3] <- f0_guess
    lowervec2[3] <- f0_guess
    uppervec2[3] <- f0_guess
  }
  if (!ifneut) {
    res1 <- nloptr::nloptr(
      c(s_guess, alpha_guess, f0_guess),
      my_likelihood,
      lb = lowervec,
      ub = uppervec,
      opts = list(
        algorithm = "NLOPT_GN_MLSL_LDS",
        local_opts = list(algorithm = "NLOPT_LN_SBPLX"),
        xtol_rel = 0.1,
        maxeval = 500
      )
    )
    mod <- nloptr::nloptr(
      res1$sol,
      my_likelihood,
      lb = lowervec2,
      ub = uppervec2,
      opts = list(
        algorithm = "NLOPT_LN_SBPLX",
        xtol_rel = prec,
        maxeval = maxiter
      )
    )
    output$s <- mod$sol[1]
    output$alpha <- exp(mod$sol[2])
    output$f0 <- 1 / (1 + exp(-mod$sol[3]))
    output$LL <- mod$obj
  }
  s_guess <- 0
  alpha_guess <- log(max(c(
    minalpha,
    get_rough_alpha_guess(l,
                                    nvec, bvec, tvec, s = s_guess, minalpha, maxalpha)
  )))
  lowervec[1] <- 0
  uppervec[1] <- 0
  lowervec2[1] <- 0
  uppervec2[1] <- 0
  lowervec[3] <- f0_guess
  uppervec[3] <- f0_guess
  lowervec2[3] <- f0_guess
  uppervec2[3] <- f0_guess
  res1 <- nloptr::nloptr(
    c(s_guess, alpha_guess, f0_guess),
    my_likelihood,
    lb = rep(-Inf, 3),
    ub = rep(Inf, 3),
    opts = list(
      algorithm = "NLOPT_GN_MLSL_LDS",
      local_opts = list(algorithm = "NLOPT_LN_SBPLX"),
      xtol_rel = 0.1,
      maxeval = 500
    )
  )
  mod <-
    nloptr::nloptr(
      res1$sol,
      my_likelihood,
      lb = lowervec2,
      ub = uppervec2,
      opts = list(
        algorithm = "NLOPT_LN_SBPLX",
        xtol_rel = prec,
        maxeval = maxiter
      )
    )
  output$s.0 <- mod$sol[1]
  output$alpha.0 <- exp(mod$sol[2])
  output$f0.0 <- 1 / (1 + exp(-mod$sol[3]))
  output$LL.0 <- mod$obj
  return(output)
}
