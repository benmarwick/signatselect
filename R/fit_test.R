#' frequency increment test (FIT) from Feder et al 2014
#'
#' Computes the frequency increment test (FIT) from Feder, A. F., Kryazhimskiy, S., & 
#' Plotkin, J. B. (2014). Identifying signatures of selection in genetic time 
#' series. Genetics, 196(2), 509-522. <https://doi.org/10.1534/genetics.113.158220>.  
#' The frequency increment test (FIT) rejects neutrality if the distribution of 
#' normalized variant frequency increments exhibits a mean that deviates 
#' significantly from zero. 
#' 
#' @param time an integer vector of time coordinates for the time-series 
#' @param v an integer vector of variant frequencies, i.e. counts of this variant at each time coordinate
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#'
#' @return A two-column data frame with the FIT statistic and its p-value
#' @export

fit <- function(time,
                v,
                alternative = "two.sided"){
  # number of time-steps
  q = length(time)
  # Create vector to hold fitness increments
  Y = rep(0,(q-1))
  # v is variant frequencies
  # Get t values from df
  t = time
  # Rescale increments according to definition
  for (i in c(2:q)) { # R indexes from 1 rather than 0
    Y[i-1] = (v[i] - v[i - 1])/sqrt(2 * v[i-1] * (1 - v[i-1]) * (t[i] - t[i-1]))
  }
  # Get t statistic from the rescaled fitness increments
  fit_stat = as.numeric(t.test(Y)$statistic)
  # Calculate the p-value for the test statistic:
  fit_p = t.test(Y,  alternative = c(alternative))$p.value
  
  # get output a a two-column data frame
  output <- data.frame(fit_stat = fit_stat, 
                       fit_p = fit_p)
  return(output)
}
