get_rough_alpha_guess <-
function(lgth, nvec, bvec, tvec, s, minalpha, maxalpha){
  # Create a vector of frequencies for each time point
  nuvec <- bvec/nvec
  # Initial x is different depending on the selection coefficient
  if (s == 0) {
    init_x <- mean(nuvec)
  } else {
    init_x <- nuvec[1]
  }
  # Get the estimated tradjectory of the logistic
  gvec <- sapply(tvec[-1], function(x){get_g(init_x,s,x)})
  # Scale the squared differences between the predicted and empirical frequencies by the odds of the predicted frequencies
  tmp_vec <- mapply(function(x,y){(x-y)^2/y/(1-y)},nuvec[-1],gvec)
  # Function to sum and count the predicted population size
  tmpf <- function(m,x){
    if (is.nan(x)) {
      return(m)
    } else {
      return(list(s=m$s+x,c=m$c+1))
    }
  }
  tmpval <- Reduce(tmpf,tmp_vec,list(s=0,c=0))
  # Extract the sum of the estimates
  sum_inv_n <- tmpval[[1]]
  # Extract the counts of the estimates
  cnt <- tmpval[[2]]
  # If the sum is 0 return half of the maximum population size
  if (sum_inv_n == 0) {
    return(maxalpha/2)
  # Otherwise return the larger of 2 times the minimum size or the count divided by the sum of the estimates
  } else {
    return(max(minalpha * 2, cnt/sum_inv_n))
  }
}
