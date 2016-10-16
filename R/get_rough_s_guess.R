get_rough_s_guess <-
function(lgth, tvec, nuvec) {
  # Get a vector of the time distance between each adjacent pair of points
  dt_vec <- tvec[-1]-tvec[-length(tvec)]
  # Get a vector of the change in frequency between each adjacent pair of points
  dnu_vec <- nuvec[-1]-nuvec[-length(nuvec)]
  # Get a vector of the mid-point in frequency between each adjacent pair of points
  nu_mid_vec <- mapply(function(a,b){return((a+b)/2)},nuvec[-1],nuvec[-length(nuvec)])
  # Scale the frequency changes by time changes
  dnudt_vec <- dnu_vec / dt_vec
  # Guess that the s for each pair of points equals the time scaled frequency divided by the mean odds for that pair
  s_est_vec <- mapply(function(a,b){return(a/b/(1-b))},dnudt_vec,nu_mid_vec)
  # Temporary reduction function for summing up the s estimates
  tmpf <- function(m,x){
    if (is.nan(x)) {
      return(m)
    } else {
      return(list(s=m$s+x,c=m$c+1))
    }
  }
  # Sum the s estimates
  tmpval <- Reduce(tmpf,s_est_vec,list(s=0,c=0))
  # Extract the output
  # Sum of all s estimates
  sum_s <- tmpval[[1]]
  # The number of s estimates
  cnt <- tmpval[[2]]
  # s_guess equals the average s estimate
  s_guess <- sum_s/cnt
  # If s_guess is smaller than the precision, just return 0
  if (abs(s_guess) <= def_prec) {
    return(0)
  } else {
    return(s_guess)
  }
}
