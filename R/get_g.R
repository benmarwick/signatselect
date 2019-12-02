get_g <-
function(x0, s, t) {
  if (s == 0) {
    # If the selection coefficient is 0, the logistic predicts no change
    return(x0)
  } else {
    # Otherwise, use the logistic formula to predict the new value
    return(x0/(x0+(1-x0)*exp(-s*t)))
  }
}
