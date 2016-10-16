get_muvec <-
function(lgth, x0, s, tvec) {
  # Get the predicted values from the logistic expansion
  gvec <- sapply(tvec, function(x){get_g(x0,s,x)})
  
  # For every time point, generate a function that predicts the next probability based on 
  # the predicted value from the logistic plus predicted gaussian noise
  mufun <- function(i) {
    dt <- tvec[i+1] - tvec[i]
    return(function(xprime){
		   return(gvec[i+1] + (xprime - gvec[i]) * get_M(gvec[i], s, dt))
})
  }
  return(sapply(1:(lgth-1), mufun))
}
