multivariate_gaussian_diag_logpdf <-
function(mu,cov_d,x){
  l <- length(mu)
  x1 <- x
  x2 <- cov_d
  x1 <- x1 - mu
  x2 <- x2 * x1
  exparg <- x1 %*% x2
  logdet <- Reduce(function(sumx, n){sumx - (log(n))},cov_d,0)
  return(-(exparg+l*(log(2*pi))+logdet)/2)
}
