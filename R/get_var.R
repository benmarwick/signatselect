get_var <-
function(x0, s, t) {
  y0 <- 1 - x0
  a <- 2 * x0 * y0 * t
  if (s == 0) {
    return(a)
  } else {
    b  <- x0^2 / s * (exp(s * t))
    c <- y0^2 / s * (exp(-s * t))
    d <- (y0^2 - x0^2)/s
    e <- (2 + s) * x0 * y0
    return(get_M(x0, s, t)^2 * e * (a + b - c + d))
  }
}
