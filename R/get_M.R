get_M <-
function(x0, s, t) {
  if (s == 0) {
    return(1)
  } else {
    exp(-s * t)/(x0 + (1-x0)*exp(-s*t))^2
  }
}
