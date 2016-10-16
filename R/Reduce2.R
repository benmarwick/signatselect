Reduce2 <-
function(f, a, b, c) {
  r <- c
  for (i in 1:length(a)) {
    r <- f(a, b, r)
  }
  return(r)
}
