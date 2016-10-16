find_min_interval <-
function(func, min, max, delta0, mindelta, deltafact, verbose=FALSE){
  ff <- function(x, a, b, fa, fb, delta){
    fx <- func(x)
    if (verbose){
      print(paste(a, b, fa, fb,x,fx,sep=' '))
    }
   if (is.nan(b)) {
     if (is.nan(fx)) {
       newdelta <- delta / 2
       return(ff(a+newdelta,a,NaN,fa,NaN,newdelta))
     } else {
       if (fx < fa) {
         return(ff(x+delta, a, x, fa, fx, delta))
       } else {
         newdelta <- delta / 2
         return(ff(a+newdelta,a,NaN,fa,NaN,newdelta))
       }  
     }
   } else if (b == Inf) {
     return(list(a, Inf, Inf))
   } else if (fx > fb) {
     return (list(a, x, b)) 
   } else if (x >= max) {
     return(ff(Inf, b, x, fb, fx, 0))
   } else {
     if (delta < max - min) {
       newdelta <- delta * deltafact
     } else {
       newdelta <- delta
     }
     return(ff(x+delta,b,x,fb,fx,newdelta))
   }
  }
  if (verbose) {
    print("a b c f(a) f(b) f(c) | x f(x)\n")
  }
  return(ff(min+delta0,min,NaN,func(min),NaN,delta0))
}
