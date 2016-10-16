get_log_scale_constant <-
function(x){
  return(Reduce(function(sumx, n){
    sumx - (log(n+1))
    },x,0))
}
