vec2mat_tridiag <-
function(m_diag,m_above,m_below){
  n <- length(m_diag)
  mat <- matrix(0,n,n)
  for (i in 1:n-1) {
    mat[i,i]<-m_diag[i]
    mat[i,i+1]<-m_above[i]
    mat[i+1,i]<-m_below[i]
  }
  mat[n,n]<-m_diag[n]
  return(mat)
}
