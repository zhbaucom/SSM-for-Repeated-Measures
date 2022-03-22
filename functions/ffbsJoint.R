ffbs.joint = function(y,V,W,m0,C0, timeDiff){
  
  ##### Variable initialization
  
  T = ncol(y); n = nrow(y)
  a = matrix(0, n, T); R = array(0, dim = c(n, n, T)) 
  m = matrix(0, n, T); C = array(0, dim = c(n, n, T)) ; B = array(0, dim = c(n, n, T-1))
  H = array(0, dim = c(n, n, T-1)); mm = matrix(0, n, T); CC = array(0, dim = c(n, n, T))
  x = matrix(0, n, T); llike = 0.0 
  
  ##### Forward kalman filter

  for (t in 1:T){ 
    if(t==1){
      a[,1] = m0; R[,,1] = C0 + W
    }else{
      a[,t] = m[,t-1]; R[,,t] = C[,,t-1] + timeDiff[t-1] * W
    }
    f = a[,t]
    Q = R[,,t] + V 
    A = R[,,t] %*% solve(Q)
    Av = A %*% (y[,t]-f)
    m[,t] = a[,t]+Av 
    C[,,t] = R[,,t]- A %*% R[,,t]

    
    if(t > 1){
      B[,,t-1] = C[,,t-1] %*% solve(R[,,t])
      H[,,t-1] = C[,,t-1]- B[,,t - 1] %*% tcrossprod(R[,,t],B[,,t - 1])
    }
    
    ##### Likelihood currently not supported
    
    # llike = llike + sum(dmvnorm(y[,t],f,(Q),log=TRUE), na.rm = TRUE) 
  }

  ##### Backward sampler
  
  mm[,T] = m[,T]; CC[,,T] = C[,,T]
  x[,T] = rmvnorm(1, mean = m[,T], sigma = (C[,,T]), checkSymmetry = FALSE) 

  for (t in (T-1):1){ 
    
    ##### Smoothed values currently not supported
    
    # mm[,t] = m[,t] + solve(R[,,t+1]) %*% C[,,t] %*% (mm[,t+1]-a[,t+1])
    # CC[,,t] = C[,,t] - solve(crossprod(R[,,t+1])) %*% crossprod(C[,,t])%*%(R[,,t+1]-CC[,,t+1])

    x[,t] = rmvnorm(1,mean = m[,t]+B[,,t]%*%(x[,t+1]-a[,t+1]), sigma = H[,,t], checkSymmetry = FALSE)
  } 
  
  
  
  return(list(x=x,m=m,C=C,mm=mm,CC=CC,llike=llike))
}
