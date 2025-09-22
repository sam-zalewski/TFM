extendedkalmanf_all = function(y, f, B, Q, h, R, x0, P0, alpha) {
  
  library(numDeriv)
  
  dy = nrow(y)
  T = ncol(y)
  dx = length(x0)
  I = diag(dx)
  
  ##INITIALISE##
  x.p = matrix(0, nrow = dx, ncol = T)
  P.p = array(0, c(dx, dx, T))
  x.f = matrix(0, nrow = dx, ncol = T)
  P.f = array(0, c(dx, dx, T))

  
  ##TIME 1##
  #PREDICTION#
  F = jacobian(f, x0)
  x.p[, 1] = f(x0)  
  P.p[,, 1] = F %*% P0 %*% t(F) + Q  #estimate process covariiance without additional matrix B
  
  # UPDATE
  H = jacobian(h, x.p[, 1]) 
  nu = y[, 1] - h(x.p[, 1])  #innovation: actual - expected measurement
  S = H %*% P.p[,, 1] %*% t(H) + R   #calculate innovation variance
  K = P.p[,, 1] %*% t(H) %*% solve(S)
  x.f[, 1] = x.p[, 1] + K %*% nu
  P.f[,, 1] = (I - K %*% H) %*% P.p[,, 1]
  
  ##RECURSION TIME 2:T##
  for (tt in 2:T) {  # safer to use 'tt' instead of 't'
    # PREDICTION
    F = jacobian(f, x.f[, tt - 1])
    x.p[, tt] = f(x.f[, tt - 1])  # Apply function directly
    P.p[,, tt] = F %*% P.f[,, tt - 1] %*% t(F) + Q  #covariance estimate
    H = jacobian(h, x.p[, tt])  
    nu = y[, tt] - h(x.p[, tt])  #inovation
    S = H %*% P.p[,, tt] %*% t(H) + R    #calculate innovation variance
    K = P.p[,, tt] %*% t(H) %*% solve(S) #kalman gain
    x.f[, tt] = x.p[, tt] + K %*% nu
    P.f[,, tt] = (I - K %*% H) %*% P.p[,, tt]
  }
  
return(list(
    x.p = x.p, #predicted state at each time 
    P.p = P.p, #predicted state covariance at each step
    x.f = x.f,  #filtered mean at each step
    P.f = P.f, #filtered step 

  ))

}