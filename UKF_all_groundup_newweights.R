unscentedkalmanf_all_newweights = function(y, f, B, Q, h, R, x0, P0, alpha, beta, kappa) {
  
  dx = length(x0)
  I = diag(dx)
  dy = nrow(y)
  T = ncol(y)
  N = NROW(x0)
  
  ##INITIALISE##
  N = dx
  
  lambda = alpha^2 * ( N + kappa) - N
  W0_m = (lambda )/ (N + lambda) 
  W0_c = W0_m + 1 - alpha^2 + beta
  W = rep(1 / (2 *(N +lambda)), 2 * N)
  
  Wm = as.matrix(c(W0_m, W))
  Wc = as.matrix(c(W0_c, W))
  
  
  
  x.p = matrix(0, nrow = dx, ncol = T)
  Yk = array(0, c(dy, (2 * dx + 1), T))
  Yk.p = matrix(0, nrow = dy, ncol = T)
  Xk.p = array(0, c(dx, (2 * dx + 1), T))
  P.p = array(0, c(dx, dx, T))
  x.f = matrix(0, nrow = dx, ncol = T)
  P.f = array(0, c(dx, dx, T))

  
  ##TIME 1##
  #CALCULATE SIGMAPOINTS#
  L = chol(P0)  # Cholesky decomposition matrix L

  Sigmapoints = matrix(x0, nrow = dx, ncol = (2 * dx + 1))
  
  for (i in 1:dx) {
    Sigmapoints[, i + 1] = x0 + sqrt(lambda + dx) * L[, i]
    Sigmapoints[, i + dx + 1] = x0 - sqrt(lambda + dx) * L[, i]
  }
  
  #PREDICTION#
  for (i in 1:(2 * dx + 1)) {
    Xk.p[, i, 1] = f(Sigmapoints[, i])
  }
  x.p[, 1] = Xk.p[, , 1] %*% Wm
  
  diff_matrix = Xk.p[, , 1] - matrix(x.p[, 1], nrow = dx, ncol = 2 * dx + 1)
  Wc_diag = diag(as.vector(Wc))
  #P.p[,, 1] = diff_matrix %*% Wc_diag %*% t(diff_matrix) + B(x0) %*% Q %*% t(B(x0))
  P.p[,, 1] = diff_matrix %*% Wc_diag %*% t(diff_matrix) + Q  #calculate P.p without using process matrix B
  
  
  
  for (i in 1:(2 * dx + 1)) {
    Yk[, i, 1] = h(Xk.p[, i, 1])
  }
  Yk.p[, 1] = rowSums(Yk[, , 1] * matrix(Wm, nrow = dy, ncol = 2 * dx + 1, byrow = TRUE))
  
  
  #UPDATE#
  diff_y_matrix = Yk[, , 1] - matrix(Yk.p[, 1], nrow = dy, ncol = 2 * dx + 1)
  Pyy = diff_y_matrix %*% Wc_diag %*% t(diff_y_matrix) + R
  
  
  diff_x_matrix = Xk.p[, , 1] - matrix(x.p[, 1], nrow = dx, ncol = 2 * dx + 1)
  Pxy = diff_x_matrix %*% Wc_diag %*% t(diff_y_matrix)
  
  K = Pxy %*% solve(Pyy)
  x.f[, 1] = x.p[, 1] + K %*% (y[, 1] - Yk.p[, 1])
  P.f[,, 1] = P.p[,, 1] - K %*% Pyy %*% t(K)
  
  
  for (t in 2:T) {
    # CALCULATE SIGMAPOINTS
    
    L = chol(P.f[,, t - 1])  
    
    Sigmapoints = matrix(x.f[, t - 1], nrow = dx, ncol = (2 * dx + 1))
    for (i in 1:dx) {
      Sigmapoints[, i + 1] = x.f[, t - 1] + sqrt(lambda + dx) * L[, i]
      Sigmapoints[, i + dx + 1] = x.f[, t - 1] - sqrt(lambda + dx) * L[, i]
    }
    
    # PREDICTION #
    for (i in 1:(2 * dx + 1)) {
      Xk.p[, i, t] = f(Sigmapoints[, i])
    }
    x.p[, t] = Xk.p[, , t] %*% Wm  # Mean prediction
    
    # Covariance Prediction
    diff_matrix = Xk.p[, , t] - matrix(x.p[, t], nrow = dx, ncol = 2 * dx + 1)
    Wc_diag = diag(as.vector(Wc))
    #P.p[,, t] = diff_matrix %*% Wc_diag %*% t(diff_matrix) + B(x.f[, t - 1]) %*% Q %*% t(B(x.f[, t - 1]))
    P.p[,, t] = diff_matrix %*% Wc_diag %*% t(diff_matrix) + Q 
    
    # TRANSFORM SIGMA POINTS INTO MEASUREMENT SPACE #
    for (i in 1:(2 * dx + 1)) {
      Yk[, i, t] = h(Xk.p[, i, t])
    }
    Yk.p[, t] = rowSums(Yk[, , t] * matrix(Wm, nrow = dy, ncol = 2 * dx + 1, byrow = TRUE))
    
    # UPDATE STEP #
    diff_y_matrix = Yk[, , t] - matrix(Yk.p[, t], nrow = dy, ncol = 2 * dx + 1)
    Pyy = diff_y_matrix %*% Wc_diag %*% t(diff_y_matrix) + R
    
    diff_x_matrix = Xk.p[, , t] - matrix(x.p[, t], nrow = dx, ncol = 2 * dx + 1)
    Pxy = diff_x_matrix %*% Wc_diag %*% t(diff_y_matrix)
    
    
    
    K = Pxy %*% solve(Pyy) #compute kalman gain
    x.f[, t] = x.p[, t] + K %*% (y[, t] - Yk.p[, t])
    P.f[,, t] = P.p[,, t] - K %*% Pyy %*% t(K)
    

  }
  
  return(list(
    x.p = x.p,   # Predicted state means
    P.p = P.p,   # Predicted state covariances
    x.f = x.f,   # Filtered state means
    P.f = P.f    # Filtered state covariances
  ))
}
