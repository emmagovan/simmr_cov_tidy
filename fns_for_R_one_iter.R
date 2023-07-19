
# sim_theta ---------------------------------------------------------------
sim_theta <- function(S = 100, lambda) {
  
  ## Create a loop instead to do this I think? Will need to make an array?
  mean_beta <- matrix(0, nrow = (n_covariates), ncol = K)
  for(i in 1:(n_covariates)){
    mean_beta[i,] <- lambda[lambda_index$mu_beta[i,]]
  }
  
  chol_prec_beta <- array(data = 0, dim = c(K, K, (n_covariates)))
  
  for(i in 1:(n_covariates)){
    chol_prec_beta[,,i][upper.tri(chol_prec_beta[,,i], diag = TRUE)] <-
      lambda[lambda_index$sigma_beta[i,]]
  }
  a<-array(NA, dim =c(S, K, (n_covariates)))
  thetabeta<-matrix(NA, ncol = (n_covariates) * K, nrow = S)
  
  
  
  for(i in 1:(n_covariates)){
    
    thetabeta[,(1+(i-1)*K):((i)*K)] = t(rMVNormC(S, mu = mean_beta[i,], U = chol_prec_beta[,,i]))
    
  }
  
  
  
  
  # This is a placeholder for more advanced code using chol_prec directly
  theta <- cbind(
    thetabeta,
    matrix(
      rgamma(S * n_isotopes,
             shape = lambda[lambda_index$c],
             rate = lambda[lambda_index$d]
      ),
      nrow = S,
      ncol = n_isotopes,
      byrow = TRUE
    )
  )
  
  return(theta)
}



# h -----------------------------------------------------------------------

# Log of likelihood added to prior
h <- function(theta) {
  # Create betas and sigma
  beta <- matrix(theta[1:((n_covariates) * K)], nrow = (n_covariates), byrow = TRUE)
  tau <- theta[((n_covariates) * K +1):(((n_covariates) * K)+n_isotopes)]
  f <- matrix(NA, ncol = n, nrow = K) 
  
  #Need to double check that this maths is right!!
  
  # for (k in 1:K) {
  #  f[,k] <-  (x_scaled %*% beta[,k])
  # }
  # 
  
  f = x_scaled %*% beta
  
  p <- matrix(NA, ncol = K, nrow = n)
  
  for (i in 1:n) {
    p[i, ] <- exp(f[i,1:K]) / (sum((exp(f[i,1:K]))))
  }
  
  ## Im not sure this bit needs to be looped over?
  mumat = matrix(c(rep(0, n*n_isotopes)), nrow = n, ncol = n_isotopes)
  hold <- 0
  
  for (i in 1:n) {
    for (j in 1:n_isotopes) {
      hold <- hold + sum(dnorm(y[i, j],
                               mean = sum(p[i, ] * q[, j] * (mu_s[, j] + mu_c[, j])) /
                                 sum(p[i, ] * q[, j]),
                               sd = sqrt(sum(p[i, ]^2 * q[, j]^2 * (sigma_s[, j]^2 + sigma_c[, j]^2)) /
                                           sum(p[i, ]^2 * q[, j]^2) + 1/tau[j]),
                               log = TRUE
                               
      ))
      mumat[i,j] = sqrt(sum(p[i, ]^2 * q[, j]^2 * (sigma_s[, j]^2 + sigma_c[, j]^2)) /
                          sum(p[i, ]^2 * q[, j]^2) + 1/tau[j])
      
      
    }
  }
  
  beta_sum <- 0
  #change to k
  for (i in 1:(n_covariates)){
    beta_sum = beta_sum +sum(dnorm(beta[i,], mu_beta_zero[i,], sigma_beta_zero[i,], log = TRUE))
  }
  
  return(hold + beta_sum + sum(dgamma(1/sqrt(tau), shape = c_0, rate = d_0, log = TRUE)))
}

# 
# n_sources = K 
# gammaprior = 0
# 
# 
# for(i in 1:n_isotopes){
#   gammaprior = gammaprior + c_0[i] * log(d_0[i]) - (lgamma(c_0[i])) +
#     (c_0[i] - 1) * theta[1,i+n_sources*n_covariates]-
#     d_0[i] * theta[1,i+n_sources*n_covariates];
# }


# log_q -------------------------------------------------------------------

log_q <- function(lambda, theta) {
  
  ## Create a loop instead to do this I think? Will need to make an array?
  mean_beta <- matrix(0, nrow = (n_covariates), ncol = K)
  for(i in 1:(n_covariates)){
    mean_beta[i,] <- lambda[lambda_index$mu_beta[i,]]
  }
  
  chol_prec_beta <- array(data = 0, dim = c(K, K, (n_covariates)))
  
  for(i in 1:(n_covariates)){
    chol_prec_beta[,,i][upper.tri(chol_prec_beta[,,i], diag = TRUE)] <-
      lambda[lambda_index$sigma_beta[i,]]
  }
  a<-array(NA, dim =c(S, K, (n_covariates)))
  thetabeta<-matrix(NA, ncol = (n_covariates) * K, nrow = S)
  
  shape_tau <- lambda[lambda_index$c]
  rate_tau <- lambda[lambda_index$d]
  
  # Extract alpha, beta and sigma from theta
  beta <- matrix(theta[1:((n_covariates) * K)], nrow = (n_covariates), ncol = K,  byrow = TRUE)
  tau <- theta[((n_covariates) * K +1):(((n_covariates) * K)+n_isotopes)]
  
  # This is a placeholder for more advanced code using chol_prec directly
  # prec <- crossprod(chol_prec)
  p_mat <- matrix(NA, nrow = n_covariates, ncol = K) #row for each beta
  for(l in 1:(n_covariates)){
    p_mat[l,] <- (matrix(beta[l,] - mean_beta[l,], nrow = 1) %*% t(chol_prec_beta[,,l]))
  }
  
  sum_p = 0
  a<-c(NA, NA)
  
  for(l in 1:(n_covariates)){
    sum_p = sum_p - 0.5 * K * log(2 * pi)-
      0.5 * sum(log(diag(chol_prec_beta[,,l])))-
      0.5 * matrix(p_mat[l,], nrow = 1) %*% (p_mat[l,])
    a[l] <- matrix(p_mat[l,], nrow = 1) %*% (p_mat[l,])
    
  }
  
  return(sum_p + sum(dgamma(1/sqrt(tau),
                            shape = shape_tau,
                            rate = rate_tau,
                            log = TRUE
  ))
  )
}

# This contains all the generic functions for FF VB

# Function to estimate different between joint and variational approx
h_lambda <- function(lambda, theta, y) {
  return(h(theta) - log_q(lambda, theta))
}
# h_lambda(lambda, theta[1,], y)

# Nable LB is the mean of delta_lqlt element-wise multiplied by h_lambda
nabla_LB <- function(lambda, theta, c = rep(0, length(lambda))) {
  big_delta_lqlt <- t(apply(theta, 1, delta_lqlt, lambda = lambda))
  big_h_lambda <- t(apply(theta, 1, h_lambda, lambda = lambda, y = y))
  big_h_lambda_rep <- matrix(rep(big_h_lambda, length(lambda)),
                             nrow = nrow(theta),
                             ncol = length(lambda)
  )
  big_c <- matrix(rep(c, nrow(theta)),
                  ncol = length(c),
                  nrow = nrow(theta),
                  byrow = TRUE
  )
  return(colMeans(big_delta_lqlt * (big_h_lambda_rep - big_c)))
}
# nabla_LB(lambda, theta)

# Now the control variate
control_var <- function(lambda, theta) {
  # Get delta log q
  big_delta_lqlt <- t(apply(theta, 1, delta_lqlt, lambda = lambda))
  # Get h_lambda
  big_h_lambda <- t(apply(theta, 1, h_lambda, lambda = lambda, y = y))
  # Make it bigger
  big_h_lambda_rep <- matrix(rep(big_h_lambda, length(lambda)),
                             nrow = nrow(theta),
                             ncol = length(lambda)
  )
  # Now get delta log q times h
  big_nabla_log_q_h <- big_delta_lqlt * big_h_lambda_rep
  # Return the diagonals of the covariance and scale
  return(diag(cov(big_nabla_log_q_h, big_delta_lqlt)) / apply(big_delta_lqlt, 2, var))
}
# c <- control_var(lambda, theta)

# LB estimate
LB_lambda <- function(lambda, theta) {
  mean(apply(theta, 1, h_lambda, lambda = lambda))
}
# LB_lambda(lambda, theta)

# Empirical version of derivative - comment this out if you want to create it
# properly
delta_lqlt <- function(lambda, theta, eps = 0.001) {
  k <- length(lambda)
  ans <- rep(NA, k)
  for (i in 1:k) {
    d <- rep(0, k)
    d[i] <- eps
    ans[i] <- (log_q(lambda + d, theta) - log_q(lambda - d, theta)) / (2 * max(d))
  }
  return(ans)
}
# delta_lqlt(lambda, theta[1,])

# Natural

# Run the VB function
run_VB_one_iter <- function(lambda, # Starting value of lambda
                   S = 100, # Number of samples to take
                   P = 10, # Maximum patience before you stop
                   beta_1 = 0.9, # Learning rates
                   beta_2 = 0.9, # Learning rates
                   tau = 1000, # Iteration at which learning rate starts to decrease
                   eps_0 = 0.1, # Raw learning rate multiplier
                   t_W = 50, # Time window for working out convergence
                   t = 1 #set iteration
) {
  
  # Starting
 # set.seed(12345)
  theta <- sim_theta(S, lambda)
  c <- control_var(lambda, theta)
  g_0 <- nabla_LB(lambda, theta)
  nu_0 <- g_0^2
  g_bar <- g_0
  nu_bar <- nu_0
  
  # Set up
  #t <- 1
  patience <- 0
  stop <- FALSE
  LB <- rep(NA, t_W)
  max_LB_bar <- -Inf
  
  #  while (!stop) {
  #if (t %% 10 == 0) print(t)
  print("t")
  print(t)
  # Generate new samples
 # set.seed(12345)
  theta <- sim_theta(S, lambda)
  
  # Compute g_t
  g_t <- nabla_LB(lambda, theta, c)
  print("g_t")
  print(g_t)
  
  # Compute new control variate
  c <- control_var(lambda, theta)
  print("c")
  print(c)
  
  # Update the learning rates
  nu_t <- g_t^2
  print("nu_t")
  print(nu_t)
  g_bar <- beta_1 * g_bar + (1 - beta_1) * g_t
  nu_bar <- beta_2 * nu_bar + (1 - beta_2) * nu_t
  
  # Update the learning rate
  alpha_t <- min(eps_0, eps_0 * tau / t)
  print("alpha_t")
  print(alpha_t)
  
  # Update lambda
  lambda <- lambda + alpha_t * g_bar / sqrt(nu_bar)
  print("lambda")
  print(lambda)
  
  # Compute the moving average LB if out of warm-up
  if (t <= t_W) {
    # Compute a new lower bound estimate
    LB[t] <- LB_lambda(lambda, theta)
  } else {
    LB[1:(t_W - 1)] <- LB[2:t_W]
    LB[t_W] <- LB_lambda(lambda, theta)
    LB_bar <- mean(LB)
    max_LB_bar <- max(max_LB_bar, LB_bar)
    if (LB_bar >= max_LB_bar) {
      patience <- 0
    } else {
      patience <- patience + 1
    }
  }
  
  if (patience > P) {
    print("Completed!")
    stop <- TRUE
  }
  t <- t + 1
  #}
  return(lambda)
}

rMVNormC <- function(n, mu, U){
  p <- length(mu)
  Z <- matrix(rnorm(p*n), p, n)
  # U <- chol(Omega) # By default R's chol fxn returns upper cholesky factor
  X <- backsolve(U, Z) # more efficient and stable than actually inverting
  X <- sweep(X, 1, mu, FUN=`+`)
  return(X)
}







