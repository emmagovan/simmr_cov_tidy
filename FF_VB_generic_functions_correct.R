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
run_VB <- function(lambda, # Starting value of lambda
                   S = 100, # Number of samples to take
                   P = 10, # Maximum patience before you stop
                   beta_1 = 0.9, # Learning rates
                   beta_2 = 0.9, # Learning rates
                   tau = 1000, # Iteration at which learning rate starts to decrease
                   eps_0 = 0.1, # Raw learning rate multiplier
                   t_W = 50 # Time window for working out convergence
) {

  # Starting
  theta <- sim_theta(S, lambda)
  c <- control_var(lambda, theta)
  g_0 <- nabla_LB(lambda, theta)
  nu_0 <- g_0^2
  g_bar <- g_0
  nu_bar <- nu_0

  # Set up
  t <- 1
  patience <- 0
  stop <- FALSE
  LB <- rep(NA, t_W)
  max_LB_bar <- -Inf

  while (!stop) {
    if (t %% 10 == 0) print(t)

    # Generate new samples
    theta <- sim_theta(S, lambda)

    # Compute g_t
    g_t <- nabla_LB(lambda, theta, c)

    # Compute new control variate
    c <- control_var(lambda, theta)

    # Update the learning rates
    nu_t <- g_t^2
    g_bar <- beta_1 * g_bar + (1 - beta_1) * g_t
    nu_bar <- beta_2 * nu_bar + (1 - beta_2) * nu_t

    # Update the learning rate
    alpha_t <- min(eps_0, eps_0 * tau / t)

    # Update lambda
    lambda <- lambda + alpha_t * g_bar / sqrt(nu_bar)

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
  }
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