#using everything rcpp EXCEPT sim_theta


# Load in package
library(simmr)
library(readxl)
library(tidyverse)
library(R2jags)

Rcpp::sourceCpp("run_VB.cpp")


# Extract data ------------------------------------------------------------
path <- system.file("extdata", "geese_data.xls", package = "simmr")
geese_data <- lapply(            excel_sheets(path), read_excel, path = path)

# Just use time point 1 for now
consumer <- geese_data[[1]] |> filter(Time == 1)
sources <- geese_data[[2]]
TEFs <- geese_data[[3]]
conc <- geese_data[[4]]

# Put them into the right names
n <- nrow(consumer)
n_isotopes <- 2
K <- nrow(sources)
n_sources = K
mu_s <- sources[, c(2, 3)] 
sigma_s <- sources[, c(4, 5)]
mu_c <- TEFs[, c(2, 3)]
sigma_c <- TEFs[, c(4, 5)]
q <- conc[, c(2:3)]

concentrationmeans = as.matrix(q)
sourcemeans = as.matrix(mu_s)
correctionmeans = as.matrix(mu_c)
corrsds = as.matrix(sigma_c)
sourcesds = as.matrix(sigma_s)
beta_prior = 0.001

########## SET THESE
#   x <- matrix(c(consumer$Skull), 
# ncol = 1)
# x_scaled <- cbind(matrix(rep(1, nrow(x)), ncol = 1), scale(x))
x_scaled <- matrix(c(rep(1, 9)), ncol = 1)

n_covariates <- (ncol(x_scaled))
y <- consumer |>
  select(d13C_Pl, d15N_Pl) |>
  as.matrix()
# y <- consumer |>
#   select(d13C_Pl) |>
#   as.matrix()

# Get the data into simmr
simmr_groups = simmr_load(mixtures=as.matrix(y),
                          source_names=unlist(sources[,1]),
                          source_means=as.matrix(sources[,c(2,3)]),
                          source_sds=as.matrix(sources[,c(4,5)]),
                          correction_means=as.matrix(TEFs[,c(2,3)]),
                          correction_sds=as.matrix(TEFs[,c(4,5)]),
                          concentration_means = as.matrix(conc[,2:3]))


#-------------------- FFVB variables-----------------------
#Don't actually need to edit anything here, just run it


# Variables for the FFVB
S <- 100
mu_beta_zero <- matrix(c(rep(0, K * (n_covariates))), 
                       nrow = (n_covariates), 
                       ncol = K)
#n_covariates + 1 for alpha (or beta_0) and then 1 beta for each covariate
sigma_beta_zero <- matrix(c(rep(1, K * (n_covariates))), 
                          nrow = (n_covariates), 
                          ncol = K)

n_isotopes <- ncol(mu_c)
c_0 <- c(rep(3, n_isotopes)) #Change to 0.001
d_0 <- c(rep(1, n_isotopes))
beta_lambda<-c(rep(0, K),rep(1, K * (K + 1) / 2))

lambda <- c(
  rep(beta_lambda, n_covariates),
  rep(3, n_isotopes), #shape
  rep(1, n_isotopes) #rate
)


# function to extract lambdas --------------------------------------------
lambda_extract <- function(n_covariates, K, n_isotopes){
  mat_size = K * (K+1) /2
  mu_beta = matrix(data = NA, nrow = (n_covariates), ncol = K)
  sigma_beta = matrix(data = NA, nrow = (n_covariates), ncol = mat_size)
  
  for(i in 1:(n_covariates)){
    mu_beta[i,] = ((i-1) * mat_size + (i-1) * K +1):((i-1) * mat_size + (i-1) * K + K)
    sigma_beta[i,] = ((i-1) * mat_size + (i) * K +1): ((i-1) * mat_size + (i) * K +mat_size)
  }
  
  c = (sigma_beta[n_covariates, mat_size] + 1):(sigma_beta[n_covariates, mat_size] + n_isotopes)
  d = (sigma_beta[n_covariates, mat_size] + n_isotopes + 1):(sigma_beta[n_covariates, mat_size] + 2 * n_isotopes)
  
  return(list(mu_beta = mu_beta,
              sigma_beta = sigma_beta,
              c = c,
              d = d
  ))
}
#just use here for now - then build in to r
lambda_index <- lambda_extract(n_covariates, K, n_isotopes)


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

#Needed for theta

rMVNormC <- function(n, mu, U){
  p <- length(mu)
  Z <- matrix(rnorm(p*n), p, n)
  # U <- chol(Omega) # By default R's chol fxn returns upper cholesky factor
  X <- backsolve(U, Z) # more efficient and stable than actually inverting
  X <- sweep(X, 1, mu, FUN=`+`)
  return(X)
}
theta <- sim_theta(S, lambda)
#theta <- sim_theta(S, lambda_new)
# Theta is alpha (K of these), beta_1 (K of these), beta_2 (K of these), and sigma (J of these)
# lambda is mu_alpha, chol(sigma_alpha), mu_beta_1, chol(sigma_beta_1), mu_beta_2,
# chol(sigma_beta_2), c, d


# h -----------------------------------------------------------------------
#replace with hcpp
#hcpp(K, n_isotopes, n_covariates, 0.001, x_scaled, as.matrix(q), as.matrix(mu_s), 
# as.matrix(mu_c), 
# as.matrix(sigma_c), as.matrix(sigma_s), 
# theta[1,], as.matrix(y))

# log_q -------------------------------------------------------------------
#replace with
#log_q_cpp(theta[1,], lambda_new, K, n_isotopes, 100, n_covariates)

# Algorithm ---------------------------------------------------------------


# log_q -------------------------------------------------------------------


# Function to estimate different between joint and variational approx


# Run the VB function
run_VB_one_iter_cpp_r_hybrid <- function(lambda, # Starting value of lambda
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
  set.seed(t^2)
 theta <- sim_theta(S, lambda)
  c <- control_var_cpp(lambda, theta, K, n_isotopes,
                       0.001, n_covariates, x_scaled, concentrationmeans,
                       sourcemeans, correctionmeans, corrsds,sourcesds,y)
  g_0 <- nabla_LB_cpp(lambda,  theta, 
                      K, n_isotopes, beta_prior,
                      S,  n_covariates,
                      x_scaled,
                      concentrationmeans,  sourcemeans,
                      correctionmeans,
                      corrsds,  sourcesds,  y,
                      rep(0, length(lambda)))
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
  set.seed(t^2 + 7)
  theta <- sim_theta(S, lambda)
  
  # Compute g_t
  g_t <- nabla_LB_cpp(lambda,  theta, 
                      K, n_isotopes, beta_prior,
                      S,  n_covariates,
                      x_scaled,
                      concentrationmeans,  sourcemeans,
                      correctionmeans,
                      corrsds,  sourcesds,  y,
                      c)
  print("g_t")
  print(g_t)
  
  # Compute new control variate
  c <- control_var_cpp(lambda, theta, K, n_isotopes,
                       0.001, n_covariates, x_scaled, concentrationmeans,
                       sourcemeans, correctionmeans, corrsds,sourcesds,y)
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
 # if (t <= t_W) {
    # Compute a new lower bound estimate
    LB[t] <- LB_lambda_cpp( theta,  lambda, 
                            hfn(theta, n_sources, n, n_covariates, x_scaled),  n_sources,  n_isotopes, 
                            beta_prior,
                            n_covariates,
                            x_scaled,
                            concentrationmeans,  sourcemeans,
                            correctionmeans,
                            corrsds,  sourcesds,  y)
  # } else {
  #   LB[1:(t_W - 1)] <- LB[2:t_W]
  #   LB[t_W] <- LB_lambda_cpp( theta,  lambda, 
  #                             hfn(theta, n_sources, n, n_covariates, x_scaled),  n_sources,  n_isotopes, 
  #                             beta_prior,
  #                             n_covariates,
  #                             x_scaled,
  #                             concentrationmeans,  sourcemeans,
  #                             correctionmeans,
  #                             corrsds,  sourcesds,  y)
    
    # LB_bar <- mean(LB)
    # max_LB_bar <- max(max_LB_bar, LB_bar)
    # if (LB_bar >= max_LB_bar) {
    #   patience <- 0
    # } else {
    #   patience <- patience + 1
    # }
 # }
  
  if (patience > P) {
    print("Completed!")
    stop <- TRUE
  }
  t <- t + 1
  #}
  return(lambda)
}


lambda_out = lambda

n_iter = 200

lambda_out_matrix<-matrix(NA, nrow = (n_iter +1), ncol = length(lambda))
lambda_out_matrix[1,] <- lambda

for(i in 1:n_iter){
  lambda_out_matrix[i+1,] <- run_VB_one_iter_cpp_r_hybrid(lambda_out_matrix[i,], t = i)
  
}

lambda_out = run_VB_one_iter_cpp_r_hybrid(lambda_out, t = 1)




# Run the VB function
run_VB_cpp_r_hybrid <- function(lambda, # Starting value of lambda
                                         S = 100, # Number of samples to take
                                         P = 10, # Maximum patience before you stop
                                         beta_1 = 0.9, # Learning rates
                                         beta_2 = 0.9, # Learning rates
                                         tau = 100, # Iteration at which learning rate starts to decrease
                                         eps_0 = 0.1, # Raw learning rate multiplier
                                         t_W = 50 # Time window for working out convergence
) {
  
  # Starting
  set.seed(12345)
  theta <- sim_theta(S, lambda)
  c <- control_var_cpp(lambda, theta, K, n_isotopes,
                       0.001, n_covariates, x_scaled, concentrationmeans,
                       sourcemeans, correctionmeans, corrsds,sourcesds,y)
  g_0 <- nabla_LB_cpp(lambda,  theta, 
                      K, n_isotopes, beta_prior,
                      S,  n_covariates,
                      x_scaled,
                      concentrationmeans,  sourcemeans,
                      correctionmeans,
                      corrsds,  sourcesds,  y,
                      rep(0, length(lambda)))
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
  print("t")
  print(t)
  # Generate new samples
  set.seed(12345)
  theta <- sim_theta(S, lambda)
  
  # Compute g_t
  g_t <- nabla_LB_cpp(lambda,  theta, 
                      K, n_isotopes, beta_prior,
                      S,  n_covariates,
                      x_scaled,
                      concentrationmeans,  sourcemeans,
                      correctionmeans,
                      corrsds,  sourcesds,  y,
                      c)
  print("g_t")
  print(g_t)
  
  # Compute new control variate
  c <- control_var_cpp(lambda, theta, K, n_isotopes,
                       0.001, n_covariates, x_scaled, concentrationmeans,
                       sourcemeans, correctionmeans, corrsds,sourcesds,y)
  
  # c <- control_var(lambda, theta)
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
    LB[t] <- LB_lambda_cpp( theta,  lambda, 
                            hfn(theta, n_sources, n, n_covariates, x_scaled),  n_sources,  n_isotopes, 
                            beta_prior,
                            n_covariates,
                            x_scaled,
                            concentrationmeans,  sourcemeans,
                            correctionmeans,
                            corrsds,  sourcesds,  y)
  } else {
    LB[1:(t_W - 1)] <- LB[2:t_W]
    LB[t_W] <- LB_lambda_cpp(theta,  lambda, 
                              hfn(theta, n_sources, n, n_covariates, x_scaled),  n_sources,  n_isotopes, 
                              beta_prior,
                              n_covariates,
                              x_scaled,
                              concentrationmeans,  sourcemeans,
                              correctionmeans,
                              corrsds,  sourcesds,  y)
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





lambdaout <- run_VB_cpp_r_hybrid(lambda)


n_samples = 3600

theta_out_rcpp <- sim_theta(n_samples, lambdaout)

beta_rcpp<-matrix(colMeans(theta_out_rcpp[,1:(K*(n_covariates))]), nrow = (n_covariates))
sigma_rcpp<- colMeans(theta_out_rcpp[,(n_sources*(n_covariates)+1):(n_sources*(n_covariates)+n_isotopes)])

f1_rcpp <- matrix(NA, ncol = K, nrow = n)

for(k in 1:K){
  f1_rcpp[,k] = x_scaled %*% beta_rcpp[,k]
}
p1_rcpp <- matrix(NA, ncol = K, nrow = n)

for (i in 1:n) {
  p1_rcpp[i, ] <- exp(f1_rcpp[i, 1:K]) / (sum((exp(f1_rcpp[i, 1:K]))))
}

print(p1_rcpp)





