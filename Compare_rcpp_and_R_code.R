#general number of covariates
#Specify the number at the start and then automatically generate it??

# Load in package
library(simmr)
library(readxl)
library(tidyverse)
library(R2jags)

# Source in the generic functions
source("FF_VB_generic_functions_correct.R")
Rcpp::sourceCpp("run_VB.cpp")


# Extract data ------------------------------------------------------------

path <- system.file("extdata", "geese_data.xls", package = "simmr")
geese_data <- lapply(excel_sheets(path), read_excel, path = path)

# Just use time point 1 for now
consumer <- geese_data[[1]] |> filter(Time == 1)
sources <- geese_data[[2]]
TEFs <- geese_data[[3]]
conc <- geese_data[[4]]

# Put them into the right names
n <- nrow(consumer)
n_isotopes <- 2
K <- nrow(sources)
mu_s <- sources[, c(2, 3)] 
sigma_s <- sources[, c(4, 5)]
mu_c <- TEFs[, c(2, 3)]
sigma_c <- TEFs[, c(4, 5)]
q <- conc[, c(2:3)]

#1 iso
# n <- nrow(consumer)
# n_isotopes <- 1
# K <- nrow(sources)
# mu_s <- sources[, c(2)] 
# sigma_s <- sources[, c(4)]
# mu_c <- TEFs[, c(2)]
# sigma_c <- TEFs[, c(4)]
# q <- conc[, c(2)]

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


#--------------------JAGS--------------

model_code <- "
model{
  for (i in 1:N) {
    for (j in 1:J) { 
      y[i,j] ~ dnorm(inprod(p[i,]*q[,j], s_mean[,j]+c_mean[,j]) / inprod(p[i,],q[,j]), 1/var_y[i,j])
      var_y[i,j] = inprod(pow(p[i,]*q[,j],2),s_sd[,j]^2+c_sd[,j]^2)/pow(inprod(p[i,],q[,j]),2) 
        + pow(sigma[j],2)
    }
  }
  for(i in 1:N) {
    p[i,1:K] = expf[i,]/sum(expf[i,])
    for(k in 1:K) {
      expf[i,k] = exp(f[i,k])
      f[i,k] = mu_f[i,k]
    }
  }
  for(k in 1:K) {
  for(i in 1:N) {
      mu_f[i,k] = inprod(x[i,], beta[,k])
  }
  
  
  for(a in 1:(ncov)){
   beta[a,k] ~ dnorm(0,1)
  }
  }   
  
  for(j in 1:J) { 
      sigma[j] ~ dgamma(0.1, 0.1)
      }
}
"

model_data = list(y=y,
                  s_mean=as.matrix(mu_s),
                  s_sd=as.matrix(sigma_s),
                  c_mean=as.matrix(mu_c),
                  c_sd=as.matrix(sigma_c),
                  q=as.matrix(q), 
                  N=nrow(y),
                  K=4,
                  ncov = n_covariates,
                  J=2,
                  x = x_scaled)


model_parameters <- c("beta", "sigma", "p")

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code),
  n.chains = 4, # Number of different starting positions
  n.iter = 10000, # Number of iterations
  n.burnin = 2000, # Number of iterations to remove at start
  n.thin = 5
) # Amount of thinning)


print(model_run)


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
c_0 <- c(rep(0.001, n_isotopes)) #Change to 0.0001
d_0 <- c(rep(0.001, n_isotopes))
beta_lambda<-c(rep(0, K),rep(1, K * (K + 1) / 2))
lambda <- c(
  rep(beta_lambda, n_covariates),
  rep(1, n_isotopes), #shape
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
theta <- sim_theta(S, lambda)
#theta <- sim_theta(S, lambda_new)
theta_cpp <- sim_thetacpp(S, lambda, K, n_isotopes, n_covariates)
# Theta is alpha (K of these), beta_1 (K of these), beta_2 (K of these), and sigma (J of these)
# lambda is mu_alpha, chol(sigma_alpha), mu_beta_1, chol(sigma_beta_1), mu_beta_2,
# chol(sigma_beta_2), c, d


# h -----------------------------------------------------------------------

# Log of likelihood added to prior
h <- function(theta) {
  # Create betas and sigma
  beta <- matrix(theta[1:((n_covariates) * K)], nrow = (n_covariates), byrow = TRUE)
  sigma <- theta[((n_covariates) * K +1):(((n_covariates) * K)+n_isotopes)]
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
                                           sum(p[i, ]^2 * q[, j]^2) + sigma[j]^2),
                               log = TRUE

                      ))
      mumat[i,j] = sqrt(sum(p[i, ]^2 * q[, j]^2 * (sigma_s[, j]^2 + sigma_c[, j]^2)) /
                          sum(p[i, ]^2 * q[, j]^2) + sigma[j]^2)


    }
  }
  
  beta_sum <- 0
  #change to k
  for (i in 1:(n_covariates)){
    beta_sum = beta_sum +sum(dnorm(beta[i,], mu_beta_zero[i,], sigma_beta_zero[i,], log = TRUE))
  }
  
  return(hold + beta_sum + sum(dgamma(sigma, shape = c_0, rate = d_0, log = TRUE)))
}
h(theta[1,])
hcpp(K, n_isotopes, n_covariates, 0.001, x_scaled, as.matrix(q), as.matrix(mu_s), 
     as.matrix(mu_c), 
     as.matrix(sigma_c), as.matrix(sigma_s), 
     theta[1,], as.matrix(y))
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
  
  shape_sigma <- lambda[lambda_index$c]
  rate_sigma <- lambda[lambda_index$d]
  
  # Extract alpha, beta and sigma from theta
  beta <- matrix(theta[1:((n_covariates) * K)], nrow = (n_covariates), ncol = K,  byrow = TRUE)
  sigma <- theta[((n_covariates) * K +1):(((n_covariates) * K)+n_isotopes)]
  
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
  
    return(sum_p + sum(dgamma(sigma,
                               shape = shape_sigma,
                               rate = rate_sigma,
                               log = TRUE
 ))
 )
}
log_q(lambda, theta[1,])
log_q_cpp(theta[1,], lambda, K, n_isotopes, 100, n_covariates)

lambda_new <- c(1:32)

log_q(lambda_new, theta[1,])
log_q_cpp(theta[1,], lambda_new, K, n_isotopes, 100, n_covariates)

# Algorithm ---------------------------------------------------------------

lambda_out <- run_VB(lambda)
lambda_outrcpp <- run_VB_cpp(lambda, K, n_isotopes, n_covariates, n, 0.001, 
                             as.matrix(q), as.matrix(mu_s), as.matrix(mu_c), 
                             as.matrix(sigma_c), as.matrix(sigma_s), y, x_scaled,
                             100, 10, 0.9, 0.9, 1000, 0.1, 50)


library(simmr)
simmr_in<-simmr_load(as.matrix(y), sources$Sources, as.matrix(mu_s), as.matrix(sigma_s), as.matrix(mu_c),
                     as.matrix(sigma_c), as.matrix(q))

simmr_out <- simmr_ffvb(simmr_in)



library(microbenchmark)
microbenchmark(run_VB_cpp(lambda, K, n_isotopes, n_covariates, n, 0.001, 
                          as.matrix(q), as.matrix(mu_s), as.matrix(mu_c), 
                          as.matrix(sigma_c), as.matrix(sigma_s), y, x_scaled,
                          100, 10, 0.9, 0.9, 1000, 0.1, 50), 
               model_run <- jags(
                 data = model_data,
                 parameters.to.save = model_parameters,
                 model.file = textConnection(model_code),
                 n.chains = 4, # Number of different starting positions
                 n.iter = 10000, # Number of iterations
                 n.burnin = 2000, # Number of iterations to remove at start
                 n.thin = 5
               ), times = 10L)

n_samples <- 3600

# Check results ---------------------------------
theta_out <- sim_theta(n_samples, lambda_out)
theta_out_rcpp <- sim_thetacpp(S, lambda_outrcpp, K, n_isotopes, n_covariates)

#Easy way
beta<-matrix(colMeans(theta_out[,1:(K*(n_covariates))]), nrow = (n_covariates))
sigma <- colMeans(theta_out[,(K*(n_covariates)+1):(K*(n_covariates)+n_isotopes)])

beta_rcpp<-matrix(colMeans(theta_out_rcpp[,1:(K*(n_covariates))]), nrow = (n_covariates))
sigma_rcpp<- colMeans(theta_out_rcpp[,(n_sources*(n_covariates)+1):(n_sources*(n_covariates)+n_isotopes)])

f1 <- matrix(NA, ncol = K, nrow = n)

for(k in 1:K){
  f1[,k] = x_scaled %*% beta[,k]
}

f1_rcpp <- matrix(NA, ncol = K, nrow = n)

for(k in 1:K){
  f1_rcpp[,k] = x_scaled %*% beta_rcpp[,k]
}



p1 <- matrix(NA, ncol = K, nrow = n)

for (i in 1:n) {
  p1[i, ] <- exp(f1[i, 1:K]) / (sum((exp(f1[i, 1:K]))))
}

p1_rcpp <- matrix(NA, ncol = K, nrow = n)

for (i in 1:n) {
  p1_rcpp[i, ] <- exp(f1_rcpp[i, 1:K]) / (sum((exp(f1_rcpp[i, 1:K]))))
}

