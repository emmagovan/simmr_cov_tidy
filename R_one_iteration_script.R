#R code for running one iteration at a time and printing out everything

#general number of covariates
#Specify the number at the start and then automatically generate it??

# Load in package
library(simmr)
library(readxl)
library(tidyverse)
library(R2jags)

# Source in the generic functions
source("fns_for_R_one_iter.R")




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

########## SET THESE-------------------------------------
#   x <- matrix(c(consumer$Skull), 
# ncol = 1)
# x_scaled <- cbind(matrix(rep(1, nrow(x)), ncol = 1), scale(x))

x_scaled <- matrix(c(rep(1, 9)), ncol = 1)

n_covariates <- (ncol(x_scaled))
y <- consumer |>
  select(d13C_Pl, d15N_Pl) |>
  as.matrix()



#One isotope

# y <- consumer |>
#   select(d13C_Pl) |>
#   as.matrix()



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

#---------SET LAMBDA----------------

lambda <- c(
  rep(beta_lambda, n_covariates),
  rep(1, n_isotopes), #shape
  rep(1, n_isotopes) #rate
)

# lambda<-c(1.30126, -0.975428, -0.198901, -0.179772,
#           4.42257, 0.575808, 4.34706, 2.10692, 0.686071, 4.42486,
#           0.831033, 0.88614, 0.973301, 4.4264,
#           1.68114, 2.15085, 0.171745, 0.420329)


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

# set.seed(12222)
theta <- sim_theta(S, lambda)


#Alt

# Rcpp::sourceCpp("run_VB_one_iter.cpp")
# set.seed(12333)
# theta <- sim_thetacpp(S, lambda, K, n_isotopes, n_covariates)

h(theta[1,])

log_q(lambda, theta[1,])

d_r <- delta_lqlt(lambda, theta[1,], 0.001)
print(d_r)

h_l <- h_lambda(lambda, theta[1,], y)
print(h_l)

nabla_LB_r <- nabla_LB(lambda, theta, c = rep(0, length(lambda)))
print(nabla_LB_r)


big_delta_lqlt <- t(apply(theta, 1, delta_lqlt, lambda = lambda))
print(big_delta_lqlt)


c_r = control_var(lambda, theta) 
print(c_r)

LB_lambdar <- LB_lambda(lambda, theta)
print(LB_lambdar)

lambda_out <- run_VB_one_iter(lambda, t = 1)
print(lambda_out)

n_iter = 10

lambda_out_matrix<-matrix(NA, nrow = (n_iter +1), ncol = length(lambda))
lambda_out_matrix[1,] <- lambda

for(i in 1:n_iter){
  lambda_out_matrix[i+1,] <- run_VB_one_iter(lambda_out_matrix[i,], t = i)
  
}
print(lambda_out_matrix[2,])



theta_out <- sim_theta(S, lambda_out, K, n_isotopes, n_covariates)
beta_out<-matrix(colMeans(theta_out[,1:(K*(n_covariates))]), nrow = (n_covariates))
sigma_out<- colMeans(theta_out[,(n_sources*(n_covariates)+1):(n_sources*(n_covariates)+n_isotopes)])
f1_r <- matrix(NA, ncol = K, nrow = n)
for(k in 1:K){
  f1_r[,k] = x_scaled %*% beta_out[,k]
}
p1_r <- matrix(NA, ncol = K, nrow = n)
for (i in 1:n) {
  p1_r[i, ] <- exp(f1_r[i, 1:K]) / (sum((exp(f1_r[i, 1:K]))))
}


print(p1_r)

