#Rcpp code running one iteration at a time and printing everything
library(simmr)
library(readxl)
library(tidyverse)
library(R2jags)

#Source in functions
Rcpp::sourceCpp("run_VB_one_iter.cpp")


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
n_sources = K

concentrationmeans = as.matrix(q)
sourcemeans = as.matrix(mu_s)
correctionmeans = as.matrix(mu_c)
corrsds = as.matrix(sigma_c)
sourcesds = as.matrix(sigma_s)
beta_prior = 0.001

#-----------Add whatever covariates here-----------------
x_scaled <- matrix(c(rep(1, 9)), ncol = 1)

# x <- matrix(c(consumer$Skull, consumer$Age, consumer$Sex, consumer$`Net Wt`),
#  ncol = 4)
#  x_scaled <- cbind(matrix(rep(1, nrow(x)), ncol = 1), scale(x))

 
n_covariates <- (ncol(x_scaled))
y <- consumer |>
  select(d13C_Pl, d15N_Pl) |>
  as.matrix()



##---------------------------
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


#--------SET LAMBDA--------------------
lambda <- c(
  rep(beta_lambda, n_covariates),
  rep(1, n_isotopes), #shape
  rep(1, n_isotopes) #rate
)


#lambda<-c(1.30126, -0.975428, -0.198901, -0.179772,
#4.42257, 0.575808, 4.34706, 2.10692, 0.686071, 4.42486,
#0.831033, 0.88614, 0.973301, 4.4264,
#1.68114, 2.15085, 0.171745, 0.420329)

#--------------------------------------------------------------
set.seed(12333)
theta <- sim_thetacpp(S, lambda, K, n_isotopes, n_covariates)

#Alt for identical thetas: 
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


# source("fns_for_R_one_iter.R")
# set.seed(12222)
# theta <- sim_theta(S, lambda)



hcpp(K, n_isotopes, n_covariates, 0.001, x_scaled, concentrationmeans, sourcemeans, 
     correctionmeans, 
     corrsds, sourcesds, 
     theta[1,], as.matrix(y))

log_q_cpp(theta[1,], lambda, K, n_isotopes, 100, n_covariates)

d_rcpp <- delta_lqltcpp(lambda, theta[1,], 0.001, n_sources, n_isotopes, n_covariates, S)
print(d_rcpp)

h_l_cpp <- h_lambdacpp(n_sources, n_isotopes, 0.001, n_covariates, S, 
                       concentrationmeans, sourcemeans, correctionmeans, 
                       corrsds, sourcesds, 
                       theta[1,], y, lambda, x_scaled) 
print(h_l_cpp)

nabla_LB_rcpp = nabla_LB_cpp(lambda,  theta, 
                             n_sources, n_isotopes, beta_prior,
                             S,  n_covariates,
                             x_scaled,
                             concentrationmeans,  sourcemeans,
                             correctionmeans,
                             corrsds,  sourcesds,  y,
                             rep(0, length(lambda)))

print(nabla_LB_rcpp)

big_delta_lqlt <- t(apply(theta, 1, delta_lqltcpp, lambda = lambda, eps = 0.001, n_sources = n_sources,
                          n_tracers = n_isotopes, n_covariates = n_covariates , S = S))
print(big_delta_lqlt)



c_rcpp <- control_var_cpp(lambda, theta, K, n_isotopes,
                          0.001, n_covariates, x_scaled, concentrationmeans,
                          sourcemeans, correctionmeans, corrsds,sourcesds,y)
print(c_rcpp)

LBlambda<- LB_lambda_cpp( theta,  lambda, 
                          hfn(theta, n_sources, n, n_covariates, x_scaled),  n_sources,  n_isotopes, 
                          beta_prior,
                          n_covariates,
                          x_scaled,
                          concentrationmeans,  sourcemeans,
                          correctionmeans,
                          corrsds,  sourcesds,  y)
print(LBlambda)


lambda_outrcpp <- run_VB_cpp_one_iter(lambda, K, n_isotopes, n_covariates, n, 0.001, 
                             as.matrix(q), as.matrix(mu_s), as.matrix(mu_c), 
                             as.matrix(sigma_c), as.matrix(sigma_s), y, x_scaled,
                             100, 10, 0.9, 0.9, 100, 0.1, 50, 1)



#Set number of iterations
n_iter = 100
lambda_outrcpp_matrix<-matrix(NA, nrow = (n_iter+1), ncol = length(lambda))
lambda_outrcpp_matrix[1,] <- lambda


#Run multiple times
for(i in 1:n_iter){
  lambda_outrcpp_matrix[i+1,] <- run_VB_cpp_one_iter(lambda_outrcpp_matrix[i,], K, n_isotopes, n_covariates, n, 0.001, 
                                          as.matrix(q), as.matrix(mu_s), as.matrix(mu_c), 
                                          as.matrix(sigma_c), as.matrix(sigma_s), y, x_scaled,
                                          100, 10, 0.9, 0.9, 100, 0.1, 50, t = i)
  
}

lambda_outrcpp <- lambda_outrcpp_matrix[1,]



theta_out_rcpp <- sim_thetacpp(S, lambda_outrcpp, K, n_isotopes, n_covariates)
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
print(model_run)

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




print(p1_rcpp)
print(p1_r)




#For this example these should be approx 0.6, 0.01, 0.2, 0.2