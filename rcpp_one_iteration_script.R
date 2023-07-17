#Rcpp code running one iteration at a time and printing everything
library(simmr)
library(readxl)
library(tidyverse)
library(R2jags)
Rcpp::sourceCpp("run_VB_one_iter.cpp")

lambda<-c(1.60233, -1.31974, 0.20779, -0.299813,
          2.06332, -0.556013, 1.94068, 2.08669,
          -0.306925, 1.91809, 2.25527, 1.89869,
          2.16079, 2.26828, 0.862199, 1.82158,
          2.28694, 0.00612707)

lambda <- c(1.32008, -0.609668, 1.0523, 1.16631,
            1.83509, -0.145871, 1.5787, 0.624907,
            0.0455274, 1.94142, 2.1015, 1.44756,
            0.0118001, 1.91335, -0.551781, 0.178901,
            -0.70396, -0.761399, 1.81927, -0.0731803,
            1.96346, 1.9464, 0.813659, 2.02584,
            0.570628, 0.101592, 2.03861, 2.06261,
            0.994666, 0.770858, 1.02269, 0.87536,
            1.91784, 0.00197615, 1.93863, -0.054107,
            -0.00383162, 1.95647, 2.2411, 1.95962,
            0.639015, 1.99183, 1.09652, 1.04422,
            -0.287044, 1.2924, 1.92963, 1.92513, 1.89454,
            -0.146175, -0.0316691, 1.94749, 2.10041,
            1.71755, 1.92093, 2.02576, 0.94695,
            -0.937233, -0.962299, -1.04214, 1.97058,
            1.31386, 1.98335, 1.86831, 1.31255,
            1.91689, 0.0486656, 0.801795, 0.231705,
            1.98232, 0.300367, 1.82334, 2.24998, 0.000279973)

lambda <- c(0.82089, -1.44535, 0.394004, 0.380435,
             2.12727, 2.13987, 1.71719, 0.319057,
             -0.392343, 1.84948, 1.20698, 1.71934,
             1.15237, 2.20957, 0.663171, 2.20121,
             2.29304, 0.00184574)

lambda<- c(1.76074, -0.620311, -0.333896, 0.58943, 2.6487, -1.37606, 2.68549,
           1.88536, -1.19752, 2.25904, 2.88603, -0.736874, 2.46533, 2.83938,
           2.06689, 2.53046, 1.70085, 0.0123056)


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

#-----------Add whatever covariates here-----------------
x_scaled <- matrix(c(rep(1, 9)), ncol = 1)

x <- matrix(c(consumer$Skull, consumer$Age, consumer$Sex, consumer$`Net Wt`),
 ncol = 4)
 x_scaled <- cbind(matrix(rep(1, nrow(x)), ncol = 1), scale(x))

 
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

#Need to initialise this somewhere
#if using loop
lambda_outrcpp <- c(
  rep(beta_lambda, n_covariates),
  rep(1, n_isotopes), #shape
  rep(1, n_isotopes) #rate
)
#For just one run
lambda <- c(
  rep(beta_lambda, n_covariates),
  rep(1, n_isotopes), #shape
  rep(1, n_isotopes) #rate
)

#--------------------------------------------------------------

e<- function(lambda_outrcpp, t){
  lambda_outrcpp <- run_VB_cpp_one_iter(lambda_outrcpp, K, n_isotopes, n_covariates, n, 0.001, 
                                        as.matrix(q), as.matrix(mu_s), as.matrix(mu_c), 
                                        as.matrix(sigma_c), as.matrix(sigma_s), y, x_scaled,
                                        100, 10, 0.9, 0.9, 1000, 0.1, 50, t)
  # newList <- list("lambda" = lambda_outrcpp, 
  #                 "theta" = sim_thetacpp(1, lambda_outrcpp, n_sources, n_isotopes, n_covariates))
  print(lambda_outrcpp)
}


for(i in 1:40){
  e(lambda_outrcpp, 1)
  
}

#Alternatively
Rcpp::sourceCpp("run_VB.cpp")

#This just prints everything as it runs - but doesn't run into the issues of doing 1 iteration at a time
#e.g. no stopping rule

#Set number of iterations
n_iter = 100
lambda_outrcpp_matrix<-matrix(NA, nrow = n_iter, ncol = length(lambda))

#Run once
lambda_outrcpp<- run_VB_cpp(lambda, K, n_isotopes, n_covariates, n, 0.001, 
                                        as.matrix(q), as.matrix(mu_s), as.matrix(mu_c), 
                                        as.matrix(sigma_c), as.matrix(sigma_s), y, x_scaled,
                                        100, 10, 0.9, 0.9, 100, 0.1, 50)

#Altering tau (the threshold after which the learning rate is decreased) seems to affect the result
#too low and it doesn't converge to the right answer


#Run multiple times
for(i in 1:n_iter){
lambda_outrcpp_matrix[i,] <- run_VB_cpp(lambda, K, n_isotopes, n_covariates, n, 0.001, 
                             as.matrix(q), as.matrix(mu_s), as.matrix(mu_c), 
                             as.matrix(sigma_c), as.matrix(sigma_s), y, x_scaled,
                             100, 10, 0.9, 0.9, 10, 0.1, 50)

}

lambda_outrcpp <- lambda_outrcpp_matrix[1,]


#just run all this to check output - see if its giving similar answers to JAGS

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
#For this example these should be approx 0.6, 0.01, 0.2, 0.2