#Comparing all the generic functions

#Get data --------------
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

lambda2 <- c(1:32)
n_sources <- K
n_tracers = n_isotopes

#-------------------
#This works
d_rcpp <- delta_lqltcpp(lambda, theta[1,], 0.01, n_sources, n_tracers, n_covariates, S)

d_r <- delta_lqlt(lambda, theta[1,], 0.01)

#--------h lambda
#ALso agrees
concentrationmeans = as.matrix(q)
sourcemeans = as.matrix(mu_s)
correctionmeans = as.matrix(mu_c)
corrsds = as.matrix(sigma_c)
sourcesds = as.matrix(sigma_s)
beta_prior = 0.001


h_l_cpp <- h_lambdacpp(n_sources, n_isotopes, 0.001, n_covariates, S, 
                     concentrationmeans, sourcemeans, correctionmeans, 
                     corrsds, sourcesds, 
                     theta[1,], y, lambda, x_scaled) 

h_l <- h_lambda(lambda, theta[1,], y)
  
  
  
#---------------

 
nabla_LB_rcpp = nabla_LB_cpp(lambda,  theta, 
                            n_sources, n_tracers, beta_prior,
                            S,  n_covariates,
                            x_scaled,
                            concentrationmeans,  sourcemeans,
                            correctionmeans,
                            corrsds,  sourcesds,  y,
                            rep(0, length(lambda)))


nabla_LB_r <- nabla_LB(lambda, theta, c = rep(0, length(lambda)))

#These are not identical - biggest difference between any entry is 0.0007473353
#Not sure if this is causing any issue
#Probably rounding errors? Maybe

#------------------
#control variate
c_rcpp <- control_var_cpp(lambda, theta, n_sources, n_tracers,
                          beta_prior,n_covariates, x_scaled, concentrationmeans,
                          sourcemeans, correctionmeans, corrsds,sourcesds,y)

  
c_r = control_var(lambda, theta)  
  
  #Same thing here - biggest difference is  0.0001922434


#-------------------
LBlambda<- LB_lambda_cpp( theta,  lambda, 
                          hfn(theta, n_sources, n, n_covariates, x_scaled),  n_sources,  n_isotopes, 
beta_prior,
 n_covariates,
  x_scaled,
  concentrationmeans,  sourcemeans,
  correctionmeans,
  corrsds,  sourcesds,  y)


LB_lambdar <- LB_lambda(lambda, theta) 


#Works

#-----------------
#okay these all agree
#Changing the eps seemed to do something
#But i have no idea if the p values are right or not
#Which is a problem.

lambda <- c(0.955376849,-1.555770607,-0.686713113,-0.806404492, 
            2.232814868,1.187747933,1.542799752,1.802486717,
            0.173387376,1.698359957,0.623842862,0.013428788,
            1.878521064,2.382494327,2.417136436,0.006784888,
            0.373499402,2.155241805)

lambda_outrcpp <- run_VB_cpp(lambda, K, n_isotopes, n_covariates, n, 0.001, 
                             as.matrix(q), as.matrix(mu_s), as.matrix(mu_c), 
                             as.matrix(sigma_c), as.matrix(sigma_s), y, x_scaled,
                             100, 10, 0.9, 0.9, 1000, 0.1, 50)







  