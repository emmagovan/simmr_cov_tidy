#Lake script
#compare VB and JAGS for Whale and Lake Examples

#source data
lake_consumer <- read.csv("lake_consumer.csv")
lake_disc <- read.csv("lake_discrimination.csv")
lake_sources <- read.csv("lake_sources.csv")


# Put them into the right names
n <- nrow(lake_consumer)
n_isotopes <- 2
K <- (lake_sources$Source) %>% unique() %>% length()
mu_s <-  matrix(c(lake_sources %>% group_by(Source) %>% 
                    summarise(Md13C = mean(d13C), 
                              Md15N = mean(d15N))))
sigma_s <- lake_sources[, c(3, 5)]
mu_c <- lake_disc[, c(2, 4)]
sigma_c <- lake_disc[, c(3, 5)]
q <- matrix(c(rep(1/(n_isotopes*K), n_isotopes *K)), nrow = K)#conc[, c(2:3)]
n_sources = K

simmr_groups = simmr_load(mixtures=as.matrix(y),
                          source_names=unlist(whale_sources[,1]),
                          source_means=as.matrix(whale_sources[,c(2,4)]),
                          source_sds=as.matrix(whale_sources[,c(3,5)]),
                          correction_means=as.matrix(whale_disc[,c(2,4)]),
                          correction_sds=as.matrix(whale_disc[,c(3,5)]),
                          concentration_means = q)

library(simmr)
plot(simmr_groups)

#-----------Add whatever covariates here-----------------
x_scaled <- matrix(c(rep(1, n)), ncol = 1)

#x <- matrix(c(consumer$Skull, consumer$Age, consumer$Sex, consumer$`Net Wt`),
#            ncol = 4)
#x_scaled <- cbind(matrix(rep(1, nrow(x)), ncol = 1), scale(x))


n_covariates <- (ncol(x_scaled))

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

#-------RUN VB-------------------------------------------------------

lambda_outrcpp<- run_VB_cpp(lambda, K, n_isotopes, n_covariates, n, 0.001, 
                            as.matrix(q), as.matrix(mu_s), as.matrix(mu_c), 
                            as.matrix(sigma_c), as.matrix(sigma_s), y, x_scaled,
                            100, 10, 0.9, 0.9, 1000, 0.1, 50)



#--------------------JAGS----------------------
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
      sigma[j] ~ dgamma(0.001, 0.001)
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
                  K=K,
                  ncov = n_covariates,
                  J=n_isotopes,
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

#-----------COMPARISON-----------------------


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
print(colMeans(model_run$BUGSoutput$sims.list$p))









