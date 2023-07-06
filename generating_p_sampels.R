# Code for generating samples of p and putting those in a matrix
# so that it can be plotted like the other plots

#So what we get out is
#lambda = K(mu) + K*(K+1)/2(sigma) for each covariate, then c and d (one for each isotope)

# beta_lambda<-c(rep(0, K),rep(1, K * (K + 1) / 2))
# lambda <- c(
#   rep(beta_lambda, n_covariates),
#   rep(1, n_isotopes), #shape
#   rep(1, n_isotopes) #rate
# )

#Use sim_theta to generate beta for each covariate
# It outputs 3600 * (K*n_cov)
#Then p is 
# f = x_scaled %*% beta
# 
# p <- matrix(NA, ncol = K, nrow = n)
# 
# for (i in 1:n) {
#   p[i, ] <- exp(f[i,1:K]) / (sum((exp(f[i,1:K]))))
# }

#But we want instead to do from 1:3600
# So we would have
p_sample = array(NA, dim = c(K, 3600, n))

n_isotopes = 2
n_covariates = 2
S = 3600
K = 4
n = 9
x_scaled = matrix(c(rep(c(1,0,1), 9)), ncol = 3)

#theta_out_rcpp <- sim_thetacpp(S, lambda_outrcpp, n_sources, n_isotopes, n_covariates)

beta = array(simmr_out$output$'1'$theta[,1:(n_covariates * K)], dim = c(S, n_covariates, K))

f <- array(NA, dim = c(n,K, S)) 

for(s in 1:S){
  f[,,s] = x_scaled %*% beta[s,,]
}

p_sample = array(NA, dim = c(n, S, K))



for(s in 1:S){
for (i in 1:n) {
  p_sample[i,s, ] <- exp(f[i,1:K, s]) / (sum((exp(f[i,1:K, s]))))
}
}

hist(p_sample[1,,1])
hist(p_sample[1,,2])
hist(p_sample[1,,3])
hist(p_sample[1,,4])


#I think this works - and could set up a plot function ? Not sure if necessary?
#Could do basic histogram plot and then leave it at that?
#Or could try those pieglyph plots

#Would get colMeans for each p

p_means = matrix(NA, nrow = n, ncol = K)

for(i in 1:n){
  p_means[i,] = colMeans(p_sample[i,,])
}

df <- reshape2::melt(p_means)
colnames(df) <- c("Num", "Source", "Proportion")

#Now figure out pielyph package thingy
library(PieGlyph)
library(ggplot2)

y = data.frame("iso1" = c(rnorm(9, 5, 2)), "iso2" = c(rnorm(9, -12, 1)), "A" = c(p_means[,1]), "B" = c(p_means[,2]),
               "C" = c(p_means[,3]), "D" = c(p_means[,4]))

ggplot(y, aes(x = iso1, y = iso2)) + geom_pie_glyph(slices = c("A", "B", "C", "D"))


#Can probably(?) just pop this in to the simmr_input code - so the new simmr_output function
#Would prob just have histogram/density/isospace with this as the new option for isospace
#Not sure if we give the option of a histogram? We'd end up with n * K plots - or else we make them
#specify which individual they want? So then it would just be K plots?






