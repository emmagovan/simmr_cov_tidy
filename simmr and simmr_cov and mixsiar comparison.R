#Compare simmr and simmr_cov (when the only cov is an intercept)
library(readxl)
library(tidyverse)
library(checkmate)

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
y <- consumer |>
  select(d13C_Pl, d15N_Pl) |>
  as.matrix()

library(simmr)
test <- simmr_load(as.matrix(y), sources$Sources, as.matrix(mu_s), as.matrix(sigma_s), as.matrix(mu_c),
                   as.matrix(sigma_c), as.matrix(q))
plot(test)


test2<-simmr_load_cov(as.matrix(y), ~c(rep(1,9)), sources$Sources, as.matrix(mu_s), as.matrix(sigma_s),
                      intercept = FALSE,
                      as.matrix(mu_c), as.matrix(sigma_c), as.matrix(q))


original_out = simmr_ffvb(test)


cov_out = simmr_ffvb_cov(test2)



hist(cov_out$output$'1'$BUGSoutput$sims.list$p[1,,])
hist(cov_out$output$'1'$BUGSoutput$sims.list$p[2,,])
hist(cov_out$output$'1'$BUGSoutput$sims.list$p[3,,])
hist(cov_out$output$'1'$BUGSoutput$sims.list$p[4,,])

colnames(x_scaled) = "A"

colnames(y) = c("d13C", "d15N")
write.csv(cbind(y, (x_scaled)), "y.csv")



#Compare simmr_cov and mixsiar
library(MixSIAR)

mix.filename <- system.file("extdata", "geese_consumer.csv", package = "MixSIAR")

mix <- load_mix_data(filename=mix.filename,
                     iso_names=c("d13C","d15N"),
                     factors="Group",
                     fac_random=FALSE,
                     fac_nested=FALSE,
                     cont_effects=NULL)

source.filename <- system.file("extdata", "geese_sources.csv", package = "MixSIAR")

source <- load_source_data(filename=source.filename,
                           source_factors=NULL,
                           conc_dep=TRUE,
                           data_type="means",
                           mix)


discr.filename <- system.file("extdata", "geese_discrimination.csv", package = "MixSIAR")

discr <- load_discr_data(filename=discr.filename, mix)


plot_data(filename="isospace_plot", plot_save_pdf=TRUE, plot_save_png=FALSE, mix,source,discr)


plot_prior(alpha.prior=1,source)

model_filename <- "MixSIAR_model.txt"
resid_err <- TRUE
process_err <- FALSE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)




jags.1 <- run_model(run="test", mix, source, discr, model_filename,
                    alpha.prior = 1, resid_err, process_err)


test2<-simmr_load_cov(matrix(c(geese_data[[1]]$d13C_Pl, geese_data[[1]]$d15N_Pl), ncol = 2), 
                      ~c(rep(1,251)), sources$Sources, as.matrix(mu_s), as.matrix(sigma_s),
                      intercept = FALSE,
                      as.matrix(mu_c), as.matrix(sigma_c), as.matrix(q),
                      group = c(rep(1,20), rep(2,30), rep(3,30), rep(4, 15), rep(5, 87), rep(6, 30),
                      rep(7, 10), rep(8, 29)))

simmr_ffvb_cov(test2)

library(microbenchmark)
microbenchmark(simmr_ffvb_cov(test2),
               run_model(run="test", mix, source, discr, model_filename,
                         alpha.prior = 1, resid_err, process_err),
               times = 5L
               )

output_JAGS(jags.1, mix, source)




