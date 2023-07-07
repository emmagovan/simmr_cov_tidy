#Test simmr_cov functions
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

test2<-simmr_load_cov(y, ~consumer$Sex, sources$Sources, mu_s, sigma_s, intercept = TRUE)

Rcpp::source("run_VB.cpp")
#simmr_out = simmr_ffvb_cov(simmr_1)



simmr_1 <- with(
  geese_data_day1,
  simmr_load_cov(
    mixtures = mixtures,
    ~matrix(data = c(rgamma(9,2,1)), ncol = 1),
    source_names = source_names,
    source_means = source_means,
    source_sds = source_sds,
    correction_means = correction_means,
    correction_sds = correction_sds,
    concentration_means = concentration_means,
    intercept = TRUE
  )
)

simmr_out = simmr_ffvb_cov(simmr_1,  ffvb_control = list(
  n_output = 3600,
  S = 100,
  P = 9,
  beta_1 = 0.9,
  beta_2 = 0.9,
  tau = 1000,
  eps_0 = 0.05,
  t_W = 50
))

