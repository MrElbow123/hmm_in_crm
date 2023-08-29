library(Rcpp)
library(HDInterval)
sourceCpp('src/cat_mult_fw_cpp.cpp')

#load("data/nonverbal.rda")
#load("data/nonverbal_cov.rda")
load("data/data.rda")
load("data/data_cov.rda")
load("data/emiss_mu_int_bar.rda")
source('R/mHMM.R')
source('R/utility_func_mHMM.R')  # dif_matrix
source('R/int_to_prob.R')
source('R/forward_prob_cpp.R')  # cat_mult_fw_r_to_cpp
source('R/logl_mnl.R')         # llmnl_int
source('R/mnl_hess.R')
source('R/mnl_RW_once.R')

#sourceCpp('src/RcppExports.cpp')
library(MCMCpack)
library(mvtnorm)
m <- 3
n_dep <- 1
q_emiss <- c(2)

# specifying starting values
start_TM <- diag(.8, m)
start_TM[lower.tri(start_TM) | upper.tri(start_TM)] <- .2
start_EM <- list(
                 matrix(c(0.3, 0.4, 0.3,
                          0.3, 0.4, 0.3), byrow = TRUE, nrow = m,
                        ncol = q_emiss[1])
                )

n_subj <- 39
n_time <- 233
xt1 <- list(matrix(data_cov$sum_mess,ncol=n_subj,nrow=n_time))
xt2 <- list(matrix(data_cov$sum_mess2,ncol=n_subj,nrow=n_time))
xt3 <- list(matrix(data_cov$live_type,ncol=n_subj,nrow=n_time))
xt4 <- list(matrix(data_cov$user_type,ncol=n_subj,nrow=n_time))
xt5 <- list(matrix(data_cov$time_zone,ncol=n_subj,nrow=n_time))
xt6 <- list(matrix(data_cov$time_gap,ncol=n_subj,nrow=n_time))
xt7 <- list(matrix(data_cov$last_time,ncol=n_subj,nrow=n_time))
xt8 <- list(matrix(data_cov$host_number,ncol=n_subj,nrow=n_time))
xt9 <- list(matrix(data_cov$position,ncol=n_subj,nrow=n_time))
xt <- c(xt1,xt2,xt3,xt4,xt5,xt6,xt7,xt8,xt9)

xx <- rep( list(matrix( c(rep(1, n_subj),rowMeans(t(xt1[[1]])),rowMeans(t(xt2[[1]])),rowMeans(t(xt3[[1]])),rowMeans(t(xt4[[1]])),
                          rowMeans(t(xt5[[1]])),rowMeans(t(xt6[[1]])),rowMeans(t(xt7[[1]])),rowMeans(t(xt8[[1]])),rowMeans(t(xt9[[1]]))
                          ), ncol = 10, nrow = n_subj)), 1+n_dep)

out_2st_c <- mHMM(s_data = data, xx = xx, xt = xt, emiss_mu_int_bar = emiss_mu_int_bar,gen = list(m = m, n_dep = n_dep, q_emiss = q_emiss),start_val = c(list(start_TM), start_EM),mcmc = list(J = 50000, burn_in = 40000))
