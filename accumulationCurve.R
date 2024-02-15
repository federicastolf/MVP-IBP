
library(Rcpp) 
library(RcppArmadillo)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(bigMVP)

rm(list=ls())
source("rfuncts/MVPIBP_functions.R")
Rcpp::sourceCpp("rcppfuncts/MVPIBPh.cpp")

#-----# load data #-----#
load("Data/X_fungi.Rdata")
load("data/fungi_binary.Rdata")
cov_data = as.data.frame(X_fungi)
Nrep = 50

#--------# observed p^*_n #---------#

## site 1 and 2
fungi1 = fungi[cov_data$site_id1 == 1,]
cov1 = cov_data[cov_data$site_id1 == 1,]
nf1 = nrow(fungi1)
fungi2 = fungi[cov_data$site_id1 == 0 & cov_data$site_id3 == 0  & 
                 cov_data$site_id4 == 0,]
cov2 = cov_data[cov_data$site_id1 == 0 & cov_data$site_id3 == 0  & 
                  cov_data$site_id4 == 0,]
nf2 = nrow(fungi2)
pstar1_obs = compute_pstar(fungi1)
pstar2_obs = compute_pstar(fungi2)

#---------------# test and train set #-------------------#

set.seed(854)
# site 1
Nobs_test1 = round(nf1*0.25)
sampled_rows1 = sample(nf1, Nobs_test1)
row_indices1 = which(cov_data$site_id1 == 1)[sampled_rows1]
X_train1 = X_fungi[-row_indices1,]
Y_train1 = fungi[-row_indices1, ]
# site2
Nobs_test2 = round(nf2*0.25)
sampled_rows2 = sample(nf2, Nobs_test2)
row_indices2 = which(cov_data$site_id1 == 0 & cov_data$site_id3 == 0
                     & cov_data$site_id4 == 0)[sampled_rows2]
X_train2 = X_fungi[-row_indices2,]
Y_train2 = fungi[-row_indices2, ]

##########################################################################
####----------------#### out of sample prediction #####--------------#####
##########################################################################

####--------------------------#### fit bigMVPh ####-------------------------####

# NB in bigMVP the input model matrix X doesn't include the intercept
q1 = ncol(X_train1) -1
q2 = ncol(X_train2)-1
eta01 = rep(0, q1+1)
eta02 = rep(0, q2+1)
nu0 = 0.001
gamma01 = q1+1
gamma02 = q2+1
Lambda01 = diag(q1+1)
Lambda02 = diag(q1+1)
eta0_rho = 0
nu0_rho = 0.01
gamma0_rho = 0.01
lambda0_sq_rho = 0.01
max_it = 100
epsilon = 0.0001
m = 15
nmcmc = 200
burnin = 50
truncP = 200
Y1b = cbind(Y_train1, matrix(0, nrow(Y_train1), truncP))
Y2b= cbind(Y_train2, matrix(0, nrow(Y_train2), truncP))
set.seed(434)
bigMVPfit1 = bigMVPh(Y1b, X_train1[,-1], eta01, nu0, gamma01, Lambda01, 
                     eta0_rho, nu0_rho, gamma0_rho, lambda0_sq_rho, max_it, 
                     epsilon, m, nmcmc, burnin)
bigMVPfit2 = bigMVPh(Y2b, X_train2[,-1], eta02, nu0, gamma02, Lambda02, 
                     eta0_rho, nu0_rho, gamma0_rho, lambda0_sq_rho, max_it, 
                     epsilon, m, nmcmc, burnin)


#------# predictions for p*_n #-------#
mseed = 354
# site 1
pstar1_bM = pred_pstarMVP(bigMVPfit1, Nrep, mseed, cov1, fungi1)
pstar1q_bM = apply(pstar1_bM, 2, quantile, probs=c(0.05,0.95))
pstar1m_bM = apply(pstar1_bM, 2, mean)
# site 2
pstar2_bM = pred_pstarMVP(bigMVPfit2, Nrep, mseed, cov2, fungi2)
pstar2q_bM = apply(pstar2_bM, 2, quantile, probs=c(0.05,0.95))
pstar2m_bM = apply(pstar2_bM, 2, mean)


###-----------------### fit hierarchical MVP-IBP model ####------------------###
a_alpha = 17.5
b_alpha = 0.25
eps_MH = 0.2
prior_var = 15

mvpIBPfit1 = MVP_IBPh(Y_train1, X_train1, eta0_rho, nu0_rho, gamma0_rho, lambda0_sq_rho, 
                      a_alpha, b_alpha,prior_var,max_it, epsilon, m, nmcmc, burnin, 
                      eps_MH, truncP)
mvpIBPfit2 = MVP_IBPh(Y_train2, X_train2, eta0_rho, nu0_rho, gamma0_rho, lambda0_sq_rho, 
                      a_alpha, b_alpha,prior_var,max_it, epsilon, m, nmcmc, burnin, 
                      eps_MH, truncP)  

#------# predictions for p*_n #-------#

# site 1
pstar1_MVPIBP = pred_pstarMVP(mvpIBPfit1, Nrep, mseed, cov1, fungi1)
pstar1o_q = apply(pstar1_MVPIBP, 2, quantile, probs=c(0.05,0.95))
pstar1o_m = apply(pstar1_MVPIBP, 2, mean)

# site 2
pstar2_MVPIBP = pred_pstarMVP(mvpIBPfit2, Nrep, mseed, cov2, fungi2)
pstar2o_q = apply(pstar2_MVPIBP, 2, quantile, probs=c(0.05,0.95))
pstar2o_m = apply(pstar2_MVPIBP, 2, mean)


###-------------------------###  IBP ####------------------------------####

#---------# fit IBP #-------#
Y_train1IBP = cbind(fungi1[-sampled_rows1, ])
Y_train2IBP = cbind(fungi2[-sampled_rows2, ])
MCMC_IBP = 2000
burnin_IBP = 1000
param_IBP = list(truncP = truncP, a_alpha = a_alpha, b_alpha = b_alpha)
IBP1 = Beta_Binom_IBP(Y_train1IBP, param_IBP, MCMC_IBP, mseed)
IBP2 = Beta_Binom_IBP(Y_train2IBP, param_IBP, MCMC_IBP, mseed)

#-------# predictions #-------#
pstar1_IBP = pred_pstarIBP(IBP1, Nrep, mseed, burnin_IBP, MCMC_IBP, fungi1)
pstar1IBP_q = apply(pstar1_IBP, 2, quantile, probs=c(0.05,0.95))
pstar1IBP_m = apply(pstar1_IBP, 2, mean)

## site 2
pstar2_IBP = pred_pstarIBP(IBP2, Nrep, mseed, burnin_IBP, MCMC_IBP, fungi2)
pstar2IBP_q = apply(pstar2_IBP, 2, quantile, probs=c(0.05,0.95))
pstar2IBP_m = apply(pstar2_IBP, 2, mean)


###################################################
####-------------##### PLOT ####-------------##### 


data_site1 = cbind.data.frame(c(1:nf1),  t(pstar1o_q), pstar1o_m, t(pstar1IBP_q),
                              pstar1IBP_m, t(pstar1q_bM), pstar1m_bM, pstar1_obs, 
                              rep(nf1 - Nobs_test1, nf1), rep("Site1", nf1))
colnames(data_site1) = c("n", "lb_MV", "ub_MV", "median_MV", "lb_IB", "ub_IB",
                         "median_IB", "lb_bM", "ub_bM", "median_bM","obs",
                         "test_tresh","site")
data_site2 = cbind.data.frame(c(1:nf2),  t(pstar2o_q), pstar2o_m, t(pstar2IBP_q),
                              pstar2IBP_m, t(pstar2q_bM), pstar2m_bM, pstar2_obs, 
                              rep(nf2 - Nobs_test2, nf2), rep("Site2", nf2))
colnames(data_site2) = c("n", "lb_MV", "ub_MV", "median_MV", "lb_IB", "ub_IB",
                         "median_IB","lb_bM", "ub_bM", "median_bM", "obs", 
                         "test_tresh","site")
dataplot = rbind.data.frame(data_site1, data_site2)

shadealpha = 0.4
linesize = 1
szp = 2
stp = 1
f1 = ggplot(data = dataplot) +
  geom_ribbon(aes(x = n, ymax = ub_IB, ymin = lb_IB), alpha = shadealpha,  fill = "steelblue") +
  geom_line(aes(n, ub_IB), linetype= "dotted", color = "steelblue") +
  geom_line(aes(n, lb_IB), linetype = "dotted", color = "steelblue") +
  geom_line(aes(n, median_IB, color = "forestgreen"), linewidth = linesize) +
  geom_ribbon(aes(x = n, ymax = ub_bM, ymin = lb_bM), alpha = shadealpha, fill = "forestgreen") +
  geom_line(aes(n, lb_bM), linetype= "dotted", color = "forestgreen") +
  geom_line(aes(n, ub_bM), linetype = "dotted", color = "forestgreen") +
  geom_line(aes(n, median_bM, color = "steelblue"), linewidth = linesize) +
  geom_ribbon(aes(x = n, ymax = ub_MV, ymin = lb_MV), alpha = shadealpha, fill = "red") +
  geom_line(aes(n, lb_MV), linetype= "dotted", color = "red") +
  geom_line(aes(n, ub_MV), linetype = "dotted", color = "red") +
  geom_line(aes(n, median_MV, color = "darkred"), linewidth = linesize) +
  geom_point(aes(n, obs), color = 'black', shape = 19,
             size = stp, stroke = stp, alpha=0.5) +
  scale_color_manual(values=c("darkred","steelblue", "forestgreen"),
                     labels = c("MVP-IBP", "IBP", "bigMVP")) +
  geom_vline(aes(xintercept = test_tresh), linetype= "longdash", linewidth = 0.8) +
  labs(y = TeX("$p^*_n$"), x= "n") +
  theme_light() + 
  facet_wrap(~ site, scales = "free") +
  theme(text = element_text(size = 16), legend.title=element_blank(),
        legend.text = element_text(size=16), 
        legend.position = "top", strip.text = element_text(size=16))

f1

# ggsave(filename = "pstar_fungi1.png", plot=f1,  width = 12, height = 7)
