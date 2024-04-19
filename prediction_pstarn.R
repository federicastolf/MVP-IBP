
library(Rcpp) 
library(RcppArmadillo)
library(bigMVP)

rm(list=ls())
source("rfuncts/MVPIBP_functions.R")
Rcpp::sourceCpp("rcppfuncts/MVPIBPh.cpp")

#-----# load data #-----#
load("Data/X_fungi.Rdata")
load("data/fungi_binary.Rdata")
cov_data = as.data.frame(X_fungi)
Nrep = 100

#--------# observed pstar #---------#

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
bigMVPfit1 = bigMVPh(Y1b, X_train1[,-1], eta01, nu0, gamma01, Lambda01, 
                     eta0_rho, nu0_rho, gamma0_rho, lambda0_sq_rho, max_it, 
                     epsilon, m, nmcmc, burnin)
bigMVPfit2 = bigMVPh(Y2b, X_train2[,-1], eta02, nu0, gamma02, Lambda02, 
                     eta0_rho, nu0_rho, gamma0_rho, lambda0_sq_rho, max_it, 
                     epsilon, m, nmcmc, burnin)


#------# predictions for p*_n #-------#
pstar1_bM = pred_pstarMVP(bigMVPfit1, Nrep, cov1, fungi1)[,(nf1-Nobs_test1):nf1]
pstar2_bM = pred_pstarMVP(bigMVPfit2, Nrep, cov2, fungi2)[,(nf2-Nobs_test2):nf2]


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
pstar1_MVPIBP = pred_pstarMVP(mvpIBPfit1, Nrep, cov1, fungi1)[,(nf1-Nobs_test1):nf1]
pstar2_MVPIBP = pred_pstarMVP(mvpIBPfit2, Nrep, cov2, fungi2)[,(nf2-Nobs_test2):nf2]

###-------------------------###  IBP ####------------------------------####

#---------# fit IBP #-------#
Y_train1IBP = cbind(fungi1[-sampled_rows1, ])
Y_train2IBP = cbind(fungi2[-sampled_rows2, ])
MCMC_IBP = 2000
burnin_IBP = 1000
param_IBP = list(truncP = truncP, a_alpha = a_alpha, b_alpha = b_alpha)
mseed = 454
IBP1 = Beta_Binom_IBP(Y_train1IBP, param_IBP, MCMC_IBP, mseed)
IBP2 = Beta_Binom_IBP(Y_train2IBP, param_IBP, MCMC_IBP, mseed)

#-------# predictions #-------#
pstar1_IBP = pred_pstarIBP(IBP1, Nrep, burnin_IBP, MCMC_IBP,
                           fungi1)[,(nf1-Nobs_test1):nf1]
pstar2_IBP = pred_pstarIBP(IBP2, Nrep, burnin_IBP, MCMC_IBP, 
                           fungi2)[,(nf2-Nobs_test2):nf2]



##############################################################################
####----------------------##### compute errors ####----------------------##### 

# out-of-sample true values
ptrue1_test = pstar1_obs[(nf1-Nobs_test1):nf1]
ptrue2_test = pstar2_obs[(nf2-Nobs_test2):nf2]

# MSE for site 1 and 2
e1_IBP = e1_MVP = e1_MVPIBP = rep(0, length(ptrue1_test))
e2_IBP = e2_MVP = e2_MVPIBP = rep(0, length(ptrue1_test))

# site 1
for(i in 1:length(ptrue1_test)){
  e1_IBP[i] = mean((pstar1_IBP[,i]-ptrue1_test[i])^2)
  e1_MVP[i] = mean((pstar1_bM[,i]-ptrue1_test[i])^2)
  e1_MVPIBP[i] = mean((pstar1_MVPIBP[,i]-ptrue1_test[i])^2)
}

# site 2
for(i in 1:length(ptrue2_test)){
  e2_IBP[i] = mean((pstar2_IBP[,i]-ptrue2_test[i])^2)
  e2_MVP[i] = mean((pstar2_bM[,i]-ptrue2_test[i])^2)
  e2_MVPIBP[i] = mean((pstar2_MVPIBP[,i]-ptrue2_test[i])^2)
}

