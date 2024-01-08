rm(list = ls())
library(Rcpp)
library(RcppArmadillo)

source("rfuncts/MVPIBP_functions.R")
source("rfuncts/gibbs_factorMVPIBP.R")
sourceCpp("rcppfuncts/MVPIBPh_nocov.cpp")
sourceCpp("rcppfuncts/bigMVP_nocov.cpp")

#--------------------------# Simulate data #---------------------------#
n = 50 # number of samples
alpha_data = seq(5, 100, by = 5)
my_seed = 4694
p = 300
type = c("factor", "ts", "block")
all_data = simulation_data(alpha_data, n, p, type, my_seed)


###### ------------- ###### Gibbs type simulations ######----------------######
###############################################################################

## -------------## parameters ##----------------##
Niter = 2000
burnin = 500
set.seed(my_seed)
seeds_g = sample.int(5000, length(all_data))

MSEpi_fMVPIBP = MSEpi_IBP = rep(0, length(all_data))
MSEps_fMVPIBP = MSEps_IBP = rep(0, length(all_data))
truncP = 50
Npred = 100

##---------------## Run gibbs ##---------------##

for(i in 1:length(all_data)){
  d1 = all_data[[i]]$data 
  
  # MVP-IBP factor model
  param_MVPIBP = list(truncP = truncP, # truncation level for MVP-IBP
                      # parameter of gamma prior for alpha
                      a_alpha = mean(rowSums(d1)), b_alpha = 1, 
                      eps_MH = 2.5, # epsilon for MH step
                      gamma = 10, # put the expected number of factors
                      # CUSP parameters
                      a_theta = 2, b_theta = 1.5, theta_inf = 0.05, 
                      # adaptation
                      start_adapt = 500, Kmax = p, a0_ad = -1, a1_ad = -5*10^(-4))
  fit_gibbs = gibbs_MVPIBP(d1, param_MVPIBP, seeds_g[i], Niter)
  
  # IBP
  param_IBP = list(truncP = truncP, # set P for MVP-IBP
                   # parameter of gamma prior for alpha
                   a_alpha = mean(rowSums(d1)), b_alpha = 1)
  fit_IBP = Beta_Binom_IBP(d1, param_IBP, Niter, seeds_g[i])
  
  # real values
  pi_0 = c(all_data[[i]]$pi_active, rep(0, param_MVPIBP$truncP))
  MSEpi_IBP1 = MSEpi_MI = rep(0, Niter-burnin)
  
  #--------# compute MSE  pi #--------#
  for(s in 1:(Niter-burnin)){
    # MSE for pi IBPP
    MSEpi_IBP1[s] = mean((fit_IBP$prob[,burnin+s] - pi_0)^2)
    # MSE for pi MVP-IBP
    MSEpi_MI[s] = mean((pnorm(fit_gibbs$beta_tilde[,burnin+s]) - pi_0)^2)
  }
  MSEpi_fMVPIBP[i] = mean(MSEpi_MI)
  MSEpi_IBP[i] = mean(MSEpi_IBP1)
  
  #--------# compute MSE  pstar #--------#
  ep1_MVP = ep1_IBP = rep(0, Npred)
  Pd1 = NCOL(d1)
  nd1 = nrow(d1)
  pstar_true = all_data[[i]]$pstar
  
  # MVP-IBP
  betap = fit_gibbs$beta[,(Niter-Npred+1):Niter]
  Lambdap = fit_gibbs$Lambda[(Niter-Npred+1):Niter]
  etap = fit_gibbs$eta[(Niter-Npred+1):Niter]
  for(n in 1:Npred){
    datas = matrix(0, nd1, Pd1)
    for(k in 1:nd1){
      for(j in 1:Pd1){
        z = betap[j,n] + sum(Lambdap[[n]][j,]*etap[[n]][k,]) + rnorm(1)
        datas[k,j] = as.numeric(z>0)
      }
    }
    pst = sum(apply(datas != 0, 2, any))
    ep1_MVP[n] = (pst - pstar_true)^2
  }
  MSEps_fMVPIBP[i] = mean(ep1_MVP)
  
  # IBP
  ep1_IBP  = rep(0, Npred)
  pi_hat = fit_IBP$prob[,(Niter-Npred+1):Niter]
  for(n in 1:Npred){
    datasim1 = matrix(0, nd1, Pd1)
    for(k in 1:nd1){
      datasim1[k,]=rbinom(Pd1, 1, pi_hat)
    }
    pst = sum(apply(datasim1 != 0, 2, any))
    ep1_IBP[n] = (pst - pstar_true)^2
  }
  MSEps_IBP[i] = mean(ep1_IBP)
  
  cat("simulation ",i,"out of",length(all_data),"\n")
}

v1 = seq(1, length(all_data), by = 3)
v2 = seq(2, length(all_data), by = 3)
v3 = seq(3, length(all_data), by = 3)

#---# results #---#

c(median(MSEpi_fMVPIBP[v1]), median(MSEpi_IBP[v1]), median(MSEpi_fMVPIBP[v2]),
  median(MSEpi_IBP[v2]), median(MSEpi_fMVPIBP[v3]), median(MSEpi_IBP[v3]))
c(summary(MSEpi_fMVPIBP[v1])[5]-summary(MSEpi_fMVPIBP[v1])[2],
  summary(MSEpi_IBP[v1])[5]-summary(MSEpi_IBP[v1])[2],
  summary(MSEpi_fMVPIBP[v2])[5]-summary(MSEpi_fMVPIBP[v2])[2],
  summary(MSEpi_IBP[v2])[5]-summary(MSEpi_IBP[v2])[2],
  summary(MSEpi_fMVPIBP[v3])[5]-summary(MSEpi_fMVPIBP[v3])[2],
  summary(MSEpi_IBP[v3])[5]-summary(MSEpi_IBP[v3])[2])

c(median(MSEps_fMVPIBP[v1]), median(MSEps_IBP[v1]), median(MSEps_fMVPIBP[v2]),
  median(MSEps_IBP[v2]), median(MSEps_fMVPIBP[v3]), median(MSEps_IBP[v3]))
c(summary(MSEps_fMVPIBP[v1])[5]-summary(MSEps_fMVPIBP[v1])[2],
  summary(MSEps_IBP[v1])[5]-summary(MSEps_IBP[v1])[2],
  summary(MSEps_fMVPIBP[v2])[5]-summary(MSEps_fMVPIBP[v2])[2],
  summary(MSEps_IBP[v2])[5]-summary(MSEps_IBP[v2])[2],
  summary(MSEps_fMVPIBP[v3])[5]-summary(MSEps_fMVPIBP[v3])[2],
  summary(MSEps_IBP[v3])[5]-summary(MSEps_IBP[v3])[2])

### ------------- ### Two-stage alghoritm type simulations ###--------------####
###############################################################################

##--------------------## parameters ##---------------------##
max_it = 100
epsilon = 0.00001
m = 15
nmcmc = 150
burnin_emp = 50
eta0_rho = 0 
nu0_rho = 0.01 
gamma0_rho = 0.01 # a_omega
lambda0_sq_rho = 0.01 # b_omega
eps_MH_ts = 0.2
b_alpha = 1
prior_var = 15
truncP_ts = 50

bigM = vector("list", length(all_data))
MIBP = vector("list", length(all_data))

##---------------## run bigMVP and MVP-IBP ##-------------##

for(i in 1:length(all_data)){
  d1 = all_data[[i]]$data 
  # MVP-IBPh
  a_alpha = mean(rowSums(d1))
  MIBP[[i]] = MVP_IBPh(d1, eta0_rho, nu0_rho, gamma0_rho, lambda0_sq_rho, 
                       a_alpha, b_alpha, max_it, epsilon, m, nmcmc, burnin_emp, 
                       eps_MH_ts, truncP_ts)
  # bigMVP
  bigM[[i]] = bigMVPh(d1, eta0_rho, nu0_rho, gamma0_rho, lambda0_sq_rho, max_it, 
                      epsilon, m, nmcmc, burnin_emp, prior_var, truncP_ts)
  
  cat("simulation ",i,"out of",length(all_data),"\n")
}

##---------------## compute errors ##---------------------##

Nsim = length(all_data)
npi_bMVP = npi_IMVP = rep(NA, Nsim)
nsigma_bMVP = nsigma_IMVP = rep(NA, Nsim)
truncP_ts = 50
for(i in 1:Nsim){
  # real values
  pi_0 = c(all_data[[i]]$pi_active, rep(0, truncP_ts))
  Sigma0_act = all_data[[i]]$Sigma_active
  Sigma0 = rbind(cbind(Sigma0_act, matrix(0, nrow(Sigma0_act), truncP_ts)), 
                 matrix(0, truncP_ts, nrow(Sigma0_act) + truncP_ts))
  diag(Sigma0) = 1
  # estimated values
  bigM_fit = bigM[[i]]
  MIBP_fit = MIBP[[i]]

  #### compute errors
  # bigMVP
  npi_bMVP[i] = mean((pnorm(bigM_fit$coefficients) - pi_0)^2)
  nsigma_bMVP[i] = mean((bigM_fit$post_mean - Sigma0)^2)
  # MVP-IBP
  npi_IMVP[i] = mean((pnorm(MIBP_fit$coefficients) - pi_0)^2)
  nsigma_IMVP[i] = mean((MIBP_fit$post_mean - Sigma0)^2)
}

#---# results #---#

c(median(npi_IMVP[v1]), median(npi_bMVP[v1]),
  median(npi_IMVP[v2]), median(npi_bMVP[v2]),
  median(npi_IMVP[v3]), median(npi_bMVP[v3]))
c(summary(npi_IMVP[v1])[5] - summary(npi_IMVP[v1])[2],
  summary(npi_bMVP[v1])[5] - summary(npi_bMVP[v1])[2],
  summary(npi_IMVP[v2])[5] - summary(npi_IMVP[v2])[2],
  summary(npi_bMVP[v2])[5] - summary(npi_IMVP[v2])[2],
  summary(npi_IMVP[v3])[5] - summary(npi_IMVP[v3])[2],
  summary(npi_bMVP[v3])[5] - summary(npi_bMVP[v3])[2])

c(median(nsigma_IMVP[v1]), median(nsigma_bMVP[v1]),
  median(nsigma_IMVP[v2]), median(nsigma_bMVP[v2]),
  median(nsigma_IMVP[v3]), median(nsigma_bMVP[v3]))
c(summary(nsigma_IMVP[v1])[5] - summary(nsigma_IMVP[v1])[2],
  summary(nsigma_bMVP[v1])[5] - summary(nsigma_bMVP[v1])[2],
  summary(nsigma_IMVP[v2])[5] - summary(nsigma_IMVP[v2])[2],
  summary(nsigma_bMVP[v2])[5] - summary(nsigma_IMVP[v2])[2],
  summary(nsigma_IMVP[v3])[5] - summary(nsigma_IMVP[v3])[2],
  summary(nsigma_bMVP[v3])[5] - summary(nsigma_bMVP[v3])[2])