rm(list = ls())
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(gridExtra)

source("rfuncts/MVPIBP_functions.R")
source("rfuncts/gibbs_factorMVPIBP.R")
sourceCpp("rcppfuncts/MVPIBPh_nocov.cpp")
sourceCpp("rcppfuncts/bigMVP_nocov.cpp")

#--------------------------# Simulate data #---------------------------#
n = 80 # number of samples
alpha_data =seq(2, 40, by = 2)
my_seed = 4694
p = 300
type = c("factor", "ts", "common")
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
MSEsigma_IBP = MSEsigma_fMVPIBP = rep(0, length(all_data))
truncP = 0

##---------------## Run gibbs ##---------------##

for(i in 1:length(all_data)){
  d1 = all_data[[i]]$data_all
  
  # MVP-IBP factor model
  param_MVPIBP = list(truncP = truncP,
                      # parameter of gamma prior for alpha
                      a_alpha = mean(rowSums(d1)), b_alpha = 1, 
                      eps_MH = 2.5, # epsilon for MH step
                      gamma = 10,
                      # CUSP parameters
                      a_theta = 2, b_theta = 1.5, theta_inf = 0.05, 
                      # adaptation CUSP
                      start_adapt = 500, Kmax = p, a0_ad = -1, a1_ad = -5*10^(-4))
  fit_gibbs = gibbs_MVPIBP(d1, param_MVPIBP, seeds_g[i], Niter)
  
  # IBP
  param_IBP = list(truncP = truncP, 
                   # parameter of gamma prior for alpha
                   a_alpha = mean(rowSums(d1)), b_alpha = 1)
  fit_IBP = Beta_Binom_IBP(d1, param_IBP, Niter, seeds_g[i])
  
  # real values
  pi_0 = all_data[[i]]$pi
  Sigma0 = all_data[[i]]$Sigma
  MSEpi_IBP1 = MSEpi_MI = MSEsig_MI = rep(0, Niter-burnin)
  
  #--------# compute MSE  pi #--------#
  for(s in 1:(Niter-burnin)){
    MSEpi_IBP1[s] = mean((fit_IBP$prob[,burnin+s] - pi_0)^2)
    MSEpi_MI[s] = mean((pnorm(fit_gibbs$beta_tilde[,burnin+s]) - pi_0)^2)
  }
  MSEpi_fMVPIBP[i] = mean(MSEpi_MI)
  MSEpi_IBP[i] = mean(MSEpi_IBP1)
  
  #--------# compute MSE  Sigma #--------#
  for(s in 1:(Niter-burnin)){
    Lambda = fit_gibbs$Lambda[[burnin+s]]
    sq_D = sqrt(rowSums(Lambda^2) + 1)
    sigma_gibbs = diag(1/sq_D) %*% (Lambda %*% t(Lambda) +
                                      diag(NROW(Lambda))) %*% diag(1/sq_D)
    MSEsig_MI[s] = mean((sigma_gibbs  - Sigma0)^2)
  }
  MSEsigma_fMVPIBP[i] = mean(MSEsig_MI)
  MSEsigma_IBP[i] = mean((Sigma0 - diag(nrow(Sigma0)))^2)
  
  #--------# compute MSE  pstar #--------#
  MSEpstar_IBP1 = MSEpstar_MI = rep(0, Niter-burnin)
  pst_exp = sum(1-(1-pi_0)^n)
  for(s in 1:(Niter-burnin)){
    # MSE for MVP-IBP
    pstar_MVP = sum(1-(1-pnorm(fit_gibbs$beta_tilde[,burnin+s]))^n)
    MSEpstar_MI[s] = (pstar_MVP - pst_exp)^2
    # MSE for IBP
    pstar_IBP = sum(1-(1-fit_IBP$prob[,burnin+s])^n)
    MSEpstar_IBP1[s] = (pstar_IBP - pst_exp)^2
  }
  MSEps_fMVPIBP[i] = mean(MSEpstar_MI)
  MSEps_IBP[i] = mean(MSEpstar_IBP1)
  
  cat("simulation ",i,"out of",length(all_data),"\n")
}


### ------------- ### Two-stage alghoritm type simulations ###--------------####
###############################################################################

##--------------------## parameters ##---------------------##
max_it = 100
epsilon = 0.00001
m = 15
nmcmc = 200
burnin = 50
eta0_rho = 0 
nu0_rho = 0.01 
gamma0_rho = 0.01 # a_omega
lambda0_sq_rho = 0.01 # b_omega
eps_MH_ts = 0.2
b_alpha = 1
prior_var = 15
truncP_ts = 0

bigM = vector("list", length(all_data))
MIBP = vector("list", length(all_data))

##---------------## run bigMVP and MVP-IBP ##-------------##

for(i in 1:length(all_data)){
  d1 = all_data[[i]]$data_all
  # MVP-IBPh
  a_alpha = mean(rowSums(d1))
  MIBP[[i]] = MVP_IBPh(d1, eta0_rho, nu0_rho, gamma0_rho, lambda0_sq_rho, 
                       a_alpha, b_alpha, max_it, epsilon, m, nmcmc, burnin, 
                       eps_MH_ts, truncP_ts)
  # bigMVP
  bigM[[i]] = bigMVPh(d1, eta0_rho, nu0_rho, gamma0_rho, lambda0_sq_rho, max_it,
                      epsilon, m, nmcmc, burnin, prior_var, truncP_ts)
  cat("simulation ",i,"out of",length(all_data),"\n")
}

##---------------## compute errors ##---------------------##

Nsim = length(all_data)
npi_bMVP = npi_IMVP = nps_bMVP = nps_IMVP  = rep(NA, Nsim)
nsigma_bMVP = nsigma_IMVP = rep(NA, Nsim)
for(i in 1:Nsim){
  # real values
  pi_0 = all_data[[i]]$pi
  Sigma0 = all_data[[i]]$Sigma
  pst_exp = sum(1-(1-pi_0)^n)
  # estimated values
  bigM_fit = bigM[[i]]
  MIBP_fit = MIBP[[i]]
  
  #### compute errors
  # bigMVP
  npi_bMVP[i] = mean((pnorm(bigM_fit$coefficients) - pi_0)^2)
  pstar = sum(1-(1-pnorm(bigM_fit$coefficients))^n)
  nps_bMVP[i] = (pstar - pst_exp)^2
  nsigma_bMVP[i] = mean((bigM_fit$post_mean - Sigma0)^2)
  # MVP-IBP
  pstarM = sum(1-(1-pnorm(MIBP_fit$coefficients))^n)
  nps_IMVP[i] = (pstarM - pst_exp)^2
  npi_IMVP[i] = mean((pnorm(MIBP_fit$coefficients) - pi_0)^2)
  nsigma_IMVP[i] = mean((MIBP_fit$post_mean - Sigma0)^2)
}


### ----------------------- ### Plots results ###--------------------------####
###############################################################################

#----# construct data for plot #----#
# MSE data
MSEdata = cbind.data.frame(c(MSEpi_fMVPIBP, MSEps_fMVPIBP, MSEsigma_fMVPIBP, 
                             MSEpi_IBP, MSEps_IBP, MSEsigma_IBP), 
                           rep(c(rep("pi",60), rep("ps",60), rep("sigma",60)),2),
                           c(rep("MVP-IBP", 180), rep("IBP", 180)),
                           rep(c("factor", "tobit", "common"), 1200))
colnames(MSEdata) = c("MSE", "type","model", "scenario")
MSEdata$model = as.factor(MSEdata$model)
MSEdata$scenario = as.factor(MSEdata$scenario)
MSEdata$type = as.factor(MSEdata$type)

# Frobenius data
Fdata = cbind.data.frame(c(npi_IMVP,  nps_IMVP, nsigma_IMVP, npi_bMVP, 
                          nps_bMVP, nsigma_bMVP),
                         rep(c(rep("pi",60), rep("ps",60), rep("sigma",60)),2),
                         c(rep("MVP-IBP", 180), rep("IBP", 180)),
                         rep(c("factor", "tobit", "common"), 1200))
colnames(Fdata) = c("F1", "type","model", "scenario")
Fdata$model = as.factor(Fdata$model)
Fdata$scenario = as.factor(Fdata$scenario)
Fdata$type = as.factor(Fdata$type)


nl = c("pi"="pi", "ps"="Delta[10]", "sigma"="Sigma")
gibbsMSE = ggplot(MSEdata, aes(x = scenario, y = MSE, fill = model))+
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(values = c("steelblue", "red")) +
  facet_wrap(~ type, scales = "free", labeller = as_labeller(nl, default = label_parsed)) +
  xlab("") +
  theme_light() +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size=8),legend.key.size = unit(1,"line"),
        legend.box.spacing = unit(0.1,"line"),
        strip.text = element_text(size = 10, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"))
gibbsMSE

Fplot = ggplot(Fdata, aes(x = scenario, y = F1, fill = model))+
  geom_boxplot(alpha=0.7) +
  scale_fill_manual(values = c("forestgreen", "red")) +
  facet_wrap(~ type, scales = "free", labeller = as_labeller(nl, default = label_parsed)) +
  xlab("") +
  ylab("Frobenius error")+
  theme_light() +
  theme(legend.position = "top", legend.title = element_blank(),
        legend.text = element_text(size=8),legend.key.size = unit(1,"line"),
        legend.box.spacing = unit(0.1,"line"),
        strip.text = element_text(size = 10, colour = "black"),
        strip.background = element_rect(fill = "gray82"),
        panel.grid.major = element_line(size = 0.3, colour = "gray93"),
        panel.grid.minor = element_line(size = 0.15, colour = "gray93"))
Fplot

simp = grid.arrange(gibbsMSE, Fplot, nrow=2)
# ggsave(filename = "MSEgibbs.png", plot=simp,  width = 5.5, height = 5)
