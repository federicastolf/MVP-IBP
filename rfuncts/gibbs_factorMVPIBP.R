library(truncnorm) # for rtruncnorm (truncated normal)
library(mvnfast) # for fast rmvn
library(Rcpp) 
library(RcppArmadillo)
Rcpp::sourceCpp("rcppfuncts/mat_mult.cpp") # for speed-up some matrix operation


gibbs_MVPIBP = function(data, param, my_seed, N_sampl){

  #------# set parameters #-------#
  set.seed(my_seed)
  n = dim(data)[1]
  data = cbind(data, matrix(0, n, param$truncP))
  p = dim(data)[2]
  u = runif(N_sampl)   # uniform rvs for adaptation
  btilde_sd = sqrt(2*log(p)) # sd for beta_tilde prior
  
  #--------# Variables to update #-------#
  z = vector("list", N_sampl) # probit latent variable z
  beta = beta_tilde = matrix(NA, p, N_sampl) # random intercept
  eta = vector("list", N_sampl) # factor scores matrix
  Lambda =  vector("list", N_sampl)   # factor loadings matrix
  K = rep(NA, N_sampl) # number of factors
  Kstar = rep(NA, N_sampl) # number of active factors
  theta_inv = vector("list", N_sampl)   # column-specific loadings precisions
  w_post = vector("list", N_sampl)   # stick-breaking weights for CUSP
  s_post = vector("list", N_sampl)   # augmented data for CUSP
  alpha = rep(NA, N_sampl) # alpha parameter of IBP-MVP
  
  # Initialization
  z[[1]] = matrix(0, n, p)
  beta_tilde[,1] = rnorm(p)
  Kstar[1] = floor(p/2)
  K[1] = Kstar[1] + 1
  Lambda[[1]] = matrix(rnorm(p*K[1]), p, K[1])
  eta[[1]] = matrix(rnorm(n*K[1]), n, K[1])
  theta_inv[[1]] = rep(1, K[1])
  w_post[[1]] = rep(1/K[1], K[1])
  s_post[[1]] = rep(K[1], K[1])
  alpha[1] = rgamma(1, param$a_alpha, param$b_alpha)
  # initilaization MH
  alpha_MH = alpha[1]
  
  # compute beta
  sq_D = sqrt(rowSums(Lambda[[1]]^2) + 1) 
  beta[,1] = beta_tilde[,1]/sq_D
  
  # mean for beta_tilde prior
  mu_btilde = sqrt(1 + btilde_sd^2)*qnorm(alpha[1]/(alpha[1] + p)) 
  mu_beta_vec = rep(mu_btilde, p)
  
  #------------------# Gibbs #-----------------#
  t0 = proc.time()
  for(t in 2:N_sampl){
    eta[[t]] = matrix(NA, n, K[t-1])
    z[[t]] =  matrix(NA, n, p)
    Lambda[[t]] = matrix(NA, p, K[t-1])
    theta_inv[[t]] = rep(NA, K[t-1])
    w_post[[t]] = rep(NA,K[t-1])
    s_post[[t]] = rep(NA,K[t-1])
    
    #---# 1) sample zeta #---#
    for(i in 1:n){
      for(j in 1:p){
        mu_z = beta[j,t-1] + Lambda[[t-1]][j,] %*% eta[[t-1]][i,]
        if(data[i,j]==0){
          z[[t]][i,j] =  rtruncnorm(1, mean = mu_z, sd = 1, a = -Inf, b = 0)
        }
        else{
          z[[t]][i,j] =  rtruncnorm(1, mean = mu_z, sd = 1, a = 0, b = Inf)
        }
      }
    }
    
    #---# 2) sample eta #---#
    V_eta = chol2inv(chol(eigenMapMatMult(t(Lambda[[t-1]]), Lambda[[t-1]]) 
                          + diag(K[t-1])))
    for(i in 1:n){
      mu_eta = V_eta %*% (t(Lambda[[t-1]]) %*% (z[[t]][i,] - beta[,t-1]))
      eta[[t]][i,] = mvnfast::rmvn(1, mu_eta, V_eta)
    }
    
    #---# 3) sample Lambda #---#
    V_j = chol2inv(chol(diag(theta_inv[[t-1]], nrow = K[t-1]) + 
                          t(eta[[t]]) %*% eta[[t]]))
    mult_mu = V_j %*% t(eta[[t]]) 
    for (j in 1:p){
      mu_j = mult_mu %*% (z[[t]][,j]- rep(beta[j,t-1], n))
      Lambda[[t]][j,] = mvnfast::rmvn(1, mu_j, V_j)
    }
    
    #---# 4) sample beta tilde #---#
    # vec of dimension p with element of D
    sq_D = sqrt(rowSums(Lambda[[t]]^2) + 1) 
    # compute variance of beta_tilde
    inv_tmp =  inverse_SHM(Lambda[[t]], c = sq_D^2*n*btilde_sd^2 + 1)
    Psi = eigenMapMatMult(Lambda[[t]], t(Lambda[[t]])) + diag(p) 
    inv_psi2 = eigenMapMatMult(inv_tmp, Psi)
    V_btilde = t(t(inv_psi2*(btilde_sd^2*sq_D))/sq_D) 
    # compute mean of beta_tilde
    invPsi = inverse_SHM(Lambda[[t]]) # inverse of Psi
    invSigma = t(t(invPsi*sq_D)*sq_D)  # inverse of Sigma
    sumz = colSums(z[[t]])
    mb_tmp = t(t(mu_beta_vec)*rep(btilde_sd^(-2), p)) +
      eigenMapMatMult(invSigma, sumz)
    mean_bt = eigenMapMatMult(V_btilde, mb_tmp)
    beta_tilde[,t] = mvnfast::rmvn(1, mean_bt, V_btilde)
    
    ## update beta
    beta[,t] = beta_tilde[,t]/sq_D
    
    #---# 5) update alpha through MH step #---#
    logp = logpost_alpha(alpha_MH, beta_tilde[,t], param$a_alpha, param$b_alpha)
    #logp = logpost_alpha_log(alpha_MH, beta_tilde[,t], param$a_alpha, param$b_alpha)
    alpha_MH_new = rnorm(1, alpha_MH, param$eps_MH)
    logp_new = logpost_alpha(alpha_MH_new, beta_tilde[,t],
                             param$a_alpha, param$b_alpha)
    # logp_new = logpost_alpha_log(alpha_MH_new, beta_tilde[,t],
    #                          param$a_alpha, param$b_alpha)
    a_acc = min(1, exp(logp_new - logp))
    if (runif(1) < a_acc){
      logp = logp_new
      alpha_MH = alpha_MH_new
    }
    alpha[t] = alpha_MH
    #alpha[t] = exp(alpha_MH)
    
    # mean for beta_tilde prior
    mu_btilde = sqrt(1 + btilde_sd^2)*qnorm(alpha[t]/(alpha[t] + p)) 
    mu_beta_vec = rep(mu_btilde, p)
    
    #---# 6) sample s_augmented #---#
    lhd_spike = rep(0, K[t-1])
    lhd_slab = rep(0, K[t-1])
    for (h in 1:K[t-1]){
      lhd_spike[h] = exp(sum(log(dnorm(Lambda[[t]][,h], mean = 0, 
                                       sd = param$theta_inf^(1/2), 
                                       log = FALSE))))
      lhd_slab[h] = mvnfast::dmvt(X = Lambda[[t]][,h], mu = rep(0,p), 
                    sigma = (param$b_theta/param$a_theta)*diag(p),
                    df = 2*param$a_theta)
      prob_h = w_post[[t-1]]*c(rep(lhd_spike[h], h), rep(lhd_slab[h], K[t-1] - h))
      if (sum(prob_h) == 0){
        prob_h = c(rep(0, K[t-1] - 1), 1)
      }
      else{
        prob_h = prob_h/sum(prob_h)
      }
      s_post[[t]][h] = c(1:K[t-1])%*%rmultinom(n = 1, size = 1, prob = prob_h)
    }
    
    #---# 7) sample v and update omega #---#
    v = rep(NA, K[t-1])
    for (h in 1:(K[t-1] - 1)){
      v[h] = rbeta(1, shape1 = 1 + sum(s_post[[t]] == h), 
                   shape2 = param$gamma + sum(s_post[[t]]>h))
    }
    v[K[t-1]] = 1
    w_post[[t]][1] = v[1]
    for (h in 2:K[t-1]){
      w_post[[t]][h] = v[h]*prod(1-v[1:(h-1)])  
    }
    
    #---# 8) sample theta^{-1} #---#
    for (h in 1:K[t-1]){
      if (s_post[[t]][h]<=h){
        theta_inv[[t]][h] = param$theta_inf^(-1)
      }
      else{
        theta_inv[[t]][h] = rgamma(n = 1, shape = param$a_theta + 0.5*p, 
                                   rate = param$b_theta + 
                                     0.5*t(Lambda[[t]][,h]) %*% Lambda[[t]][,h])
      }
    }
    
    #---# update K[t] #---#
    active = which(s_post[[t]] > c(1:K[t-1]))
    Kstar[t] = length(active)
    K[t] = K[t-1]
    

    # adpatation for K
    if (t >= param$start_adapt & u[t] <= exp(param$a0_ad + param$a1_ad*t)){
      if (Kstar[t] < K[t-1] - 1){
        # set truncation to Kstar[t] and subset all variables, keeping only active columns
        K[t] = Kstar[t]+1
        eta[[t]] = cbind(eta[[t]][,active], rnorm(n))
        theta_inv[[t]] = c(theta_inv[[t]][active], param$theta_inf^(-1))
        w_post[[t]] = c(w_post[[t]][active], 1 - sum(w_post[[t]][active]))
        Lambda[[t]] = cbind(Lambda[[t]][,active], 
                            rnorm(p, mean = 0, sd = sqrt(param$theta_inf)))
      } else if (K[t-1] < param$Kmax) {
        # increase truncation by 1 and extend all variables, sampling from the prior/model
        K[t] = K[t-1]+1
        eta[[t]] = cbind(eta[[t]], rnorm(n))
        v[K[t-1]] = rbeta(1, shape1 = 1, shape2 = param$gamma)
        v = c(v,1)
        w_post[[t]] = rep(NA, K[t])
        w_post[[t]][1] = v[1]
        for (h in 2:K[t]){
          w_post[[t]][h] = v[h]*prod(1-v[1:(h-1)])  
        }
        theta_inv[[t]] = c(theta_inv[[t]], param$theta_inf^(-1))
        Lambda[[t]] = cbind(Lambda[[t]],
                            rnorm(p, mean = 0, sd = sqrt(param$theta_inf)))
      } 
    }
  }
  
  runtime = proc.time()-t0
  output = list("data" = data, "my_seed" = my_seed, "N_sampl" = N_sampl, 
                "param" = param, "runtime" = runtime, "Lambda" = Lambda,  
                "eta" = eta, "beta_tilde" = beta_tilde, "z" = z, "beta" = beta,
                "K" = K, "Kstar" = Kstar, "theta_inv" = theta_inv, 
                "w_post" = w_post, "s_post" = s_post, "alpha" = alpha)
  return(output)
}



