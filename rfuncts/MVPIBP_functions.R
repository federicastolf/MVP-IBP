require(Matrix)
require(mvnfast)

######################## functions GIBBS MVP-IBP #######################
########################################################################

inverse_SHM = function(M, c = 1){
  "----------------------------------------------------------------------------
  compute efficiently the inverse of a matrix with structure MM^T + cI 
  exploiting the Sherman–Morrison–Woodbury formula
  ----------------------------------------------------------------------------"
  psi0 = diag(1/c, NROW(M))
  for (i in 1:NCOL(M)){
    lambdak = M[,i]
    wk = eigenMapMatMult(psi0, lambdak)
    psi_k = psi0 - (eigenMapMatMult(wk, t(wk)))/(1+ sum(lambdak*wk))
    psi0 = psi_k
  }
  return(psi_k)
}

compute_invR = function(Lambda){
  sd_psi = sqrt(rowSums(Lambda^2) + 1) 
  inv_psi = inverse_SHM(Lambda)
  invR = t(t(inv_psi*sd_psi)*sd_psi)
  invR
}

# Log likelihood for alpha for MH
log_lik_alpha = function(alpha, beta){
  p = length(beta)
  tau2 = 2*log(p)
  mu = sqrt(1 + tau2)*qnorm(alpha/(alpha + p))
  if (alpha > 0) {
    logl = dnorm(beta, mu, sqrt(tau2), log = TRUE)
    return(sum(logl))
  } else {
    return(-Inf)
  }
}

# Log-posterior for MH - log parametrization
logpost_alpha_log = function(alpha, beta, a_alpha, b_alpha){
  log_likelihood = log_lik_alpha(alpha, beta)
  return(log_likelihood + dgamma(exp(alpha), a_alpha, b_alpha, log = TRUE) + alpha)
}

# Log-posterior for MH 
logpost_alpha = function(alpha, beta, a_alpha, b_alpha){
  log_likelihood = log_lik_alpha(alpha, beta)
  return(log_likelihood + dgamma(alpha, a_alpha, b_alpha, log = TRUE))
}

################################ simulation #################################
#############################################################################

synthetic_data = function(n, p, alpha, my_seed, type, k_factor = 10){
  "----------------------------------------------------------------------------
  simulate binary data matrix in different scenarios: 
  - factor (MVP-IBP factor model)
  - common (MVP-IBP where Sigma has an equi-correlation structure)
  - ts (MVP-IBP where the latent variables follow a t-student distribution)
  ----------------------------------------------------------------------------"
  set.seed(my_seed)
  b_sd = sqrt(2*log(p))
  b_mu = sqrt(1+b_sd^2)*qnorm(alpha/(alpha + p))
  beta = rnorm(p, b_mu, b_sd)
  data_MVPIBP = matrix(0, nrow = n, ncol = p)
  z = rep(NA, p)

  
  # factor 
  if(type=="factor"){
    Lambda = matrix(rnorm(p*k_factor, mean = 1), p, k_factor)
    beta_tilde = beta*sqrt(rowSums(Lambda^2) +1)
    for (i in 1:n){
      eps = rnorm(p)
      eta = rnorm(k_factor)
      z = beta_tilde + Lambda %*% eta + eps
      data_MVPIBP[i,] = as.numeric(z>0)
    }
    sq_D = sqrt(rowSums(Lambda^2) + 1)
    Sigma =  diag(1/sq_D) %*% (Lambda %*% t(Lambda) + diag(p)) %*% diag(1/sq_D)
    beta = beta_tilde
  }

  #common
  if(type=="common"){
    rhot = runif(1, 0, 0.8)
    Sigma = matrix(rhot, p, p)
    diag(Sigma) = 1
    for(i in 1:n){
      eps = mvnfast::rmvn(1, mu=rep(0,p), sigma = Sigma)
      z = beta + eps
      data_MVPIBP[i,] = as.numeric(z>0)
    }
  }
  
  # t-student
  if(type=="ts"){
    Lambda = matrix(rnorm(p*k_factor), p, k_factor)
    beta_tilde = beta*sqrt(rowSums(Lambda^2) +1)
    for (i in 1:n){
      eps = rt(p, df = 10)
      eta = rnorm(k_factor)
      z = beta_tilde + Lambda %*% eta + eps
      data_MVPIBP[i,] = as.numeric(z>0)
    }
    sq_D = sqrt(rowSums(Lambda^2) + 1)
    Sigma =  diag(1/sq_D) %*% (Lambda %*% t(Lambda) + diag(p)) %*% diag(1/sq_D)
    beta = beta_tilde
  }
  active = apply(data_MVPIBP != 0, 2, any)
  pstar = sum(active)
  data_active = data_MVPIBP[,active]
  pi = pnorm(beta)
  pi_active =  pnorm(beta)[active]
  Sigma_active = Sigma[active, active]
  return(list("data_all" = data_MVPIBP, "Sigma" = Sigma, "beta" = beta, 
              "pi" = pi, "Sigma_active" = Sigma_active, "pi_active" = pi_active,
              "active"= active, "pstar" = pstar, "data" = data_active,
              "type" = type))
}

simulation_data = function(alpha, n, p, type){
  Ndata = length(alpha)*length(type)
  seeds_D = sample.int(5000, Ndata)
  dataALLsynth = vector("list", Ndata)
  for(i in 1:length(type)){
    for(k in 1:length(alpha)){
      dataALLsynth[[(k-1)*length(type) + i]] = synthetic_data(n, p, alpha[k],
                                                 seeds_D[(k-1)*length(type) + i],
                                                type[i])
    }
  }
  dataALLsynth
}


############################ Gibbs for IBP ################################
###########################################################################


Beta_Binom_IBP = function(data, param, Niter, my_seed){
  "----------------------------------------------------------------------------
  Gibbs for truncated Beta-Binomial rapresentation of the IBP
  ----------------------------------------------------------------------------"
  set.seed(my_seed)
  n = dim(data)[1]
  data = cbind(data, matrix(0, n, param$truncP))
  p = dim(data)[2]
  m_j = colSums(data) # occurence for each species
  # posterior draws of the logarithmic occurrence probabilities within each habitats
  PROB = matrix(0, p, Niter)  
  alpha = rep(0, Niter) # posterior draws of the alpha parameters
  # initialize alpha
  alpha[1] = rgamma(1, shape = param$a_alpha, rate = param$b_alpha)  
  # Gibbs sampling steps
  for (iter in 2:Niter){
    # full conditional for probabilities
    for (j in 1:p){
      PROB[j, iter] = log(rbeta(1, alpha[iter-1]/p + m_j[j], n - m_j[j] + 1))
    }
    # full conditional for alphas
    alpha[iter] = rgamma(1, param$a_alpha, param$b_alpha - sum(PROB[iter])/p)
  }
  expP = exp(PROB)
  return(list("prob" = expP, "alpha" = alpha))
}


############################ plot ############################################
###############################################################################

build_phylo = function(covariate, fit_MVPIBP, tax_matrix, Xmodel){
  "----------------------------------------------------------------------------
  Return a phyloseq object to make Krona wheels plot
  - covariate: covariate from w
  - fit_MVP_IBP: output from a hierarchical MVP-IBP model
  - tax_matrix: data frame with taxonoical tree classification of data 
  - Xmodel: model matrix covariate used for fitting fit_MVPIBP
  ----------------------------------------------------------------------------"
  
  p = ncol(fit_MVPIBP$Y)
  n = nrow(fit_MVPIBP$Y)
  
  # posterior mean of marginal species occurence
  pi = pnorm(Xmodel %*% fit_MVPIBP$coefficients)[,1:p]
  
  # otu_matrix
  otu_mat = as.data.frame(t(pi))
  rownames(otu_mat) = paste0("otu", 1:p)
  colnames(otu_mat) = paste0("s", 1:n)
  otu_mat = as.matrix(otu_mat)
  
  # construct samples_df
  covariate_fungi = cbind.data.frame(c(1:n), covariate_fungi)
  colnames(covariate_fungi)[1] = "sample"
  covariate_fungi = cbind.data.frame(covariate_fungi, rep("a",NROW(covariate_fungi)))
  colnames(covariate_fungi)[4] = "no_cov"
  samples_df = covariate_fungi %>% 
    tibble::column_to_rownames("sample") 
  rownames(samples_df) = paste0("s", 1:n)
  colnames(samples_df)[1]="site_id"
  
  # obtain phyloseq object
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE) 
  samples = sample_data(samples_df)
  TAX = tax_table(tax_matrix)
  fungi_phylo = phyloseq(OTU, TAX, samples)
  return(fungi_phylo)
}


compute_pstar = function(data){
  # compute pstar for all rows of data
  n = nrow(data)
  pstar = rep(0, n)
  pstar[1] = sum(data[1,])
  for(i in 2:n){
    datai = data[c(1:i),]
    pstar[i] = sum(apply(datai != 0, 2, any))
  }
  pstar
}


pred_pstarMVP = function(fit_MVPIBP, Nrep, covariate, Y){
  nf = nrow(covariate)
  X = as.matrix(covariate)
  p = ncol(Y)
  beta_est = fit_MVPIBP$coefficients[,1:p]
  q = nrow(beta_est)
  pstar_out = matrix(0, Nrep, nf)
  for(n in 1:Nrep){
    beta_sample = matrix(NA, q, NCOL(beta_est))
    for(j in 1:NCOL(beta_est)){
      sigma_est = matrix(fit_MVPIBP$cov_mats[,j], q, q)
      muj = beta_est[,j]
      beta_sample[,j] = mvnfast::rmvn(1, muj, sigma_est)
    }
    pi_hat = pnorm(X %*% beta_sample)
    pi_mean = apply(pi_hat,2,mean)
    datasim1 = matrix(0, nrow(pi_hat), length(pi_mean))
    for(i in 1:nrow(pi_hat)){
      datasim1[i,]=rbinom(length(pi_mean), 1, pi_mean)
    }
    pstar_out[n,] = compute_pstar(datasim1)
  }
  pstar_out
}


pred_pstarIBP = function(IBP_fit, Nrep, burnin, MCMC, data){
  p = ncol(data)
  pi_hat = apply(IBP_fit[[1]][,burnin:MCMC],1,mean)[1:p]
  pstar_IBP = matrix(0, Nrep, nrow(data))
  for(n in 1:Nrep){
    datasim1 = matrix(0,nrow(data), length(pi_hat))
    for(i in 1:nrow(data)){
      datasim1[i,]=rbinom(length(pi_hat),1,pi_hat)
    }
    pstar_IBP[n,] = compute_pstar(datasim1)
  }
  pstar_IBP
}
