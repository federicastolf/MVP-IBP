#include "RcppArmadillo.h"
#include "Rcpp.h"

// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;


List GLRcpp(int n, double a, double b){
  const arma::vec& i = arma::linspace(1, n-1, n-1);
  arma::mat J(n,n,arma::fill::zeros);
  const arma::vec& d = i / (arma::sqrt(4 * arma::square(i) - 1));
  
  // Setting off-diagonal elements
  J.diag(1) = d;
  J.diag(-1) = d;
  // Initialize matrix and vector for eigenvalues/vectors
  arma::vec L;
  arma::mat V;
  arma::eig_sym(L, V, J);
  
  arma::vec w = 2 * arma::vectorise(arma::square(V.row(0)));
  arma::vec x = .5 * ((b - a) * L + a + b);
  w = -.5 * (a - b) * w;
  return List::create(
    Named("x") = x,
    Named("w") = w);
}


arma::vec mat_diag(arma::mat A, arma::mat B, arma::mat C){
  int n = A.n_rows;
  int p = A.n_cols;
  arma::vec d(n, arma::fill::zeros);
  arma::mat D = A*B;
  for(int i=0; i<n; ++i){
    d(i) = sum(D(i, arma::span(0,p-1)).t()%C(arma::span(0,p-1),i));
  }
  return(d);
}

arma::mat gen_mvnrnd(arma::mat M, arma::mat Covs){
  int p = M.n_rows;
  int q = M.n_cols;
  arma::mat X(p,q);
  for(int j=0; j<q; ++j)
  {
    arma::vec m = M(arma::span(0,p-1), j);
    arma::mat S = reshape(Covs(arma::span(0, p*p - 1), j), p, p);
    X(arma::span(0,p-1),j) = arma::mvnrnd(m, S);
  }
  return X;
}

arma::mat vec_log_post_beta_d1_d2(arma::vec y, arma::mat X, arma::vec beta, arma::vec mean_b, arma::mat inv_covb){
  int n = y.size();
  int q = X.n_cols;
  arma::vec qr = 2*y - arma::ones(n);
  arma::vec lambda(n, arma::fill::zeros);
  for(int i=0; i<n; ++i){
    lambda(i) = qr(i)*arma::normpdf(qr(i)*sum(X(i,arma::span(0,q-1))*beta))/arma::normcdf(qr(i)*sum(X(i,arma::span(0,q-1))*beta));
  }
  arma::vec del_log_posterior(q, arma::fill::zeros);
  for(int j=0; j<q; ++j){
    del_log_posterior(j) = sum(lambda.t()*X(arma::span(0,n-1), j));
  }
  arma::vec score = del_log_posterior - inv_covb*(beta - mean_b);
  
  arma::mat hessian(q, q, arma::fill::zeros);
  for(int i=0; i<n; ++i){
    hessian = hessian + lambda(i)*((sum(X(i,arma::span(0,q-1))*beta)) + lambda(i))*X(i, arma::span(0, q-1)).t()*X(i, arma::span(0,q-1));
  }
  hessian = hessian + inv_covb;
  arma::mat result(q, q+1);
  result(arma::span(0, q-1), 0) = score;
  result(arma::span(0, q -1), arma::span(1, q)) = -hessian;
  return result;
}



arma::vec vec_log_post_beta_laplace(arma::vec y, arma::mat X, arma::vec mean_b, arma::mat inv_covb, int max_it, double epsilon){
  
  int n = y.size();
  int q = X.n_cols;
  arma::vec b(q, arma::fill::zeros);
  arma::mat H(q, q, arma::fill::zeros);
  for(int i = 1; i<=max_it; ++i){
    arma::mat res = vec_log_post_beta_d1_d2(y, X, b, mean_b, inv_covb);
    arma::vec u = res(arma::span(0, q-1), 0);
    H = res(arma::span(0, q-1), arma::span(1, q));
    arma::vec err = -inv(H)*u;
    double norm_err = norm(err, 2);
    if(norm_err<epsilon){
      break;
    }
    b = b - inv(H)*u;
  }

  double vec_size = q + q*q + 3*n;
  arma::vec result(vec_size);
  result(arma::span(0, q-1)) = b;
  arma::mat I = arma::diagmat(arma::ones(q));
  arma::mat H_inv = solve(-H, I);
  arma::vec cov_vec = arma::vectorise(H_inv);
  result(arma::span(q, q-1 + pow(q, 2))) = cov_vec;
  arma::vec s = mat_diag(X, H_inv, X.t());
  result(arma::span(q-1 + pow(q, 2) + 1, q-1 + pow(q, 2) + n)) = arma::ones(n) + s;
  arma::vec qr = 2*y - arma::ones(n);
  arma::vec m = X*b/pow(arma::ones(n) + s, 0.5);
  result(arma::span(q-1 + pow(q, 2)+ n + 1, q-1 + pow(q, 2) + 2*n)) = qr;
  result(arma::span(q-1 + pow(q, 2)+ 2*n + 1, q-1 + pow(q, 2) + 3*n)) = qr%m;
  return result;
}



arma::mat marginal_probit(arma::mat Y, arma::mat X, double alpha_hat, double prior_var, double epsilon, int max_it){
  int n = Y.n_rows;
  int p = Y.n_cols;
  int q = X.n_cols;
  arma::vec mean_b(q, arma::fill::zeros); 
  // compute tau_p and mu_p
  double taup = 2*log(p);
  double mup = sqrt(1+taup)*R::qnorm(alpha_hat / (alpha_hat + p), 0.0, 1.0, 1, 0);
  mean_b(0) = mup;
  arma::vec var_b(q);
  var_b(0) = taup;
  var_b.subvec(1, q-1).fill(prior_var);
  arma::mat inv_covb = diagmat(1/var_b);

  double vec_size = q + q*q + 3*n;
  arma::mat all_res(vec_size, p);
  for(int j=0; j<p; ++j){
    all_res(arma::span(0, vec_size - 1), j) = vec_log_post_beta_laplace(Y(arma::span(0,n-1), j), X, mean_b, inv_covb, max_it, epsilon);   
  }
  return all_res;
}


double cox_bvn0(double a, double b, double rho){
  double l1;
  if(rho>0){
    double l = 0.0;
    if(a>=0 || b>=0){
      double c = std::max(a,b);
      double d = std::min(a,b);
      double phi_c = arma::normcdf(-c);
      double mu_c = arma::normpdf(c)/phi_c;
      double xi = (rho*mu_c - d)/pow(1 - rho*rho, 0.5);
      l = phi_c*arma::normcdf(xi) - 1 + arma::normcdf(c) + arma::normcdf(d);
    } else if(a<0 && b<0){
      double c = std::max(-a,-b);
      double d = std::min(-a,-b);
      double phi_c = arma::normcdf(-c);
      double mu_c = arma::normpdf(c)/phi_c;
      double xi = (rho*mu_c - d)/pow(1 - rho*rho, 0.5);
      l = phi_c*arma::normcdf(xi);
    }
    l1 = l;
  } else {
    double l;
    l = arma::normcdf(a) - cox_bvn0(a, -b, -rho);
    l1 = l;
  }
  return l1;
}


arma::vec cox_bvn(arma::vec a, arma::vec b, arma::vec rhos){
  int n = rhos.size();
  arma::vec l(n, arma::fill::zeros);
  for(int i=0; i<n; ++i){
    l(i) = log(cox_bvn0(a(i), b(i), rhos(i)));
  }
  return l;
}

double post_cor_new_h(double rho, arma::vec q1, arma::vec q2, arma::vec m1, arma::vec m2, arma::vec s1, arma::vec s2, double alpha1, double alpha2_sq){
  int n = m1.size();
  arma::vec rho_all = q1(arma::span(0,n-1))%q2(arma::span(0, n-1))%(rho*arma::ones(n))/(pow(s1, 0.5)%pow(s2, 0.5));
  arma::vec l(n);
  l = cox_bvn(m1(arma::span(0,n-1)), m2(arma::span(0,n-1)), rho_all);   
  double prior_prob = -log(1-pow(rho, 2)) + R::dnorm(0.5*log(1+rho) - 0.5*log(1-rho), alpha1, pow(alpha2_sq, 0.5), true);
  return sum(l) + prior_prob;
}



List marginal_pairwise_probit_h(arma::mat params, int m, arma::mat Y, arma::mat X, double alpha1, double alpha2_sq){
  arma::wall_clock timer;
  timer.tic();
  int n = Y.n_rows;
  int p = Y.n_cols;
  int q = X.n_cols;
  List GL_res = GLRcpp(m, -1, 1);
  arma::vec rhos = GL_res[0];
  arma::vec wts = GL_res[1];
  arma::mat beta = params(arma::span(0, q-1), arma::span(0, p-1));
  arma::mat S = params(arma::span(q-1 + pow(q, 2) + 1, q-1 + pow(q, 2) + n), arma::span(0, p-1));
  arma::mat Q = params(arma::span(q-1 + pow(q, 2)+ n + 1, q-1 + pow(q, 2) + 2*n), arma::span(0, p-1));
  arma::mat M = params(arma::span(q-1 + pow(q, 2)+ 2*n + 1, q-1 + pow(q, 2) + 3*n), arma::span(0, p-1));
  arma::mat corr_mean = arma::mat(p, p);
  arma::mat corr_var = arma::mat(p, p);
  corr_mean.diag() = arma::ones(p);
  corr_var.diag() = arma::zeros(p);
  for(int j=0; j<(p-1); ++j){
    for(int k=(j+1); k<p; ++k){
      arma::vec d(m);
      for(int l=0; l<m; ++l){
        d(l) = post_cor_new_h(rhos(l), Q(arma::span(0,n-1),j), Q(arma::span(0,n-1),k), M(arma::span(0,n-1),j), M(arma::span(0,n-1),k), S(arma::span(0,n-1),j), S(arma::span(0,n-1),k), alpha1, alpha2_sq);
      }
      arma::uvec id = find_finite(d);
      double m_const = median(d(id));
      arma::vec d1 = exp(-m_const*arma::ones(id.size()) + d(id));
      corr_mean(j,k) = corr_mean(k,j) = sum(wts(id)%d1%rhos(id))/sum(wts(id)%d1);
      corr_var(j,k) = corr_var(k,j) = sum(wts(id)%d1%pow(rhos(id), 2))/sum(wts(id)%d1) - pow(corr_mean(j,k),2);
    }
  }
  List result;
  double runtime = timer.toc();
  result["post_mean"] = corr_mean;
  result["post_var"] = corr_var;
  result["runtime"] = runtime;
  return result;
}


// log lik for alpha for MH
double log_lik_alpha(double alpha, arma::vec beta){
  int p = beta.size();
  double tau2 = 2 * log(p);
  double mu = sqrt(1 + tau2) * R::qnorm(alpha / (alpha + p), 0.0, 1.0, 1, 0);
  if (alpha > 0){
    arma::vec logl = log_normpdf(beta, mu, sqrt(tau2));
   return sum(logl);
  }
  else{
    return R_NegInf;
  }
  
}

// Log-posterior for MH
double logpost_alpha(double alpha, arma::vec beta, double a_alpha, double b_alpha) {
  double log_likelihood = log_lik_alpha(alpha, beta);
  return log_likelihood + R::dgamma(exp(alpha), a_alpha, b_alpha, true) + alpha;
}



arma::vec one_stage_sampling(arma::mat Y, arma::mat X, double a_alpha, double b_alpha, double prior_var, int max_it, double epsilon, double nmcmc, double burnin, int m, double eps_MH, int truncP){
  int n = Y.n_rows;
  Y = arma::join_rows(Y, arma::zeros<arma::mat>(n, truncP));
  int p = Y.n_cols;
  int q = X.n_cols;

  arma::vec alpha_beta_samples((nmcmc - burnin));
  double alpha_beta = arma::randg(1, distr_param(a_alpha,b_alpha))(0);
  for(int ii =1; ii<=nmcmc; ++ii){
    // estimate of beta
    arma::mat beta_res = marginal_probit(Y, X, alpha_beta, prior_var, epsilon, max_it);
    arma::mat beta_mean = beta_res(arma::span(0, q-1), arma::span(0, p-1));    // Draw samples of beta
        arma::mat cov_mats = beta_res(arma::span(q, q-1 + pow(q, 2)), arma::span(0, p-1));
    arma::mat b_samples = gen_mvnrnd(beta_mean, cov_mats);
    arma::vec beta0 = b_samples.row(0).t();
    // MH step
    double logp = logpost_alpha(alpha_beta, beta0, a_alpha, b_alpha);
    double alpha_beta_new = arma::randn() * eps_MH + alpha_beta;
    double logp_new = logpost_alpha(alpha_beta_new, beta0, a_alpha, b_alpha);
    double a_acc = std::min(1.0, exp(logp_new - logp));
    if (R::runif(0.0, 1.0) < a_acc){
      logp = logp_new;
      alpha_beta = alpha_beta_new;
    }
    if(ii> burnin)
    {
      alpha_beta_samples(ii - burnin -1) = alpha_beta;
    }
    
    if(ii%100 == 0)
    {
      Rprintf("Iteration = %d \n", ii);
    }
  }
  return alpha_beta_samples;
}


arma::vec two_stage_sampling(arma::mat Y, arma::mat X, double eta0_rho, double nu0_rho, double gamma0_rho, double lambda0_rho, double a_alpha, double b_alpha, double prior_var, int max_it, double epsilon, double nmcmc, double burnin, int m, double eps_MH){
  int p = Y.n_cols;
  int q = X.n_cols;

  arma::vec eta_rho_samples((nmcmc - burnin));
  arma::vec omega_sq_rho_samples((nmcmc - burnin));
  arma::vec alpha_beta_samples((nmcmc - burnin));
  double eta_rho = 0;
  double omega_sq_rho = 1;
  // MH initializiation 
  double alpha_beta = arma::randg(1, distr_param(a_alpha,b_alpha))(0);
  for(int ii =1; ii<=nmcmc; ++ii){
    // estimate of beta
    arma::mat beta_res = marginal_probit(Y, X, alpha_beta, prior_var, epsilon, max_it);
    arma::mat beta_mean = beta_res(arma::span(0, q-1), arma::span(0, p-1));
    arma::mat cov_mats = beta_res(arma::span(q, q-1 + pow(q, 2)), arma::span(0, p-1));
    // Draw samples of beta
    arma::mat b_samples = gen_mvnrnd(beta_mean, cov_mats);
    arma::vec beta0 = b_samples.row(0).t();
    // MH step
    double logp = logpost_alpha(alpha_beta, beta0, a_alpha, b_alpha);
    double alpha_beta_new = arma::randn() * eps_MH + alpha_beta;
    double logp_new = logpost_alpha(alpha_beta_new, beta0, a_alpha, b_alpha);
    double a_acc = std::min(1.0, exp(logp_new - logp));
    if (R::runif(0.0, 1.0) < a_acc){
      logp = logp_new;
      alpha_beta = alpha_beta_new;
    }
    // estimate of Sigma
    List rho_res = marginal_pairwise_probit_h(beta_res, m, Y, X, eta_rho, omega_sq_rho);
    arma::mat rho_mean = rho_res[0];
    arma::mat rho_var = rho_res[1];
    arma::uvec alt_lower_indices = trimatl_ind(size(rho_mean), -1);
    arma::vec rho1 = rho_mean(alt_lower_indices);
    arma::vec rho2 = rho_var(alt_lower_indices);
    // Draw samples
    arma::vec rho_samp = rho1 + pow(rho2, 0.5)%randn(0.5*p*(p-1));
    arma::vec gamma_samp = 0.5*log(rho_samp+arma::ones(0.5*p*(p-1))) - 0.5*log(arma::ones(0.5*p*(p-1))-rho_samp);
    arma::uvec id = arma::find_finite(gamma_samp);
    double gamma_mean = arma::mean(gamma_samp(id));
    int p_l = id.n_elem;
  
    // Posterior quantities
    double nu_rho_q = nu0_rho + 0.5*p_l*(p_l-1);
    double gamma_rho_q = gamma0_rho + 0.25*p_l*(p_l-1);
    double gamma_sq_sum = (0.5*p_l*(p_l-1) - 1)*arma::var(gamma_samp(id));
    double gamma_post_var = lambda0_rho + 0.5*gamma_sq_sum + (nu0_rho*0.5*p_l*(p_l-1)/(0.5*nu0_rho+ 0.5*p_l*(p_l-1)))*pow(eta0_rho - gamma_mean,2);
    double omega_sq_rho = 1/arma::randg<double>(distr_param(gamma_rho_q, 1/gamma_post_var));
    double eta_rho = (nu0_rho*eta0_rho + 0.5*p_l*(p_l-1)*gamma_mean)/(nu0_rho+ 0.5*p_l*(p_l-1)) + (pow(omega_sq_rho, 0.5)/nu_rho_q)*arma::randn();
    if(ii> burnin)
    {
      eta_rho_samples(ii - burnin -1) = eta_rho;
      omega_sq_rho_samples(ii - burnin - 1) = omega_sq_rho;
      alpha_beta_samples(ii - burnin -1) = alpha_beta;
    }
    
    if(ii%20 == 0)
    {
      Rprintf("Iteration = %d \n", ii);
    }
  }
  double alpha_hat = mean(alpha_beta_samples);
  double eta_rho_hat = mean(eta_rho_samples);
  double omega_sq_rho_hat = mean(omega_sq_rho_samples);
  arma::vec par = {alpha_hat, eta_rho_hat, omega_sq_rho_hat};
  return par;
}

// [[Rcpp::export]]

List TRACEh(arma::mat Y, arma::mat X, double eta0_rho, double nu0_rho, double gamma0_rho, double lambda0_sq_rho, double a_alpha, double b_alpha, double prior_var, int max_it, double epsilon, int m, int nmcmc, int burnin, double eps_MH, int truncP){
  arma::mat data = Y;
  int n = Y.n_rows;
  Y = arma::join_rows(Y, arma::zeros<arma::mat>(n, truncP));
  int p = Y.n_cols;
  int q = X.n_cols;
  arma::wall_clock timer;
  timer.tic();
  arma::vec emp_bayes_estimates = two_stage_sampling(Y, X, eta0_rho, nu0_rho, gamma0_rho, lambda0_sq_rho, a_alpha, b_alpha, prior_var, max_it, epsilon, nmcmc, burnin, m, eps_MH);
  double alpha_hat = exp(emp_bayes_estimates(0));
  double alpha1 = emp_bayes_estimates(1);
  double alpha2_sq = emp_bayes_estimates(2);
  arma::mat first_stage = marginal_probit(Y, X, alpha_hat, prior_var, epsilon, max_it);
  List second_stage = marginal_pairwise_probit_h(first_stage, m, Y, X, alpha1, alpha2_sq);
  List result;
  result["coefficients"] = first_stage(arma::span(0, q-1), arma::span(0, p-1));
  result["cov_mats"] = first_stage(arma::span(q, q-1 + pow(q, 2)), arma::span(0, p-1));
  
  result["post_mean"] = second_stage[0];
  result["post_var"] = second_stage[1];
  double runtime = timer.toc();
  result["runtime"] = runtime;
   result["Y"] = data;
  return result;
}

