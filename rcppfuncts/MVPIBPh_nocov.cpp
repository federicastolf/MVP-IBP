#include "RcppArmadillo.h"
#include "Rcpp.h"
#include "RcppEigen.h"

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]


using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
    Eigen::MatrixXd C = A * B;

    return Rcpp::wrap(C);
}


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
  // Only need first row...
  arma::vec w = 2 * arma::vectorise(arma::square(V.row(0)));
  arma::vec x = .5 * ((b - a) * L + a + b);
  w = -.5 * (a - b) * w;

  return List::create(
    Named("x") = x,
    Named("w") = w);
  
}


// it calculates the CDF of a bivariate normal distribution
// with correlation coefficient rho for given values a and b.

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
  l = cox_bvn(m1(arma::span(0,n-1)), m2(arma::span(0,n-1)), rho_all);   //+ lkj_marginal(rho, nu_lkj - 1 + 0.5*q, nu_lkj - 1 + 0.5*q)*arma::ones(n);
  double prior_prob = -log(1-pow(rho, 2)) + R::dnorm(0.5*log(1+rho) - 0.5*log(1-rho), alpha1, pow(alpha2_sq, 0.5), true);
  return sum(l) + prior_prob;//+ lkj_marginal(rho, nu_lkj - 1 + 0.5*q, nu_lkj - 1 + 0.5*q);
}



List marginal_pairwise_probit_h(arma::mat params, int m, arma::mat Y, double alpha1, double alpha2_sq){
  arma::wall_clock timer;
  timer.tic();
  int n = Y.n_rows;
  int q = Y.n_cols;
  List GL_res = GLRcpp(m, -1, 1);
  arma::vec rhos = GL_res[0];
  arma::vec wts = GL_res[1];
  arma::mat beta = params(0, arma::span(0, q-1)); //beta_hat (1xq)
  //List cov_mats = params[1];
  arma::mat S = params(arma::span(2, 2 + n-1), arma::span(0, q-1));
  arma::mat Q = params(arma::span(2+ n, 2+ 2*n-1), arma::span(0, q-1));
  arma::mat M = params(arma::span(2+ 2*n , 2 + 3*n-1), arma::span(0, q-1));
  arma::mat corr_mean = arma::mat(q, q);
  arma::mat corr_var = arma::mat(q, q);
  corr_mean.diag() = arma::ones(q);
  corr_var.diag() = arma::zeros(q);
  for(int j=0; j<(q-1); ++j){
    for(int k=(j+1); k<q; ++k){
      arma::vec d(m);
      for(int l=0; l<m; ++l){
        d(l) = post_cor_new_h(rhos(l), Q(arma::span(0,n-1),j), Q(arma::span(0,n-1),k), M(arma::span(0,n-1),j), M(arma::span(0,n-1),k), S(arma::span(0,n-1),j), S(arma::span(0,n-1),k), alpha1, alpha2_sq);
      }
      arma::uvec id = arma::find_finite(d);
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


arma::vec vec_log_post_beta_d1_d2_h(arma::vec y, double beta, double prior_var, double prior_mean){
  int n = y.size();
  arma::vec q = 2*y - arma::ones(n); //1xn
  arma::vec lambda(n, arma::fill::zeros);
  for(int i=0; i<n; ++i){
    lambda(i) = q(i)*arma::normpdf(q(i)*beta)/arma::normcdf(q(i)*beta);
  }
  double del_log_posterior = sum(lambda); //derivative log posterior
  double score = del_log_posterior - (beta - prior_mean)/prior_var; // first derivative, try!
  
  double hessian=0;
  for(int i=0; i<n; ++i){
    hessian = hessian + lambda(i)*(beta + lambda(i)); 
  }
  hessian = hessian + 1/prior_var;
  arma::vec result(2);
  result(0) = score;
  result(1) = -hessian;
  return result;
}


arma::vec vec_log_post_beta_laplace_h(arma::vec y, double alpha, int p, int max_it, double epsilon){
  
  //new
  int n = y.size();
  double prior_var = 2*log(p);
  double prior_mean =  sqrt(1 + prior_var) * R::qnorm(alpha / (alpha + p), 0.0, 1.0, 1, 0);
  double b = prior_mean; // you could inizialize with mu_p o 0?
  double H = 0;
  int it_n = 0;
  for(int i = 1; i<=max_it; ++i){
    it_n=i;
    arma::vec res = vec_log_post_beta_d1_d2_h(y, b, prior_var, prior_mean);
    double u = res(0);
    H = res(1);
    double err = -u/H; // update of beta in NR
    if(err<epsilon){
      break;
      //Rprintf("Iteration = %d \n", i);
    }
    b = b - u/H;
    //if(i== max_it)
    //{
      //Rprintf("Iteration max");
    //}
  }
  double vec_size = 2 + 3*n +1;
  arma::vec result(vec_size);
  result(0) = b;
  double H_inv = -1/H; 
  result(1) = H_inv;
  result(arma::span(2, 2 + n-1)) = arma::ones(n) + arma::ones(n)*H_inv; //check
  arma::vec q = 2*y - arma::ones(n);
  arma::vec m = b/pow(arma::ones(n) + arma::ones(n)*H_inv, 0.5);
  result(arma::span(2+ n, 2+ 2*n - 1)) = q;
  result(arma::span(2 + 2*n, 2+ 3*n - 1)) = q%m; //check
  result(2+ 3*n)=it_n;
  return result;
  // puoi ritornarti il num di it
}

arma::mat marginal_probit(arma::mat Y, double alpha, double epsilon, int max_it){
  int n = Y.n_rows;
  int p = Y.n_cols;
  double vec_size = 2 + 3*n +1; // +3n is elements useful for second stage
  arma::mat all_res(vec_size, p);

  for(int j=0; j<p; ++j){
    all_res(arma::span(0, vec_size - 1), j) = vec_log_post_beta_laplace_h(Y(arma::span(0,n-1), j), alpha, p, max_it, epsilon); 
  }
  return all_res;
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


arma::vec two_stage_sampling(arma::mat Y, double eta0_rho, double nu0_rho, double gamma0_rho, double lambda0_rho, double a_alpha, double b_alpha, int max_it, double epsilon, double nmcmc, double burnin, int m, double eps_MH){
  int p = Y.n_cols; 
  arma::vec eta_rho_samples((nmcmc - burnin)); //don't understand
  arma::vec omega_sq_rho_samples((nmcmc - burnin)); // you need that!
  arma::vec alpha_beta_samples((nmcmc - burnin));
  double eta_rho = 0;
  double omega_sq_rho = 1;
  // MH initializiation
  double alpha_beta = arma::randg(1)(0);

 for(int ii =1; ii<=nmcmc; ++ii){
    
    //first stage
    arma::mat beta_res = marginal_probit(Y, alpha_beta, epsilon, max_it);
    arma::vec beta_mean = beta_res.row(0).t();
    arma::vec var_beta = beta_res.row(1).t();
    arma::vec epsilon = arma::randn<arma::vec>(p);
    arma::vec b_samples = beta_mean + var_beta % epsilon;
    // from b_samples sample alpha trhough a metropolis

    // MH step
    double logp = logpost_alpha(alpha_beta, b_samples, a_alpha, b_alpha);
    double alpha_beta_new = arma::randn() * eps_MH + alpha_beta;
    double logp_new = logpost_alpha(alpha_beta_new, b_samples, a_alpha, b_alpha);
    double a_acc = std::min(1.0, exp(logp_new - logp));
    if (R::runif(0.0, 1.0) < a_acc){
      logp = logp_new;
      alpha_beta = alpha_beta_new;
    }
    // second stage
    List rho_res = marginal_pairwise_probit_h(beta_res, m, Y, eta_rho, omega_sq_rho);
    // draw sigma_jk (6)
    arma::mat rho_mean = rho_res[0]; //rho_jk mean
    arma::mat rho_var = rho_res[1]; //rho_jk variance
    arma::uvec alt_lower_indices = trimatl_ind( size(rho_mean), -1); //take the lower triangular
    arma::vec rho1 = rho_mean(alt_lower_indices);
    arma::vec rho2 = rho_var(alt_lower_indices);
    arma::vec rho_samp = rho1 + pow(rho2, 0.5)%randn(0.5*p*(p-1));
    // set gamma_jk
    arma::vec gamma_samp = 0.5*log(rho_samp+arma::ones(0.5*p*(p-1))) - 0.5*log(arma::ones(0.5*p*(p-1))-rho_samp);
    arma::uvec id = arma::find_finite(gamma_samp);
    double gamma_mean = arma::mean(gamma_samp(id));
    int q_l = id.n_elem;

    // sample omega^2 from a in inverse gamma (7)
    double nu_rho_q = nu0_rho + 0.5*q_l*(q_l-1);
    double gamma_rho_q = gamma0_rho + 0.25*q_l*(q_l-1); // WHY 0.25?
    double gamma_sq_sum = (0.5*q_l*(q_l-1) - 1)*arma::var(gamma_samp(id));
    // I don't understand this parameter
    double gamma_post_var = lambda0_rho + 0.5*gamma_sq_sum + (nu0_rho*0.5*q_l*(q_l-1)/(0.5*nu0_rho+ 0.5*q_l*(q_l-1)))*pow(eta0_rho - gamma_mean,2);
    double omega_sq_rho = 1/arma::randg<double>(distr_param(gamma_rho_q, 1/gamma_post_var));
    
    // WHAT's that?
    double eta_rho = (nu0_rho*eta0_rho + 0.5*q_l*(q_l-1)*gamma_mean)/(nu0_rho+ 0.5*q_l*(q_l-1)) + (pow(omega_sq_rho, 0.5)/nu_rho_q)*arma::randn();
    //Rprintf("eta_rho = %f \n", eta_rho);
    if(ii> burnin)
    {
      eta_rho_samples(ii - burnin -1) = eta_rho;
      omega_sq_rho_samples(ii - burnin - 1) = omega_sq_rho;
      alpha_beta_samples(ii - burnin -1) = alpha_beta;
    }
    
    //if(ii%20 == 0)
    //{
    //  Rprintf("Iteration = %d \n", ii);
    //}

 }
  double alpha_hat = mean(alpha_beta_samples);
  double eta_rho_hat = mean(eta_rho_samples);
  double omega_sq_rho_hat = mean(omega_sq_rho_samples);
  arma::vec par = {alpha_hat, eta_rho_hat, omega_sq_rho_hat};
  return par;
}


// [[Rcpp::export]]

List MVP_IBPh(arma::mat Y, double eta0_rho, double nu0_rho, double gamma0_rho, double lambda0_sq_rho, double a_alpha, double b_alpha, int max_it, double epsilon, int m, int nmcmc, int burnin, double eps_MH, int truncP){
  int n = Y.n_rows;
  Y = arma::join_rows(Y, arma::zeros<arma::mat>(n, truncP));
  int p = Y.n_cols;
  arma::wall_clock timer;
  timer.tic();

  arma::vec emp_bayes_estimates = two_stage_sampling(Y, eta0_rho, nu0_rho, gamma0_rho, lambda0_sq_rho, a_alpha, b_alpha, max_it, epsilon, nmcmc, burnin, m, eps_MH);
  double alpha_beta = exp(emp_bayes_estimates(0));
  double alpha1 = emp_bayes_estimates(1);
  double alpha2_sq = emp_bayes_estimates(2);

  arma::mat first_stage = marginal_probit(Y, alpha_beta, epsilon, max_it);
  List second_stage = marginal_pairwise_probit_h(first_stage, m, Y, alpha1, alpha2_sq);
  List result;
  result["coefficients"] = first_stage(0, arma::span(0, p-1));
  result["cov_mats"] = first_stage(1, arma::span(0, p-1));
  
  result["post_mean"] = second_stage[0];
  result["post_var"] = second_stage[1];
  double runtime = timer.toc();
  result["runtime"] = runtime;
  return result;
}