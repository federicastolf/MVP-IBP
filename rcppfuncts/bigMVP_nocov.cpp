#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;



// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


double bvn0(double h, double k, double r) {
  NumericVector w;
  NumericVector x;
  if (abs(r)<0.3){ 
    w = {0.171324492379170,   0.360761573048138,   0.467913934572690,   0.171324492379170,   0.360761573048138,   0.467913934572690};
    x = {0.067530485796848,   0.338790613533735,   0.761380813916803,   1.932469514203152,   1.661209386466265,   1.238619186083197};
  }
  else if (abs(r)<0.75){
    w = {0.047175336386512,   0.106939325995318,   0.160078328543346,   0.203167426723066,   0.233492536538355,   0.249147045813403,   0.047175336386512, 0.106939325995318,   0.160078328543346,  0.203167426723066,   0.233492536538355,  0.249147045813403};
    x = {0.018439365753281,   0.095882743629525,   0.230097325805695,   0.412682045713383,   0.632168501001820,   0.874766591488531,   1.981560634246719, 1.904117256370475,   1.769902674194305,   1.587317954286617,   1.367831498998180,   1.125233408511469};
  }
  else{
    w = {0.017614007139152,   0.040601429800387,   0.062672048334109,   0.083276741576705,   0.101930119817240,   0.118194531961518,   0.131688638449177, 0.142096109318382,   0.149172986472604,   0.152753387130726,   0.017614007139152,   0.040601429800387,   0.062672048334109,   0.083276741576705, 0.101930119817240,   0.118194531961518,   0.131688638449177,   0.142096109318382,   0.149172986472604,   0.152753387130726};
    x = {0.006871400814905,   0.036028072722086,   0.087765571748674,   0.160883028177781,   0.253668093539849,   0.363946319273485,   0.489132998049173, 0.626293911284580,   0.772214148858355,   0.923473478866503,   1.993128599185095,   1.963971927277914,   1.912234428251326,   1.839116971822219, 1.746331906460151,   1.636053680726515,   1.510867001950827,   1.373706088715420,   1.227785851141645,   1.076526521133497};
  }
  double tp = 2 * M_PI ;
  double hk = h * k;
  double hs = (h*h + k*k)/2;
  double asr = asin(r)/2;
  NumericVector sn = sin(asr*x);
  NumericVector bvn = exp((sn*hk-hs)/(1-pow(sn,2)));
  double prob = std::inner_product(bvn.begin(), bvn.end(), w.begin(), 0.0);
  prob = prob * asr / tp + 0.25 * erfc(h/sqrt(2)) * erfc(k/sqrt(2));
  return prob;
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

arma::vec bvn(arma::vec a, arma::vec b, arma::vec rhos){
  int n = rhos.size();
  arma::vec l(n, arma::fill::zeros);
  for(int i=0; i<n; ++i){
    l(i) = log(bvn0(a(i), b(i), rhos(i)));
  }
  return l;
}

arma::vec cox_bvn(arma::vec a, arma::vec b, arma::vec rhos){
  int n = rhos.size();
  arma::vec l(n, arma::fill::zeros);
  for(int i=0; i<n; ++i){
    l(i) = log(cox_bvn0(a(i), b(i), rhos(i)));
  }
  return l;
}



arma::vec vec_log_post_beta_d1_d2_h(arma::vec y, double beta, double prior_var){
  int n = y.size();
  arma::vec q = 2*y - arma::ones(n);
  arma::vec lambda(n, arma::fill::zeros);
  for(int i=0; i<n; ++i){
    lambda(i) = q(i)*arma::normpdf(q(i)*beta)/arma::normcdf(q(i)*beta);
  }
  double del_log_posterior = sum(lambda); //derivative log posterior
  double score = del_log_posterior - beta/prior_var; // first derivative, try!
  

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





arma::vec vec_log_post_beta_laplace_h(arma::vec y, double prior_var, int max_it, double epsilon){
  
  // change so that it returns a vector
  int n = y.size();
  double b = 0; // you could inizialize with mu_p o 0?
  double H = 0;
  for(int i = 1; i<=max_it; ++i){
    arma::vec res = vec_log_post_beta_d1_d2_h(y, b, prior_var);
    double u = res(0);
    H = res(1);
    double err = -u/H; // update of beta in NR
    if(err<epsilon){
      break;
    }
    b = b - u/H;
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
  return result;
}



arma::mat marginal_probit_h(arma::mat Y, double prior_var, double epsilon, int max_it){
  int n = Y.n_rows;
  int p = Y.n_cols;
  double vec_size = 2 + 3*n + 1;
  arma::mat all_res(vec_size, p);
  for(int j=0; j<p; ++j){
    all_res(arma::span(0, vec_size - 1), j) = vec_log_post_beta_laplace_h(Y(arma::span(0,n-1), j), prior_var, max_it, epsilon); 
  }
  return all_res;
}



double post_cor_new_h(double rho, arma::vec q1, arma::vec q2, arma::vec m1, arma::vec m2, arma::vec s1, arma::vec s2, double alpha1, double alpha2_sq){
  int n = m1.size();
  arma::vec rho_all = q1(arma::span(0,n-1))%q2(arma::span(0, n-1))%(rho*arma::ones(n))/(pow(s1, 0.5)%pow(s2, 0.5));
  arma::vec l(n);
  l = cox_bvn(m1(arma::span(0,n-1)), m2(arma::span(0,n-1)), rho_all);   //+ lkj_marginal(rho, nu_lkj - 1 + 0.5*q, nu_lkj - 1 + 0.5*q)*arma::ones(n);
  double prior_prob = -log(1-pow(rho, 2)) + R::dnorm(0.5*log(1+rho) - 0.5*log(1-rho), alpha1, pow(alpha2_sq, 0.5), true);
  return sum(l) + prior_prob;//+ lkj_marginal(rho, nu_lkj - 1 + 0.5*q, nu_lkj - 1 + 0.5*q);
}



List GLRcpp(int n, double a, double b){
  const arma::vec& i = arma::linspace(1, n-1, n-1);
  arma::mat J(n,n,arma::fill::zeros);
  const arma::vec& d = i / (arma::sqrt(4 * arma::square(i) - 1));
  J.diag(1) = d;
  J.diag(-1) = d;
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





List marginal_pairwise_probit_h(arma::mat params, int m, arma::mat Y, double alpha1, double alpha2_sq){
  arma::wall_clock timer;
  timer.tic();
  int n = Y.n_rows;
  int p = Y.n_cols;
  List GL_res = GLRcpp(m, -1, 1);
  arma::vec rhos = GL_res[0];
  arma::vec wts = GL_res[1];
  arma::mat beta = params(0, arma::span(0, p-1));
  //List cov_mats = params[1];
  arma::mat S = params(arma::span(2, 2 + n-1), arma::span(0, p-1));
  arma::mat Q = params(arma::span(2+ n, 2+ 2*n-1), arma::span(0, p-1));
  arma::mat M = params(arma::span(2+ 2*n , 2 + 3*n-1), arma::span(0, p-1));
  arma::mat corr_mean = arma::mat(p, p);
  arma::mat corr_var = arma::mat(p, p);
  corr_mean.diag() = arma::ones(p);
  corr_var.diag() = arma::zeros(p);
  for(int j=0; j<(p-1); ++j){
    for(int k=(j+1); k<p; ++k){
      arma::vec d(m);
      for(int l=0; l<m; ++l){
        d(l) = post_cor_new_h(rhos(l), Q(arma::span(0,n-1),j), Q(arma::span(0,n-1),k), M(arma::span(0,n-1),j), M(arma::span(0,n-1),k), S(arma::span(0,n-1),j), S(arma::span(0,n-1),k), alpha1, alpha2_sq);
        //Rprintf("check = %f", d(l));
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



arma::vec two_stage_sampling(arma::mat Y,  double eta0_rho, double nu0_rho, double gamma0_rho, double lambda0_rho, int max_it, double epsilon, double nmcmc, double burnin, int m, double prior_var){
  int p = Y.n_cols;
  arma::vec eta_rho_samples((nmcmc - burnin));
  arma::vec omega_sq_rho_samples((nmcmc - burnin));
  double eta_rho = 0;
  double omega_sq_rho = 1;

  for(int ii =1; ii<=nmcmc; ++ii){
    // estimate of beta
    arma::mat beta_res = marginal_probit_h(Y, prior_var, epsilon, max_it);
    // estimate of Sigma
    List rho_res = marginal_pairwise_probit_h(beta_res, m, Y, eta_rho, omega_sq_rho);
    arma::mat rho_mean = rho_res[0];
    arma::mat rho_var = rho_res[1];
    arma::uvec alt_lower_indices = trimatl_ind( size(rho_mean), -1);
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
    }
  }

  double eta_rho_hat = mean(eta_rho_samples);
  double omega_sq_rho_hat = mean(omega_sq_rho_samples);
  arma::vec par = {eta_rho_hat, omega_sq_rho_hat};
  return par;
}




// [[Rcpp::export]]

List bigMVPh(arma::mat Y, double eta0_rho, double nu0_rho, double gamma0_rho, double lambda0_sq_rho, int max_it, double epsilon, int m, int nmcmc, int burnin, double prior_var, int truncP){
  int n = Y.n_rows;
  Y = arma::join_rows(Y, arma::zeros<arma::mat>(n, truncP));
  int p = Y.n_cols;
  arma::wall_clock timer;
  timer.tic();
  arma::vec emp_bayes_estimates = two_stage_sampling(Y, eta0_rho, nu0_rho, gamma0_rho, lambda0_sq_rho, max_it, epsilon, nmcmc, burnin, m, prior_var);
  double alpha1 = emp_bayes_estimates(0);
  double alpha2_sq = emp_bayes_estimates(1);
  //Rprintf("done %d ", p);
  
  arma::mat first_stage = marginal_probit_h(Y, prior_var, epsilon, max_it); //same as non-hierarchical
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
