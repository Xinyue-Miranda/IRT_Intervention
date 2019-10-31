#include <RcppArmadillo.h>
#include <cmath>
#include <truncnorm.h>
#include <mvnorm.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec update_b(int N, int K, double sigma_sqr, mat theta, cube M, mat lambda, mat Z){
  vec b(K);
  for (int k = 0; k < K; k++) {
    double mean = sigma_sqr * sum(theta * M.slice(k) * lambda.row(k).t() - 
                                  Z.col(k)) / (sigma_sqr*N+1);
    double variance = std::sqrt(float(sigma_sqr/(sigma_sqr*N+1)));
    b[k] = R::rnorm(mean, variance);
  }
  return b;
}

// [[Rcpp::export]]
mat update_lambda(int K, mat theta, cube M, vec b, mat Sigma, mat Z, int d){
  mat lambda(K,d);
  mat I;
  for(int k = 0; k < K; k++){
    mat V_lambda  = inv_sympd( M.slice(k) * theta.t() * theta * M.slice(k) + inv(Sigma));
    mat mu_lambda = V_lambda *  M.slice(k) * theta.t() * (b(k) + Z.col(k));
    lambda.row(k) = rmvnorm(1, mu_lambda, V_lambda);
  }
  return lambda;
}

// [[Rcpp::export]]
mat update_lambda_star(int K, mat lambda, cube M, int d){
  mat lambda_star(K,d);
  for(int k = 0; k < K; k++) {
    lambda_star.row(k) = lambda.row(k) * M.slice(k);
  }
  return lambda_star;
}

// [[Rcpp::export]]
mat update_theta(int N, mat lambda_star, vec b, mat Z, cube inv_v, vec ind, 
                 mat theta_start, int d){
  mat theta(N,d);
  for(int n = 0; n < N; n++){
    mat V_theta = inv_sympd(lambda_star.t() * lambda_star + inv_v.slice(n));
    auto it = std::find (std::begin(ind), std::end(ind), n+1);
    if (it == std::end(ind)) {
      vec mu_theta = V_theta * lambda_star.t() * (b + trans(Z.row(n)));
      theta.row(n) = rmvnorm(1, mu_theta, V_theta);
    } else {
      theta.row(n) = theta_start.row(n);
    }
  }
  return theta;
}

// [[Rcpp::export]]
mat update_Z(int N, int K, mat lambda, cube M, mat theta, vec b, mat Y){
  mat Z(N,K);
  for (int i = 0; i < N; i++) {
    for (int k = 0; k < K; k++) {
      vec mu_z = lambda.row(k) * M.slice(k) * theta.row(i).t();
      if(Y(i,k) == 1){
        Z(i,k) = r_truncnorm(mu_z[0], 1, 0, 100000);
      }
      else{Z(i,k) = r_truncnorm(mu_z[0], 1, -100000, 0);}
    }
  }
  return Z;
}
