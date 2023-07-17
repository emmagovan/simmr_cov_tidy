#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

static double const log2pi = std::log(2.0 * M_PI);


void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

//dnorm 
// [[Rcpp::export]]
arma::vec dmvnrm_arma_fast(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = true) { 
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}

//calculate crossproduct of 2 matrices
// [[Rcpp::export]]
NumericMatrix crossprod(NumericMatrix X){
  NumericMatrix ans(X.nrow(), X.ncol());
  
  for(int i = 0; i<X.ncol(); i++){
    for(int j=0; j<X.ncol(); j++){
      for(int n =0; n<X.nrow(); n++){
        
        ans(i,j) += X(n,i) * X(n,j);
      }
    }}
  return(ans);
}

//matrix multiplication
// [[Rcpp::export]]
NumericMatrix matmult(NumericMatrix x, NumericMatrix y) {
  NumericMatrix ans(x.nrow(), y.ncol());
  
  for(int i = 0; i<x.nrow(); i++){
    for(int j=0; j<y.ncol(); j++){
      for(int k =0; k<y.nrow(); k++){
        
        ans(i,j) += x(i,k) * y(k,j);
      }
    }}
  
  return ans;
}

//rnorm sample
// [[Rcpp::export]]
NumericMatrix rMVNormCpp(const double n,
                         const arma::vec mu,
                         const NumericMatrix U) {
  
  
  // Dimension of MVN
  int p = mu.size();
  
  // Simulate iid standard normals
  arma::mat Z(p, n);
  Z.imbue(norm_rand);
  
  // Now backsolve and add back on the means
  arma::mat X = solve(as<arma::mat>(U), Z);
  for ( int i = 0; i < n; ++i ) {
    X.col(i) += mu;
  }
  
  return Rcpp::wrap(X.t());
}


//solves matrix
// [[Rcpp::export]]
NumericMatrix solvearma(const NumericMatrix X) {
  
  arma::mat b = arma::eye(X.nrow(), X.ncol());
  
  
  // Now backsolve and add back on the means
  arma::mat ans = solve(as<arma::mat>(X), b);
  
  
  return Rcpp::wrap(ans.t());
}

//simulates S samples of theta - so several rnorms + rgamma
//[[Rcpp::export]]
NumericMatrix sim_thetacpp(int S, NumericVector lambda, int n_sources,
                           int n_tracers, int n_cov){
  NumericMatrix theta(S, (n_cov*n_sources + n_tracers));
  NumericMatrix mean_beta((n_cov), n_sources);
  int mat_size = n_sources * (n_sources+1) /2;
  
  for(int i=0; i<n_cov; i++){
    for(int k=0; k<n_sources;k++){
      mean_beta(i,k) = lambda(i * mat_size + i * n_sources + k);
    }
  }
  
  NumericMatrix sig_beta(n_cov, mat_size);
  
  for(int m = 0; m<mat_size; m++){
    for(int i =0; i<n_cov; i++){
      sig_beta(i,m) = lambda(i* mat_size + (i+1) * n_sources + m);
      
    }
  }
  
  NumericVector count(n_cov);
  
  for(int i =0; i<n_cov; i++){
    count(i) = 0;
  }
  
  arma::cube chol_prec(n_sources, n_sources, n_cov);
  
  for(int j = 0; j< n_sources; j++){
    for(int i = 0; i<n_sources; i++){
      for(int m = 0; m<n_cov; m++){
        if (i <= j){
          count(m) +=1;
          chol_prec(i,j,m) = sig_beta(m, count(m)-1);
          
          
        }
        
        else{
          chol_prec(i,j,m) = 0;
        }
      }
    }
  }
  
  
  
  
  arma::mat theta_arma(S, (n_cov*n_sources + n_tracers)); //Want to go from chol_prec array to
  // A matrix of thetas generated using rMVNormCpp
  
  for(int i=0; i<n_cov; i++){
    theta_arma.submat(0, (i)*n_sources, S-1, (i+1)*n_sources - 1) = as<arma::mat>(rMVNormCpp(S, mean_beta(i,_), Rcpp::wrap(chol_prec.slice(i))));
  }
  
  theta = Rcpp::wrap(theta_arma);
  
  
  
  for(int i = 0; i<n_tracers; i++){
    theta(_,i+n_sources*n_cov) = (Rcpp::rgamma(S,  lambda(n_cov * mat_size + n_cov *n_sources +i),
                                  lambda(n_cov * mat_size + n_cov *n_sources +i + n_tracers)));
  }
  
  
  return theta;
}