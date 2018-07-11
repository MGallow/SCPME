#ifndef ADMM_H
#define ADMM_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::List ADMMc(const arma::mat &S, const arma::mat &A, const arma::mat &B, const arma::mat &C, const arma::mat &initOmega, const arma::mat &initZ2, const arma::mat &initY, const double lam, const double tau = 10, double rho = 2, const double mu = 10, const double tau_inc = 2, const double tau_dec = 2, std::string crit = "ADMM", const double tol_abs = 1e-4, const double tol_rel = 1e-4, const int maxit = 1e4);

#endif