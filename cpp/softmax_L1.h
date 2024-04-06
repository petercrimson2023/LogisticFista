#ifndef SOFTMAX_L1_H
#define SOFTMAX_L1_H

#include <RcppArmadillo.h>

using namespace arma;

mat softmax(const arma::mat& X)
{
  mat expX = arma::exp(X-repmat(max(X,1),1,X.n_cols));
  mat expXsum = arma::sum(expX, 1);
  mat expXsum_rep = repmat(expXsum, 1, X.n_cols);
  mat P = expX / expXsum_rep;
  return P;
}


double f(const mat& Beta,const mat& X,const mat& y_one_hot, int n,int p, int k)
{
  mat beta_resize = reshape(Beta, p, k);
  mat score = X * beta_resize;
  mat P = softmax(score);
  double epsilon = 1e-10;
  mat P_safe = P + epsilon;
  double loss = -arma::accu(y_one_hot % arma::log(P_safe)) / n;
  return(loss);
}


mat gradf(const mat& Beta,const mat& X,const mat& y_one_hot, int n,int p, int k)
{
  mat beta_resize = reshape(Beta, p, k);
  mat score = X * beta_resize;
  mat P = softmax(score);
  mat grad_matrix = X.t() * (P - y_one_hot) / n;
  return(reshape(grad_matrix,p*k,1));
}


double g(const mat& Beta, double lambda)
{
  double norm_Beta = norm(Beta, 1);
  double gvalue = lambda * norm_Beta;
  return(gvalue);
}

mat gradg(const mat& Beta,double tau, double lambda)
{
  mat result =  abs(Beta)-tau*lambda;
  result = max(result,zeros<mat>(result.n_rows, result.n_cols));
  mat gradg_value = result % sign(Beta);
  return(gradg_value);
}

mat p_y(double& lambda,double& L,mat& theta,const mat& X,const mat& y_one_hot,int& n,int& p,int& k)
{
  mat u = theta - 1/L * gradf(theta,X,y_one_hot,n,p,k);
  mat v = gradg(u,1/L,lambda);
  return(v);
}

double Q(mat& theta1,mat& theta2,const mat& X,const mat& y_one_hot,double lambda, double& L,int& n,int& p,int& k)
{
  double arg = f(theta2,X,y_one_hot,n,p,k)+dot((theta1-theta2).t(),gradf(theta2,X,y_one_hot,n,p,k)) + L/2 * dot((theta1-theta2).t(),theta1-theta2)+g(theta1,lambda);
  return(arg);
}


#endif // SOFTMAX_L1_H
