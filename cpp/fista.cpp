#include <RcppArmadillo.h>
#include "softmax_L1.h"

using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
Rcpp::List fista(double lambda, double L_init, const arma::mat& theta0, const arma::mat& X, const arma::mat& y_one_hot,
                 int max_iter=10000, double eps = 1e-6, double eita = 1.2, bool loss_compute =true,
                 int n=1, int p=1, int k=1)
{
  // Initialization
  double L_old = L_init;
  mat gama = theta0;
  mat p_l_gama = theta0;
  mat theta_old = theta0;
  mat theta_new = theta0;
  double t_old = 1;
  std::vector<double> loss_list = {100};
  double L_bar = 0;
  double fvalue = 0;
  double qvalue = -1;
  bool smallest_ik_condition = true;
  double L_new = 0;
  double tk = 0;

  // Iteration
  bool condition = true;
  int times = 1;
  int ik = 1;



  while (condition) {
    ik = 1;
    smallest_ik_condition = true;

    while (smallest_ik_condition) {
      L_bar = pow(eita, ik) * L_old;

      p_l_gama = p_y(lambda, L_bar, gama, X, y_one_hot, n, p, k);
      fvalue = f(p_l_gama, X, y_one_hot, n, p, k) + g(p_l_gama, lambda);
      qvalue = Q(p_l_gama, gama, X, y_one_hot, lambda, L_bar, n, p, k);

      if (fvalue <= qvalue) {
        smallest_ik_condition = false;
      } else {
        ik++;
      }
    }

    loss_list.push_back(fvalue);

    L_new = L_bar;
    theta_new = p_y(lambda, L_new, gama, X, y_one_hot, n, p, k);
    tk = (1 + sqrt(1 + 4 * t_old * t_old)) / 2;
    gama = theta_new + (t_old - 1) / tk * (theta_new - theta_old);
    t_old = tk;
    L_old = L_new;

    //bool end_loop_condition = false;

    double beta_delta = arma::max(arma::max(arma::abs(theta_new - theta_old)));
    //double loss_delta = std::abs(loss_list[times] - loss_list[times - 1]);

    // if (!loss_compute) {
    //   end_loop_condition = (times > max_iter) || (bool)(arma::max(arma::abs(theta_new - theta_old)) < eps);
    // } else {
    //   end_loop_condition = (times > max_iter) || (bool)(arma::max(arma::abs(theta_new - theta_old)) < eps) || (bool)(std::abs(loss_list[times] - loss_list[times - 1]) < 1e-3);
    // }

    //if ((times > max_iter) || beta_delta < eps || loss_delta < 1e-3) {
    if ((times > max_iter) || beta_delta < eps ) {
      condition = false;
    } else {
      times++;
      theta_old = theta_new;
    }
  }

  return Rcpp::List::create(Rcpp::Named("theta") = reshape(theta_new, p, k),
                            Rcpp::Named("loss") = loss_list,
                            Rcpp::Named("iter_times") = times);
}
