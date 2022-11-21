#include "cpp11/doubles.hpp"
#include "cpp11/integers.hpp"
#include "Rmath.h"
#include <limits>
#include "local_rng.h"
using namespace cpp11;



double mean_local(doubles vec) {
  double mean = 0.0;
  for (double elem : vec) {
    mean += elem;
  }
  mean /= vec.size();
  return mean;
}


[[cpp11::register]] doubles UpdatePhiTau_cpp(
    doubles calendar_ages,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2) {

  local_rng rng_state;
  int nclus = calendar_ages.size();
  writable::doubles thetadiff(calendar_ages);
  double thetabar;
  double s;
  double mu_phi_new;
  double lambda_new, nu1_new, nu2_new;
  writable::doubles phi_tau(2);

  thetabar = mean_local(calendar_ages);
  for (int i = 0; i < thetadiff.size(); ++i) {
    thetadiff[i] = pow(thetadiff[i] - thetabar, 2);
  }
  s = mean_local(thetadiff);

  // Update parameters according to conjugate prior
  nu1_new = nu1 + nclus / 2.0;
  nu2_new = nu2;
  nu2_new += 0.5 * nclus * (s + lambda * pow(thetabar - mu_phi, 2)/(lambda + nclus));

  lambda_new = lambda + nclus;
  mu_phi_new = (lambda * mu_phi + nclus * thetabar) / (lambda + nclus);

  phi_tau[1] = Rf_rgamma(nu1_new, 1.0/nu2_new);
  phi_tau[0] = Rf_rnorm(mu_phi_new, 1.0/sqrt(lambda_new * phi_tau[1]));

  return phi_tau;
}

