#include "cpp11/doubles.hpp"
#include "cpp11/integers.hpp"
#include "Rmath.h"
#include <limits>
#include <vector>
#include "local_rng.h"
using namespace cpp11;



double mean_local(std::vector<double> vec) {
  double mean = 0.0;
  for (double elem : vec) {
    mean += elem;
  }
  mean /= vec.size();
  return mean;
}


void UpdatePhiTau(
    std::vector<double> calendar_ages,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2,
    double* phi,
    double* tau) {

  int nclus = calendar_ages.size();
  std::vector<double> thetadiff(nclus);
  double thetabar;
  double s;
  double mu_phi_new;
  double lambda_new, nu1_new, nu2_new;

  thetabar = mean_local(calendar_ages);
  for (int i = 0; i < nclus; ++i) {
    thetadiff[i] = pow(calendar_ages[i] - thetabar, 2);
  }
  s = mean_local(thetadiff);

  // Update parameters according to conjugate prior
  nu1_new = nu1 + nclus / 2.;
  nu2_new = nu2;
  nu2_new += 0.5 * nclus * (s + lambda * pow(thetabar - mu_phi, 2)/(lambda + nclus));

  lambda_new = lambda + nclus;
  mu_phi_new = (lambda * mu_phi + nclus * thetabar) / (lambda + nclus);

  *tau = Rf_rgamma(nu1_new, 1./nu2_new);
  *phi = Rf_rnorm(mu_phi_new, 1./sqrt(lambda_new * *tau));
}


[[cpp11::register]] std::vector<double> UpdatePhiTau_cpp(
    doubles calendar_ages,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2) {

  local_rng rng_state;

  int nclus = calendar_ages.size();
  std::vector<double> cluster_calendar_ages(nclus);
  std::vector<double> phi_tau(2);
  writable::doubles retvalue(2);

  for (int i = 0; i < nclus; ++i) {
    cluster_calendar_ages[i] = calendar_ages[i];
  }

  UpdatePhiTau(cluster_calendar_ages, mu_phi, lambda, nu1, nu2, &phi_tau[0], &phi_tau[1]);

  retvalue[0] = phi_tau[0];
  retvalue[1] = phi_tau[1];

  return phi_tau;
}

