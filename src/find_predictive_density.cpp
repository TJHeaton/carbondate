#include "cpp11.hpp"
#include <array>
#include <vector>
#include "Rmath.h"
using namespace cpp11;

double LogMarginalNormalGamma(
    double calendar_age,
    double lambda,
    double nu1,
    double nu2,
    double mu_phi);


std::vector<double> MixtureDensity_cpp(
    doubles calendar_ages,
    doubles weight,
    doubles phi,
    doubles sd) {

  int n = calendar_ages.size();
  int nclust = weight.size();
  std::vector<double> density(n, 0);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < nclust; j++) {
      density[i] += weight[j] * Rf_dnorm4(calendar_ages[i], phi[j], sd[j], 0);
    }
  }
  return density;
}


[[cpp11::register]] std::vector<double> FindPredictiveDensityWalker_cpp(
    doubles calendar_ages,
    doubles weight,
    doubles phi,
    doubles tau,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2) {

  int n = calendar_ages.size();
  int nclust = weight.size();
  std::vector<double> density(n, 0);
  double sum_weight, logmarg;

  for (int i = 0; i < n; i++) {
    sum_weight = 0.;
    for (int j = 0; j < nclust; j++) {
      density[i] += weight[j] * Rf_dnorm4(calendar_ages[i], phi[j], 1. / sqrt(tau[j]), 0);
      sum_weight += weight[j];
    }
    logmarg = LogMarginalNormalGamma(calendar_ages[i], lambda, nu1, nu2, mu_phi);
    density[i] += (1. - sum_weight) * exp(logmarg);
  }
  return density;
}


[[cpp11::register]] std::vector<double> FindPredictiveDensityPolyaUrn_cpp(
    doubles calendar_ages,
    integers cluster_identifiers,
    doubles phi,
    doubles tau,
    double alpha,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2) {

  int n = calendar_ages.size();
  int nobs = cluster_identifiers.size();
  int nclust = phi.size();
  std::vector<double> density(n, 0.);
  std::vector<int> observations_per_clust(nclust, 0);
  double logmarg;

  for (int i = 0; i < nobs; i++) {
    observations_per_clust[cluster_identifiers[i] - 1] += 1;
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < nclust; j++) {
      density[i] += observations_per_clust[j]
        * Rf_dnorm4(calendar_ages[i], phi[j], 1. / sqrt(tau[j]), 0)
        / (nobs + alpha);
    }
    logmarg = LogMarginalNormalGamma(calendar_ages[i], lambda, nu1, nu2, mu_phi);
    density[i] += alpha * exp(logmarg)/(nobs + alpha);
  }
  return density;
}
