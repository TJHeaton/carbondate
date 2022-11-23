#include <limits>
#include <algorithm>
#include <vector>
#include "cpp11.hpp"
#include "Rmath.h"
#include "local_rng.h"
using namespace cpp11;


void UpdatePhiTau(
    std::vector<double> calendar_ages,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2,
    double* phi,
    double* tau);

double find_a(integers cluster_identifiers, int clust_num, double brprod, std::vector<double> u);

double find_b(integers cluster_identifiers, int clust_num, std::vector<double> u, std::vector<double> v);

double min_value(std::vector<double> vec);

int sample_int(int n, std::vector<double> prob, bool one_based);


[[cpp11::register]] list WalkerUpdateWeights_cpp(
    doubles weightprev,
    doubles vprev,
    integers cluster_identifiers,
    int n_clust,
    double alpha) {

  local_rng rng_state;
  int n = cluster_identifiers.size();  // Number of observations
  std::vector<double> u(n);            // auxilliary variables
  std::vector<double> v(n_clust);
  std::vector<double> weight;
  weight.reserve(2*n_clust);
  double compvar;
  int clust_num = 0;
  double brprod = 1.;
  double sum_weight = 0.;
  double new_weight;

  using namespace cpp11::literals;
  cpp11::writable::list retlist;

  for (int k = 0; k < n_clust; k++) v[k] = vprev[k];

  //////////////////////////////////////////////////////////////////////////////
  // Create auxiliary u variables
  for (int k = 0; k < n; k++) {
    u[k] = Rf_runif(0., weightprev[cluster_identifiers[k] - 1]);
  }
  compvar = 1. - min_value(u);

  //////////////////////////////////////////////////////////////////////////////
  // Iteratively update the weights until we have all we need

  while (sum_weight < compvar) {
    clust_num++;
    if (clust_num <= n_clust) {

      double a = find_a(cluster_identifiers, clust_num, brprod, u);
      double b = find_b(cluster_identifiers, clust_num, u, v);
      double A = pow(1. - a, alpha);
      double B = A - pow(1. - b, alpha);

      v[clust_num - 1] = 1. - pow(A - B * Rf_runif(0., 1.), 1. / alpha);
    } else {
      v.push_back(Rf_rbeta(1., alpha));
    }
    new_weight = brprod * v[clust_num - 1];
    sum_weight += new_weight;
    weight.push_back(new_weight);
    brprod *= (1. - v[clust_num - 1]);
  }
  n_clust = clust_num;
  v.resize(n_clust);

  retlist.push_back({"u"_nm = u});
  retlist.push_back({"v"_nm = v});
  retlist.push_back({"weight"_nm = weight});
  retlist.push_back({"n_clust"_nm = n_clust});
  return retlist;
}



[[cpp11::register]] list WalkerUpdateClusterPhiTau_cpp(
    int n_clust,
    doubles calendar_ages,
    integers cluster_identifiers,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2) {

  local_rng rng_state;
  std::vector<double> phi(n_clust);
  std::vector<double> tau(n_clust);
  int n = cluster_identifiers.size();
  std::vector<double> cluster_calendar_ages(0);
  cluster_calendar_ages.reserve(n);

  using namespace cpp11::literals;
  cpp11::writable::list retlist;

  for (int c = 1; c <= n_clust; c++) {
    // Find out which observations belong in this cluster
    for (int j = 0; j < n; j++) {
      if (cluster_identifiers[j] == c) cluster_calendar_ages.push_back(calendar_ages[j]);
    }
    if (cluster_calendar_ages.size() == 0) {
      // No observations in this cluster, so sample from the prior
      tau[c-1] = Rf_rgamma(nu1, 1./nu2);
      phi[c-1] = Rf_rnorm(mu_phi, 1./sqrt(lambda * tau[c-1]));
    } else {
      // There are some observations, so update the phi and tau for this
      // cluster (conjugate according to NormalGamma prior)
      UpdatePhiTau(cluster_calendar_ages, mu_phi, lambda, nu1, nu2, &phi[c-1], &tau[c-1]);
      cluster_calendar_ages.resize(0);
    }
  }

  retlist.push_back({"phi"_nm = phi});
  retlist.push_back({"tau"_nm = tau});
  return retlist;
}


[[cpp11::register]] std::vector<int> WalkerUpdateClusterIdentifiers_cpp(
    doubles calendar_ages,
    doubles u,
    doubles weight,
    doubles phi,
    doubles tau) {

  local_rng rng_state;
  int n = calendar_ages.size();
  int nclust = weight.size();
  std::vector<int> cluster_identifiers(n);
  std::vector<int> poss_cluster_ids(0);
  poss_cluster_ids.reserve(nclust);
  std::vector<double> dens(0);
  dens.reserve(nclust);

  for (int i = 0; i < n; i++) {
    for (int c = 1; c <= nclust; c++) {
      if (weight[c-1] > u[i]) {
        poss_cluster_ids.push_back(c);
        dens.push_back(Rf_dnorm4(phi[c-1], calendar_ages[i], sqrt(1. / tau[c-1]), 0));
      }
    }
    cluster_identifiers[i] = poss_cluster_ids[sample_int(poss_cluster_ids.size(), dens, false)];
    poss_cluster_ids.resize(0);
    dens.resize(0);
  }

  return cluster_identifiers;
}
