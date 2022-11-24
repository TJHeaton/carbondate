#include <vector>
#include "cpp11.hpp"
#include "Rmath.h"
using namespace cpp11;

int SampleInt(int n, std::vector<double> prob, bool one_based);


double Mean(const std::vector<double>& vec) {
  double mean = 0.0;
  for (double elem : vec) {
    mean += elem;
  }
  mean /= vec.size();
  return mean;
}


void UpdatePhiTau(
    const std::vector<double>& calendar_ages,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2,
    double& phi,
    double& tau);


double FindNewV(
    integers& cluster_ids,
    int cluster_id,
    double brprod,
    double alpha,
    const std::vector<double>& u,
    const std::vector<double>& v) {

  int n = cluster_ids.size();
  int m = v.size();
  writable::doubles prodv(m + 1);
  bool prodv_set = false;  // Flag for whether we've calculated prodv yet
  double b_sub = 0.0;
  double b_sub_max;
  int index;
  double a = 0., b = 1., A, B;

 // We have to work out a and b and sample a new v_(cluster_id) by inverse cdf
  for (int k = 0; k < n; ++k) {
    if (cluster_ids[k] > cluster_id) {
      if (!prodv_set) {
        // Vector to aid denominator for b: ith entry is Prod_{l < i, l != j} (1- v_l)
        // In fact strictly here  l == j included, but we multiply by (1 - v_j) below
        prodv[0] = 1.0 - v[0];
        for (int k = 1; k < m; ++k) prodv[k] = prodv[k-1] * (1.0 - v[k]);
        prodv_set = true;
      }
      index = cluster_ids[k] - 1;
      b_sub_max = u[k] / (v[index] * prodv[index - 1]);
      if (b_sub_max > b_sub) b_sub = b_sub_max;
    } else if ((cluster_ids[k] == cluster_id) && (u[k] > a)) {
      // Set a to the max value of u[k] for the observations in the current cluster
      a = u[k];
    }
  }
  a /= brprod;
  b -= b_sub * (1.0 - v[cluster_id - 1]);

  A = pow(1. - a, alpha);
  B = A - pow(1. - b, alpha);

  return 1. - pow(A - B * Rf_runif(0., 1.), 1. / alpha);
}


// Iteratively updates the weights until we have all we need
// This can change the number of weights (i.e. number of clusters)
void WalkerUpdateWeights(
    integers& cluster_ids,
    const std::vector<double>& u,
    int n,
    int current_n_clust,
    double min_u,
    double alpha,
    std::vector<double>& v,
    std::vector<double>& weight,
    int& n_clust) {

  double compvar;
  int clust_num = 0;
  double brprod = 1.;
  double sum_weight = 0.;
  double new_weight;

  compvar = 1. - min_u;
  while (sum_weight < compvar) {
    clust_num++;
    if (clust_num <= current_n_clust) {
      v[clust_num - 1] = FindNewV(cluster_ids, clust_num, brprod, alpha, u, v);
    } else {
      // v_(clust_num) just from prior beta
      v.push_back(Rf_rbeta(1., alpha));
    }
    new_weight = brprod * v[clust_num - 1];
    sum_weight += new_weight;
    weight.push_back(new_weight);
    brprod *= (1. - v[clust_num - 1]);
  }
  n_clust = clust_num;
  v.resize(n_clust);
}


void UpdatePhiTau(
    const std::vector<double>& calendar_ages,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2,
    double& phi,
    double& tau) {

  int nclus = calendar_ages.size();
  std::vector<double> thetadiff(nclus);
  double thetabar;
  double s;
  double mu_phi_new;
  double lambda_new, nu1_new, nu2_new;

  thetabar = Mean(calendar_ages);
  for (int i = 0; i < nclus; ++i) {
    thetadiff[i] = pow(calendar_ages[i] - thetabar, 2);
  }
  s = Mean(thetadiff);

  // Update parameters according to conjugate prior
  nu1_new = nu1 + nclus / 2.;
  nu2_new = nu2;
  nu2_new += 0.5 * nclus * (s + lambda * pow(thetabar - mu_phi, 2)/(lambda + nclus));

  lambda_new = lambda + nclus;
  mu_phi_new = (lambda * mu_phi + nclus * thetabar) / (lambda + nclus);

  tau = Rf_rgamma(nu1_new, 1./nu2_new);
  phi = Rf_rnorm(mu_phi_new, 1./sqrt(lambda_new * tau));
}


void WalkerUpdateClusterPhiTau(
    int n_clust,
    const doubles& calendar_ages,
    const integers& cluster_identifiers,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2,
    std::vector<double>& phi,
    std::vector<double>& tau) {

  int n = cluster_identifiers.size();
  std::vector<double> cluster_calendar_ages(0);
  cluster_calendar_ages.reserve(n);

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
      UpdatePhiTau(cluster_calendar_ages, mu_phi, lambda, nu1, nu2, phi[c-1], tau[c-1]);
      cluster_calendar_ages.resize(0);
    }
  }
}


void WalkerUpdateClusterIdentifiers(
    doubles &calendar_ages,
    const std::vector<double>& u,
    const std::vector<double>& weight,
    const std::vector<double>& phi,
    const std::vector<double>& tau,
    std::vector<int>& cluster_ids) {

  int n = calendar_ages.size();
  int nclust = weight.size();
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
    cluster_ids[i] = poss_cluster_ids[SampleInt(poss_cluster_ids.size(), dens, false)];
    poss_cluster_ids.resize(0);
    dens.resize(0);
  }
}