#include <vector>
#include "cpp11.hpp"
#include "Rmath.h"
using namespace cpp11;

int SampleInt(int, std::vector<double>, bool);

double LogMarginalNormalGamma(double, double, double, double, double);

void UpdatePhiTau(
    const std::vector<double>&, double, double, double, double, double&, double&);


void CreateNewPhiTau(
    double calendar_age,
    double lambda,
    double nu1,
    double nu2,
    double mu_phi,
    double &phi,
    double &tau) {

  nu1 += 0.5;
  nu2 += lambda * pow(calendar_age - mu_phi, 2) / (2. * (lambda + 1.));
  mu_phi = (lambda * mu_phi + calendar_age) / (lambda + 1);
  lambda++;

  tau = Rf_rgamma(nu1, 1. / nu2);
  phi = Rf_rnorm(mu_phi, 1. / sqrt(lambda * tau));
}


void PolyaUrnUpdateClusterIds(
    const doubles& calendar_ages,
    double alpha,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2,
    std::vector<int>& cluster_ids,
    std::vector<double>& phi,
    std::vector<double>& tau,
    std::vector<int>& observations_per_cluster) {

  int n = calendar_ages.size();       // Number of observations
  int n_clust = phi.size();           // Number of clusters
  int cluster_id, new_cluster_id;
  std::vector<double> cprob(n_clust + 1); // Probability of sampling each cluster ID
  double logmarg;
  double phi_new, tau_new;           // New phi and tau values when we sample a new cluster
  int n_non_empty_clust = n_clust;   // Keeps track of how many clusters we have each iteration
  std::vector<int> cluster_id_map;   // Maps the old cluster ids to the new shifted cluster ids
  cpp11::writable::list retlist;

  using namespace cpp11::literals;

  // Reserve memory to increase vectors if we introduce new clusters
  cprob.reserve(2 * n_clust);
  phi.reserve(2 * n_clust);
  tau.reserve(2 * n_clust);

  // First initialise the number of observations in each cluster
  for (int i = 0; i < n; i++) observations_per_cluster[cluster_ids[i] - 1]++;

  for (int i = 0; i < n; i++) {
    cluster_id = cluster_ids[i];

    // Updating this cluster ID, so remove from count
    observations_per_cluster[cluster_id - 1]--;

    if (observations_per_cluster[cluster_id - 1] == 0) n_non_empty_clust--;

    // Calculate the probability for sampling each new cluster_id for this observation
    // (including a new cluster ID)
    for (int c = 1; c <= n_clust; c++) {
      if (observations_per_cluster[c-1] == 0) {
        cprob[c - 1] = 0.;
      } else {
        cprob[c - 1] = Rf_dnorm4(calendar_ages[i], phi[c - 1], 1./sqrt(tau[c - 1]), 0);
        cprob[c - 1] *= observations_per_cluster[c - 1];
      }
      logmarg = LogMarginalNormalGamma(calendar_ages[i], lambda, nu1, nu2, mu_phi);
      cprob[n_clust] = exp(logmarg) * alpha;
    }
    // Sample cluster ID for the new cluster
    new_cluster_id = SampleInt(n_clust + 1, cprob, 1);

    // If we've sampled a new cluster, add new phi and tau and update other variables to reflect
    if (new_cluster_id == n_clust + 1) {
      CreateNewPhiTau(calendar_ages[i], lambda, nu1, nu2, mu_phi, phi_new, tau_new);
      phi.push_back(phi_new);
      tau.push_back(tau_new);
      cprob.push_back(0.);
      observations_per_cluster.push_back(1);
      n_clust++;
      n_non_empty_clust++;
    } else {
      observations_per_cluster[new_cluster_id - 1]++;
    }
    cluster_ids[i] = new_cluster_id;
  }

  // Create a map of old cluster labelling to new cluster labelling
  // Shift phi and tau to remove values where there are no observations
  cluster_id_map.resize(n_clust + 1);
  int newc = 1;
  for (int c = 1; c <= n_clust; c++) {
    if (observations_per_cluster[c-1] > 0) {
      cluster_id_map[c] = newc;
      phi[newc - 1] = phi[c - 1];
      tau[newc - 1] = tau[c - 1];
      observations_per_cluster[newc - 1] = observations_per_cluster[c-1];
      newc++;
    }
  }
  phi.resize(n_non_empty_clust);
  tau.resize(n_non_empty_clust);
  observations_per_cluster.resize(n_non_empty_clust);

  // Change the cluster ID labelling so that there are no skipped values
  for (int i = 0; i < n; i++) cluster_ids[i] = cluster_id_map[cluster_ids[i]];
}


void PolyaUrnUpdateClusterPhiTau(
    const doubles& calendar_ages,
    const std::vector<int>& cluster_ids,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2,
    std::vector<double>& phi,
    std::vector<double>& tau) {

  int n_clust = phi.size();
  int n = cluster_ids.size();
  std::vector<double> cluster_calendar_ages(0);
  cluster_calendar_ages.reserve(n);

  for (int c = 1; c <= n_clust; c++) {
    // Find out which observations belong in this cluster
    for (int j = 0; j < n; j++) {
      if (cluster_ids[j] == c) cluster_calendar_ages.push_back(calendar_ages[j]);
    }
    UpdatePhiTau(cluster_calendar_ages, mu_phi, lambda, nu1, nu2, phi[c-1], tau[c-1]);
    cluster_calendar_ages.resize(0);
  }
}

// Works out the likelihood of a particular Dir(alpha) given a partition cluster_identifiers
double PolyaUrnAlphaLogLikelihood(
    const std::vector<int>& observations_per_cluster, double alpha, double n) {
  int n_clust = observations_per_cluster.size();
  double log_likelihood = n_clust * log(alpha);
  int c_shift;

  for (int i = 0; i < n_clust; i++) {
    // Use max of value or 1 since 0! is 1
    c_shift = observations_per_cluster[i] <= 1 ? 1 : observations_per_cluster[i] - 1;
    for (int j = 1; j <= c_shift; j++) log_likelihood += log(j);
  }
  for (int i = 0; i < n; i++) log_likelihood -= log(alpha + i);
  return log_likelihood;
}


double PolyaUrnUpdateAlpha(
    int n,
    const std::vector<int>& observations_per_cluster,
    double current_alpha,
    double alpha_shape,
    double alpha_rate){

  double alpha = -1.;         // Updated value of alpha
  double prop_sd = 1.;        // Standard deviation for sampling proposed value of alpha
  double log_prior_rate, log_likelihood_rate, log_proposal_rate, hr;

  // Sample alpha from truncated normal distribution
  while (alpha <= 0.) {
    alpha = Rf_rnorm(current_alpha, prop_sd);
  }

  log_prior_rate = Rf_dgamma(alpha, alpha_shape, 1./alpha_rate, 1);
  log_prior_rate -= Rf_dgamma(current_alpha, alpha_shape, 1./alpha_rate, 1);
  log_likelihood_rate = PolyaUrnAlphaLogLikelihood(observations_per_cluster, alpha, n);
  log_likelihood_rate -= PolyaUrnAlphaLogLikelihood(observations_per_cluster, current_alpha, n);
  // Adjust for non-symmetric truncated normal proposal
  log_proposal_rate = Rf_pnorm5(current_alpha, 0., 1., 1, 1) - Rf_pnorm5(alpha, 0., 1., 1, 1);
  hr = exp(log_prior_rate + log_likelihood_rate + log_proposal_rate);
  // Accept or reject new alpha
  if (Rf_runif(0., 1.) < hr) {
    return alpha;
  }
  return current_alpha;
}


[[cpp11::register]] double PolyaUrnUpdateAlpha_test(
    int n,
    integers nci,
    double current_alpha,
    double alpha_shape,
    double alpha_rate){

  double alpha = -1.;         // Updated value of alpha
  double prop_sd = 1.;        // Standard deviation for sampling proposed value of alpha
  double log_prior_rate, log_likelihood_rate, log_proposal_rate, hr;
  std::vector<int> observations_per_cluster(nci.begin(), nci.end());

  // Sample alpha from truncated normal distribution
  while (alpha <= 0.) {
    alpha = Rf_rnorm(current_alpha, prop_sd);
  }

  log_prior_rate = Rf_dgamma(alpha, alpha_shape, 1./alpha_rate, 1);
  log_prior_rate -= Rf_dgamma(current_alpha, alpha_shape, 1./alpha_rate, 1);
  log_likelihood_rate = PolyaUrnAlphaLogLikelihood(observations_per_cluster, alpha, n);
  log_likelihood_rate -= PolyaUrnAlphaLogLikelihood(observations_per_cluster, current_alpha, n);
  // Adjust for non-symmetric truncated normal proposal
  log_proposal_rate = Rf_pnorm5(current_alpha, 0., 1., 1, 1) - Rf_pnorm5(alpha, 0., 1., 1, 1);
  hr = exp(log_prior_rate + log_likelihood_rate + log_proposal_rate);
  // Accept or reject new alpha
  if (Rf_runif(0., 1.) < hr) {
    return alpha;
  }
  return current_alpha;
}
