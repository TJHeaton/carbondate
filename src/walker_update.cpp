#include <vector>
#include "cpp11.hpp"
#include "Rmath.h"
#include "local_rng.h"
using namespace cpp11;

void WalkerUpdateWeights(
    integers& cluster_ids,
    const std::vector<double>& u,
    int n,
    int current_n_clust,
    double min_u,
    double alpha,
    std::vector<double>& v,
    std::vector<double>& weight,
    int& n_clust);

void WalkerUpdateClusterPhiTau(
    int n_clust,
    const doubles& calendar_ages,
    const integers& cluster_identifiers,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2,
    std::vector<double>& phi,
    std::vector<double>& tau);

void WalkerUpdateClusterIdentifiers(
    doubles &calendar_ages,
    const std::vector<double>& u,
    const std::vector<double>& weight,
    const std::vector<double>& phi,
    const std::vector<double>& tau,
    std::vector<int>& cluster_ids);

double AlphaLogLikelihood(double n_clust, double alpha, double n);


// Performs one iteration of Walker DP update on (w, phi, tau)
[[cpp11::register]] list DPWalkerUpdate(
  doubles calendar_ages,         // calendar each for each observation
  doubles current_weight,        // The weight per cluster,
  doubles current_v,
  integers current_cluster_ids,  // The cluster each observation belongs to
  int current_n_clust,           // The current number of clusters
  double alpha,
  double mu_phi,
  double lambda,
  double nu1,
  double nu2) {

  local_rng rng_state;           // Ensures RNG follows R and R follows after
  int n = calendar_ages.size();  // Number of observations
  std::vector<double> u(n);      // Auxiliary variables
  std::vector<double> weight;    // Updated weights
  weight.reserve(2*current_n_clust);
  std::vector<double> v(current_v.begin(), current_v.end());
  std::vector<int> cluster_ids(n);  // Updated
  std::vector<double> phi, tau;  // Updated cluster means and precisions
  int n_clust;                   // Updated number of clusters
  double min_u = 1.;             // Minimum value of auxiliary variables

  using namespace cpp11::literals;
  cpp11::writable::list retlist;

  // Create auxiliary variables u
  for (int k = 0; k < n; k++) {
    u[k] = Rf_runif(0., current_weight[current_cluster_ids[k] - 1]);
    if (u[k] < min_u) min_u = u[k];   // update minimum value of u
  }

  WalkerUpdateWeights(current_cluster_ids, u, n, current_n_clust, min_u, alpha, v, weight, n_clust);

  // Update the cluster means and precisions, introducing new ones for those without observations
  // First set the vectors to the correct size now n_clust has been updated
  phi.resize(n_clust);
  tau.resize(n_clust);
  WalkerUpdateClusterPhiTau(
    n_clust, calendar_ages, current_cluster_ids, mu_phi, lambda, nu1, nu2, phi, tau);

  // Now update the cluster id for each observation
  WalkerUpdateClusterIdentifiers(calendar_ages, u, weight, phi, tau, cluster_ids);

  // Return the updated weights, cluster ids, means and precisions
  retlist.push_back({"weight"_nm = weight});
  retlist.push_back({"v"_nm = v});
  retlist.push_back({"cluster_ids"_nm = cluster_ids});
  retlist.push_back({"phi"_nm = phi});
  retlist.push_back({"tau"_nm = tau});
  retlist.push_back({"n_clust"_nm = n_clust});
  return retlist;

}


// Updates the Dirichlet process parameter alpha via Metropolis-Hastings
// with truncated proposal distribution
[[cpp11::register]] double WalkerUpdateAlpha(
    integers cluster_ids,   // The cluster each observation belongs to
    double current_alpha,
    double alpha_shape,
    double alpha_rate) {

  local_rng rng_state;  // Ensures RNG follows R and R follows after
  int n = cluster_ids.size();
  double alpha;         // Updated value of alpha
  int n_clust = 0, cluster_id;
  std::vector<int> observations_per_cluster(2*n, 0);  // Allocate plenty of space
  double log_prior_rate, log_likelihood_rate, log_proposal_rate, hr;

  // Sample new alpha from truncated normal distribution
  while (true) {
    alpha = Rf_rnorm(current_alpha, 1.);
    if (alpha > 0.) break;
  }

  // Find number of distinct populated clusters
  for (int i = 0; i < n; i++) {
    cluster_id = cluster_ids[i];
    if (observations_per_cluster[cluster_id - 1] == 0) n_clust++;
    observations_per_cluster[cluster_id - 1]++;
  }

  log_prior_rate = Rf_dgamma(alpha, alpha_shape, 1./alpha_rate, 1);
  log_prior_rate -= Rf_dgamma(current_alpha, alpha_shape, 1./alpha_rate, 1);
  log_likelihood_rate = AlphaLogLikelihood(n_clust, alpha, n);
  log_likelihood_rate -= AlphaLogLikelihood(n_clust, current_alpha, n);
  // Adjust for non-symmetric truncated normal proposal
  log_proposal_rate = Rf_pnorm5(current_alpha, 0., 1., 1, 1) - Rf_pnorm5(alpha, 0., 1., 1, 1);
  hr = exp(log_prior_rate + log_likelihood_rate + log_proposal_rate);
  // Accept or reject new alpha
  if (Rf_runif(0., 1.) < hr) {
    return alpha;
  }
  return current_alpha;
  }


// Updates mu_phi via Gibbs sampling based upon current (phi, tau) values
[[cpp11::register]] double UpdateMuPhi(
    doubles phi,   // means of clusters
    doubles tau,   // precisions of clusters
    double lambda,
    double A,      // prior mean of mu_phi
    double B) {    // prior precision of mu_phi

  local_rng rng_state;  // Ensures RNG follows R and R follows after
  double mu_phi;
  double posterior_mean, posterior_precision;
  int n = phi.size();
  double sum_tau = 0.;
  double sum_tau_mult_phi = 0.;

  for (int j = 0; j < n; j++) {
    sum_tau += tau[j];
    sum_tau_mult_phi += tau[j] * phi[j];
  }

  posterior_mean = (A * B + lambda * sum_tau_mult_phi) / (B + lambda * sum_tau);
  posterior_precision = B + lambda * sum_tau;
  mu_phi = Rf_rnorm(posterior_mean, 1. / sqrt(posterior_precision));
  return mu_phi;
}
