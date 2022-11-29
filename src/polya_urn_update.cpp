#include <vector>
#include "cpp11.hpp"
#include "Rmath.h"
#include "local_rng.h"
using namespace cpp11;

int SampleInt(int n, std::vector<double> prob, bool one_based);


double LogMarginalNormalGamma(
    double calendar_age,
    double lambda,
    double nu1,
    double nu2,
    double mu_phi) {

  double logden, margprec, margdf;

  margprec = (nu1 * lambda) / (nu2 * (lambda + 1.));
  margdf = 2. * nu1;

  logden = lgamma((margdf + 1.) / 2.) - lgamma(margdf / 2.);
  logden += 0.5 * (log(margprec) - log(margdf) - log(M_PI));
  logden -= ((margdf + 1) / 2) * log(1 + margprec * pow(calendar_age - mu_phi, 2) / margdf);

  return logden;
}


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


[[cpp11::register]] list PolyaUrnUpdateClusterIdentifier(
    doubles calendar_ages,         // calendar each for each observation
    integers current_cluster_ids,  // The cluster each observation belongs to
    doubles current_phi,
    doubles current_tau,
    double alpha,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2) {

  local_rng rng_state;                // Ensures RNG follows R and R follows after
  int n = calendar_ages.size();       // Number of observations
  int n_clust = current_phi.size();   // Number of clusters
  std::vector<int> cluster_ids(current_cluster_ids.begin(), current_cluster_ids.end());
  std::vector<double> phi(current_phi.begin(), current_phi.end());
  std::vector<double> tau(current_tau.begin(), current_tau.end());
  std::vector<int> observations_per_cluster(n_clust);
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

  // Change the cluster ID labelling so that there are no skipped values
  for (int i = 0; i < n; i++) cluster_ids[i] = cluster_id_map[cluster_ids[i]];

  // Return the updated cluster ids, means and precisions
  retlist.push_back({"cluster_ids"_nm = cluster_ids});
  retlist.push_back({"phi"_nm = phi});
  retlist.push_back({"tau"_nm = tau});
  retlist.push_back({"observations_per_cluster"_nm = observations_per_cluster});
  return retlist;
}
