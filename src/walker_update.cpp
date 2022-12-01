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


[[cpp11::register]] list DPWalkerUpdate_cpp(
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
