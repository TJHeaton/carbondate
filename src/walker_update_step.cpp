#include <vector>
#include "cpp11.hpp"
#include "Rmath.h"
#include "local_rng.h"
using namespace cpp11;


void WalkerUpdateWeights(
    integers&, const std::vector<double>&, int, double, double,
    std::vector<double>&, std::vector<double>&, int&);

void WalkerUpdateClusterPhiTau(
    int, const doubles&, const integers&, double, double, double, double,
    std::vector<double>&, std::vector<double>&);

void WalkerUpdateClusterIdentifiers(
    doubles&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&,
    const std::vector<double>&, std::vector<int>&);

double WalkerUpdateAlpha(
    const std::vector<int>&, double, double, double);

double UpdateMuPhi(
    const std::vector<double>&, const std::vector<double>&, double, double, double);

std::vector<double> UpdateCalendarAges(
        int, const doubles&, double, double, const std::vector<int>&, const std::vector<double>&,
        const std::vector<double>&, const doubles&, const doubles&, int, const doubles&, const doubles&);

// Performs one iteration of Walker DP update on all model parameters.
// Updates the following parameters:
// - calendar_ages
// - weights
// - v
// - phi, tau
// - cluster_ids
// - n_clust
// - alpha
// - mu_phi
// Uses the following hyperparameters which are unchanged
// - alpha_shape, alpha_rate, lambda, nu1, nu2, A, B
[[cpp11::register]] list WalkerUpdateStep(
    doubles current_calendar_ages, // current calendar age each for each observation
    doubles current_weight,        // weight per cluster
    doubles current_v,
    integers current_cluster_ids,  // cluster each observation belongs to
    int current_n_clust,           // current number of clusters
    double current_alpha,          // current value of the DPMM concentration parameter
    double current_mu_phi,         // current value of the overall cluster centering
    double alpha_shape,
    double alpha_rate,
    double lambda,
    double nu1,
    double nu2,
    double A,
    double B,
    double w,
    double m,
    doubles c14_determinations,
    doubles c14_sigmas,
    int calcurve_yr_index_offset,
    doubles mucalallyr,
    doubles sigcalallyr) {

  local_rng rng_state;                   // Ensures RNG follows R and R follows after
  int n = current_calendar_ages.size();  // Number of observations
  std::vector<double> u(n);              // Auxiliary variables
  std::vector<double> weight;            // Updated weights
  //weight.reserve(5*current_n_clust);     // Reserve extra space for new clusters
  std::vector<double> v(current_v.begin(), current_v.end());  // Updated v
  std::vector<int> cluster_ids(n);      // Updated cluster_ids
  std::vector<double> phi, tau;         // Updated cluster means and precisions
  std::vector<double> calendar_ages;    // Updated calendar_ages
  int n_clust;                          // Updated number of clusters
  double mu_phi, alpha;                 // Updated value of mu_phi and alpha
  double min_u = 1.;                    // Minimum value of auxiliary variables

  using namespace cpp11::literals;
  cpp11::writable::list retlist;

  // Create auxiliary variables u
  for (int k = 0; k < n; k++) {
    if (current_cluster_ids[k] > current_weight.size()) printf("error 1 \n");
    u[k] = Rf_runif(0., current_weight[current_cluster_ids[k] - 1]);
    if (u[k] < min_u) min_u = u[k];   // update minimum value of u
  }

  WalkerUpdateWeights(
    current_cluster_ids, u, current_n_clust, min_u, current_alpha, v, weight, n_clust);

  // Update the cluster means and precisions, introducing new ones for those without observations
  // First set the vectors to the correct size now n_clust has been updated
  phi.resize(n_clust);
  tau.resize(n_clust);
  WalkerUpdateClusterPhiTau(
    n_clust, current_calendar_ages, current_cluster_ids, current_mu_phi, lambda, nu1, nu2, phi, tau);

  // Now update the cluster id for each observation
  WalkerUpdateClusterIdentifiers(current_calendar_ages, u, weight, phi, tau, cluster_ids);

  alpha = WalkerUpdateAlpha(cluster_ids, current_alpha, alpha_shape, alpha_rate);
  mu_phi = UpdateMuPhi(phi, tau, lambda, A, B);

  calendar_ages = UpdateCalendarAges(
    n,
    current_calendar_ages,
    w,
    m,
    cluster_ids,
    phi,
    tau,
    c14_determinations,
    c14_sigmas,
    calcurve_yr_index_offset,
    mucalallyr,
    sigcalallyr);

  cpp11::writable::doubles weight_out(weight.begin(), weight.end());
  cpp11::writable::doubles v_out(v.begin(), v.end());
  cpp11::writable::doubles phi_out(phi.begin(), phi.end());
  cpp11::writable::doubles tau_out(tau.begin(), tau.end());
  cpp11::writable::doubles calendar_ages_out(calendar_ages.begin(), calendar_ages.end());
  cpp11::writable::integers cluster_ids_out(cluster_ids.begin(), cluster_ids.end());

  // Return the updated parameters
  retlist.push_back({"weight"_nm = weight_out});
  retlist.push_back({"v"_nm = v_out});
  retlist.push_back({"cluster_ids"_nm = cluster_ids_out});
  retlist.push_back({"phi"_nm = phi_out});
  retlist.push_back({"tau"_nm = tau_out});
  retlist.push_back({"n_clust"_nm = n_clust});
  retlist.push_back({"alpha"_nm = alpha});
  retlist.push_back({"mu_phi"_nm = mu_phi});
  retlist.push_back({"calendar_ages"_nm = calendar_ages_out});
  return retlist;
}
