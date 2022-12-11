#include <vector>
#include "cpp11.hpp"
#include "Rmath.h"
#include "local_rng.h"
using namespace cpp11;

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
    std::vector<int>& observations_per_cluster);

double UpdateMuPhi(
    const std::vector<double>&, const std::vector<double>&, double, double, double);

std::vector<double> UpdateCalendarAges(
    int, doubles, double, double, std::vector<int>, std::vector<double>, std::vector<double>,
    doubles, doubles, doubles, doubles);

void UpdatePhiTau(
    const std::vector<double>&, double, double, double, double, double&, double&);

void PolyaUrnUpdateClusterPhiTau(
    const doubles&, const std::vector<int>&, double, double, double, double, std::vector<double>&,
    std::vector<double>&);

double PolyaUrnUpdateAlpha(int, const std::vector<int>&, double, double, double);


// Performs one iteration of Polya Urn DP update on all model parameters.
// Updates the following parameters:
// - calendar_ages
// - phi, tau
// - cluster_ids
// - n_clust
// - alpha
// - mu_phi
// Uses the following hyperparameters which are unchanged
// - alpha_shape, alpha_rate, lambda, nu1, nu2, A, B
[[cpp11::register]] list PolyaUrnUpdateStep(
    doubles current_calendar_ages, // current calendar age each for each observation
    integers current_cluster_ids,  // cluster each observation belongs to
    doubles current_phi,           // Vector of cluster means
    doubles current_tau,           // Vector of cluster precisions
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
    doubles mucalallyr,
    doubles sigcalallyr) {

  local_rng rng_state;                        // Ensures RNG follows R and R follows after
  int n = current_calendar_ages.size();       // Number of observations
  int n_clust = current_phi.size();           // Number of clusters
  std::vector<int> cluster_ids(current_cluster_ids.begin(), current_cluster_ids.end());
  std::vector<double> phi(current_phi.begin(), current_phi.end());
  std::vector<double> tau(current_tau.begin(), current_tau.end());
  std::vector<int> observations_per_cluster(n_clust);
  std::vector<double> calendar_ages;          // Updated calendar_ages
  double mu_phi, alpha;                       // Updated value of mu_phi and alpha

  using namespace cpp11::literals;
  cpp11::writable::list retlist;

  PolyaUrnUpdateClusterIds(
    current_calendar_ages,
    current_alpha,
    current_mu_phi,
    lambda,
    nu1,
    nu2,
    cluster_ids,
    phi,
    tau,
    observations_per_cluster);

  PolyaUrnUpdateClusterPhiTau(current_calendar_ages, cluster_ids, current_mu_phi, lambda, nu1, nu2, phi, tau);
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
    mucalallyr,
    sigcalallyr);

  alpha = PolyaUrnUpdateAlpha(n, observations_per_cluster, current_alpha, alpha_shape, alpha_rate);

  // Return the updated parameters
  retlist.push_back({"cluster_ids"_nm = cluster_ids});
  retlist.push_back({"phi"_nm = phi});
  retlist.push_back({"tau"_nm = tau});
  retlist.push_back({"alpha"_nm = alpha});
  retlist.push_back({"mu_phi"_nm = mu_phi});
  retlist.push_back({"calendar_ages"_nm = calendar_ages});
  retlist.push_back({"observations_per_cluster"_nm = observations_per_cluster});
  return retlist;
  }
