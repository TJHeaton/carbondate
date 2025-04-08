#include <algorithm>
#include <vector>
#include "Rmath.h"
#include "cpp11.hpp"
#include "local_rng.h"
using namespace cpp11;

double LogMarginalNormalGamma(
    double calendar_age,
    double lambda,
    double nu1,
    double nu2,
    double mu_phi);

std::vector<int> GetSampleIds(int start_index, int finish_index, int size);

// Finds the quantiles at a given edge width away from the start and end of the distribution.
// Note this partially sorts the vector provided as an argument
void EdgeQuantiles(
    std::vector<double>& vec,
    double edge_width,
    double& lower_quantile,
    double& upper_quantile) {

  double indl, indu, gl, gu;
  int jl, ju;

  // The indices we use are defined by type = 7 for the R quantile function
  indl = edge_width * ((double) vec.size() - 1.) + 1.;
  jl = std::floor(indl);
  indu = (1. - edge_width) * ((double) vec.size() - 1.) + 1.;
  ju = std::floor(indu);

  // Rather than sorting the entire vector, just sort partially to find the elements we'll look up
  std::nth_element(vec.begin(), vec.begin() + jl - 1, vec.end());
  std::nth_element(vec.begin() + jl, vec.begin() + jl, vec.end());
  std::nth_element(vec.begin() + jl + 1, vec.begin() + ju - 1, vec.end());
  std::nth_element(vec.begin() + ju, vec.begin() + ju, vec.end());

  // quantiles found using the formula for type = 7 in the R quantile function
  gl = indl - jl;
  lower_quantile = (1. - gl) * vec[jl - 1] + gl * vec[jl];

  gu = indu - ju;
  upper_quantile = (1. - gu) * vec[ju - 1] + gu * vec[ju];
}

// Pass a set of means, sds and weights and it returns the density of the corresponding mixture of
// normals.
std::vector<double> MixtureDensity_cpp(
    const doubles& calendar_ages,
    const doubles& weight,
    const doubles& phi,
    const doubles& sd) {

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

double WalkerDensityForCalendarAge(
    double calendar_age,
    const doubles& weight,
    const doubles& phi,
    const doubles& tau,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2) {

  int nclust = weight.size();
  double density = 0, sum_weight = 0., logmarg;

  for (int j = 0; j < nclust; j++) {
    density += weight[j] * Rf_dnorm4(calendar_age, phi[j], 1. / sqrt(tau[j]), 0);
    sum_weight += weight[j];
  }
  // The predictive density for a new observation is a scaled t-distribution
  logmarg = LogMarginalNormalGamma(calendar_age, lambda, nu1, nu2, mu_phi);
  density += (1. - sum_weight) * exp(logmarg);

  return density;
}

double PolyaUrnDensityForCalendarAge(
    double calendar_age,
    const integers& observations_per_cluster,
    const doubles& phi,
    const doubles& tau,
    double alpha,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2,
    double n_obs) {

  int nclust = phi.size();
  double density = 0, logmarg;

  for (int j = 0; j < nclust; j++) {
    density += observations_per_cluster[j] * Rf_dnorm4(calendar_age, phi[j], 1. / sqrt(tau[j]), 0);
  }
  // The predictive density for a new observation is a scaled t-distribution
  logmarg = LogMarginalNormalGamma(calendar_age, lambda, nu1, nu2, mu_phi);
  density += alpha * exp(logmarg);
  density /= n_obs + alpha;

  return density;
}

[[cpp11::register]] data_frame FindPredictiveDensityAndCIWalker(
    doubles calendar_ages,
    list weights,
    list phis,
    list taus,
    doubles mu_phis,
    double lambda,
    double nu1,
    double nu2,
    int n_posterior_samples,
    double quantile_edge_width,
    int n_burn,
    int n_end) {

  local_rng rng_state;// Ensures RNG follows R and R follows after
  int n = calendar_ages.size();
  // A vector of vectors to represent the density matrix
  std::vector<std::vector<double>> density_samples(n, std::vector<double>(n_posterior_samples));
  std::vector<double> mean_density(n, 0);
  std::vector<double> ci_lower(n), ci_upper(n);
  int s; //current sample id
  std::vector<int> sample_ids;

  sample_ids = GetSampleIds(n_burn, n_end - 1, n_posterior_samples);

  for (int j = 0; j < n_posterior_samples; j++) {
    s = sample_ids[j];
    doubles weight = weights[s];
    doubles phi = phis[s];
    doubles tau = taus[s];
    for (int i = 0; i < n; i++) {
      density_samples[i][j] = WalkerDensityForCalendarAge(
        calendar_ages[i], weight, phi, tau, mu_phis[s], lambda, nu1, nu2
      );
      mean_density[i] += density_samples[i][j];
    }
  }

  for (int i = 0; i < n; i++) {
    mean_density[i] /= n_posterior_samples;
    EdgeQuantiles(density_samples[i], quantile_edge_width, ci_lower[i], ci_upper[i]);
  }

  writable::data_frame retdata({
    "calendar_age_BP"_nm = calendar_ages,
    "density_mean"_nm = mean_density,
    "density_ci_lower"_nm = ci_lower,
    "density_ci_upper"_nm = ci_upper,
  });

  return std::move(retdata);
}

// Finds mean predictive density and confidence intervals, sampling over a given number of points
// in the provided data
[[cpp11::register]] data_frame FindPredictiveDensityandCIPolyaUrn(
    doubles calendar_ages,
    list observations_per_clusters,
    list phis,
    list taus,
    doubles alphas,
    doubles mu_phis,
    int n_obs,
    double lambda,
    double nu1,
    double nu2,
    int n_posterior_samples,
    double quantile_edge_width,
    int n_burn,
    int n_end) {

  local_rng rng_state;// Ensures RNG follows R and R follows after
  int n = calendar_ages.size();
  // A vector of vectors to represent the density matrix
  std::vector<std::vector<double>> density_samples(n, std::vector<double>(n_posterior_samples));
  std::vector<double> mean_density(n, 0);
  std::vector<double> ci_lower(n), ci_upper(n);
  int s; //current sample id
  std::vector<int> sample_ids;

  sample_ids = GetSampleIds(n_burn, n_end - 1, n_posterior_samples);

  for (int j = 0; j < n_posterior_samples; j++) {
    s = sample_ids[j];
    integers observations_per_cluster = observations_per_clusters[s];
    doubles phi = phis[s];
    doubles tau = taus[s];
    for (int i = 0; i < n; i++) {
      density_samples[i][j] = PolyaUrnDensityForCalendarAge(
        calendar_ages[i], observations_per_cluster, phi, tau, alphas[s], mu_phis[s], lambda, nu1,
        nu2, n_obs
      );
      mean_density[i] += density_samples[i][j];
    }
  }

  for (int i = 0; i < n; i++) {
    mean_density[i] /= n_posterior_samples;
    EdgeQuantiles(density_samples[i], quantile_edge_width, ci_lower[i], ci_upper[i]);
  }

  writable::data_frame retdata({
      "calendar_age_BP"_nm = calendar_ages,
      "density_mean"_nm = mean_density,
      "density_ci_lower"_nm = ci_lower,
      "density_ci_upper"_nm = ci_upper,
  });

  return std::move(retdata);
}

// Finds predictive density at a single data point
[[cpp11::register]] doubles FindInstantPredictiveDensityWalker(
    doubles calendar_ages,
    doubles weight,
    doubles phi,
    doubles tau,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2) {

  int n = calendar_ages.size();
  writable::doubles instant_density(n);

  for (int i = 0; i < n; i++) {
    instant_density[i] = WalkerDensityForCalendarAge(
      calendar_ages[i], weight, phi, tau, mu_phi, lambda, nu1, nu2
    );
  }
  return instant_density;
}


// Finds predictive density at a single data point
[[cpp11::register]] doubles FindInstantPredictiveDensityPolyaUrn(
    doubles calendar_ages,
    integers observations_per_cluster,
    doubles phi,
    doubles tau,
    double alpha,
    double mu_phi,
    double n_obs,
    double lambda,
    double nu1,
    double nu2) {

  int n = calendar_ages.size();
  writable::doubles instant_density(n);

  for (int i = 0; i < n; i++) {
    instant_density[i] = PolyaUrnDensityForCalendarAge(
      calendar_ages[i], observations_per_cluster, phi, tau, alpha, mu_phi, lambda, nu1, nu2, n_obs);
  }
  return instant_density;
}

// Finds mean predictive density over a set of points, using every point
[[cpp11::register]] doubles FindPredictiveDensityWalker(
    doubles calendar_ages,
    list weights,
    list phis,
    list taus,
    doubles mu_phis,
    double lambda,
    double nu1,
    double nu2) {

  int n_points = calendar_ages.size();
  int n_out = weights.size();
  writable::doubles mean_density(n_points);
  std::vector<int> sample_ids;

  for (int i = 0; i < n_points; i++) {
    mean_density[i] = 0.;
    for (int j = 0; j < n_out; j++) {
      mean_density[i] += WalkerDensityForCalendarAge(
        calendar_ages[i], weights[j], phis[j], taus[j], mu_phis[j], lambda, nu1, nu2
      );
    }
    mean_density[i] /= n_out;
  }

  return mean_density;
}
