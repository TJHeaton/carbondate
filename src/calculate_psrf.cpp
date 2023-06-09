#include <algorithm>
#include <vector>
#include "cpp11.hpp"
using namespace cpp11;


double mean(doubles vec, int N) {

  double mean = 0;

  for (int i = 0; i < N; i++) mean += vec[i] / N;

  return mean;
}


double var(doubles vec, int N, double mean) {

  double var = 0;

  for (int i = 0; i < N; i++) var += pow(vec[i] - mean, 2) / N;

  return var;
}


[[cpp11::register]] doubles CalculatePsrf(list thetas, int start_sample, int end_sample) {

  int M = thetas.size();
  int n_points = end_sample - start_sample + 1;
  double between_chain_variance, within_chain_variance, pooled_variance;
  double overall_mean;
  writable::doubles R(n_points);

  std::vector<double> theta_mean(M);
  std::vector<double> theta_variance(M);

  for (int N = start_sample; N <= end_sample; N++) {
    between_chain_variance = overall_mean = within_chain_variance = 0;
    for (int m = 0; m < M; m++) {
      theta_mean[m] = mean(thetas[m], N);
      within_chain_variance += var(thetas[m], N, theta_mean[m]) / M;
      overall_mean += theta_mean[m] / M;
    }
    for (int m = 0; m < M; m++) {
      between_chain_variance += N / (M - 1) * pow(theta_mean[m] - overall_mean, 2);
    }
    pooled_variance = (N - 1.0) / N * within_chain_variance;
    pooled_variance += (M + 1.0) / (M * N) * between_chain_variance;

    R[N - start_sample] = pooled_variance / within_chain_variance;
  }

  return R;
}
