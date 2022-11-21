#include <limits>
#include <algorithm>
#include <vector>
#include "cpp11.hpp"
#include "Rmath.h"
#include "local_rng.h"
using namespace cpp11;


doubles UpdatePhiTau_cpp(
    doubles calendar_ages,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2);

double find_a(integers cluster_identifiers, int clust_num, double brprod, std::vector<double> u);

double find_b(integers cluster_identifiers, int clust_num, std::vector<double> u, std::vector<double> v);

double min_value(std::vector<double> vec);


[[cpp11::register]] list DPWalkerUpdate_cpp(
    doubles calendar_ages,
    doubles weightprev,
    doubles vprev,
    integers cluster_identifiers,
    doubles phi,
    doubles tau,
    int n_clust,
    double alpha,
    double mu_phi,
    double lambda,
    double nu1,
    double nu2) {

  int n = calendar_ages.size();
  using namespace cpp11::literals;
  local_rng rng_state;
  std::vector<double> u(n);  // auxilliary variables
  double compvar;
  cpp11::writable::list retlist;
  std::vector<double> weight;
  weight.reserve(2*n_clust);
  int clust_num = 0;
  double brprod = 1.;
  double sum_weight = 0.;
  double new_weight;

  std::vector<double> v(n_clust);
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

