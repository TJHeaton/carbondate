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


double min_value(doubles vec) {
  auto pos = std::min_element(vec.begin(), vec.end());
  return *pos;
}


[[cpp11::register]] integers which_equal(integers vec, int i) {
  std::vector<int> indices;
  indices.reserve(20);
  for (int k = 0; k < vec.size(); ++k) {
    if (vec[k] == i) {
      indices.push_back(k + 1);
    }
  }
  return integers(indices);
}


[[cpp11::register]] integers which_mt(doubles vec, double i) {
  std::vector<int> indices;
  indices.reserve(20);
  for (int k = 0; k < vec.size(); ++k) {
    if (vec[k] > i) {
      indices.push_back(k + 1);
    }
  }
  return integers(indices);
}

[[cpp11::register]] integers which_mt_int(integers vec, int i) {
  std::vector<int> indices;
  indices.reserve(20);
  for (int k = 0; k < vec.size(); ++k) {
    if (vec[k] > i) {
      indices.push_back(k + 1);
    }
  }
  return integers(indices);
}


double find_a(
    integers delta,
    int clust_num,
    double brprod,
    doubles u) {
  double a = 0;

  for (int k = 0; k < delta.size(); ++k) {
    if ((delta[k] == clust_num) && (u[k] > a)) {
      a = u[k];
    }
  }
  return a/brprod;
}



double find_b(
    integers delta,
    int clust_num,
    doubles u,
    doubles v) {

  int n = delta.size();
  int m = v.size();
  writable::doubles prodv(m + 1);
  double b = 1.0;
  double b_sub = 0.0;
  double b_sub_comp;
  int index;

  if (std::any_of(delta.cbegin(), delta.cend(), [=](int i){ return i > clust_num; })) {
    prodv[0] = 1.0 - v[0];
    for (int k = 1; k < m; ++k) {
      prodv[k] = prodv[k-1] * (1.0 - v[k]);
      prodv[k-1] /= (1.0 - v[clust_num - 1]);
    }
    prodv[m-1] /= (1.0 - v[clust_num - 1]);

    for (int k = 0; k < n; ++k) {
      if (delta[k] > clust_num) {
        index = delta[k] - 1;
        b_sub_comp = u[k] / (v[index] * prodv[index - 1]);
        if (b_sub_comp > b_sub) {
          b_sub = b_sub_comp;
        }
      }
    }
    b -= b_sub;
  }

  return b;
}



[[cpp11::register]] double update_v_j(
    integers delta,
    int clust_num,
    int n_clust,
    double brprod,
    doubles u,
    doubles v,
    double alpha) {

  local_rng rng_state;
  double v_j;
  double a;
  double b;
  double A;
  double B;

  if (clust_num <= n_clust) {

    a = find_a(delta, clust_num, brprod, u);
    b = find_b(delta, clust_num, u, v);

    A = pow(1. - a, alpha);
    B = A - pow(1. - b, alpha);
    v_j = 1. - pow(A - B * Rf_runif(0., 1.), 1. / alpha);

  } else {

    v_j = Rf_rbeta(1., alpha);
  }

  return v_j;
}



[[cpp11::register]] list DPWalkerUpdate_cpp(
    doubles calendar_ages,
    doubles weight,
    doubles v,
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
  writable::doubles u(n);
  double umin;
  cpp11::writable::list retlist;

  for (int k = 0; k < n; k++) {
    u[k] = Rf_runif(0, weight[cluster_identifiers[k] - 1]);
  }
  umin = min_value(u);

  retlist.push_back({"u"_nm = u});
  retlist.push_back({"umin"_nm = umin});
  return retlist;
}

