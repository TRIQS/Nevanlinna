#include <triqs/test_tools/gfs.hpp>
#include <triqs_Nevanlinna/solver_core.hpp>
#include "test_common.hpp"

using namespace triqs_Nevanlinna;
using nda::clef::placeholder;
using namespace std::complex_literals;

TEST(NevanlinnaSolver, TriqsGFData) { // NOLINT
  Nevanlinna_parameters_t p;
  solver_core a(p);
  const double eta  = 0.1;
  const double mu   = 1.0;
  const double beta = 10.0;
  const int n_iw    = 100;
  const int n_omega = 100;

  placeholder<0> w_;

  triqs::gfs::gf<triqs::mesh::imfreq> gf(triqs::mesh::imfreq(beta, triqs::mesh::Fermion, n_iw), {1, 1});
  triqs::mesh::refreq re_grid(-5.0, 5.0, n_omega);
  triqs::gfs::gf<triqs::mesh::refreq> re_data(re_grid, {1, 1});
  gf(w_) << 1.0 / (w_ - mu);
  re_data(w_) << 1.0 / (w_ + eta * 1i - mu);
  a.solve(gf);
  triqs::gfs::gf<triqs::mesh::refreq> result = a.evaluate(re_grid, eta);
  ASSERT_TRUE(array_are_close(re_data.data(), result.data(), 1e-10));
}

TEST(CaratheodorySolver, TriqsGFData) { // NOLINT
  Nevanlinna_parameters_t p{CARATHEODORY, 128};
  solver_core a(p);
  const double eta  = 0.1;
  const double mu   = 1.0;
  const double beta = 10.0;
  const int n_iw    = 100;
  const int n_omega = 100;

  placeholder<0> w_;

  triqs::gfs::gf<triqs::mesh::imfreq> gf(triqs::mesh::imfreq(beta, triqs::mesh::Fermion, n_iw), {2, 2});
  triqs::mesh::refreq re_grid(-5.0, 5.0, n_omega);
  triqs::gfs::gf<triqs::mesh::refreq> re_data(re_grid, {2, 2});
  nda::matrix<double> F(2, 2);
  F = {{-mu, 0.5}, {0.5, -mu}};
  gf(w_) << 1.0 / (w_ - F);
  re_data(w_) << 1.0 / (w_ + eta * 1i - F);
  a.solve(gf);
  triqs::gfs::gf<triqs::mesh::refreq> result = a.evaluate(re_grid, eta);
  ASSERT_TRUE(array_are_close(re_data.data(), result.data(), 1e-10));
}

MAKE_MAIN_MPI