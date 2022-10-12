#include <triqs/test_tools/gfs.hpp>
#include "nevanlinna/solver.hpp"

using namespace nevanlinna;
using nda::clef::placeholder;
using namespace std::complex_literals;

TEST(NevanlinnaSolver, TriqsGFData) { // NOLINT
  nevanlinna_parameters_t p;
  solver a(p);
  const double eta = 0.1;
  const double mu = 1.0;
  const double beta = 10.0;
  const int n_iw = 100;
  const int n_omega = 100;

  placeholder<0> w_;

  triqs::gfs::gf<triqs::mesh::imfreq> gf(triqs::mesh::imfreq(beta, triqs::mesh::Fermion, n_iw), {1, 1});
  triqs::mesh::refreq re_grid(-5.0, 5.0, n_omega);
  triqs::gfs::gf<triqs::mesh::refreq> re_data(re_grid, {1, 1});
  gf(w_) << 1.0/(w_ - mu);
  re_data(w_) << 1.0/(w_ + eta*1i - mu);
  a.solve(gf);
  triqs::gfs::gf<triqs::mesh::refreq> result = a.evaluate(re_grid, eta);
  ASSERT_TRUE(array_are_close(re_data.data(), result.data(), 1e-10));
}

