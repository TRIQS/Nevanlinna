#include <cmath>
#include "./nevanlinna.hpp"
#include "nevanlinna_error.hpp"

namespace nevanlinna {

  void solver::solve(const triqs::gfs::gf<triqs::mesh::imfreq,triqs::gfs::scalar_valued>& g_iw) {
    nda::array_view<std::complex<double>, 1> view = g_iw.data();
    nda::array<std::complex<double>, 1> mesh(view.size());
    std::copy(g_iw.mesh().begin(), g_iw.mesh().end(), mesh.begin());
    solve(mesh, view);
  }
//  template<typename T>
  void solver::solve(const nda::array<std::complex<double>, 1> & mesh, const nda::array<std::complex<double>, 1> & data) {
    assert(mesh.shape(0) == data.shape(0));
    size_t M = mesh.shape(0);
    _phis.resize(M);
    _abcds.resize(M);
    _mesh.resize(M);
    _phis[0] = complex_t(data(0));
    for (int k = 0; k < M; k++) {
      _abcds[k] = matrix_t::Identity(2, 2);
      _mesh[k] = mesh(k);
    }
    matrix_t prod(2, 2);
    for (int j = 0; j < M - 1; j++) {
      for (int k = j; k < M; k++) {
        prod(0,0) = (_mesh[k] - _mesh[j]) / (_mesh[k] - std::conj(_mesh[j]));
        prod(0,1) = _phis[j];
        prod(1,0) = std::conj(_phis[j]) * (_mesh[k] - _mesh[j]) / (_mesh[k] - std::conj(_mesh[j]));
        prod(1,1) = complex_t{1., 0.};
        _abcds[k] *= prod;
      }
      _phis[j + 1] = (- _abcds[j + 1](1, 1) * complex_t(data(j + 1)) + _abcds[j + 1](0, 1)) /
         (_abcds[j + 1](1, 0) * complex_t(data(j + 1)) - _abcds[j + 1](0, 0));
    }
    return;
  }

  nda::array<double, 1> solver::evaluate(const nda::array<double, 1> &grid, double eta) const {
     nda::array<std::complex<double>, 1> complex_grid(grid.shape());
     complex_grid = grid + eta*1i;
     return evaluate(complex_grid);
  }

  nda::array<double, 1> solver::evaluate(const nda::array<std::complex<double>, 1> &grid) const {
    size_t M = _phis.size();
    if(M == 0) {
      throw nevanlinna_error("Empty continuation data. Please run solve(...) first.");
    }
    complex_t I {0., 1.};
    complex_t One {1., 0.};
    nda::array<double, 1> results(grid.shape());
    matrix_t prod(2, 2);
    for (int i = 0; i< grid.shape(0); ++i) {
      matrix_t result = matrix_t::Identity(2, 2);
      complex_t z = complex_t(grid(i));
      for (int j = 0; j < M; j++) {
        prod(0,0) =  (z - _mesh[j]) / (z - std::conj(_mesh[j]));
        prod(0,1) = _phis[j];
        prod(1,0) = std::conj(_phis[j])* ((z - _mesh[j]) / (z - std::conj(_mesh[j])));
        prod(1,1) = complex_t{1., 0.};
        result *= prod;
      }
      complex_t param {0., 0.}; //theta_{M+1}, choose to be constant function 0 here
      complex_t theta = (result(0, 0) * param + result(0, 1)) / (result(1, 0) * param + result(1, 1));
      results(i) = (1/M_PI * (I * (One + theta) / (One - theta)).imag().get_d()); //inverse Mobius transform from theta to NG
    }
    return results;
  }

  solver::solver(nevanlinna_parameters_t const & p) {
    mpf_set_default_prec(p.precision);
  }

} // namespace nevanlinna
