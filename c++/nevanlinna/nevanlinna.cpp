#include <cmath>
#include "./nevanlinna.hpp"

namespace nevanlinna {

//  template<typename T>
  void solver::solve(const nda::array<std::complex<double>, 1> & mesh, const nda::array<std::complex<double>, 1> & data) {
    assert(mesh.shape(0) == data.shape(0));
    size_t M = mesh.shape(0);
    _phis.resize(M);
    _abcds.resize(M);
    _phis[0] = data(0);
    for (int k = 0; k < M; k++) _abcds[k] = matrix_t::Identity(2, 2);
    for (int j = 0; j < M - 1; j++) {
      for (int k = j; k < M; k++) {
        matrix_t prod(2, 2);
        prod(0, 0) = (mesh(k) - mesh(j)) / (mesh(k) - std::conj(mesh(j)));
        prod(0, 1) = _phis[j];
        prod(1, 0) = std::conj(_phis[j]) * (complex_t(mesh(k) - mesh(j)) / complex_t(mesh(k) - std::conj(mesh(j))));
        prod(1, 1) = complex_t{1., 0.};
        _abcds[k] *= prod;
      }
      _phis[j + 1] = (- _abcds[j + 1](1, 1) * complex_t(data(j + 1)) + _abcds[j + 1](0, 1)) /
         (_abcds[j + 1](1, 0) * complex_t(data(j + 1)) - _abcds[j + 1](0, 0));
    }
    return;
  }

  nda::array<double, 1> solver::evaluate(const nda::array<double, 1> &grid) {
    return nda::array<double, 1>(grid.shape());
  }

  solver::solver(nevanlinna_parameters_t const & p) {
    mpf_set_default_prec(p.precision);
  }

} // namespace nevanlinna
