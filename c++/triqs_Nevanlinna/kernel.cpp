//
// Created by iskakoff on 10/12/22.
//

#include "kernel.hpp"


namespace triqs_Nevanlinna {

  void kernel::solve(const nda::array<std::complex<double>, 1> & mesh, const nda::array<std::complex<double>, 1> & data) {
    assert(mesh.shape(0) == data.shape(0));
    if(std::any_of(mesh.begin(), mesh.end(), [](const std::complex<double> & v) {return v.imag()<0;})) {
      throw Nevanlinna_negative_grid_error("Data should be defined on the positive Matsubara frequencies.");
    }
    size_t M = mesh.shape(0);
    _phis.resize(M);
    _abcds.resize(M);
    _mesh.resize(M);
    auto mdata = mobius_trasformation(data);
    _phis[0] = mdata[0];
    for (int k = 0; k < M; k++) {
      _abcds[k] = matrix_t::Identity(2, 2);
      _mesh[k] = mesh(k);
    }
    matrix_t prod(2, 2);
    for (int j = 0; j < M - 1; j++) {
      for (int k = j; k < M; k++) {
        prod <<  (_mesh[k] - _mesh[j]) / (_mesh[k] - std::conj(_mesh[j])),  _phis[j],
           std::conj(_phis[j]) * (_mesh[k] - _mesh[j]) / (_mesh[k] - std::conj(_mesh[j])), complex_t{1., 0.};
        _abcds[k] *= prod;
      }
      _phis[j + 1] = (- _abcds[j + 1](1, 1) * mdata[j+1] + _abcds[j + 1](0, 1)) /
         (_abcds[j + 1](1, 0) * mdata[j+1] - _abcds[j + 1](0, 0));
    }
    return;
  }

  nda::array<double, 1> kernel::evaluate(const nda::array<double, 1> &grid, double eta) const {
    nda::array<std::complex<double>, 1> complex_grid(grid.shape());
    complex_grid = grid + eta*1i;
    nda::array<std::complex<double>, 1> G_w = evaluate(complex_grid);
    nda::array<double, 1> A_w(G_w.shape());
    std::transform(G_w.begin(), G_w.end(), A_w.begin(), [] (const std::complex<double> & v) {return -v.imag()/M_PI;} );
    return A_w;
  }

  nda::array<std::complex<double>, 1> kernel::evaluate(const nda::array<std::complex<double>, 1> &grid) const {
    size_t M = _phis.size();
    if(M == 0) {
      throw Nevanlinna_error("Empty continuation data. Please run solve(...) first.");
    }
    complex_t I {0., 1.};
    complex_t One {1., 0.};
    nda::array<std::complex<double>, 1> results(grid.shape());
    matrix_t prod(2, 2);
    for (int i = 0; i< grid.shape(0); ++i) {
      matrix_t result = matrix_t::Identity(2, 2);
      complex_t z = complex_t(grid(i));
      for (int j = 0; j < M; j++) {
        prod <<  (z - _mesh[j]) / (z - std::conj(_mesh[j])), _phis[j] ,
           std::conj(_phis[j])* ((z - _mesh[j]) / (z - std::conj(_mesh[j]))), complex_t{1., 0.};
        result *= prod;
      }
      complex_t param {0., 0.}; //theta_{M+1}, choose to be constant function 0 here
      complex_t theta = (result(0, 0) * param + result(0, 1)) / (result(1, 0) * param + result(1, 1));
      complex_t value = I * (One + theta) / (One - theta);
      results(i) = -std::complex<double>(value.real().get_d(), value.imag().get_d()); //inverse Mobius transform from theta to NG
    }
    return results;
  }

  std::vector<kernel::complex_t> kernel::mobius_trasformation(const nda::array<std::complex<double>, 1> & data) const {
    std::vector<complex_t> mdata(data.shape(0));
    std::transform(data.begin(), data.end(), mdata.begin(), [](const std::complex<double>&d) {return complex_t(-d - 1i) / complex_t(-d + 1i);});
    return mdata;
  }



}
