#include "Nevanlinna_factorization.hpp"

using namespace std::complex_literals;

namespace triqs_Nevanlinna {

  nda::vector<complex_mpt> Nevanlinna_factorization::mobius_trasformation(nda::vector_const_view<std::complex<double>> data) const {
    nda::vector<complex_mpt> mdata(data.size());
    std::transform(data.begin(), data.end(), mdata.begin(),
                   [](const std::complex<double> &d) { return complex_mpt(-d - 1i) / complex_mpt(-d + 1i); });
    return mdata;
  }

  void Nevanlinna_factorization::build(nda::vector_const_view<std::complex<double>> mesh, nda::vector_const_view<std::complex<double>> data) {
    assert(mesh.size() == data.size());
    if (std::any_of(mesh.begin(), mesh.end(), [](const std::complex<double> &v) { return v.real() != 0.0 or v.imag() < 0; })) {
      throw Nevanlinna_negative_grid_error("Data should be defined on the positive Matsubara frequencies.");
    }
    size_t M = mesh.size();
    _phis.resize(M);
    _abcds.resize(M);
    _mesh.resize(M);
    // invalidate previous real frequency grid
    _grid.resize(0);
    _data    = mobius_trasformation(data);
    _phis[0] = _data[0];
    for (int k = 0; k < M; k++) {
      _abcds[k] = matrix_cplx_mpt::Identity(2, 2);
      _mesh[k]  = mesh(k);
    }
    auto prod = matrix_cplx_mpt(2, 2);
    for (int j = 0; j < M - 1; j++) {
      for (int k = j; k < M; k++) {
        prod << (_mesh[k] - _mesh[j]) / (_mesh[k] - std::conj(_mesh[j])), _phis[j],
           std::conj(_phis[j]) * (_mesh[k] - _mesh[j]) / (_mesh[k] - std::conj(_mesh[j])), complex_mpt{1., 0.};
        _abcds[k] *= prod;
      }
      _phis[j + 1] = (-_abcds[j + 1](1, 1) * _data[j + 1] + _abcds[j + 1](0, 1)) / (_abcds[j + 1](1, 0) * _data[j + 1] - _abcds[j + 1](0, 0));
    }
    return;
  }

  nda::vector<double> Nevanlinna_factorization::evaluate(nda::vector_const_view<double> grid, double eta) {
    auto complex_grid                     = make_regular(grid + eta * 1i);
    nda::vector<std::complex<double>> G_w = evaluate(complex_grid, nda::vector_const_view<std::complex<double>>());
    nda::vector<double> A_w(G_w.shape());
    std::transform(G_w.begin(), G_w.end(), A_w.begin(), [](const std::complex<double> &v) { return -v.imag() / M_PI; });
    return A_w;
  }

  nda::vector<std::complex<double>> Nevanlinna_factorization::evaluate(nda::vector_const_view<std::complex<double>> grid,
                                                                       nda::vector_const_view<std::complex<double>> theta_M_1) {
    size_t M = _phis.size();
    if (M == 0) { throw Nevanlinna_uninitialized_error("Empty continuation data. Please run solve(...) first."); }
    if (theta_M_1.size() != grid.size() && theta_M_1.size() != 0) {
      throw Nevanlinna_error("theta_{M+1} should either have a value at every frequency point or be empty.");
    }
    // check if grid has not been changed
    if (grid.size() == _grid.size()
        && std::equal(grid.begin(), grid.end(), _grid.begin(), [](const std::complex<double> &w1, const complex_mpt &w2) {
             return std::abs(w1.real() - w2.real().convert_to<double>()) < 1e-9 && std::abs(w1.imag() - w2.imag().convert_to<double>()) < 1e-9;
           })) {
      return evaluate_for_theta(grid, theta_M_1);
    }
    _coeffs = std::vector<matrix_cplx_mpt>(grid.size());
    _grid.resize(grid.size());
    std::transform(grid.begin(), grid.end(), _grid.begin(), [](const std::complex<double> &w) { return complex_mpt{w.real(), w.imag()}; });
    auto I       = complex_mpt{0., 1.};
    auto One     = complex_mpt{1., 0.};
    auto results = nda::vector<std::complex<double>>(grid.size());
    auto prod    = matrix_cplx_mpt(2, 2);
    for (int i = 0; i < grid.size(); ++i) {
      matrix_cplx_mpt result = matrix_cplx_mpt::Identity(2, 2);
      auto z                 = _grid[i];
      for (int j = 0; j < M; j++) {
        prod << (z - _mesh[j]) / (z - std::conj(_mesh[j])), _phis[j], std::conj(_phis[j]) * ((z - _mesh[j]) / (z - std::conj(_mesh[j]))), One;
        result *= prod;
      }
      _coeffs[i] = result;
      auto param = complex_mpt{0., 0.}; //theta_{M+1}, choose to be constant function 0 here
      auto theta = complex_mpt(result(0, 0) * param + result(0, 1)) / (result(1, 0) * param + result(1, 1));
      auto value = complex_mpt(I * (One + theta) / (One - theta));
      results(i) =
         -std::complex<double>(value.real().convert_to<double>(), value.imag().convert_to<double>()); //inverse Mobius transform from theta to NG
    }
    return results;
  }

  nda::vector<std::complex<double>> Nevanlinna_factorization::evaluate_for_theta(nda::vector_const_view<std::complex<double>> grid,
                                                                                 nda::vector_const_view<std::complex<double>> theta_M_1) const {
    size_t M = _phis.size();
    if (M == 0) { throw Nevanlinna_uninitialized_error("Empty continuation data. Please run solve(...) first."); }
    if (theta_M_1.size() != grid.size() || theta_M_1.size() != 0) {
      throw Nevanlinna_error("theta_{M+1} should either have a value at every frequency point or be empty.");
    }
    auto I   = complex_mpt{0., 1.};
    auto One = complex_mpt{1., 0.};
    nda::vector<std::complex<double>> results(grid.shape());
    for (int i = 0; i < _grid.size(); ++i) {
      auto result         = _coeffs[i];
      auto theta_M_plus_1 = theta_M_1.size() == 0 ? complex_mpt{0., 0.} : complex_mpt{theta_M_1(i).real(), theta_M_1(i).imag()};

      auto theta = (result(0, 0) * theta_M_plus_1 + result(0, 1)) / (result(1, 0) * theta_M_plus_1 + result(1, 1));
      auto value = I * (One + theta) / (One - theta);
      //inverse Mobius transform from theta to NG
      results(i) = -std::complex<double>(value.real().convert_to<double>(), value.imag().convert_to<double>());
    }
    return results;
  }

} // namespace triqs_Nevanlinna
