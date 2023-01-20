#include "Caratheodory_kernel.hpp"
#include "Nevanlinna_error.hpp"

namespace triqs_Nevanlinna {
  void Caratheodory_kernel::init(nda::vector_const_view<std::complex<double>> mesh, nda::array_const_view<std::complex<double>, 3> data) {
    _Ws.clear();
    _mesh.clear();
    _sqrt_one.clear();
    _sqrt_two.clear();
    _dim = data.shape()[1];
    _mesh.resize(mesh.shape()[0]);
    _Ws.resize(mesh.shape()[0]);
    _sqrt_one.resize(mesh.shape()[0]);
    _sqrt_two.resize(mesh.shape()[0]);
    complex_mpt One {1., 0.};
    complex_mpt I {0., 1.};
    matrix_cplx_mpt id (_dim, _dim);
    matrix_cplx_mpt val (_dim, _dim);
    id.setIdentity();
    for (int iw = mesh.shape()[0] - 1, w = 0; iw >= 0 ; --iw, ++w) {
      _mesh[w] = (complex_mpt(mesh(iw)) - I) / (complex_mpt(mesh(iw)) + I);
      _Ws[w].resize(_dim, _dim);
      for(int i =0; i< _dim; ++i) {
        for (int j = 0; j < _dim; ++j) {
          val(i, j) = complex_mpt(data(iw, i, j)) ;
        }
      }
      val = (id - I * val) * (id + I * val).inverse();
      _Ws[w] = val;
      /* ca_real norm_ave = _Ws[iw].norm()/_dim/_dim;
        if (norm_ave > 1) {_Ws[iw] /= norm_ave; std::cout << "hey ";} */
    }
    std::cout << std::endl;
    for (int i = mesh.size() - 1; i > 0; i--) {
      auto &zi = _mesh[i];
      auto &Wi = _Ws[i];
      bool is_Schur_1 = true;
      auto sqrt_one_i = sqrt_m(id - Wi * Wi.adjoint(), is_Schur_1);
      auto sqrt_two_i = sqrt_m(id - Wi.adjoint() * Wi, is_Schur_1);
      // See Eq. 8 PhysRevB.104.165111
      for (int j = i - 1; j >= 0; j--) {
        auto &zj = _mesh[j];
        auto &Wj = _Ws[j];
        auto y_ij = complex_mpt{std::abs(zi), 0.} * (zi - zj) / zi / (One - std::conj(zi) * zj);
        _Ws[j] = sqrt_one_i.inverse() * (Wj - Wi) * (id - Wi.adjoint() * Wj).inverse().eval() * sqrt_two_i / y_ij;
//        _Ws[j] = sqrt_one_i.inverse() * (id - Wi* Wj.inverse() ) * (Wj.inverse() - Wi.adjoint()).inverse() * sqrt_two_i / y_ij;
      }
      _sqrt_one[i] = sqrt_one_i;
      _sqrt_two[i] = sqrt_two_i.inverse(); // original
                                          //_sqrt_two[i] = sqrt_two_i;
    }
    bool is_Schur = true;
    _sqrt_one[0] = sqrt_m(id - _Ws[0] * _Ws[0].adjoint(), is_Schur);
    _sqrt_two[0] = sqrt_m(id - _Ws[0].adjoint() * _Ws[0], is_Schur).inverse();
  }
  nda::array<std::complex<double>, 3> Caratheodory_kernel::evaluate(nda::vector_const_view<std::complex<double>> grid) const {
    if(_dim == 0) {
      throw Nevanlinna_uninitialized_error("Empty continuation data. Please run solve(...) first.");
    }
    std::vector<matrix_cplx_mpt> Vs(_mesh.size()); //intermediate Vs (for calculating Psis)
    std::vector<matrix_cplx_mpt> Fs(_mesh.size()); //intermediate Psis (Schur class functions)
    auto One = complex_mpt{1., 0.};
    auto I = complex_mpt{0., 1.};
    auto id = matrix_cplx_mpt::Identity(_dim, _dim);
    nda::array<std::complex<double>, 3> results(grid.shape()[0], _dim, _dim);
    id.setIdentity();
    for (int i = 0; i < grid.shape()[0]; i++) {
      auto z = (complex_mpt(grid(i)) - I) / (complex_mpt(grid(i)) + I);
      auto& z0 = _mesh[0];
      auto& W0 = _Ws[0];
      Vs[0] = complex_mpt{std::abs(z0), 0.} * (z0 - z) / z0 / (One - std::conj(z0) * z) * id;
      Fs[0] = (id + Vs[0] * W0.adjoint()).inverse() * (Vs[0] + W0);
      for (int j = 1; j < _mesh.size(); j++) {
        auto& zj = _mesh[j];
        auto& Wj = _Ws[j];
        // See Eq. 9 PhysRevB.104.165111
        Vs[j] = complex_mpt{std::abs(zj), 0.} * (zj - z) / zj / (One - std::conj(zj) * z) * _sqrt_one[j] * Fs[j - 1] * _sqrt_two[j];
        // See Eq. 10 PhysRevB.104.165111
        Fs[j] = (id + Vs[j] * Wj.adjoint()).inverse() * (Vs[j] + Wj);
      }
      // See Eq. 11 PhysRevB.104.165111
      auto val = matrix_cplx_mpt(-I * (id + Fs[_mesh.size() - 1]).inverse() * (id - Fs[_mesh.size() - 1]));
      for (int n = 0; n < _dim; ++n) {
        for (int m = 0; m < _dim; ++m) {
          results(i, n, m) =std::complex<double>(val(n, m).real().convert_to<double>(), val(n, m).imag().convert_to<double>());
        }
      }
    }
    return results;
  }

}
