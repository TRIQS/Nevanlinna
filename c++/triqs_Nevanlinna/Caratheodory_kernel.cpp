#include "Caratheodory_kernel.hpp"
#include "Nevanlinna_error.hpp"

namespace triqs_Nevanlinna {
  void Caratheodory_kernel::init(const nda::array<std::complex<double>, 1> &mesh, const nda::array<std::complex<double>, 3> &data) {
    _Ws.clear();
    _mesh.clear();
    _sqrt_one.clear();
    _sqrt_two.clear();
    _dim = data.shape()[1];
    _mesh.resize(mesh.shape()[0]);
    _Ws.resize(mesh.shape()[0]);
    _sqrt_one.resize(mesh.shape()[0]);
    _sqrt_two.resize(mesh.shape()[0]);
    complex_t One {1., 0.};
    complex_t I {0., 1.};
    matrix_t id (_dim, _dim);
    matrix_t val (_dim, _dim);
    id.setIdentity();
    for (int iw = mesh.shape()[0] - 1, w = 0; iw >= 0 ; --iw, ++w) {
      _mesh[w] = (complex_t(mesh(iw)) - I) / (complex_t(mesh(iw)) + I);
      _Ws[w].resize(_dim, _dim);
      for(int i =0; i< _dim; ++i) {
        for (int j = 0; j < _dim; ++j) {
          val(i, j) = complex_t(data(iw, i, j)) ;
        }
      }
      val = (id - I * val) * (id + I * val).inverse();
      _Ws[w] = val;
      /* ca_real norm_ave = _Ws[iw].norm()/_dim/_dim;
        if (norm_ave > 1) {_Ws[iw] /= norm_ave; std::cout << "hey ";} */
    }
    std::cout << std::endl;
    for (int i = mesh.shape()[0] - 1; i > 0; i--) {
      complex_t zi = _mesh[i];
      matrix_t Wi = _Ws[i];
      bool is_Schur_1 = true;
      matrix_t sqrt_one_i = sqrt_m(id - Wi * Wi.adjoint(), is_Schur_1);
      matrix_t sqrt_two_i = sqrt_m(id - Wi.adjoint() * Wi, is_Schur_1);
      // See Eq. 8 PhysRevB.104.165111
      for (int j = i - 1; j >= 0; j--) {
        complex_t zj = _mesh[j];
        matrix_t Wj = _Ws[j];
        complex_t y_ij = complex_t{std::abs(zi), 0.} * (zi - zj) / zi / (One - std::conj(zi) * zj);
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
  nda::array<std::complex<double>, 3> Caratheodory_kernel::evaluate(const nda::array<std::complex<double>, 1> &grid) const {
    if(_dim == 0) {
      throw Nevanlinna_uninitialized_error("Empty continuation data. Please run solve(...) first.");
    }
    std::vector<matrix_t> Vs(_mesh.size()); //intermediate Vs (for calculating Psis)
    std::vector<matrix_t> Fs(_mesh.size()); //intermediate Psis (Schur class functions)
    complex_t I {0., 1.};
    complex_t One {1., 0.};
    matrix_t id (_dim, _dim);
    nda::array<std::complex<double>, 3> results(grid.shape()[0], _dim, _dim);
    id.setIdentity();
    for (int i = 0; i < grid.shape()[0]; i++) {
      complex_t z = (complex_t(grid(i)) - I) / (complex_t(grid(i)) + I);
      complex_t z0 = _mesh[0];
      matrix_t W0 = _Ws[0];
      Vs[0] = complex_t{std::abs(z0), 0.} * (z0 - z) / z0 / (One - std::conj(z0) * z) * id;
      Fs[0] = (id + Vs[0] * W0.adjoint()).inverse() * (Vs[0] + W0);
      for (int j = 1; j < _mesh.size(); j++) {
        complex_t zj = _mesh[j];
        matrix_t Wj = _Ws[j];
        // See Eq. 9 PhysRevB.104.165111
        Vs[j] = complex_t{std::abs(zj), 0.} * (zj - z) / zj / (One - std::conj(zj) * z) * _sqrt_one[j] * Fs[j - 1] * _sqrt_two[j];
        // See Eq. 10 PhysRevB.104.165111
        Fs[j] = (id + Vs[j] * Wj.adjoint()).inverse() * (Vs[j] + Wj);
      }
      matrix_t val = Fs[_mesh.size() - 1];
      // See Eq. 11 PhysRevB.104.165111
      val = -I * (id + val).inverse() * (id - val);
      for (int n = 0; n < _dim; ++n) {
        for (int m = 0; m < _dim; ++m) {
          results(i, n, m) =std::complex<double>(val(n, m).real().convert_to<double>(), val(n, m).imag().convert_to<double>());
        }
      }
    }
    return results;
  }

}