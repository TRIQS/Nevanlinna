#include "Caratheodory_kernel.hpp"
#include "Nevanlinna_error.hpp"

namespace triqs_Nevanlinna {
  void Caratheodory_kernel::init(nda::vector_const_view<std::complex<double>> mesh, nda::array_const_view<std::complex<double>, 3> data) {
    _dim      = data.shape()[1];
    size_t nw = std::count_if(mesh.begin(), mesh.end(), [](const std::complex<double> &w) { return w.imag() > 0; });
    _mesh.resize(nw);
    _Ws.resize(nw);
    _sqrt_one.resize(nw);
    _sqrt_two.resize(nw);
    auto id  = matrix_cplx_mpt::Identity(_dim, _dim);
    auto val = matrix_cplx_mpt(_dim, _dim);
    for (int iw = mesh.shape()[0] - 1, w = 0; iw >= 0; --iw) {
      if (mesh(iw).imag() < 0.0) { continue; }
      _mesh[w] = (complex_mpt(mesh(iw)) - I) / (complex_mpt(mesh(iw)) + I);
      _Ws[w].resize(_dim, _dim);
      for (int i = 0; i < _dim; ++i) {
        for (int j = 0; j < _dim; ++j) { val(i, j) = complex_mpt(data(iw, i, j)); }
      }
      val    = (id - I * val) * (id + I * val).inverse();
      _Ws[w] = val;
      ++w;
    }
    _data = _Ws;
    for (int i = _mesh.size() - 1; i > 0; i--) {
      auto &zi        = _mesh[i];
      auto &Wi        = _Ws[i];
      bool is_Schur_1 = true;
      auto sqrt_one_i = sqrt_m(id - Wi * Wi.adjoint(), is_Schur_1);
      auto sqrt_two_i = sqrt_m(id - Wi.adjoint() * Wi, is_Schur_1);
      // See Eq. 8 PhysRevB.104.165111
      for (int j = i - 1; j >= 0; j--) {
        auto &zj  = _mesh[j];
        auto &Wj  = _Ws[j];
        auto y_ij = complex_mpt{std::abs(zi), 0.} * (zi - zj) / zi / (One - std::conj(zi) * zj);
        _Ws[j]    = sqrt_one_i.inverse() * (Wj - Wi) * (id - Wi.adjoint() * Wj).inverse().eval() * sqrt_two_i / y_ij;
        //        _Ws[j] = sqrt_one_i.inverse() * (id - Wi* Wj.inverse() ) * (Wj.inverse() - Wi.adjoint()).inverse() * sqrt_two_i / y_ij;
      }
      _sqrt_one[i] = sqrt_one_i;
      _sqrt_two[i] = sqrt_two_i.inverse(); // original
                                           //_sqrt_two[i] = sqrt_two_i;
    }
    bool is_Schur = true;
    _sqrt_one[0]  = sqrt_m(id - _Ws[0] * _Ws[0].adjoint(), is_Schur);
    _sqrt_two[0]  = sqrt_m(id - _Ws[0].adjoint() * _Ws[0], is_Schur).inverse();
  }

  nda::array<std::complex<double>, 3> Caratheodory_kernel::evaluate(nda::vector_const_view<std::complex<double>> grid) {
    if (_dim == 0) { throw Nevanlinna_uninitialized_error("Empty continuation data. Please run solve(...) first."); }
    std::vector<matrix_cplx_mpt> Vs(_mesh.size()); //intermediate Vs (for calculating Psis)
    std::vector<matrix_cplx_mpt> Fs(_mesh.size()); //intermediate Psis (Schur class functions)
    auto id  = matrix_cplx_mpt::Identity(_dim, _dim);
    nda::array<std::complex<double>, 3> results(grid.shape()[0], _dim, _dim);
    for (int i = 0; i < grid.shape()[0]; i++) {
      auto z   = (complex_mpt(grid(i)) - I) / (complex_mpt(grid(i)) + I);
      auto &z0 = _mesh[0];
      auto &W0 = _Ws[0];
      Vs[0]    = complex_mpt{std::abs(z0), 0.} * (z0 - z) / z0 / (One - std::conj(z0) * z) * id;
      Fs[0]    = (id + Vs[0] * W0.adjoint()).inverse() * (Vs[0] + W0);
      for (int j = 1; j < _mesh.size(); j++) {
        auto &zj = _mesh[j];
        auto &Wj = _Ws[j];
        // See Eq. 9 PhysRevB.104.165111
        Vs[j] = complex_mpt{std::abs(zj), 0.} * (zj - z) / zj / (One - std::conj(zj) * z) * _sqrt_one[j] * Fs[j - 1] * _sqrt_two[j];
        // See Eq. 10 PhysRevB.104.165111
        Fs[j] = (id + Vs[j] * Wj.adjoint()).inverse() * (Vs[j] + Wj);
      }
      // See Eq. 11 PhysRevB.104.165111
      auto val = matrix_cplx_mpt(-I * (id + Fs[_mesh.size() - 1]).inverse() * (id - Fs[_mesh.size() - 1]));
      for (int n = 0; n < _dim; ++n) {
        for (int m = 0; m < _dim; ++m) {
          results(i, n, m) = std::complex<double>(val(n, m).real().convert_to<double>(), val(n, m).imag().convert_to<double>());
        }
      }
    }
    return results;
  }

  nda::vector<double> Caratheodory_kernel::get_Pick_eigenvalues() const {
    auto nw = _data.shape(0);
    if (nw == 0) { return nda::vector<double>(); }
    auto N    = _data(0).cols();
    auto Pick = Eigen::MatrixXcd(nw * N, nw * N);
    auto id   = matrix_cplx_mpt::Identity(N, N);
    for (int i = 0; i < nw; i++) {
      for (int j = 0; j < nw; j++) {
        auto val = (id - _data(i).adjoint() * _data(j)) / (One - std::conj(_mesh(i)) * _mesh(j));
        for (int n = 0; n < N; ++n) {
          for (int m = 0; m < N; ++m) {
            Pick.block(i * N, j * N, N, N)(n, m) = std::complex<double>(val(n, m).real().convert_to<double>(), val(n, m).imag().convert_to<double>());
          }
        }
      }
    }
    auto eigenvalues      = Pick.eigenvalues();
    auto pick_eigenvalues = nda::vector<double>(nw * N);
    std::transform(eigenvalues.begin(), eigenvalues.end(), pick_eigenvalues.begin(), [](const std::complex<double> &r) { return r.real(); });
    return pick_eigenvalues;
  }

} // namespace triqs_Nevanlinna
