
#include "Nevanlinna_kernel.hpp"

namespace triqs_Nevanlinna {

  void Nevanlinna_kernel::init(const nda::array<std::complex<double>, 1> &mesh, const nda::array<std::complex<double>, 3> &data) {
    _factorizations.clear();
    size_t N  = data.shape()[1];
    size_t nw = std::count_if(mesh.begin(), mesh.end(), [](const std::complex<double> &w) { return w.imag() > 0; });
    nda::array<std::complex<double>, 1> data_in(nw);
    nda::array<std::complex<double>, 1> mesh_in(nw);
    for (size_t n = 0; n < N; ++n) {
      for (size_t iw = 0, iww = 0; iw < mesh.shape()[0]; ++iw) {
        if (mesh(iw).imag() < 0) continue;
        data_in(iww) = data(iw, n, n);
        mesh_in(iww) = mesh(iw);
        ++iww;
      }
      Nevanlinna_factorization f;
      f.build(mesh_in, data_in);
      _factorizations.push_back(f);
    }
  }
  nda::array<std::complex<double>, 3> Nevanlinna_kernel::evaluate(const nda::array<std::complex<double>, 1> &grid) const {
    nda::array<std::complex<double>, 3> results(grid.shape()[0], size(), size());
    for (size_t n = 0; n < size(); ++n) {
      nda::array<std::complex<double>, 1> data = _factorizations[n].evaluate(grid);
      for (size_t iw = 0; iw < grid.shape()[0]; ++iw) { results(iw, n, n) = data(iw); }
    }
    return results;
  }
}
