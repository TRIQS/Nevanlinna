#include "Nevanlinna_kernel.hpp"

namespace triqs_Nevanlinna {

  void Nevanlinna_kernel::init(nda::vector_const_view<std::complex<double>> mesh, nda::array_const_view<std::complex<double>, 3> data) {
    _factorizations.clear();
    size_t N  = data.shape()[1];
    size_t nw = std::count_if(mesh.begin(), mesh.end(), [](const std::complex<double> &w) { return w.imag() > 0; });
    nda::vector<std::complex<double>> data_in(nw);
    nda::vector<std::complex<double>> mesh_in(nw);
    _N_im_freq = nw;
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

  nda::array<std::complex<double>, 3> Nevanlinna_kernel::evaluate(nda::vector_const_view<std::complex<double>> grid) {
    return evaluate(grid, nda::array_const_view<std::complex<double>, 3>());
  }

  nda::array<std::complex<double>, 3> Nevanlinna_kernel::evaluate(nda::vector_const_view<std::complex<double>> grid,
                                                                  nda::array_const_view<std::complex<double>, 3> theta) {
    nda::array<std::complex<double>, 3> results(grid.shape()[0], size(), size());
    nda::vector<std::complex<double>> theta_(theta.shape()[0]);
    for (size_t n = 0; n < size(); ++n) {
      for (size_t iwm = 0; iwm < theta.shape()[0]; ++iwm) { theta_(iwm) = theta(iwm, n, n); }
      nda::vector<std::complex<double>> data = _factorizations[n].evaluate(grid, theta_);
      for (size_t iw = 0; iw < grid.shape()[0]; ++iw) { results(iw, n, n) = data(iw); }
    }
    return results;
  }
  nda::vector<double> Nevanlinna_kernel::pick_eigenvalues() const {
    if (size() == 0) { return {}; }
    nda::vector<double> pick_eigenvalues(_N_im_freq * size());

    for (size_t n = 0; n < size(); ++n) {
      nda::vector<double> eigenvalues = _factorizations[n].get_pick_eigenvalues();
      std::copy(eigenvalues.begin(), eigenvalues.end(), pick_eigenvalues.begin() + _N_im_freq * n);
    }
    std::sort(pick_eigenvalues.begin(), pick_eigenvalues.end());
    return pick_eigenvalues;
  }
} // namespace triqs_Nevanlinna
