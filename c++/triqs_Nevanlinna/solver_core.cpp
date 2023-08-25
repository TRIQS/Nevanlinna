#include <cmath>
#include "solver_core.hpp"
#include "kernel_factory.hpp"

namespace triqs_Nevanlinna {

  void solver_core::solve(triqs::gfs::gf_const_view<triqs::mesh::imfreq> g_iw) {
    auto mesh = nda::vector<std::complex<double>>(g_iw.mesh().size());
    for (auto i : range(g_iw.mesh().size())) mesh[i] = g_iw.mesh()[i];
    _kernel->init(mesh, g_iw.data());
  }

  triqs::gfs::gf<triqs::mesh::refreq> solver_core::evaluate(const triqs::mesh::refreq &grid, double eta) {
    return evaluate(grid, eta, nda::array_const_view<std::complex<double>, 3>());
  }

  [[nodiscard]] triqs::gfs::gf<triqs::mesh::refreq> solver_core::evaluate(const triqs::mesh::refreq &grid, double eta,
                                                                          nda::array_const_view<std::complex<double>, 3> theta) {
    triqs::gfs::gf<triqs::mesh::refreq> result(grid, {static_cast<long>(_kernel->size()), static_cast<long>(_kernel->size())});
    size_t nw = grid.size();
    nda::vector<std::complex<double>> mesh(nw);
    for (auto i : range(nw)) mesh[i] = grid[i] + 1i * eta;
    nda::array<std::complex<double>, 3> data = _kernel->evaluate(mesh, theta);
    for (size_t n = 0; n < _kernel->size(); ++n) {
      for (size_t m = 0; m < _kernel->size(); ++m) {
        for (int w = 0; w < nw; ++w) { result.data()(w, n, m) = data(w, n, m); }
      }
    }
    return result;
  }

  solver_core::solver_core(Nevanlinna_parameters_t const &p) : _kernel(kernel_factory::get_kernel(p)) {}

} // namespace triqs_Nevanlinna
