#include <cmath>
#include "solver.hpp"
#include "kernel_factory.hpp"

namespace triqs_Nevanlinna {

  void solver::solve(triqs::gfs::gf_const_view<triqs::mesh::imfreq> g_iw) {
    nda::vector<std::complex<double>> mesh(g_iw.mesh().size());
    std::transform(g_iw.mesh().begin(), g_iw.mesh().end(), mesh.begin(), [](const triqs::mesh::imfreq::domain_pt_t &w) {return w;});
    _kernel->init(mesh, g_iw.data());
  }

  triqs::gfs::gf<triqs::mesh::refreq> solver::evaluate(const triqs::mesh::refreq &grid, double eta) const {
    triqs::gfs::gf<triqs::mesh::refreq> result(grid, {static_cast<long>(_kernel->size()), static_cast<long>(_kernel->size())});
    size_t nw = grid.size();
    nda::vector<std::complex<double>> mesh(nw);
    std::transform(grid.begin(), grid.end(), mesh.begin(), [&eta](double v) { return v + 1i * eta; });
    nda::array<std::complex<double>, 3> data = _kernel->evaluate(mesh);
    for (size_t n = 0; n < _kernel->size(); ++n) {
      for (size_t m = 0; m < _kernel->size(); ++m) {
        for (int w = 0; w < nw; ++w) {
          result.data()(w, n, m) = data(w, n, m);
        }
      }
    }
    return result;
  }

  solver::solver(Nevanlinna_parameters_t const & p) : _kernel(kernel_factory::get_kernel(p.kernel)) {}

} // namespace triqs_Nevanlinna
