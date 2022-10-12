#include <cmath>
#include <gmpxx.h>
#include "./solver.hpp"
#include "nevanlinna_error.hpp"

namespace nevanlinna {

  void solver::solve(const triqs::gfs::gf<triqs::mesh::imfreq>& g_iw) {
    _kernels.clear();
    size_t nw = g_iw.mesh().positive_only() ? g_iw.mesh().size() : g_iw.mesh().size()/2;
    size_t N = g_iw.data_shape()[1];
    nda::array<std::complex<double>, 1> data(nw);
    nda::array<std::complex<double>, 1> mesh(nw);
    for(int orb = 0; orb < N; ++orb) {
      kernel k;
      int i = 0, j = 0;
      for (auto pt = g_iw.mesh().begin(); pt != g_iw.mesh().end(); ++pt, ++j) {
        if (pt.to_point().imag() < 0) continue;
        data(i) = g_iw.data()(j, orb, orb);
        mesh(i) = pt.to_point();
        ++i;
      }
      k.solve(mesh, data);
      _kernels.push_back(k);
    }
  }

  triqs::gfs::gf<triqs::mesh::refreq> solver::evaluate(const triqs::mesh::refreq & grid, double eta) const {
    triqs::gfs::gf<triqs::mesh::refreq> result(grid, {static_cast<long>(_kernels.size()), static_cast<long>(_kernels.size())});
    size_t nw = grid.size();
    nda::array<std::complex<double>, 1> mesh(nw);
    std::transform(grid.begin(), grid.end(), mesh.begin(), [&eta] (double v) {return v + 1i*eta;});
    for(size_t n = 0; n< _kernels.size() ; ++n) {
      nda::array<std::complex<double>, 1> data = _kernels[n].evaluate(mesh);
      for(int w = 0; w < nw; ++w) {result.data()(w, n, n) = data(w);}
    }
    return result;
  }

  solver::solver(nevanlinna_parameters_t const & p) {
    mpf_set_default_prec(p.precision);
  }

} // namespace nevanlinna
