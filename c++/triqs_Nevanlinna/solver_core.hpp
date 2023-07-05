#pragma once
#include <memory>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <nda/nda.hpp>
#include <h5/h5.hpp>
#include "Nevanlinna_parameters_t.hpp"
#include "kernel.hpp"

namespace triqs_Nevanlinna {

  /**
   * Nevanlinna analytical continuation solver_core for TRIQS GFs
   *
   * @note Perform analytical continuation for the diagonal part of the matrix-values TRIQS Green's function
   * @include triqs_Nevanlinna/Nevanlinna.hpp
   */
  class solver_core {

    public:
    ~solver_core() = default;

    CPP2PY_ARG_AS_DICT
    solver_core(Nevanlinna_parameters_t const &p);

    // Copy/Move construction
    solver_core(solver_core const &) = delete;
    solver_core(solver_core &&)      = default;

    /// Copy/Move assignment
    solver_core &operator=(solver_core const &) = delete;
    solver_core &operator=(solver_core &&)      = default;

    /**
     * Construct a Nevanlinna factorization for matrix-valued Matsubara frequency Green's function
     *
     * @param g_iw - matrix-valued Matsubara frequency Green's function
     */
    void solve(triqs::gfs::gf_const_view<triqs::mesh::imfreq> g_iw);

    /**
     * Evaluate diagonal part of the real-frequency Green's function on a chosen grid
     * based on the precomputed Nevanlinna factorization
     *
     * @param grid - real frequency grid
     * @param eta - Lorentzian broadening
     * @return Real-frequency matrix-valued TRIQS Green's function on a chosen grid.
     */
    [[nodiscard]] triqs::gfs::gf<triqs::mesh::refreq> evaluate(const triqs::mesh::refreq &grid, double eta);

    [[nodiscard]] triqs::gfs::gf<triqs::mesh::refreq> evaluate(const triqs::mesh::refreq &grid, double eta,
                                                               nda::array_const_view<std::complex<double>, 3> theta);

    [[nodiscard]] nda::vector<double> get_Pick_eigenvalues() const { return _kernel->get_Pick_eigenvalues(); };

    [[nodiscard]] size_t size() const { return _kernel->size(); };

    private:
    // vector of Nevanlinna factorization kernels for multi-orbital factorization
    std::unique_ptr<kernel> _kernel{};
  };
} // namespace triqs_Nevanlinna
