#pragma once
#include <memory>
#include <triqs/gfs.hpp>
#include <triqs/mesh.hpp>
#include <triqs/utility/macros.hpp>
#include <nda/nda.hpp>
#include <h5/h5.hpp>
#include "Nevanlinna_parameters_t.hpp"
#include "kernel.hpp"

namespace triqs_Nevanlinna {

  /**
   * @brief Nevanlinna analytical continuation solver for TRIQS Green's functions.
   *
   * @details Performs analytical continuation for the diagonal part of the matrix-values 
   * TRIQS Green's function.
   */
  class solver_core {

    public:
    ~solver_core() = default;

    /**
     * @brief Construct the solver.
     * @param p Construction parameters (kernel choice and multiprecision precision).
     */
    solver_core(Nevanlinna_parameters_t const &p);

    // Copy/Move construction
    solver_core(solver_core const &) = delete;
    solver_core(solver_core &&)      = default;

    /// Copy/Move assignment
    solver_core &operator=(solver_core const &) = delete;
    solver_core &operator=(solver_core &&)      = default;

    /**
     * @brief Perform a Nevanlinna factorization for a matrix-valued Matsubara frequency Green's function.
     *
     * @param g_iw Matrix-valued Matsubara frequency Green's function.
     */
    void solve(triqs::gfs::gf_const_view<triqs::mesh::imfreq> g_iw);

    /**
     * @brief Evaluate diagonal part of the real-frequency Green's function on a chosen grid.
     *
     * @details Uses the precomputed Nevanlinna factorization.
     *
     * @param grid Real frequency grid.
     * @param eta Lorentzian broadening.
     * @return Real-frequency matrix-valued TRIQS Green's function on a chosen grid.
     */
    [[nodiscard]] triqs::gfs::gf<triqs::mesh::refreq> evaluate(const triqs::mesh::refreq &grid, double eta);

    /**
     * @brief Evaluate the real-frequency Green's function on a chosen grid using Hardy-function optimization.
     *
     * @param grid Real frequency grid.
     * @param eta Lorentzian broadening.
     * @param theta Hardy-function basis coefficients used to optimize the spectral function.
     * @return Real-frequency matrix-valued TRIQS Green's function on a chosen grid.
     */
    [[nodiscard]] triqs::gfs::gf<triqs::mesh::refreq> evaluate(const triqs::mesh::refreq &grid, double eta,
                                                               nda::array_const_view<std::complex<double>, 3> theta);

    /// Eigenvalues of the Pick matrix (non-negative eigenvalues indicate the data is continuable).
    [[nodiscard]] C2PY_PROPERTY_GET(Pick_eigenvalues) nda::vector<double> get_Pick_eigenvalues() const { return _kernel->get_Pick_eigenvalues(); };

    /// Number of orbitals (matrix dimension) of the continued Green's function.
    [[nodiscard]] C2PY_PROPERTY_GET(size) size_t size() const { return _kernel->size(); };

    private:
    // vector of Nevanlinna factorization kernels for multi-orbital factorization
    std::unique_ptr<kernel> _kernel{};
  };
} // namespace triqs_Nevanlinna
