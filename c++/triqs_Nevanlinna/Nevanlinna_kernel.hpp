#ifndef TRIQS_NEVANLINNA_NEVANLINNA_KERNEL_HPP
#define TRIQS_NEVANLINNA_NEVANLINNA_KERNEL_HPP

#include <complex>

#include <nda/nda.hpp>
#include <nda/mpi.hpp>
#include <triqs/utility/macros.hpp>

#include "kernel.hpp"
#include "Nevanlinna_factorization.hpp"

namespace triqs_Nevanlinna {

  /**
   * @brief Diagonal Nevanlinna continuation kernel.
   *
   * @details Builds one independent Nevanlinna factorization per orbital; off-diagonal elements of the
   * input Green's function are ignored. Supports Hardy-function (theta) optimization in ``evaluate()``.
   */
  class Nevanlinna_kernel : public kernel {

    public:
    /**
     * @brief Construct the diagonal Nevanlinna kernel.
     *
     * @param precision Number of decimal digits of internal multiprecision arithmetic (only honored with MPFR support).
     */
    Nevanlinna_kernel(int precision = mp_digits) : kernel(precision) {
      if (!mpi::communicator().rank())
        std::cerr << "This is Nevanlinna analytical continuation. All off-diagonal elements will be ignored." << std::endl;
    }

    /**
     * @brief Initialize the diagonal Nevanlinna continuation from Matsubara-frequency input data.
     *
     * @param mesh Positive Matsubara frequencies, given as complex values.
     * @param data Matrix-valued Green's function data on those frequencies.
     */
    void init(nda::vector_const_view<std::complex<double>> mesh, nda::array_const_view<std::complex<double>, 3> data) override;

    /**
     * @brief Evaluate the diagonal real-frequency Green's function on a chosen grid.
     *
     * @param grid Complex real-frequency grid (real frequency plus Lorentzian broadening).
     * @return Matrix-valued real-frequency Green's function on the grid.
     */
    [[nodiscard]] nda::array<std::complex<double>, 3> evaluate(nda::vector_const_view<std::complex<double>> grid) override;

    /**
     * @brief Evaluate the diagonal real-frequency Green's function using Hardy-function optimization.
     *
     * @param grid Complex real-frequency grid (real frequency plus Lorentzian broadening).
     * @param theta Hardy-function basis coefficients used to optimize the spectral function.
     * @return Matrix-valued real-frequency Green's function on the grid.
     */
    [[nodiscard]] nda::array<std::complex<double>, 3> evaluate(nda::vector_const_view<std::complex<double>> grid,
                                                               nda::array_const_view<std::complex<double>, 3> theta) override;

    /// Number of orbitals (matrix dimension) handled by the kernel.
    [[nodiscard]] C2PY_PROPERTY_GET(size) size_t size() const override { return _factorizations.size(); }

    /// Eigenvalues of the Pick matrix; non-negative eigenvalues indicate the data is continuable (Nevanlinna).
    [[nodiscard]] C2PY_PROPERTY_GET(Pick_eigenvalues) nda::vector<double> get_Pick_eigenvalues() const override;

    private:
    size_t _N_im_freq{};
    std::vector<Nevanlinna_factorization> _factorizations{};
  };

} // namespace triqs_Nevanlinna
#endif //TRIQS_NEVANLINNA_NEVANLINNA_KERNEL_HPP
