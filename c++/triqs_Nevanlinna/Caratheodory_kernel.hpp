#ifndef TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP
#define TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP

#include "types.hpp"
#include "kernel.hpp"

#include <triqs/utility/macros.hpp>

namespace triqs_Nevanlinna {

  /**
   * @brief Full matrix-valued Caratheodory continuation kernel (PhysRevB.104.165111).
   *
   * @details Continues the complete matrix-valued Green's function, including off-diagonal elements.
   * The Hardy-function (theta) optimization path is not implemented for this kernel.
   */
  class Caratheodory_kernel : public kernel {
    static constexpr double tol = 1e-12;

    public:
    /**
     * @brief Construct the full matrix-valued Caratheodory kernel.
     *
     * @param precision Number of decimal digits of internal multiprecision arithmetic (only honored with MPFR support).
     */
    Caratheodory_kernel(int precision = mp_digits) : kernel(precision) {}

    /**
     * @brief Build the full matrix-valued Caratheodory continuation from Matsubara-frequency input data.
     *
     * @param mesh Positive Matsubara frequencies, given as complex values.
     * @param data Matrix-valued Green's function data on those frequencies.
     */
    void init(nda::vector_const_view<std::complex<double>> mesh, nda::array_const_view<std::complex<double>, 3> data) override;

    /**
     * @brief Evaluate the full matrix-valued real-frequency Green's function on a chosen grid.
     *
     * @param grid Complex real-frequency grid (real frequency plus Lorentzian broadening).
     * @return Matrix-valued real-frequency Green's function on the grid.
     */
    [[nodiscard]] nda::array<std::complex<double>, 3> evaluate(nda::vector_const_view<std::complex<double>> grid) override;

    /**
     * @brief Evaluate the full matrix-valued real-frequency Green's function on a chosen grid.
     * 
     * @details Theta optimization is not supported here. A non-empty ``theta`` is ignored and plain ``evaluate(grid)`` is used.
     * 
     * @param grid Complex real-frequency grid (real frequency plus Lorentzian broadening).
     * @param theta Hardy-function basis coefficients used to optimize the spectral function (ignored).
     * @return Matrix-valued real-frequency Green's function on the grid.
     */
    [[nodiscard]] nda::array<std::complex<double>, 3> evaluate(nda::vector_const_view<std::complex<double>> grid,
                                                               nda::array_const_view<std::complex<double>, 3> theta) override {
      if (theta.shape()[0] != 0) {
        std::cerr << "Continuation poles optimization has not been implemented in matrix-valued continuation yet." << std::endl;
      }
      return evaluate(grid);
    };

    /// Number of orbitals (matrix dimension) handled by the kernel.
    [[nodiscard]] C2PY_PROPERTY_GET(size) size_t size() const override { return _dim; }

    /// Eigenvalues of the Pick matrix; non-negative eigenvalues indicate the data is continuable (Nevanlinna).
    [[nodiscard]] C2PY_PROPERTY_GET(Pick_eigenvalues) nda::vector<double> get_Pick_eigenvalues() const override;

    private:
    int _dim = 0;
    nda::vector<complex_mpt> _mesh{};
    nda::vector<matrix_cplx_mpt> _data{}; //W_is
    nda::vector<matrix_cplx_mpt> _Ws{};   //W_is
    // See Eq. 6 PhysRevB.104.165111
    nda::vector<matrix_cplx_mpt> _sqrt_one{}; //[1 - W_i * W_i^dagger]^0.5
    nda::vector<matrix_cplx_mpt> _sqrt_two{}; //[1 - W_i^dagger * W_i]^-0.5
    nda::vector<double> _Pick_eigenvalues{};

    /**
     * @brief Calculate the Hermitian square root of matrix \f$ M \f$.
     *
     * @param M Matrix whose square root is computed.
     * @param is_Schur Set to true if the matrix is a Schur matrix.
     * @return Square root of the matrix.
     */
    matrix_cplx_mpt sqrt_m(const matrix_cplx_mpt &M, bool &is_Schur) {
      Eigen::ComplexEigenSolver<matrix_cplx_mpt> ces;
      ces.compute(M);
      matrix_cplx_mpt D = ces.eigenvalues();
      is_Schur          = true;
      for (int i = 0; i < D.rows(); i++) {
        if (D(i, 0).real() < tol) { is_Schur = false; }
      }
      return ces.eigenvectors() * D.array().sqrt().matrix().asDiagonal() * ces.eigenvectors().inverse().eval();
    }
  };
} // namespace triqs_Nevanlinna
#endif //TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP
