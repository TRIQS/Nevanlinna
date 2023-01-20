#ifndef TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP
#define TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP

#include "types.hpp"
#include "kernel.hpp"

namespace triqs_Nevanlinna {

  class Caratheodory_kernel : public kernel {
    static constexpr double tol = 1e-12;

    public:
    void init(nda::vector_const_view<std::complex<double>> mesh, nda::array_const_view<std::complex<double>, 3> data) override;
    [[nodiscard]] nda::array<std::complex<double>, 3> evaluate(nda::vector_const_view<std::complex<double>> grid) const override;

    [[nodiscard]] size_t size() const override {
      return _dim;
    }

    private:
    int _dim;
    std::vector<complex_t> _mesh;
    std::vector<matrix_t> _Ws; //W_is
    // See Eq. 6 PhysRevB.104.165111
    std::vector<matrix_cplx_mpt> _sqrt_one; //[1 - W_i * W_i^dagger]^0.5
    std::vector<matrix_cplx_mpt> _sqrt_two; //[1 - W_i^dagger * W_i]^-0.5

    /**
     * Calculate Hermitian square root of matrix M.
     *
     * @param M - matrix to calculate square root
     * @param is_Schur set true if matrix is Schur-matrix
     * @return square root of matrix
     */
    matrix_cplx_mpt sqrt_m (const matrix_cplx_mpt & M, bool & is_Schur) {
      Eigen::ComplexEigenSolver<matrix_cplx_mpt> ces;
      ces.compute(M);
      matrix_cplx_mpt D = ces.eigenvalues();
      is_Schur = true;
      for (int i = 0; i < D.rows(); i++) {
        if (D(i, 0).real() < tol) {
          is_Schur = false;
        }
      }
      return ces.eigenvectors() * D.array().sqrt().matrix().asDiagonal() * ces.eigenvectors().inverse().eval();
    }

  };
}
#endif //TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP
