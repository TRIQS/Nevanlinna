#ifndef TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP
#define TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP

#include "types.hpp"
#include "kernel.hpp"

namespace triqs_Nevanlinna {

  class Caratheodory_kernel : public kernel {
    static constexpr double tol = 1e-12;

    public:
    Caratheodory_kernel(int precision = mp_digits) : kernel(precision) {}

    void init(nda::vector_const_view<std::complex<double>> mesh, nda::array_const_view<std::complex<double>, 3> data) override;
    [[nodiscard]] nda::array<std::complex<double>, 3> evaluate(nda::vector_const_view<std::complex<double>> grid) override;

    [[nodiscard]] nda::array<std::complex<double>, 3> evaluate(nda::vector_const_view<std::complex<double>> grid,
                                                               nda::array_const_view<std::complex<double>, 3> theta) override {
      if (theta.shape()[0] != 0) {
        std::cerr << "Continuation poles optimization has not been implemented in matrix-valued continuation yet." << std::endl;
      }
      return evaluate(grid);
    };

    [[nodiscard]] size_t size() const override { return _dim; }

    [[nodiscard]] nda::vector<double> get_Pick_eigenvalues() const override;

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
     * Calculate Hermitian square root of matrix M.
     *
     * @param M - matrix to calculate square root
     * @param is_Schur set true if matrix is Schur-matrix
     * @return square root of matrix
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
