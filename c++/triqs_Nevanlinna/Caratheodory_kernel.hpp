#ifndef TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP
#define TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP

#include "kernel.hpp"
#include <complex>

#include <Eigen/Dense>

#ifdef WITH_MPFR
#include <boost/multiprecision/mpfr.hpp>
#else
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

using namespace std::complex_literals;
namespace triqs_Nevanlinna {

  class Caratheodory_kernel : public kernel {
    static constexpr double tol = 1e-12;
#ifdef WITH_MPFR
    using real_t = boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<100>, boost::multiprecision::et_off>;
#else
    using real_t = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100>, boost::multiprecision::et_off>;
#endif
    using complex_t = std::complex<real_t>;
    using matrix_t  = Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic>;
    public:
    Caratheodory_kernel() = default;
    ~Caratheodory_kernel() override = default;

    // Copy/Move construction
    Caratheodory_kernel(Caratheodory_kernel const &) = default;
    Caratheodory_kernel(Caratheodory_kernel &&)      = default;

    /// Copy/Move assignment
    Caratheodory_kernel &operator=(Caratheodory_kernel const &) = default;
    Caratheodory_kernel &operator=(Caratheodory_kernel &&)      = default;

    void init(nda::vector_const_view<std::complex<double>> mesh, nda::array_const_view<std::complex<double>, 3> data) override;
    [[nodiscard]] nda::array<std::complex<double>, 3> evaluate(nda::vector_const_view<std::complex<double>> grid) const override;

    [[nodiscard]] virtual size_t size() const override {
      return _dim;
    }

    private:
    int _dim;
    std::vector<complex_t> _mesh;
    std::vector<matrix_t> _Ws; //W_is
    // See Eq. 6 PhysRevB.104.165111
    std::vector<matrix_t> _sqrt_one; //[1 - W_i * W_i^dagger]^0.5
    std::vector<matrix_t> _sqrt_two; //[1 - W_i^dagger * W_i]^-0.5

    /**
     * Calculate Hermitian square root of matrix M.
     *
     * @param M - matrix to calculate square root
     * @param is_Schur set true if matrix is Schur-matrix
     * @return square root of matrix
     */
    matrix_t sqrt_m (const matrix_t & M, bool & is_Schur) {
      Eigen::ComplexEigenSolver<matrix_t> ces;
      ces.compute(M);
      matrix_t D = ces.eigenvalues();
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
