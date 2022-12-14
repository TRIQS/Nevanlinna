#ifndef TRIQS_NEVANLINNA_NEVANLINNA_FACTORIZATION_HPP
#define TRIQS_NEVANLINNA_NEVANLINNA_FACTORIZATION_HPP

#include <nda/nda.hpp>

#include <complex>

#include <Eigen/Dense>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "Nevanlinna_error.hpp"

using namespace std::complex_literals;

namespace triqs_Nevanlinna {
  /**
   * Class to construct Nevanlinna factorization for single orbital
   */
  class Nevanlinna_factorization {
    using complex_t = std::complex<boost::multiprecision::cpp_dec_float_100>;
    using matrix_t  = Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic>;
    public:

    /**
     * Build Nevanlinna factorization for data defined on positive Matsubara mesh
     *
     * @param mesh - Matsubara frequency mesh
     * @param data - Matsubara frequency data
     */
    void build(nda::vector_const_view<std::complex<double>> mesh, nda::vector_const_view<std::complex<double>> data);

    /**
     * Evaluate real-frequency data on specified grid with selected Lorentzian broadening
     *
     * @param grid - real frequency grid
     * @param eta  - Lorentzian broadening
     * @return Nevanlinna continuation on the specified grid
     */
    [[nodiscard]] nda::vector<double> evaluate(nda::vector_const_view<double> grid, double eta = 0.05) const;

    /**
     * Evaluate Nevanlinna continuation on a complex-valued grid
     * @param grid - grid for continuation evaluation
     * @return Nevanlinna continuation on the specified grid
     */
    [[nodiscard]] nda::vector<std::complex<double>> evaluate(nda::vector_const_view<std::complex<double>> grid) const;

    private:
    std::vector<complex_t> _phis;
    std::vector<matrix_t> _abcds;
    std::vector<complex_t> _mesh;

    std::vector<complex_t> mobius_trasformation(nda::vector_const_view<std::complex<double>> data) const;

  };

}
#endif //TRIQS_NEVANLINNA_NEVANLINNA_FACTORIZATION_HPP
