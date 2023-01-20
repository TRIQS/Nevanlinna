#ifndef TRIQS_NEVANLINNA_NEVANLINNA_FACTORIZATION_HPP
#define TRIQS_NEVANLINNA_NEVANLINNA_FACTORIZATION_HPP

#include <nda/nda.hpp>

#include <complex>

#include "types.hpp"
#include "Nevanlinna_error.hpp"

namespace triqs_Nevanlinna {
  /**
   * Class to construct Nevanlinna factorization for single orbital
   */
  class Nevanlinna_factorization {
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
    std::vector<complex_mpt> _phis;
    std::vector<matrix_cplx_mpt> _abcds;
    std::vector<complex_mpt> _mesh;

    std::vector<complex_mpt> mobius_trasformation(nda::vector_const_view<std::complex<double>> data) const;

  };

}
#endif //TRIQS_NEVANLINNA_NEVANLINNA_FACTORIZATION_HPP
