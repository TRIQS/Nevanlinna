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
    [[nodiscard]] nda::vector<double> evaluate(nda::vector_const_view<double> grid, double eta = 0.05);

    /**
     * Evaluate Nevanlinna continuation on a complex-valued grid
     * @param grid - grid for continuation evaluation
     * @return Nevanlinna continuation on the specified grid
     */
    [[nodiscard]] nda::vector<std::complex<double>>
    evaluate(nda::vector_const_view<std::complex<double>> grid,
             nda::vector_const_view<std::complex<double>> theta_M_1 = nda::vector_const_view<std::complex<double>>());

    nda::vector<double> get_Pick_eigenvalues() const {
      auto M = _data.shape()[0];
      if (M == 0) { throw Nevanlinna_uninitialized_error("Empty continuation data. Please run solve(...) first."); }
      //fill the Pick matrix
      auto Pick = Eigen::MatrixXcd(M, M);
      auto I    = complex_mpt{0., 1.};
      auto one  = complex_mpt{1., 0.};
      for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
          complex_mpt freq_i = (_mesh(i) - I) / (_mesh(i) + I);
          complex_mpt freq_j = (_mesh(j) - I) / (_mesh(j) + I);
          auto val           = (one - _data(i) * std::conj(_data(j))) / (one - freq_i * std::conj(freq_j));
          Pick(i, j)         = std::complex<double>(val.real().convert_to<double>(), val.imag().convert_to<double>());
        }
      }
      auto evals            = Pick.eigenvalues();
      auto Pick_eigenvalues = nda::vector<double>(M);
      std::transform(evals.begin(), evals.end(), Pick_eigenvalues.begin(), [](const std::complex<double> &r) { return r.real(); });
      return Pick_eigenvalues;
    }

    private:
    nda::vector<complex_mpt> _phis;
    nda::vector<matrix_cplx_mpt> _abcds;
    nda::vector<complex_mpt> _mesh;
    nda::vector<complex_mpt> _data;
    nda::vector<complex_mpt> _grid;
    nda::vector<matrix_cplx_mpt> _coeffs;

    [[nodiscard]] nda::vector<complex_mpt> mobius_trasformation(nda::vector_const_view<std::complex<double>> data) const;

    [[nodiscard]] nda::vector<std::complex<double>> evaluate_for_theta(nda::vector_const_view<std::complex<double>> grid,
                                                                       nda::vector_const_view<std::complex<double>> theta_M_1) const;
  };

} // namespace triqs_Nevanlinna
#endif //TRIQS_NEVANLINNA_NEVANLINNA_FACTORIZATION_HPP
