#ifndef NEVANLINNA_KERNEL_HPP
#define NEVANLINNA_KERNEL_HPP
#include <nda/nda.hpp>

#include <complex>

#include <gmpxx.h>
#include <Eigen/Dense>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "Nevanlinna_error.hpp"


using namespace std::complex_literals;

namespace triqs_Nevanlinna {
  class kernel {
    using complex_t = std::complex<boost::multiprecision::cpp_dec_float_100>;
    using matrix_t  = Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic>;

    public:
    kernel() = default;

    // Copy/Move construction
    kernel(kernel const &) = default;
    kernel(kernel &&)      = default;

    /// Copy/Move assignment
    kernel &operator=(kernel const &) = default;
    kernel &operator=(kernel &&)      = default;

    void solve(const nda::array<std::complex<double>, 1> &mesh, const nda::array<std::complex<double>, 1> &data);
    nda::array<double, 1> evaluate(const nda::array<double, 1> &grid, double eta = 0.05) const;
    nda::array<std::complex<double>, 1> evaluate(const nda::array<std::complex<double>, 1> &grid) const;

    private:
    std::vector<complex_t> _phis;
    std::vector<matrix_t> _abcds;
    std::vector<complex_t> _mesh;

    std::vector<complex_t> mobius_trasformation(const nda::array<std::complex<double>, 1> &data) const;
  };

} // namespace nevanlinna
#endif //NEVANLINNA_KERNEL_HPP
