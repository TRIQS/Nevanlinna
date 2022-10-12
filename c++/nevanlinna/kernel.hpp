#ifndef NEVANLINNA_KERNEL_HPP
#define NEVANLINNA_KERNEL_HPP
#include <nda/nda.hpp>

#include <complex>

#include <gmpxx.h>

#include <Eigen/Dense>

#include "nevanlinna_error.hpp"

using namespace std::complex_literals;

namespace nevanlinna {
  class kernel {
    using complex_t = std::complex<mpf_class>;
    using matrix_t = Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic>;
    public:

    kernel() = default;

    // Copy/Move construction
    kernel(kernel const &) = default;
    kernel(kernel &&)      = default;

    /// Copy/Move assignment
    kernel &operator=(kernel const &) = default;
    kernel &operator=(kernel &&) = default;

    void solve(const nda::array<std::complex<double>, 1> & mesh, const nda::array<std::complex<double>, 1> & data);
    nda::array<double, 1> evaluate(const nda::array<double, 1> & grid, double eta = 0.05) const;
    nda::array<std::complex<double>, 1> evaluate(const nda::array<std::complex<double>, 1> & grid) const;

    private:
    std::vector<complex_t> _phis;
    std::vector<matrix_t> _abcds;
    std::vector<complex_t> _mesh;

    std::vector<complex_t> mobius_trasformation(const nda::array<std::complex<double>, 1> & data) const;
  };

}
#endif //NEVANLINNA_KERNEL_HPP
