#ifndef TRIQS_NEVANLINNA_KERNEL_HPP
#define TRIQS_NEVANLINNA_KERNEL_HPP
#include <complex>

#include <nda/nda.hpp>
#include <Eigen/Dense>
#include <boost/multiprecision/cpp_dec_float.hpp>


#include "Nevanlinna_factorization.hpp"

namespace triqs_Nevanlinna {
  class kernel {
    public:
    kernel() = default;
    virtual ~kernel() = default;

    // Copy/Move construction
    kernel(kernel const &) = default;
    kernel(kernel &&)      = default;

    /// Copy/Move assignment
    kernel &operator=(kernel const &) = default;
    kernel &operator=(kernel &&)      = default;

    virtual void init(const nda::array<std::complex<double>, 1> &mesh, const nda::array<std::complex<double>, 3> &data) = 0;
    virtual nda::array<std::complex<double>, 3> evaluate(const nda::array<std::complex<double>, 1> &grid) const = 0;

    virtual size_t size() const = 0;
  };

} // namespace nevanlinna
#endif //TRIQS_NEVANLINNA_KERNEL_HPP
