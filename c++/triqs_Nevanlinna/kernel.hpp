#ifndef TRIQS_NEVANLINNA_KERNEL_HPP
#define TRIQS_NEVANLINNA_KERNEL_HPP
#include <complex>

#include <nda/nda.hpp>
#include "Nevanlinna_error.hpp"
#include "types.hpp"

namespace triqs_Nevanlinna {
  class kernel {
    public:
    virtual ~kernel() = default;
    kernel([[maybe_unused]] int precision = mp_digits) {
#ifdef WITH_MPFR
      boost::multiprecision::mpfr_float::default_precision(precision);
#else
      std::cerr << "Precision cannot be selected for boost multiprecision implementation. To set precision rebuild code with MPFR library support."
                << std::endl;
#endif
    };

    // Copy/Move construction
    kernel(kernel const &) = default;
    kernel(kernel &&)      = default;

    /// Copy/Move assignment
    kernel &operator=(kernel const &) = default;
    kernel &operator=(kernel &&)      = default;

    virtual void init(nda::vector_const_view<std::complex<double>> /*mesh*/, nda::array_const_view<std::complex<double>, 3> /*data*/) {
      throw Nevanlinna_error("This method is not implemented in the base class.");
    };
    [[nodiscard]] virtual nda::array<std::complex<double>, 3> evaluate(nda::vector_const_view<std::complex<double>> /*grid*/) {
      throw Nevanlinna_error("This method is not implemented in the base class.");
    };
    [[nodiscard]] virtual nda::array<std::complex<double>, 3> evaluate(nda::vector_const_view<std::complex<double>> /*grid*/,
                                                                       nda::array_const_view<std::complex<double>, 3> /*theta*/) {
      throw Nevanlinna_error("This method is not implemented in the base class.");
    };

    [[nodiscard]] virtual size_t size() const { return 0; };
    [[nodiscard]] virtual nda::vector<double> pick_eigenvalues() const {
      throw Nevanlinna_error("This method is not implemented in the base class.");
    };
  };

} // namespace triqs_Nevanlinna
#endif //TRIQS_NEVANLINNA_KERNEL_HPP
