#ifndef TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP
#define TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP

#include "kernel.hpp"

namespace triqs_Nevanlinna {
  class Caratheodory_kernel : public kernel {

    public:
    void init(nda::vector_const_view<std::complex<double>> mesh, nda::array_const_view<std::complex<double>, 3> data) override;
    [[nodiscard]] nda::array<std::complex<double>, 3> evaluate(nda::vector_const_view<std::complex<double>> grid) const override;

    [[nodiscard]] size_t size() const override {
      return 0;
    }

  };
}
#endif //TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP
