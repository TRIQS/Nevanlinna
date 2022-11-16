#ifndef TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP
#define TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP

#include "kernel.hpp"

namespace triqs_Nevanlinna {
  class Caratheodory_kernel : public kernel {

    public:
    virtual ~Caratheodory_kernel() = default;
    void init(const nda::array<std::complex<double>, 1> &mesh, const nda::array<std::complex<double>, 3> &data) override;
    [[nodiscard]] nda::array<std::complex<double>, 3> evaluate(const nda::array<std::complex<double>, 1> &grid) const override;

    [[nodiscard]] virtual size_t size() const override {
      return 0;
    }

  };
}
#endif //TRIQS_NEVANLINNA_CARATHEODORY_KERNEL_HPP
