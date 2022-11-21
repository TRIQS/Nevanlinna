#ifndef TRIQS_NEVANLINNA_NEVANLINNA_KERNEL_HPP
#define TRIQS_NEVANLINNA_NEVANLINNA_KERNEL_HPP

#include <complex>

#include <nda/nda.hpp>
#include <Eigen/Dense>


#include "kernel.hpp"
#include "Nevanlinna_factorization.hpp"

namespace triqs_Nevanlinna {

  class Nevanlinna_kernel : public kernel {

    public:
    Nevanlinna_kernel() : kernel() {
      std::cerr<<"This is Nevanlinna analytical continuation. All off-diagonal elements will be ignored."<<std::endl;
    }
    ~Nevanlinna_kernel() override = default;

    // Copy/Move construction
    Nevanlinna_kernel(Nevanlinna_kernel const &) = default;
    Nevanlinna_kernel(Nevanlinna_kernel &&)      = default;

    /// Copy/Move assignment
    Nevanlinna_kernel &operator=(Nevanlinna_kernel const &) = default;
    Nevanlinna_kernel &operator=(Nevanlinna_kernel &&)      = default;

    void init(const nda::array<std::complex<double>, 1> &mesh, const nda::array<std::complex<double>, 3> &data) override;
    [[nodiscard]] nda::array<std::complex<double>, 3> evaluate(const nda::array<std::complex<double>, 1> &grid) const override;

    [[nodiscard]] virtual size_t size() const override {
      return _factorizations.size();
    }

    private:
    std::vector<Nevanlinna_factorization> _factorizations;

};

}
#endif //TRIQS_NEVANLINNA_NEVANLINNA_KERNEL_HPP
