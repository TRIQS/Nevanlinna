#ifndef TRIQS_NEVANLINNA_NEVANLINNA_KERNEL_HPP
#define TRIQS_NEVANLINNA_NEVANLINNA_KERNEL_HPP

#include <complex>

#include <nda/nda.hpp>

#include "kernel.hpp"
#include "Nevanlinna_factorization.hpp"

namespace triqs_Nevanlinna {

  class Nevanlinna_kernel : public kernel {

    public:
    Nevanlinna_kernel(int precision = mp_digits) : kernel(precision) {
      std::cerr << "This is Nevanlinna analytical continuation. All off-diagonal elements will be ignored." << std::endl;
    }

    void init(nda::vector_const_view<std::complex<double>> mesh, nda::array_const_view<std::complex<double>, 3> data) override;
    [[nodiscard]] nda::array<std::complex<double>, 3> evaluate(nda::vector_const_view<std::complex<double>> grid) override;
    [[nodiscard]] nda::array<std::complex<double>, 3> evaluate(nda::vector_const_view<std::complex<double>> grid,
                                                               nda::array_const_view<std::complex<double>, 3> theta) override;

    [[nodiscard]] size_t size() const override { return _factorizations.size(); }

    [[nodiscard]] nda::vector<double> pick_eigenvalues() const override;

    private:
    size_t _N_im_freq{};
    std::vector<Nevanlinna_factorization> _factorizations;
  };

} // namespace triqs_Nevanlinna
#endif //TRIQS_NEVANLINNA_NEVANLINNA_KERNEL_HPP
