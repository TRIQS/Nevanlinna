#include "Caratheodory_kernel.hpp"


namespace triqs_Nevanlinna {
  void Caratheodory_kernel::init(nda::vector_const_view<std::complex<double>> mesh, nda::array_const_view<std::complex<double>, 3> data) {
    throw std::runtime_error("Kernel is not implemented yet.");
  }
  nda::array<std::complex<double>, 3> Caratheodory_kernel::evaluate(nda::vector_const_view<std::complex<double>> grid) const {
    throw std::runtime_error("Kernel is not implemented yet.");
  }

}
