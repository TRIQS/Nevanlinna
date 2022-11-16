#include "Caratheodory_kernel.hpp"


namespace triqs_Nevanlinna {
  void Caratheodory_kernel::init(const nda::array<std::complex<double>, 1> &mesh, const nda::array<std::complex<double>, 3> &data) {
    throw std::runtime_error("Kernel is not implemented yet.");
  }
  nda::array<std::complex<double>, 3> Caratheodory_kernel::evaluate(const nda::array<std::complex<double>, 1> &grid) const {
    throw std::runtime_error("Kernel is not implemented yet.");
  }

}