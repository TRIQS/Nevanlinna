#ifndef TRIQS_NEVANLINNA_KERNEL_FACTORY_HPP
#define TRIQS_NEVANLINNA_KERNEL_FACTORY_HPP

#include "Nevanlinna_parameters_t.hpp"
#include "kernel.hpp"
#include "Nevanlinna_kernel.hpp"
#include "Caratheodory_kernel.hpp"

namespace triqs_Nevanlinna{
class kernel_factory {
  public:
  static std::unique_ptr<kernel> get_kernel(kernels k) {
    switch(k) {
      case(NEVANLINNA):
        return std::make_unique<Nevanlinna_kernel>();
        break;
      case(CARATHEODORY):
        return std::make_unique<Caratheodory_kernel>();
        break;
      default:
        throw unimplemented_kernel_error("Kernel is not implemented");
    }
  }
};

}
#endif //TRIQS_NEVANLINNA_KERNEL_FACTORY_HPP
