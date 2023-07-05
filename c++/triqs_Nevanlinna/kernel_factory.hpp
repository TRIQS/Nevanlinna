#ifndef TRIQS_NEVANLINNA_KERNEL_FACTORY_HPP
#define TRIQS_NEVANLINNA_KERNEL_FACTORY_HPP

#include "Nevanlinna_parameters_t.hpp"
#include "kernel.hpp"
#include "Nevanlinna_kernel.hpp"
#include "Caratheodory_kernel.hpp"

namespace triqs_Nevanlinna {
  class kernel_factory {
    public:
    static std::unique_ptr<kernel> get_kernel(Nevanlinna_parameters_t const &p) {
      switch (p.kernel) {
        case (NEVANLINNA): return std::make_unique<Nevanlinna_kernel>(p.precision); break;
        case (CARATHEODORY): return std::make_unique<Caratheodory_kernel>(p.precision); break;
        default: throw unimplemented_kernel_error("Kernel is not implemented");
      }
    }
  };

} // namespace triqs_Nevanlinna
#endif //TRIQS_NEVANLINNA_KERNEL_FACTORY_HPP
