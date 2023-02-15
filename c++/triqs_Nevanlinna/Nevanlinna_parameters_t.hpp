#pragma once

namespace triqs_Nevanlinna {

  enum kernels { NEVANLINNA, CARATHEODORY };

  struct Nevanlinna_parameters_t {
    kernels kernel = NEVANLINNA;
    int precision  = 100;
  };

} // namespace triqs_Nevanlinna
