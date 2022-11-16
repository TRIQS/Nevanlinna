#pragma once

namespace triqs_Nevanlinna {

  using kernels = int;
  static const kernels NEVANLINNA=1;
  static const kernels CARATHEODORY=2;

  struct Nevanlinna_parameters_t {
    kernels kernel = NEVANLINNA;
  };

}
