#pragma once

namespace triqs_Nevanlinna {

  /// Continuation algorithm selector.
  enum kernels { NEVANLINNA, CARATHEODORY };

  /// Construction parameters for the Nevanlinna solver_core.
  struct Nevanlinna_parameters_t {
    /// Continuation kernel to use.
    kernels kernel = NEVANLINNA;

    /// Number of decimal digits of internal multiprecision arithmetic (only honored with MPFR support).
    int precision = 100;
  };

} // namespace triqs_Nevanlinna
