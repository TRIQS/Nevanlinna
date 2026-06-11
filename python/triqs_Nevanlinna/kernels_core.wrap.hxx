#include <c2py/c2py.hpp>

#ifndef C2PY_HXX_DECLARATION_kernels_core_GUARDS
#define C2PY_HXX_DECLARATION_kernels_core_GUARDS
template <> constexpr bool c2py::is_wrapped<triqs_Nevanlinna::Nevanlinna_kernel>       = true;
template <> inline constexpr auto c2py::tp_name<triqs_Nevanlinna::Nevanlinna_kernel>   = "triqs_Nevanlinna.kernels_core.NevanlinnaKernel";
template <> constexpr bool c2py::is_wrapped<triqs_Nevanlinna::Caratheodory_kernel>     = true;
template <> inline constexpr auto c2py::tp_name<triqs_Nevanlinna::Caratheodory_kernel> = "triqs_Nevanlinna.kernels_core.CaratheodoryKernel";
#endif