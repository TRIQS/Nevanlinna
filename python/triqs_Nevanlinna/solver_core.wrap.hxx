#include <c2py/c2py.hpp>

#ifndef C2PY_HXX_DECLARATION_solver_core_GUARDS
#define C2PY_HXX_DECLARATION_solver_core_GUARDS
template <> constexpr bool c2py::is_wrapped<triqs_Nevanlinna::Nevanlinna_parameters_t>     = true;
template <> inline constexpr auto c2py::tp_name<triqs_Nevanlinna::Nevanlinna_parameters_t> = "triqs_Nevanlinna.solver_core.NevanlinnaParametersT";
template <> constexpr bool c2py::is_wrapped<triqs_Nevanlinna::solver_core>                 = true;
template <> inline constexpr auto c2py::tp_name<triqs_Nevanlinna::solver_core>             = "triqs_Nevanlinna.solver_core.SolverCore";
template <> constexpr bool c2py::is_wrapped<triqs_Nevanlinna::kernels>                     = true;
template <>
const std::map<triqs_Nevanlinna::kernels, str_t> c2py::enum_to_string<triqs_Nevanlinna::kernels> = {
   {triqs_Nevanlinna::kernels::NEVANLINNA, "NEVANLINNA"},
   {triqs_Nevanlinna::kernels::CARATHEODORY, "CARATHEODORY"}};
#endif