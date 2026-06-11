
// C.f. https://numpy.org/doc/1.21/reference/c-api/array.html#importing-the-api
#define PY_ARRAY_UNIQUE_SYMBOL _cpp2py_ARRAY_API
#ifndef CLAIR_C2PY_WRAP_GEN
#ifdef __clang__
// #pragma clang diagnostic ignored "-W#warnings"
#endif
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#pragma GCC diagnostic ignored "-Wcpp"
#endif

#define C2PY_VERSION_MAJOR 0
#define C2PY_VERSION_MINOR 1

#include <c2py/c2py.hpp>

using c2py::operator""_a;

// ==================== enums =====================

// ==================== module classes =====================

// --------- class _c2py_cls_0 -----------
using _c2py_cls_0                                            = triqs_Nevanlinna::Nevanlinna_kernel;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_0>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_0> = "triqs_Nevanlinna.kernels_core.NevanlinnaKernel";
static const auto _c2py_init_0 = c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_0, int>("precision"_a = triqs_Nevanlinna::mp_digits)};
template <> constexpr initproc c2py::tp_init<_c2py_cls_0> = c2py::pyfkw_constructor<_c2py_init_0>;
template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_0> = _c2py_init_0.doc(R"DOC(
Construct the diagonal Nevanlinna kernel.

Parameters
----------
precision : {par_0}
   Number of decimal digits of internal multiprecision arithmetic (only honored with MPFR support).
)DOC",
                                                                    {{c2py::python_typename<int>()}});
// evaluate
static auto const _c2py_fun_0 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_0 &self,
                                            nda::basic_array_view<const std::complex<double>, 1, nda::C_stride_layout, 'V', nda::default_accessor,
                                                                  nda::borrowed<nda::mem::AddressSpace::Host>>
                                               grid) -> decltype(auto) { return self.evaluate(grid); },
                                         "self", "grid"),
                           c2py::cmethod([](_c2py_cls_0 &self,
                                            nda::basic_array_view<const std::complex<double>, 1, nda::C_stride_layout, 'V', nda::default_accessor,
                                                                  nda::borrowed<nda::mem::AddressSpace::Host>>
                                               grid,
                                            nda::basic_array_view<const std::complex<double>, 3, nda::C_stride_layout, 'A', nda::default_accessor,
                                                                  nda::borrowed<nda::mem::AddressSpace::Host>>
                                               theta) -> decltype(auto) { return self.evaluate(grid, theta); },
                                         "self", "grid", "theta")};

// get_Pick_eigenvalues
static auto const _c2py_fun_1 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_0 const &self) -> decltype(auto) { return self.get_Pick_eigenvalues(); }, "self")};

// init
static auto const _c2py_fun_2 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_0 &self,
                                            nda::basic_array_view<const std::complex<double>, 1, nda::C_stride_layout, 'V', nda::default_accessor,
                                                                  nda::borrowed<nda::mem::AddressSpace::Host>>
                                               mesh,
                                            nda::basic_array_view<const std::complex<double>, 3, nda::C_stride_layout, 'A', nda::default_accessor,
                                                                  nda::borrowed<nda::mem::AddressSpace::Host>>
                                               data) -> decltype(auto) { return self.init(mesh, data); },
                                         "self", "mesh", "data")};

static const auto _c2py_doc_0 = _c2py_fun_0.doc(
   R"DOC(
[1] Evaluate the diagonal real-frequency Green's function on a chosen grid.

------

[2] Evaluate the diagonal real-frequency Green's function using Hardy-function optimization.

------

Parameters
----------
grid : {par_0}
   Complex real-frequency grid (real frequency plus Lorentzian broadening).
theta : {par_1}
   Hardy-function basis coefficients used to optimize the spectral function.

Returns
-------
{ret_0}
   Matrix-valued real-frequency Green's function on the grid.
)DOC",
   {{c2py::python_typename<nda::basic_array_view<const std::complex<double>, 1, nda::C_stride_layout, 'V', nda::default_accessor,
                                                 nda::borrowed<nda::mem::AddressSpace::Host>>>()},
    {c2py::python_typename<nda::basic_array_view<const std::complex<double>, 3, nda::C_stride_layout, 'A', nda::default_accessor,
                                                 nda::borrowed<nda::mem::AddressSpace::Host>>>()}},
   {c2py::python_typename<
      nda::basic_array<std::complex<double>, 3, nda::C_layout, 'A', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>()});
static const auto _c2py_doc_1 = _c2py_fun_1.doc(R"DOC()DOC");
static const auto _c2py_doc_2 =
   _c2py_fun_2.doc(R"DOC(
Initialize the diagonal Nevanlinna continuation from Matsubara-frequency input data.

Parameters
----------
mesh : {par_0}
   Positive Matsubara frequencies, given as complex values.
data : {par_1}
   Matrix-valued Green's function data on those frequencies.
)DOC",
                   {{c2py::python_typename<nda::basic_array_view<const std::complex<double>, 1, nda::C_stride_layout, 'V', nda::default_accessor,
                                                                 nda::borrowed<nda::mem::AddressSpace::Host>>>()},
                    {c2py::python_typename<nda::basic_array_view<const std::complex<double>, 3, nda::C_stride_layout, 'A', nda::default_accessor,
                                                                 nda::borrowed<nda::mem::AddressSpace::Host>>>()}});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_0>[] = {
   {"evaluate", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {"get_Pick_eigenvalues", (PyCFunction)c2py::pyfkw<_c2py_fun_1>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_1.c_str()},
   {"init", (PyCFunction)c2py::pyfkw<_c2py_fun_2>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_2.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

static constexpr auto prop_doc_0 = R"DOC(Eigenvalues of the Pick matrix; non-negative eigenvalues indicate the data is continuable (Nevanlinna).)DOC";

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_0>[] = {

   {"Pick_eigenvalues", c2py::getter_from_method<c2py::castmc<>(&triqs_Nevanlinna::Nevanlinna_kernel::get_Pick_eigenvalues)>, nullptr, prop_doc_0,
    nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <> PyMappingMethods c2py::tp_as_mapping<_c2py_cls_0> = {c2py::tpxx_size<_c2py_cls_0>, nullptr, nullptr};

template <>
const std::string c2py::tp_doc<_c2py_cls_0> = R"DOC(Diagonal Nevanlinna continuation kernel.

Builds one independent Nevanlinna factorization per orbital; off-diagonal elements of the
input Green's function are ignored. Supports Hardy-function (theta) optimization in ``evaluate()``.)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_0>;
// --------- class _c2py_cls_1 -----------
using _c2py_cls_1                                            = triqs_Nevanlinna::Caratheodory_kernel;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_1>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_1> = "triqs_Nevanlinna.kernels_core.CaratheodoryKernel";
static const auto _c2py_init_1 = c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_1, int>("precision"_a = triqs_Nevanlinna::mp_digits)};
template <> constexpr initproc c2py::tp_init<_c2py_cls_1> = c2py::pyfkw_constructor<_c2py_init_1>;
template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_1> = _c2py_init_1.doc(R"DOC(
Construct the full matrix-valued Caratheodory kernel.

Parameters
----------
precision : {par_0}
   Number of decimal digits of internal multiprecision arithmetic (only honored with MPFR support).
)DOC",
                                                                    {{c2py::python_typename<int>()}});
// evaluate
static auto const _c2py_fun_3 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_1 &self,
                                            nda::basic_array_view<const std::complex<double>, 1, nda::C_stride_layout, 'V', nda::default_accessor,
                                                                  nda::borrowed<nda::mem::AddressSpace::Host>>
                                               grid) -> decltype(auto) { return self.evaluate(grid); },
                                         "self", "grid"),
                           c2py::cmethod([](_c2py_cls_1 &self,
                                            nda::basic_array_view<const std::complex<double>, 1, nda::C_stride_layout, 'V', nda::default_accessor,
                                                                  nda::borrowed<nda::mem::AddressSpace::Host>>
                                               grid,
                                            nda::basic_array_view<const std::complex<double>, 3, nda::C_stride_layout, 'A', nda::default_accessor,
                                                                  nda::borrowed<nda::mem::AddressSpace::Host>>
                                               theta) -> decltype(auto) { return self.evaluate(grid, theta); },
                                         "self", "grid", "theta")};

// get_Pick_eigenvalues
static auto const _c2py_fun_4 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_1 const &self) -> decltype(auto) { return self.get_Pick_eigenvalues(); }, "self")};

// init
static auto const _c2py_fun_5 =
   c2py::dispatcher_f_kw_t{c2py::cmethod([](_c2py_cls_1 &self,
                                            nda::basic_array_view<const std::complex<double>, 1, nda::C_stride_layout, 'V', nda::default_accessor,
                                                                  nda::borrowed<nda::mem::AddressSpace::Host>>
                                               mesh,
                                            nda::basic_array_view<const std::complex<double>, 3, nda::C_stride_layout, 'A', nda::default_accessor,
                                                                  nda::borrowed<nda::mem::AddressSpace::Host>>
                                               data) -> decltype(auto) { return self.init(mesh, data); },
                                         "self", "mesh", "data")};

static const auto _c2py_doc_3 = _c2py_fun_3.doc(
   R"DOC(
[1] Evaluate the full matrix-valued real-frequency Green's function on a chosen grid.

------

[2] Evaluate the full matrix-valued real-frequency Green's function on a chosen grid.

Theta optimization is not supported here. A non-empty ``theta`` is ignored and plain ``evaluate(grid)`` is used.

------

Parameters
----------
grid : {par_0}
   Complex real-frequency grid (real frequency plus Lorentzian broadening).
theta : {par_1}
   Hardy-function basis coefficients used to optimize the spectral function (ignored).

Returns
-------
{ret_0}
   Matrix-valued real-frequency Green's function on the grid.
)DOC",
   {{c2py::python_typename<nda::basic_array_view<const std::complex<double>, 1, nda::C_stride_layout, 'V', nda::default_accessor,
                                                 nda::borrowed<nda::mem::AddressSpace::Host>>>()},
    {c2py::python_typename<nda::basic_array_view<const std::complex<double>, 3, nda::C_stride_layout, 'A', nda::default_accessor,
                                                 nda::borrowed<nda::mem::AddressSpace::Host>>>()}},
   {c2py::python_typename<
      nda::basic_array<std::complex<double>, 3, nda::C_layout, 'A', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>>>()});
static const auto _c2py_doc_4 = _c2py_fun_4.doc(R"DOC()DOC");
static const auto _c2py_doc_5 =
   _c2py_fun_5.doc(R"DOC(
Build the full matrix-valued Caratheodory continuation from Matsubara-frequency input data.

Parameters
----------
mesh : {par_0}
   Positive Matsubara frequencies, given as complex values.
data : {par_1}
   Matrix-valued Green's function data on those frequencies.
)DOC",
                   {{c2py::python_typename<nda::basic_array_view<const std::complex<double>, 1, nda::C_stride_layout, 'V', nda::default_accessor,
                                                                 nda::borrowed<nda::mem::AddressSpace::Host>>>()},
                    {c2py::python_typename<nda::basic_array_view<const std::complex<double>, 3, nda::C_stride_layout, 'A', nda::default_accessor,
                                                                 nda::borrowed<nda::mem::AddressSpace::Host>>>()}});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_1>[] = {
   {"evaluate", (PyCFunction)c2py::pyfkw<_c2py_fun_3>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_3.c_str()},
   {"get_Pick_eigenvalues", (PyCFunction)c2py::pyfkw<_c2py_fun_4>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_4.c_str()},
   {"init", (PyCFunction)c2py::pyfkw<_c2py_fun_5>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_5.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

static constexpr auto prop_doc_1 = R"DOC(Eigenvalues of the Pick matrix; non-negative eigenvalues indicate the data is continuable (Nevanlinna).)DOC";

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_1>[] = {

   {"Pick_eigenvalues", c2py::getter_from_method<c2py::castmc<>(&triqs_Nevanlinna::Caratheodory_kernel::get_Pick_eigenvalues)>, nullptr, prop_doc_1,
    nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <> PyMappingMethods c2py::tp_as_mapping<_c2py_cls_1> = {c2py::tpxx_size<_c2py_cls_1>, nullptr, nullptr};

template <>
const std::string c2py::tp_doc<_c2py_cls_1> = R"DOC(Full matrix-valued Caratheodory continuation kernel (PhysRevB.104.165111).

Continues the complete matrix-valued Green's function, including off-diagonal elements.
The Hardy-function (theta) optimization path is not implemented for this kernel.)DOC"
   + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_1>;

// ==================== module functions ====================

//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {
   PyModuleDef_HEAD_INIT,
   "kernels_core",                                                                              /* name of module */
   R"RAWDOC(Matrix-valued analytic-continuation kernels (Nevanlinna and Caratheodory).)RAWDOC", /* module documentation, may be NULL */
   -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
   module_methods,
   NULL,
   NULL,
   NULL,
   NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_kernels_core() {

  if (not c2py::check_python_version("kernels_core")) return NULL;

  // import numpy iff 'numpy/arrayobject.h' included
#ifdef Py_ARRAYOBJECT_H
  import_array();
#endif

  PyObject *m;

  if (PyType_Ready(&c2py::wrap_pytype<c2py::py_range>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_0>) < 0) return NULL;
  if (PyType_Ready(&c2py::wrap_pytype<_c2py_cls_1>) < 0) return NULL;

  m = PyModule_Create(&module_def);
  if (m == NULL) return NULL;

  auto &conv_table = *c2py::conv_table_sptr.get();

  conv_table[std::type_index(typeid(c2py::py_range)).name()] = &c2py::wrap_pytype<c2py::py_range>;
#define _add_type(T, N) c2py::add_type_object_to_main<T>(N, m, conv_table)
  _add_type(_c2py_cls_0, "NevanlinnaKernel");
  _add_type(_c2py_cls_1, "CaratheodoryKernel");
#undef _add_type

  return m;
}
#endif
// CLAIR_WRAP_GEN
