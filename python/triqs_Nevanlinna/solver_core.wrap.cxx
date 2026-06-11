
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

template <> constexpr bool c2py::is_wrapped<triqs_Nevanlinna::kernels> = true;
template <>
const std::map<triqs_Nevanlinna::kernels, str_t> c2py::enum_to_string<triqs_Nevanlinna::kernels> = {
   {triqs_Nevanlinna::kernels::NEVANLINNA, "NEVANLINNA"},
   {triqs_Nevanlinna::kernels::CARATHEODORY, "CARATHEODORY"}};

// ==================== module classes =====================

// --------- class _c2py_cls_0 -----------
using _c2py_cls_0                                            = triqs_Nevanlinna::Nevanlinna_parameters_t;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_0>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_0> = "triqs_Nevanlinna.solver_core.NevanlinnaParametersT";

static int synth_constructor_0(PyObject *self, PyObject *args, PyObject *kwargs) {
  if (args and PyTuple_Check(args) and (PyTuple_Size(args) > 0)) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing triqs_Nevanlinna::Nevanlinna_parameters_t.\nNo positional arguments allowed. Use keywords arguments"));
    return -1;
  }
  c2py::pydict_extractor de{kwargs};
  try {
    ((c2py::wrap<_c2py_cls_0> *)self)->_c = new _c2py_cls_0{};
  } catch (std::exception const &e) {
    PyErr_SetString(PyExc_RuntimeError,
                    ("Error in constructing triqs_Nevanlinna::Nevanlinna_parameters_t from a Python dict.\n   "s + e.what()).c_str());
    return -1;
  }
  auto &self_c = *(((c2py::wrap<_c2py_cls_0> *)self)->_c);
  de("kernel", self_c.kernel, true);
  de("precision", self_c.precision, true);
  return de.check();
}

template <> constexpr initproc c2py::tp_init<_c2py_cls_0> = synth_constructor_0;

template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_0> =
   c2py::replace_tags(R"DOC(Synthesized constructor with the following keyword arguments:

Parameters
----------
kernel : {par_0}, default=NEVANLINNA

precision : {par_1}, default=100

)DOC",
                      "par", {c2py::python_typename<triqs_Nevanlinna::kernels>(), c2py::python_typename<int>()});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_0>[] = {

   {nullptr, nullptr, 0, nullptr} // Sentinel
};

constexpr auto _c2py_doc_member_0 = R"DOC(Continuation kernel to use.)DOC";
constexpr auto _c2py_doc_member_1 = R"DOC(Number of decimal digits of internal multiprecision arithmetic (only honored with MPFR support).)DOC";
static PyObject *prop_get_dict_0(PyObject *self, void *) {
  auto &self_c = *(((c2py::wrap<_c2py_cls_0> *)self)->_c);
  c2py::pydict dic;
  dic["kernel"]    = self_c.kernel;
  dic["precision"] = self_c.precision;
  return dic.new_ref();
}

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_0>[] = {
   c2py::getsetdef_from_member<&_c2py_cls_0::kernel, _c2py_cls_0>("kernel", _c2py_doc_member_0),
   c2py::getsetdef_from_member<&_c2py_cls_0::precision, _c2py_cls_0>("precision", _c2py_doc_member_1),
   {"__dict__", (getter)prop_get_dict_0, nullptr, "", nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <>
const std::string c2py::tp_doc<_c2py_cls_0> =
   R"DOC(Construction parameters for the Nevanlinna solver_core.)DOC" + std::string{"\n\n----------\n\n"} + c2py::tp_ctor_doc<_c2py_cls_0>;
// --------- class _c2py_cls_1 -----------
using _c2py_cls_1                                            = triqs_Nevanlinna::solver_core;
template <> constexpr bool c2py::is_wrapped<_c2py_cls_1>     = true;
template <> inline constexpr auto c2py::tp_name<_c2py_cls_1> = "triqs_Nevanlinna.solver_core.SolverCore";
static const auto _c2py_init_0 = c2py::dispatcher_c_kw_t{c2py::c_constructor<_c2py_cls_1, const triqs_Nevanlinna::Nevanlinna_parameters_t &>("p")};
template <> constexpr initproc c2py::tp_init<_c2py_cls_1> = c2py::pyfkw_constructor<_c2py_init_0>;
template <>
const std::string c2py::tp_ctor_doc<_c2py_cls_1> = _c2py_init_0.doc(R"DOC(
Construct the solver.

Parameters
----------
p : {par_0}
   Construction parameters (kernel choice and multiprecision precision).
)DOC",
                                                                    {{c2py::python_typename<const triqs_Nevanlinna::Nevanlinna_parameters_t &>()}});
// evaluate
static auto const _c2py_fun_0 = c2py::dispatcher_f_kw_t{
   c2py::cmethod([](_c2py_cls_1 &self, const triqs::mesh::refreq &grid, double eta) -> decltype(auto) { return self.evaluate(grid, eta); }, "self",
                 "grid", "eta"),
   c2py::cmethod([](_c2py_cls_1 &self, const triqs::mesh::refreq &grid, double eta,
                    nda::basic_array_view<const std::complex<double>, 3, nda::C_stride_layout, 'A', nda::default_accessor,
                                          nda::borrowed<nda::mem::AddressSpace::Host>>
                       theta) -> decltype(auto) { return self.evaluate(grid, eta, theta); },
                 "self", "grid", "eta", "theta")};

// solve
static auto const _c2py_fun_1 = c2py::dispatcher_f_kw_t{c2py::cmethod(
   [](_c2py_cls_1 &self, triqs::gfs::gf_const_view<triqs::mesh::imfreq> g_iw) -> decltype(auto) { return self.solve(g_iw); }, "self", "g_iw")};

static const auto _c2py_doc_0 =
   _c2py_fun_0.doc(R"DOC(
[1] Evaluate diagonal part of the real-frequency Green's function on a chosen grid.

Uses the precomputed Nevanlinna factorization.

------

[2] Evaluate the real-frequency Green's function on a chosen grid using Hardy-function optimization.

------

Parameters
----------
grid : {par_0}
   Real frequency grid.
eta : {par_1}
   Lorentzian broadening.
theta : {par_2}
   Hardy-function basis coefficients used to optimize the spectral function.

Returns
-------
{ret_0}
   Real-frequency matrix-valued TRIQS Green's function on a chosen grid.
)DOC",
                   {{c2py::python_typename<const triqs::mesh::refreq &>()},
                    {c2py::python_typename<double>()},
                    {c2py::python_typename<nda::basic_array_view<const std::complex<double>, 3, nda::C_stride_layout, 'A', nda::default_accessor,
                                                                 nda::borrowed<nda::mem::AddressSpace::Host>>>()}},
                   {c2py::python_typename<triqs::gfs::gf<triqs::mesh::refreq>>()});
static const auto _c2py_doc_1 = _c2py_fun_1.doc(R"DOC(
Perform a Nevanlinna factorization for a matrix-valued Matsubara frequency Green's function.

Parameters
----------
g_iw : {par_0}
   Matrix-valued Matsubara frequency Green's function.
)DOC",
                                                {{c2py::python_typename<triqs::gfs::gf_const_view<triqs::mesh::imfreq>>()}});

// ----- Method table ----
template <>
PyMethodDef c2py::tp_methods<_c2py_cls_1>[] = {
   {"evaluate", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {"solve", (PyCFunction)c2py::pyfkw<_c2py_fun_1>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_1.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

static constexpr auto prop_doc_0 = R"DOC(Eigenvalues of the Pick matrix (non-negative eigenvalues indicate the data is continuable).)DOC";

// ----- Member and property table ----

template <>
constinit PyGetSetDef c2py::tp_getset<_c2py_cls_1>[] = {

   {"Pick_eigenvalues", c2py::getter_from_method<c2py::castmc<>(&triqs_Nevanlinna::solver_core::get_Pick_eigenvalues)>, nullptr, prop_doc_0, nullptr},
   {nullptr, nullptr, nullptr, nullptr, nullptr}};

template <> PyMappingMethods c2py::tp_as_mapping<_c2py_cls_1> = {c2py::tpxx_size<_c2py_cls_1>, nullptr, nullptr};

template <>
const std::string c2py::tp_doc<_c2py_cls_1> = R"DOC(Nevanlinna analytical continuation solver for TRIQS Green's functions.

Performs analytical continuation for the diagonal part of the matrix-values 
TRIQS Green's function.)DOC"
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
   "solver_core",                                                                                                /* name of module */
   R"RAWDOC(Nevanlinna analytic continuation of fermionic Green's functions to the real-frequency axis.)RAWDOC", /* module documentation, may be NULL */
   -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
   module_methods,
   NULL,
   NULL,
   NULL,
   NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_solver_core() {

  if (not c2py::check_python_version("solver_core")) return NULL;

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
  _add_type(_c2py_cls_0, "NevanlinnaParametersT");
  _add_type(_c2py_cls_1, "SolverCore");
#undef _add_type

  return m;
}
#endif
// CLAIR_WRAP_GEN
