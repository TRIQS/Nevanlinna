# Generated automatically using the command :
# c++2py ../../c++/triqs_Nevanlinna/solver_core.hpp -p --members_read_only -N triqs_Nevanlinna -a triqs_Nevanlinna -m solver_core -o solver_core --moduledoc="The Nevanlinna python module" -C triqs --cxxflags="-std=c++20" --only="solver_core kernels" --converter nda_py
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "solver_core", doc = r"The Nevanlinna python module", app_name = "triqs_Nevanlinna")

# Imports
module.add_imports(*['triqs.gf', 'triqs.gf.meshes'])

# Add here all includes
module.add_include("triqs_Nevanlinna/solver_core.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <nda_py/cpp2py_converters.hpp>
#include <triqs/cpp2py_converters/gf.hpp>

using namespace triqs_Nevanlinna;
""")

module.add_enum("kernels", ['NEVANLINNA', 'CARATHEODORY'], "triqs_Nevanlinna", doc = r"""""")

# The class solver_core
c = class_(
        py_type = "SolverCore",  # name of the python class
        c_type = "triqs_Nevanlinna::solver_core",   # name of the C++ class
        doc = r"""Nevanlinna analytical continuation solver_core for TRIQS GFs""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(**triqs_Nevanlinna::Nevanlinna_parameters_t)""", doc = r"""



+----------------+---------------------------+------------+---------------+
| Parameter Name | Type                      | Default    | Documentation |
+================+===========================+============+===============+
| kernel         | triqs_Nevanlinna::kernels | NEVANLINNA |               |
+----------------+---------------------------+------------+---------------+
| precision      | int                       | 100        |               |
+----------------+---------------------------+------------+---------------+
""")

c.add_method("""void solve (triqs::gfs::gf_view<triqs::mesh::imfreq> g_iw)""",
             doc = r"""Construct a Nevanlinna factorization for matrix-valued Matsubara frequency Green's function

Parameters
----------
g_iw
     - matrix-valued Matsubara frequency Green's function""")

c.add_method("""triqs::gfs::gf<triqs::mesh::refreq> evaluate (triqs::mesh::refreq grid, double eta)""",
             doc = r"""Evaluate diagonal part of the real-frequency Green's function on a chosen grid
 based on the precomputed Nevanlinna factorization

Parameters
----------
grid
     - real frequency grid

eta
     - Lorentzian broadening

Returns
-------
out
     Real-frequency matrix-valued TRIQS Green's function on a chosen grid.""")

c.add_method("""triqs::gfs::gf<triqs::mesh::refreq> evaluate (triqs::mesh::refreq grid, double eta, nda::array_view<std::complex<double>, 3> theta)""",
             doc = r"""""")

c.add_property(name = "Pick_eigenvalues",
               getter = cfunction("nda::vector<double> get_Pick_eigenvalues ()"),
               doc = r"""""")

c.add_property(name = "size",
               getter = cfunction("size_t size ()"),
               doc = r"""""")

module.add_class(c)


# Converter for Nevanlinna_parameters_t
c = converter_(
        c_type = "triqs_Nevanlinna::Nevanlinna_parameters_t",
        doc = r"""""",
)
c.add_member(c_name = "kernel",
             c_type = "triqs_Nevanlinna::kernels",
             initializer = """ NEVANLINNA """,
             doc = r"""""")

c.add_member(c_name = "precision",
             c_type = "int",
             initializer = """ 100 """,
             doc = r"""""")

module.add_converter(c)


module.generate_code()
