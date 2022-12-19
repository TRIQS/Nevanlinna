# Generated automatically using the command :
# c++2py ../../c++/triqs_Nevanlinna/solver.hpp -p --members_read_only -N triqs_Nevanlinna -a triqs_Nevanlinna -m solver -o solver --moduledoc="The Nevanlinna python module" -C triqs --cxxflags="-std=c++20" --only="solver" --converter nda_py
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "solver", doc = r"The Nevanlinna python module", app_name = "triqs_Nevanlinna")

# Imports
module.add_imports(*['triqs.gf', 'triqs.gf.meshes'])

# Add here all includes
module.add_include("triqs_Nevanlinna/solver.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <triqs/cpp2py_converters/gf.hpp>

using namespace triqs_Nevanlinna;
""")


# The class solver
c = class_(
        py_type = "Solver",  # name of the python class
        c_type = "triqs_Nevanlinna::solver",   # name of the C++ class
        doc = r"""Nevanlinna analytical continuation solver for TRIQS GFs""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(**triqs_Nevanlinna::Nevanlinna_parameters_t)""", doc = r"""



+----------------+---------------------------+------------+---------------+
| Parameter Name | Type                      | Default    | Documentation |
+================+===========================+============+===============+
| kernel         | triqs_Nevanlinna::kernels | NEVANLINNA |               |
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

module.add_converter(c)


module.generate_code()