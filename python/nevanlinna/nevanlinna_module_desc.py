# Generated automatically using the command :
# c++2py ../../c++/nevanlinna/nevanlinna.hpp -p --members_read_only -N nevanlinna -a nevanlinna -m nevanlinna_module -o nevanlinna_module --moduledoc="The nevanlinna python module" -C triqs --cxxflags="-std=c++20" --only="solver" --converter nda_py
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "nevanlinna_module", doc = r"The nevanlinna python module", app_name = "nevanlinna")

# Imports
module.add_imports(*['triqs.gf', 'triqs.gf.meshes'])

# Add here all includes
module.add_include("nevanlinna/nevanlinna.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <nda_py/cpp2py_converters.hpp>
#include <triqs/cpp2py_converters/gf.hpp>

using namespace nevanlinna;
""")


# The class solver
c = class_(
        py_type = "Solver",  # name of the python class
        c_type = "nevanlinna::solver",   # name of the C++ class
        doc = r"""Nevanlinna analytical continuation solver""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(**nevanlinna_parameters_t)""", doc = r"""



+----------------+------+---------+---------------+
| Parameter Name | Type | Default | Documentation |
+================+======+=========+===============+
| precision      | int  | --      |               |
+----------------+------+---------+---------------+
""")

c.add_method("""void solve (nda::array<std::complex<double>, 1> mesh, nda::array<std::complex<double>, 1> data)""",
             doc = r"""""")

c.add_method("""void solve (triqs::gfs::gf<triqs::mesh::imfreq, triqs::gfs::scalar_valued> g_iw)""",
             doc = r"""""")

c.add_method("""nda::array<double, 1> evaluate (nda::array<std::complex<double>, 1> grid)""",
             doc = r"""""")

c.add_method("""nda::array<double, 1> evaluate (nda::array<double, 1> grid, double eta = 0.05)""",
             doc = r"""""")

module.add_class(c)


# Converter for nevanlinna_parameters_t
c = converter_(
        c_type = "nevanlinna::nevanlinna_parameters_t",
        doc = r"""""",
)
c.add_member(c_name = "precision",
             c_type = "int",
             initializer = """  """,
             doc = r"""""")

module.add_converter(c)


module.generate_code()