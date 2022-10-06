# Generated automatically using the command :
# c++2py ../../c++/nevanlinna/nevanlinna.hpp -p --members_read_only -N nevanlinna -a nevanlinna -m nevanlinna_module -o nevanlinna_module --moduledoc="The nevanlinna python module" -C triqs --cxxflags="-std=c++20" --only="solver" --converter nda_py
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "nevanlinna_module", doc = r"The nevanlinna python module", app_name = "nevanlinna")

# Imports

# Add here all includes
module.add_include("nevanlinna/nevanlinna.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <nda_py/cpp2py_converters.hpp>

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



+----------------+--------+---------+---------------+
| Parameter Name | Type   | Default | Documentation |
+================+========+=========+===============+
| N_omega        | int    | --      |               |
+----------------+--------+---------+---------------+
| Omega_min      | double | --      |               |
+----------------+--------+---------+---------------+
| Omega_max      | double | --      |               |
+----------------+--------+---------+---------------+
""")

c.add_method("""void solve (nda::array<double, 2> input)""",
             doc = r"""""")

c.add_method("""nda::array<double, 1> evaluate (nda::array<double, 1> grid)""",
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