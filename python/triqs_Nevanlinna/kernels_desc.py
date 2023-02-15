# Generated automatically using the command :
# c++2py ../../c++/triqs_Nevanlinna/kernels.hpp -p --members_read_only -N triqs_Nevanlinna -a triqs_Nevanlinna -m kernels -o kernels --moduledoc="The Nevanlinna kernels python module" -C triqs --cxxflags="-std=c++20" --only="Nevanlinna_kernel Caratheodory_kernel  kernel kernels" --converter nda_py
from cpp2py.wrap_generator import *

# The module
module = module_(full_name = "kernels", doc = r"The Nevanlinna kernels python module", app_name = "triqs_Nevanlinna")

# Imports

# Add here all includes
module.add_include("triqs_Nevanlinna/kernels.hpp")

# Add here anything to add in the C++ code at the start, e.g. namespace using
module.add_preamble("""
#include <cpp2py/converters/complex.hpp>
#include <nda_py/cpp2py_converters.hpp>

using namespace triqs_Nevanlinna;
""")


# The class kernel
c = class_(
        py_type = "Kernel",  # name of the python class
        c_type = "triqs_Nevanlinna::kernel",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(int precison = mp_digits)""", doc = r"""""")

c.add_method("""void init (nda::vector_view<std::complex<double>> mesh, nda::array_view<std::complex<double>, 3> data)""",
             doc = r"""""")

c.add_method("""nda::array<std::complex<double>, 3> evaluate (nda::vector_view<std::complex<double>> grid)""",
             doc = r"""""")

c.add_method("""nda::array<std::complex<double>, 3> evaluate (nda::vector_view<std::complex<double>> grid, nda::array_view<std::complex<double>, 3> theta)""",
             doc = r"""""")

c.add_property(name = "size",
               getter = cfunction("size_t size ()"),
               doc = r"""""")

c.add_property(name = "pick_eigenvalues",
               getter = cfunction("nda::vector<double> pick_eigenvalues ()"),
               doc = r"""""")

module.add_class(c)

# The class Nevanlinna_kernel
c = class_(
        py_type = "NevanlinnaKernel",  # name of the python class
        c_type = "triqs_Nevanlinna::Nevanlinna_kernel",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(int precision = mp_digits)""", doc = r"""""")

c.add_method("""void init (nda::vector_view<std::complex<double>> mesh, nda::array_view<std::complex<double>, 3> data)""",
             doc = r"""""")

c.add_method("""nda::array<std::complex<double>, 3> evaluate (nda::vector_view<std::complex<double>> grid)""",
             doc = r"""""")

c.add_method("""nda::array<std::complex<double>, 3> evaluate (nda::vector_view<std::complex<double>> grid, nda::array_view<std::complex<double>, 3> theta)""",
             doc = r"""""")

c.add_property(name = "size",
               getter = cfunction("size_t size ()"),
               doc = r"""""")

c.add_property(name = "pick_eigenvalues",
               getter = cfunction("nda::vector<double> pick_eigenvalues ()"),
               doc = r"""""")

module.add_class(c)

# The class Caratheodory_kernel
c = class_(
        py_type = "CaratheodoryKernel",  # name of the python class
        c_type = "triqs_Nevanlinna::Caratheodory_kernel",   # name of the C++ class
        doc = r"""""",   # doc of the C++ class
        hdf5 = False,
)

c.add_constructor("""(int precision = mp_digits)""", doc = r"""""")

c.add_method("""void init (nda::vector_view<std::complex<double>> mesh, nda::array_view<std::complex<double>, 3> data)""",
             doc = r"""""")

c.add_method("""nda::array<std::complex<double>, 3> evaluate (nda::vector_view<std::complex<double>> grid)""",
             doc = r"""""")

c.add_method("""nda::array<std::complex<double>, 3> evaluate (nda::vector_view<std::complex<double>> grid, nda::array_view<std::complex<double>, 3> theta)""",
             doc = r"""""")

c.add_property(name = "size",
               getter = cfunction("size_t size ()"),
               doc = r"""""")

c.add_property(name = "pick_eigenvalues",
               getter = cfunction("nda::vector<double> pick_eigenvalues ()"),
               doc = r"""""")

module.add_class(c)



module.generate_code()