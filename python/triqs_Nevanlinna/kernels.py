"""
Matrix-valued analytic-continuation kernels.

Re-exports the compiled :class:`~triqs_Nevanlinna.kernels_core.NevanlinnaKernel` and
:class:`~triqs_Nevanlinna.kernels_core.CaratheodoryKernel` classes from the
:mod:`~triqs_Nevanlinna.kernels_core` c2py extension module, and imports
:mod:`triqs.utility.mpi` so that MPI is initialized before the kernels use it internally.
"""
import triqs.utility.mpi as mpi
from .kernels_core import *
