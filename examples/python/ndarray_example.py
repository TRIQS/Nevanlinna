import h5py
import numpy as np
from triqs_Nevanlinna import kernels

kernel = kernels.NevanlinnaKernel()

with h5py.File("data/input.h5", "r") as inp:
    iw = inp["mesh"][()]
    Giw = inp["data"][()]

eta = 0.1
w = np.linspace(-6, 6, 6000) + 1.j*eta

kernel.init(iw, Giw)
kernel.evaluate(w)
