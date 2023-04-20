#!/usr/bin/env python

import unittest

from triqs_Nevanlinna import Solver
from h5 import *
import numpy as np
from triqs.gf import Gf, MeshReFreq, MeshImFreq, iOmega_n, Omega, inverse

class test_solver(unittest.TestCase):

    def setUp(self):
        self.eta = 0.1
        self.w_mesh = MeshReFreq(-5.0, 5.0, 101)
        self.iw_mesh = MeshImFreq(beta=50, S='Fermion', n_iw=128)


    def test_solver(self):
        g_im = Gf(indices = [1], mesh= self.iw_mesh, name = "$G_\mathrm{imp}$")
        g_re = Gf(indices = [1], mesh = self.w_mesh)
        g_im << inverse(iOmega_n + 0.5)
        g_re << inverse(Omega + self.eta*1.j + 0.5)
        solver = Solver("kernels::NEVANLINNA", precision=128)
        solver.solve(g_im)
        g_re_solver = solver.evaluate(self.w_mesh, self.eta)
        self.assertTrue(np.allclose(g_re_solver.data, g_re.data))

    def test_matrix_solver(self):

        mu = 1.0
        g_im = Gf(indices = [1,2], mesh= self.iw_mesh, name = "$G_\mathrm{imp}$")
        g_re = Gf(indices = [1,2], mesh = self.w_mesh)
        F = np.array([[-mu, 0.5], [0.5, -mu]])
        One = np.eye(2)
        g_im << inverse(iOmega_n*One - F)
        g_re << inverse((Omega + self.eta*1.j)*One - F)
        solver = Solver(kernel="kernels::CARATHEODORY", precision=128)
        solver.solve(g_im)
        g_re_solver = solver.evaluate(self.w_mesh, self.eta)
        print(solver.Pick_eigenvalues)
        self.assertTrue(np.allclose(g_re_solver.data, g_re.data))


if __name__ == '__main__':
    unittest.main()
