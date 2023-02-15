#!/usr/bin/env python

import unittest

from triqs_Nevanlinna import Solver
from h5 import *
import numpy as np
from triqs.gf import GfImFreq, GfReFreq, MeshReFreq, iOmega_n, Omega, inverse

class test_solver(unittest.TestCase):

    def test_solver(self):
        eta = 0.1
        m = MeshReFreq(-5.0, 5.0, 101)
        g_im = GfImFreq(indices = [1], beta = 50, n_points = 128, name = "$G_\mathrm{imp}$")
        g_re = GfReFreq(indices = [1], mesh = m)
        g_im << inverse(iOmega_n + 0.5)
        g_re << inverse(Omega + eta*1.j + 0.5)
        solver = Solver()#kernel="kernels::NEVANLINNA")
        solver.solve(g_im)
        g_re_solver = solver.evaluate(m, eta)
        self.assertTrue(np.allclose(g_re_solver.data, g_re.data))

    def test_matrix_solver(self):
        eta = 0.1
        mu = 1.0
        m = MeshReFreq(-5.0, 5.0, 101)
        g_im = GfImFreq(indices = [1,2], beta = 10, n_points = 128, name = "$G_\mathrm{imp}$")
        g_re = GfReFreq(indices = [1,2], mesh = m)
        F = np.array([[-mu, 0.5], [0.5, -mu]])
        One = np.eye(2)
        g_im << inverse(iOmega_n*One - F)
        g_re << inverse((Omega + eta*1.j)*One - F)
        solver = Solver(kernel="kernels::CARATHEODORY")
        solver.solve(g_im)
        g_re_solver = solver.evaluate(m, eta)
        self.assertTrue(np.allclose(g_re_solver.data, g_re.data))


if __name__ == '__main__':
    unittest.main()
