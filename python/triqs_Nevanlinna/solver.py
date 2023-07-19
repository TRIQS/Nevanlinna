
import numpy as np
import scipy.optimize as opt
from .solver_core import SolverCore

class Solver(SolverCore):

    def __init__(self, kernel='NEVANLINNA', precision=100):
        SolverCore.__init__(self, kernel=kernel, precision=precision)

    def optimize(self, grid, eta, target, nk=15, maxiter=1000, gtol=1e-3, verbose = False):
        def get_f_k(nk, z):
            '''
            Compute Hardy functions basis

            :param nk: basis dimension
            :param z: real frequency with broadening
            :return: Basis of nk Hardy function
            '''
            f_k = np.zeros([nk, z.shape[0]], dtype=np.complex128)
            # Mobius transformation
            zz = (z - 1.j)/(z + 1.j)
            prefactor = (1.0-zz)/np.sqrt(np.pi)
            for k in range(nk):
                f_k[k, :] = prefactor * (zz ** k)
            return f_k
        def update_theta(x, f_k, theta):
            akr, aki, bkr, bki = np.split(x, 4)
            ak = akr + 1.j*aki
            bk = bkr + 1.j*bki
            nk = f_k.shape[0]
            for i in range(self.size):
                theta[:, i, i] = np.dot(ak[i * nk: (i+1)*nk], f_k) + np.dot(bk[i * nk: (i+1)*nk], np.conj(f_k))
            return theta
        z = np.array([w + eta*1.j for w in grid.values()])
        # Hardy basis functions
        f_k = get_f_k(nk, z)
        theta = np.zeros([len(grid), self.size, self.size], dtype=np.complex128)
        # Initial values for Basis coefficients
        ak = np.zeros([self.size, nk], dtype=np.complex128)
        bk = np.zeros([self.size, nk], dtype=np.complex128)
        def Nevan (x):
            update_theta(x, f_k, theta)
            Gw = self.evaluate(grid, eta, theta)
            return target(Gw)
        res = opt.minimize(Nevan, x0=np.concatenate((ak.real, ak.imag, bk.real, bk.imag)) ,
                              method = "CG", options={'maxiter': maxiter, 'gtol': gtol})
        update_theta(res.x, f_k, theta)
        if verbose:
            print(res.message)
            print(f'The optimizer performed {res.nit} functions evaluations.')
        return self.evaluate(grid, eta, theta)

