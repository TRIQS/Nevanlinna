"""High-level Nevanlinna analytic-continuation solver."""

import numpy as np
import scipy.optimize as opt
from .solver_core import SolverCore, NevanlinnaParametersT

class Solver(SolverCore):
    r"""
    Nevanlinna analytic-continuation solver.

    Continues a fermionic Matsubara-frequency Green's function to the real-frequency
    axis. Call :meth:`solve` with the Matsubara data, then :meth:`evaluate` or
    :meth:`optimize` to obtain the real-frequency Green's function.

    Parameters
    ----------
    kernel : {'NEVANLINNA', 'CARATHEODORY'}, optional
        Continuation kernel. ``'NEVANLINNA'`` continues each diagonal orbital
        independently; ``'CARATHEODORY'`` performs a full matrix-valued continuation.
        Default ``'NEVANLINNA'``.
    precision : int, optional
        Number of decimal digits of internal multiprecision arithmetic (only honored
        when built with MPFR support). Default ``100``.
    """

    def __init__(self, kernel='NEVANLINNA', precision=100):
        SolverCore.__init__(self, NevanlinnaParametersT(kernel=kernel, precision=precision))

    def optimize(self, grid, eta, target, nk=15, maxiter=1000, gtol=1e-3, verbose = False):
        r"""
        Optimize the spectral function over a Hardy-function basis.

        Adds a Hardy-function basis to the Nevanlinna interpolant and minimizes
        ``target`` over its coefficients, smoothing the real-frequency Green's function
        while preserving the input Matsubara data. Requires a prior call to
        :meth:`solve`.

        Parameters
        ----------
        grid : triqs.gfs.MeshReFreq
            Real-frequency mesh on which to evaluate.
        eta : float
            Lorentzian broadening added to the real frequencies.
        target : callable
            Objective to minimize, called as ``target(G_w)`` with the evaluated
            real-frequency Green's function and returning a scalar (e.g. a measure of
            negative spectral weight).
        nk : int, optional
            Number of Hardy basis functions per orbital. Default ``15``.
        maxiter : int, optional
            Maximum number of conjugate-gradient iterations. Default ``1000``.
        gtol : float, optional
            Gradient-norm tolerance for convergence. Default ``1e-3``.
        verbose : bool, optional
            If ``True``, print the optimizer status on completion. Default ``False``.

        Returns
        -------
        triqs.gfs.Gf
            Real-frequency Green's function evaluated with the optimized Hardy
            coefficients.
        """
        def get_f_k(nk, z):
            r"""
            Build the Hardy-function basis.

            Parameters
            ----------
            nk : int
                Number of basis functions.
            z : numpy.ndarray
                Complex real-frequency points (real frequency plus broadening).

            Returns
            -------
            numpy.ndarray
                Array of shape ``(nk, len(z))`` with the Hardy basis functions.
            """
            f_k = np.zeros([nk, z.shape[0]], dtype=np.complex128)
            # Mobius transformation
            zz = (z - 1.j)/(z + 1.j)
            prefactor = (1.0-zz)/np.sqrt(np.pi)
            for k in range(nk):
                f_k[k, :] = prefactor * (zz ** k)
            return f_k
        def update_theta(x, f_k, theta, nk):
            akr, aki, bkr, bki = np.split(x.reshape(self.size, 4*nk), 4, axis=1)
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
            update_theta(x, f_k, theta, nk)
            Gw = self.evaluate(grid, eta, theta)
            return target(Gw)
        res = opt.minimize(Nevan, x0=np.concatenate((ak.real, ak.imag, bk.real, bk.imag)).flatten() ,
                              method = "CG", options={'maxiter': maxiter, 'gtol': gtol})
        update_theta(res.x, f_k, theta, nk)
        if verbose:
            print(res.message)
            print(f'The optimizer performed {res.nit} functions evaluations.')
        return self.evaluate(grid, eta, theta)

