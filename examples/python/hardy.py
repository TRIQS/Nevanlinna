import numpy as np
from triqs_Nevanlinna import Solver
from triqs.gf import Gf, MeshReFreq, MeshImFreq, iOmega_n, Omega, inverse

import scipy.integrate

# setup semicircular DOS parameters
J = 1.0
mu = 0.0

# setup meshes
eta = 0.05
w_mesh = MeshReFreq(-5.0, 5.0, 1001)
iw_mesh = MeshImFreq(beta=20, S='Fermion', n_iw=256)

# setup Green's functions
g_im = Gf(indices = [1], mesh= iw_mesh, name = "$G_\mathrm{imp}$")
g_re = Gf(indices = [1], mesh = w_mesh)


for iw, w in enumerate(w_mesh.values()):
    g_re[0,0].data[iw] = ((w + mu) - np.sqrt((w + mu)**2 - 4 * J**2 + 0.0j)) / (2 * J**2)

for iw, w in enumerate(iw_mesh.values()):
    g_im[0,0].data[iw] = ((w + mu) - np.sqrt((w + mu)**2 - 4 * J**2 + 0.0j)) / (2 * J**2)

# initialize solver
solver = Solver("NEVANLINNA", precision=128)
solver.solve(g_im)

# compute analytical continuation
g_re_solver = solver.evaluate(w_mesh, eta)

# print data into files



w = np.array([x for x in w_mesh.values()])
def SmoothnessOptimization (Gw):
    lmbda = 1e-4
    diff = 0.0
    # loop over all orbitals
    for i in range(Gw.data.shape[1]):
        # extract spectral function
        Aw = -Gw.data[:,i,i].imag/np.pi
        # compute estimate for second derivative
        Aw2 = (Aw[:-2] - 2*Aw[1:-1] + Aw[2:])/(Gw.mesh.delta**2)
        # evaluate normalization criteria weight
        diff += np.abs(1 - scipy.integrate.simpson(Aw, w))**2
        # evaluate smoothness criteria weight
        diff += lmbda * scipy.integrate.simpson(Aw2**2, w[1:-1])
    return diff

# perform Hardy optimization
g_re_opt = solver.optimize(w_mesh, eta, SmoothnessOptimization, gtol=3*1e-4, maxiter=2000)

www = np.array([x for x in w_mesh.values()])
ddd = -g_re[0,0].imag.data/np.pi
ddd_n = -g_re_solver[0,0].imag.data/np.pi
ddd_o = -g_re_opt[0,0].imag.data/np.pi
np.savetxt("Semicircular.txt", np.column_stack([www, ddd, ddd_n, ddd_o]))
