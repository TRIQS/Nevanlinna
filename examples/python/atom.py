import numpy as np
from triqs_Nevanlinna import Solver
from triqs.gf import Gf, MeshReFreq, MeshImFreq, iOmega_n, Omega, inverse

# setup Hubbard atom parameters
U = 2.0
mu = U/2.0

# setup meshes
eta = 0.1
w_mesh = MeshReFreq(-5.0, 5.0, 101)
iw_mesh = MeshImFreq(beta=50, S='Fermion', n_iw=128)

# setup Green's functions
g_im = Gf(indices = [1], mesh= iw_mesh, name = "$G_\mathrm{imp}$")
g_re = Gf(indices = [1], mesh = w_mesh)

g_im << inverse(iOmega_n + mu) + inverse(iOmega_n - mu)
g_re << inverse(Omega + eta*1.j + mu) + inverse(Omega + eta*1.j - mu)

# initialize solver
solver = Solver("NEVANLINNA", precision=128)
solver.solve(g_im)

# compute analytical continuation
g_re_solver = solver.evaluate(w_mesh, eta)

# print data into files
www = np.array([x for x in w_mesh.values()])
ddd = -g_re[0,0].imag.data/np.pi
np.savetxt("Atom_ED.txt", np.column_stack([www, ddd]))

www = np.array([x for x in w_mesh.values()])
ddd = -g_re_solver[0,0].imag.data/np.pi
np.savetxt("Atom_Nevanlinna.txt", np.column_stack([www, ddd]))