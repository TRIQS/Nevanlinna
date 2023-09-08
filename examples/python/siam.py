import sys, os
import numpy as np

from triqs.atom_diag import *

from triqs.gf import *
from triqs.operators import c, c_dag, n, dagger
from itertools import product

from triqs_Nevanlinna import Solver


beta = 40.             # Inverse temperature
mu = 1.             # Chemical potential
t = 0.3
U = [2., 0.1]      # Density-density interaction
E = [ -0.5, 0.7 ]    # Bath-site energies
V = [ 0.4, 0.5 ]    # Couplings to Bath-sites
n_orb = 2
n_orb_bath = len(E)

# triqs block structure
block_names = ['up', 'dn']
gf_struct = [ (s, n_orb) for s in block_names ]
# block structure of full Hamiltonian
gf_struct_tot = [[s, n_orb + n_orb*n_orb_bath] for s in block_names]

# setup meshes
n_iw = 800
iw_mesh = MeshImFreq(beta, 'Fermion', n_iw)
n_w = 1501
window = [-7,7]
w_mesh = MeshReFreq(window=window,n_w=n_w)
eta = 0.1

# ==== Local Hamiltonian ====
h_0 = - mu* (( n('up',0) + n('dn',0) ) + ( n('up',1) + n('dn',1) ) ) 
h_0 += t* ( c_dag('up',0) * c('up', 1) + c_dag('up',1) * c('up', 0) +
           c_dag('dn',0) * c('dn', 1) + c_dag('dn',1) * c('dn', 0) )
h_int = U[0] * (n('up',0) * n('dn',0) + n('up',1) * n('dn',1)) + 0.5*U[1] * (n('up',0) * n('dn',1) + n('up',1) * n('dn',0) +
                                                                            n('up',1) * n('up',0) + n('dn',0) * n('dn',1))
h_imp = h_0 + h_int

# ==== Bath & Coupling Hamiltonian ====
h_bath, h_coup = 0, 0
for d in range(n_orb):
    for i, (E_i, V_i) in enumerate(zip(E, V)):
        for sig in ['up','dn']:
            h_bath += (E_i)* n(sig,n_orb + d*n_orb_bath + i)
            h_coup += V_i * (c_dag(sig,d) * c(sig,n_orb + d*n_orb_bath + i) + c_dag(sig,n_orb + d*n_orb_bath + i) * c(sig,d))

# ==== Total impurity hamiltonian ====
h_tot = h_imp + h_coup + h_bath



# define fundamental operators
fop_imp  = [(s,o) for s, n in gf_struct for o in range(n)]
fop_bath = [(s,o) for s, o in product(block_names, range(n_orb, n_orb + n_orb*n_orb_bath))]
fop_tot  = fop_imp + fop_bath

# ED object
ad_tot = AtomDiag(h_tot, fop_tot)

# imaginary frequency Green's function
G_iw_tot = atomic_g_iw(ad_tot, beta, gf_struct_tot, n_iw)
block_list = [G_iw_tot[bl][:n_orb, :n_orb] for bl, n_orb in gf_struct]
G_iw = BlockGf(name_list=block_names, block_list=block_list)

# real frequency
G_w_tot = atomic_g_w(ad_tot, beta, gf_struct_tot, window,n_w, eta)
block_list = [G_w_tot[bl][:n_orb, :n_orb] for bl, n_orb in gf_struct]
G_w = BlockGf(name_list=block_names, block_list=block_list)

# setup Nevanlinna kernel solver
solver = Solver(kernel="NEVANLINNA", precision=300)

# solve 
solver.solve(G_iw['up'])

# evaluate on real frequency mesh
G_w_nvla = solver.evaluate(w_mesh, eta)


# setup Nevanlinna kernel solver
solver_C = Solver(kernel="CARATHEODORY", precision=120)

# solve 
solver_C.solve(G_iw['up'])

# evaluate on real frequency mesh
G_w_cara = solver_C.evaluate(w_mesh, eta)

www = np.array([x for x in w_mesh.values()])
#np.savetxt("SIAM_ED.txt", np.column_stack([www, G_w['up'][0,0].real.data,
#                                                G_w['up'][0,0].imag.data,
#                                                G_w['up'][0,1].real.data,
#                                                G_w['up'][0,1].imag.data]))

#www = np.array([x for x in w_mesh.values()])
#np.savetxt("SIAM_Nev.txt", np.column_stack([www, G_w_nvla[0,0].real.data,
#                                                 G_w_nvla[0,0].imag.data,
#                                                 G_w_nvla[0,1].real.data,
#                                                 G_w_nvla[0,1].imag.data]))

#np.savetxt("SIAM_Car.txt", np.column_stack([www, G_w_cara[0,0].real.data,
#                                                 G_w_cara[0,0].imag.data,
#                                                 G_w_cara[0,1].real.data,
#                                                 G_w_cara[0,1].imag.data]))


np.savetxt("SIAM_ED.txt", np.column_stack([www, G_w['up'][0,0].imag.data,
                                                G_w['up'][0,1].imag.data]))

np.savetxt("SIAM_Nev.txt", np.column_stack([www, G_w_nvla[0,0].imag.data,
                                                 G_w_nvla[0,1].imag.data]))

np.savetxt("SIAM_Car.txt", np.column_stack([www, G_w_cara[0,0].imag.data,
                                                 G_w_cara[0,1].imag.data]))



