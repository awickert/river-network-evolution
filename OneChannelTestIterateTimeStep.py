import numpy as np
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve, isolve
from matplotlib import pyplot as plt

D = 35E-3 # [m]
porosity = lambda_p = 0.35 # [-]

nx = 101
h = 1. * np.ones(nx)
B = 100 * np.ones(nx)
x = np.linspace(0, 1000, nx)
dx = np.mean(np.diff(x))
eta = -1E-2 * x + np.max(x)*1E-2
eta = np.round(eta, 6) # coarse trick to rmv floating point issues
t = np.linspace(0, 10, 11) # start at 1 below, t0 is initial

A0 = 11.325 / (1 - lambda_p) * h/D

etatmp = eta.copy() # for iterating
eta_with_ghost = np.hstack((eta[1] + 0.5, eta, eta[-2] - 0.5))
deta = eta_with_ghost[2:] - eta_with_ghost[:-2]
dt = 10. 

for t in range(10):
  for i in range(5):
    # etatmp used to update coefficient: this is the nonlinearity that 
    # requires iteration
    ###################################################################
    etatmp_with_ghost = np.hstack((etatmp[1] + 0.2, etatmp, etatmp[-2] - 0.2))
    detatmp = etatmp_with_ghost[2:] - etatmp_with_ghost[:-2]
    # HAVE TO CHECK ABS TO LET UPSTREAM QS HAPPEN
    A1 = (- ( (h/D) * detatmp/(2*dx) ) - 0.0816)**.5

    # Minus for flipping eta(t) and eta(t+1)
    A = - A0 * A1
    #A = 0*A+1 # Making A linear for the moment -- none of the above matters!!!
    Adt = A * dt
    # upstream on left -- becomes on top w/ multiplication
    l1 =  Adt/dx**2
    c0 =  -2*Adt/dx**2 + 1 # +1 for eta(t+1)
    r1 =  Adt/dx**2
    r1[0]  += l1[0]
    l1[-1] += r1[-1]

    # RHS = B
    RHS = eta.copy() # (eta(t))
    RHS[0] += Adt[0] * (eta_with_ghost[2] - eta_with_ghost[0]) / (dx**2)
    RHS[-1] -= Adt[-1] * (eta_with_ghost[-1] - eta_with_ghost[-3]) / (dx**2)

    # Now populate tridiagonal
    l1 = np.roll(l1, -1)
    r1 = np.roll(r1, 1)

    diags = np.vstack((l1,c0,r1))
    offsets = np.array([-1,0,1])

    coeff_matrix = spdiags(diags, offsets, nx, nx, format='csr')
    # Eventually have to use this for iteration
    #etatmp = spsolve(coeff_matrix, RHS, use_umfpack=True)
    # round = coarse trick to rmv floating point issues
    etatmp = np.round(spsolve(coeff_matrix, RHS, use_umfpack=True), 6)
    #etatmp[1:-1] = coeff_matrix * eta[1:-1]
  eta = etatmp.copy()
  print etatmp[-1]

plt.plot(eta)
plt.show()
