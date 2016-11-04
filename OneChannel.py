import numpy as np
from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve, isolve
from matplotlib import pyplot as plt

def sediment__discharge_per_unit_width(D, h, z, x):
  """
  Compute q_s as a function of the bed elevations, grain size, and flow depth.
  This utility exists because the equation solved by the main code combines 
  MPM and Exner into a nonlinear diffusion equation and doesn't explicitly give
  transport rates
  
  Wang and Parker (2006) version of MPM (1947)
  Normal flow assumptions: depth--slope product
  """
  dz_dx = (z[2:] - z[:-2]) / (x[2:] - x[:-2])
  q_s = 7.55 * ( h[1:-1]/D * -dz_dx - 0.0816 )**1.5
  return q_s

def transport__slope(D, h, q_s):
  """
  This transport slope, or d(eta)/dx, is the slope required for a certain
  sediment discharge per unit channel width.
  
  This is a utility to create the ghost nodes for the Neumann boundary
  condition.
  
  It will be returned positive, even though it is a negative dz/dx

  S_t = -0.26 * D/h * (q_s**(2./3.) + 0.314) --> +0.26...
  """
  S_t = 0.26 * D/h * (q_s**(2./3.) + 0.314)
  return S_t

D = 35E-3 # [m]
porosity = lambda_p = 0.35 # [-]

nx = 1001
h = 1. * np.ones(nx)
B = 100 * np.ones(nx)
x = np.linspace(0, 10000, nx)
dx = np.mean(np.diff(x))
eta = -1E-2 * x + np.max(x)*1E-2
eta = np.round(eta, 6) # coarse trick to rmv floating point issues
t = np.linspace(0, 10, 11) # start at 1 below, t0 is initial

A0 = 11.325 / (1 - lambda_p) * h/D

q_s_in = 10.#0.69623693
#q_s_out = whatever it has to be to transport out as much material as it receives

S_t_in = transport__slope(D, h[0], q_s_in)
S_t_out = transport__slope(D, h[0], q_s_in*2)

dt = 10000.

for t in range(5):
  #S_t_out = -(eta[-1] - eta[-3])/(2*dx)
  etatmp = eta.copy() # for iterating
  eta_with_ghost = np.hstack((eta[1] + S_t_in*2*dx, eta, eta[-2] - S_t_out*2*dx))
  deta = eta_with_ghost[2:] - eta_with_ghost[:-2]
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
