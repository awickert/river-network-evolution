import numpy as np
from scipy.sparse import spdiags, block_diag
from scipy.sparse.linalg import spsolve, isolve
from matplotlib import pyplot as plt

plt.ion()

def sediment__discharge_per_unit_width(D, h, z, x):
  """
  Compute q_s as a function of the bed elevations, grain size, and flow depth.
  This utility exists because the equation solved by the main code combines 
  MPM and Exner into a nonlinear diffusion equation and doesn't explicitly give
  transport rates
  
  Wang and Parker (2006) version of MPM (1948)
  Normal flow assumptions: depth--slope product
  
  q_s = 7.55 * ( np.abs(h[1:-1]/D * S) - 0.0816 )**1.5
  with corrections for values below the threshold of motion
  and for sdiment motion that is backwards.
  S = -dz_dx
  """
  S = (z[2:] - z[:-2]) / (x[2:] - x[:-2]) # = -dz_dx
  q_s_inner = np.abs(h[1:-1]/D * S) - 0.0816
  q_s_inner[q_s_inner < 0] = 0 # no transport if less than threshold
  q_s = np.sign(q_s_inner) * 7.55 * q_s_inner**1.5
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
  S_t = np.sign(q_s) * 0.26 * D/h * (np.abs(q_s)**(2./3.) + 0.314)
  return S_t
  
def build_coeff_matrix_section(eta, etatmp, eta_with_ghost, D, h, dx, S_t_in, S_t_out):
  # etatmp used to update coefficient: this is the nonlinearity that 
  # requires iteration
  ###################################################################
  etatmp_with_ghost = np.hstack((etatmp[1] + S_t_in*2*dx, etatmp, etatmp[-2] - S_t_out*2*dx))
  detatmp_dx = (etatmp_with_ghost[2:] - etatmp_with_ghost[:-2]) / (2*dx)

  #A1 = (- ( (h/D) * detatmp_dx ) - 0.0816)**.5
  # HAVE TO CHECK ABS TO LET UPSTREAM QS HAPPEN
  A1_inside_inside = -h*detatmp_dx/D # - because MPM has slope down positive
  A1_inside = np.abs(h*detatmp_dx/D) - 0.0816
  A1_inside[A1_inside < 0] = 0 # no transport
  #print A1_inside
  A1 = np.sign(A1_inside) * (A1_inside)**0.5
  
  # Minus for flipping eta(t) and eta(t+1)
  A = - A0 * A1
  #A = 0*A+1 # Making A linear for the moment -- none of the above matters!!!
  Adt = A * dt
  # upstream on left -- becomes on top w/ multiplication
  l1 =  Adt/dx**2
  c0 =  -2*Adt/dx**2 + 1 # +1 for eta(t+1)
  r1 =  Adt/dx**2
  #r1[0]  += l1[0]
  #l1[-1] += r1[-1]

  # RHS = B
  RHS = eta.copy() # (eta(t))
  #RHS[0] += Adt[0] * (eta_with_ghost[2] - eta_with_ghost[0]) / (dx**2)
  #RHS[-1] -= Adt[-1] * (eta_with_ghost[-1] - eta_with_ghost[-3]) / (dx**2)

  # Now populate tridiagonal
  l1 = np.roll(l1, -1)
  r1 = np.roll(r1, 1)

  diags = np.vstack((l1,c0,r1))
  offsets = np.array([-1,0,1])

  coeff_matrix = spdiags(diags, offsets, nx, nx, format='csr')
  
  return coeff_matrix
  

D = 55E-3 # [m]
porosity = lambda_p = 0.35 # [-]

nx = 1E1
h = 2. * np.ones(nx)
B = 100 * np.ones(nx)
S = 1E-2
x = np.linspace(0, 1E3, nx)
dx = np.mean(np.diff(x))
eta = -S * x + np.max(x)*S
eta = np.round(eta, 6) # coarse trick to rmv floating point issues
t = np.linspace(0, 10, 11) # start at 1 below, t0 is initial

A0 = 11.325 / (1 - lambda_p) * h/D

#q_s_in = 0.69623693 # [m^3 s^{-1}]
q_s_in = sediment__discharge_per_unit_width(D, h, eta, x)[0]
#q_s_out = whatever it has to be to transport out as much material as it receives

S_t_in = transport__slope(D, h[0], q_s_in)
S_t_out = transport__slope(D, h[0], q_s_in*1.)

dt = 3.15E5

#print np.mean(eta)

eta__all_channels = np.vstack((eta, eta, eta))
etatmp__all = eta__all_channels.copy()


# Ignoring for now -- for iterating

# Assuming in order: so flow_from is really irrelavant; flow_to is the important part
flow_from_to = np.array([[0,2],[1,2]])

for t in range(1):

  coeff_matrix_blocks = []
  #S_t_out = -(eta[-1] - eta[-3])/(2*dx)
  # iterations
  #for i in range(3):

  ################################
  ### BACKBONE OF COEFF_MATRIX ###
  ################################

  for i in range(len(eta__all_channels)):
    # Do this when I have list of etas
    #eta = eta[i]
    # Use these to define which b.c.'s to apply in tridiagonal
    # Note that the first stream must always be upstream, and the last
    # must always be downstream, due to the matrix edges
    # For each of these that is false, use the array of adjacencies
    # (flow_from_to) to map to get inputs from proper blocks
    # This will mean using the first (or last) row of the new block, and placing
    # a point in the second-to-last (or second-to-first) column of the block
    # that it references
    at_upstream_end = (flow_from_to[:,1] != i).all()
    at_downstream_end = (flow_from_to[:,0] != i).all()
    # outside iteration
    etatmp = eta.copy() # for iterating
    eta_with_ghost = np.hstack((eta[1] + S_t_in*2*dx, eta, eta[-2] - S_t_out*2*dx))
    deta = eta_with_ghost[2:] - eta_with_ghost[:-2]
    # inside iteration
    coeff_matrix_blocks.append(build_coeff_matrix_section(eta, etatmp, \
      eta_with_ghost, D, h, dx, S_t_in, S_t_out))
  coeff_matrix = block_diag(coeff_matrix_blocks, format='csr')

  ##############################
  ### LINKS BETWEEN SEGMENTS ###
  ##############################

  # 1. List whether these are upstream-most or downstream-most ends
  at_upstream_end = []
  at_downstream_end = []
  # S for segment (or stream) index
  for S in range(len(eta__all_channels)):
    at_upstream_end.append( (flow_from_to[:,1] != S).all() )
    at_downstream_end.append( (flow_from_to[:,0] != S).all() )
  
  for S in range(len(eta__all_channels)):
    n = len(eta__all_channels[S]) # n x n matrix
    if at_upstream_end[S]:
      # Assuming symmetry (shortcut for constant flux, not generalized!):
      coeff_matrix[S*n, S*n+1] *= 2
    else:
      from_segments = list(flow_from_to[:,0][flow_from_to[:,1] == S])
      for Sf in from_segments:
        # Internal symmetry: always safe to assume
        # So look at next upstream on block it comes from
        # This may be a bit of a shortcut, though!
        coeff_matrix[S*n, (Sf+1)*n-1] += coeff_matrix[(Sf+1)*n-1, (Sf+1)*n-2]
    if at_downstream_end[S]:
      coeff_matrix[(S+1)*n-1, (S+1)*n-2] *= 2
    else:
      #to_segments = list(flow_from_to[S,1])
      # FOR TRIBUTARIES, LIST ISN'T NEEDED. This will work though, if I 
      # generate anastomosing networks later
      #for St in to_segments:
      St = flow_from_to[S,1] # to-segment
      # Internal symmetry: always safe to assume
      # So look at next upstream on block it comes from
      # This may be a bit of a shortcut, though!
      # last row, ghost column, to first row, first column: pick up start 
      # of next river as end of previous one
      coeff_matrix[(S+1)*n-1, St*n] += coeff_matrix[St*n, St*n+1]
    #  to_cell = flow_from_to[i,1]
    #  coeff_matrix[(i+1)*len(eta__all_channels[i])-1,to_cell*len(eta__all_channels)+1] \
    #               += coeff_matrix[(to_cell+1)*len(eta__all_channels[i])-1, \
    #                               (to_cell+1)*len(eta__all_channels[i])-2]

  ########################
  ### NOW IT IS FORMED ###
  ########################
  plt.imshow(coeff_matrix.todense(), interpolation='nearest'); plt.show()



    # Eventually have to use this for iteration
    #etatmp = spsolve(coeff_matrix, RHS, use_umfpack=True)
    # round = coarse trick to rmv floating point issues
    etatmp = spsolve(coeff_matrix, RHS, use_umfpack=True)
    #etatmp[1:-1] = coeff_matrix * eta[1:-1]
    #print etatmp[-1]
  #print ""
  eta = etatmp.copy()
print np.mean(eta)

#plt.plot(eta)
#plt.show()
