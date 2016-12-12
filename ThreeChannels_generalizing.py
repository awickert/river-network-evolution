#from __future__ import division
import numpy as np
from scipy.sparse import spdiags, block_diag
from scipy.sparse.linalg import spsolve, isolve
from matplotlib import pyplot as plt
import copy

class rnet(object):

  def __init__(self):
    pass

  def sediment__discharge_per_unit_width(self):
    """
    Compute q_s as a function of the bed elevations, grain size, and flow depth.
    This utility exists because the equation solved by the main code combines 
    MPM and Exner into a nonlinear diffusion equation and doesn't explicitly
    give transport rates
    
    Wang and Parker (2006) version of MPM (1948)
    Normal flow assumptions: depth--slope product
    
    q_s = 7.55 * ( np.abs(h[1:-1]/D * S) - 0.0816 )**1.5
    with corrections for values below the threshold of motion
    and for sdiment motion that is backwards.
    S = -dz_dx
    """
    z = self.eta
    q_s = []
    for Si in range(len(z)):
      _S = -(z[Si][2:] - z[Si][:-2]) / (self.x[Si][2:] - self.x[Si][:-2]) # = -dz_dx
      q_s_inner = np.abs(self.h[Si][1:-1] * _S/self.D) - 0.0816
      q_s_inner[q_s_inner < 0] = 0 # no transport if less than threshold
      q_s_i = np.sign(_S) * 7.55 * q_s_inner**1.5
      q_s.append(q_s_i)
    return q_s

  def sediment__discharge_per_unit_width_at_link(self, h, z, x):
    _S = -(float(z[1]) - float(z[0])) / (float(x[1]) - float(x[0])) # = -dz_dx
    q_s_inner = np.abs(h*_S/self.D) - 0.0816
    if q_s_inner < 0:
      q_s_inner = 0 # no transport if less than threshold
    q_s = np.sign(_S) * 7.55 * q_s_inner**1.5
    return q_s
    
  def transport__slope(self, q_s, h):
    """
    This transport slope, or d(eta)/dx, is the slope required for a certain
    sediment discharge per unit channel width.
    
    This is a utility to create the ghost nodes for the Neumann boundary
    condition.
    
    It will be returned positive, even though it is a negative dz/dx

    S_t = -0.26 * D/h * (q_s**(2./3.) + 0.314) --> +0.26...
    """
    S_t = np.sign(q_s) * 0.26 * self.D/h * (np.abs(q_s)**(2./3.) + 0.314)
    return S_t
    """
    if q_s is not None:
      pass
    else:
      q_s = self.q_s_in # default
    S_t = []
    for Si in range(len(self.x)):
      S_t.append( np.sign(q_s)[Si] * 0.26 * self.D/self.h[Si] * (np.abs(q_s[Si])**(2./3.) + 0.314) )
    return S_t
    """
    
  def build_coeff_matrix_section_core(self, Si):
    """
    Just the center of the coefficient matrix -- but including the outside 
    terms (collected in A) that require iteration
    """

    # Add boundary conditions, as needed
    #"""
    # OLD -- ALL JUST USING TRANSPORT SLOPE (UPSTREAM AND DOWN)
    # REGARDLESS OF WHETHER OR NOT IT IS THE GOOD CHOICE
    eta_iter__with_ghost = np.hstack((self.eta_iter[Si][1] + self.S_t_in[Si]*2*self.dx[Si], self.eta_iter[Si], self.eta_iter[Si][-2] - self.S_t_out[Si]*2*self.dx[Si]))
    #d_eta_iter__dx = (eta_iter__with_ghost[2:] - eta_iter__with_ghost[:-2]) / (2*self.dx[Si])
    """
    # NEW -- 
    eta_iter__with_ghost = self.eta_iter[Si]
    # Same as self.flow_from has len 0
    #if (self.headwaters_segments[:,0] == [Si]).any():
    # SAME, but more consistent with below
    if len(self.flow_from[Si]) == 0:
      eta_iter__with_ghost = np.hstack((self.eta_iter[Si][1] + self.S_t_in[Si]*2*self.dx[Si], eta_iter__with_ghost))
    else:
      # Compute equivalent slope for sum of sediment inputs from upstream
      # This way, one value will work in array -- will get more complicated if 
      # multiple grain sizes are used!
      q_s_in__forward_calc = []
      for flow_from_Si in self.flow_from[Si]:
        upstream_cell__eta = self.eta[flow_from_Si][-1]
        downstream_cell__eta = self.eta[Si][1] # for central differencing: this is eta[Si]
        # real slope -- not needed but calc for debug
        S_in = - (downstream_cell__eta - upstream_cell__eta)/(2*self.dx)
        # Flow depth: use a mean
        _h = np.mean((self.h[flow_from_Si][-1], self.h[Si][1]))
        _z = [upstream_cell__eta, downstream_cell__eta]
        _x = [self.x[flow_from_Si][-1], self.x[Si][1]]
        _q_s = self.sediment__discharge_per_unit_width_at_link(_h, _z, _x)
        q_s_in__forward_calc.append(_q_s)
      q_s_in__forward_calc_total = np.sum(q_s_in__forward_calc)
      # Use slope at cell for depth for transport slope
      S_t_eq = self.transport__slope(q_s_in__forward_calc_total, self.h[Si][0])
      eta_iter__with_ghost = np.hstack((self.eta_iter[Si][1] + S_t_eq*2*self.dx[Si], eta_iter__with_ghost))
    if len(self.flow_from[Si]) != 0:
      eta_iter__with_ghost = np.hstack((eta_iter__with_ghost, self.eta_iter[Si][-2] - self.S_t_out[Si]*2*self.dx[Si]))
    else:
      # Allow flow from multiple rivers
      q_s_out__forward_calc = []
      for flow_to_Si in self.flow_to[Si]:
        upstream_cell__eta = self.eta[Si][-2] # for central differencing: this is eta[Si]
        downstream_cell__eta = self.eta[flow_to_Si][0]
        # real slope -- not needed but calc for debug
        S_out = - (downstream_cell__eta - upstream_cell__eta)/(2*self.dx)
        # Flow depth: use a mean
        _h = np.mean((self.h[flow_to_Si][-1], self.h[Si][1]))
        _z = [upstream_cell__eta, downstream_cell__eta]
        _x = [self.x[flow_to_Si][-1], self.x[Si][1]]
        _q_s = self.sediment__discharge_per_unit_width_at_link(_h, _z, _x)
        q_s_out__forward_calc.append(_q_s)
      q_s_out__forward_calc_total = np.sum(q_s_out__forward_calc)
      # Use slope at cell for depth for transport slope
      # Not recovering input slope!
      S_t_eq = self.transport__slope(q_s_out__forward_calc_total, self.h[Si][0])
      eta_iter__with_ghost = np.hstack((self.eta_iter[Si][1] + S_t_eq*2*self.dx[Si], eta_iter__with_ghost))
    """
      
    """
      upstream_cell__eta = self.eta[Si][1] # for central differencing
      downstream_cell__eta = self.eta[self.flow_to[Si][0]][-1]
      S_in = (downstream_cell__eta - upstream_cell__eta)/(2*self.dt)
      # OLD ONE -- no updates, yet.
      eta_iter__with_ghost = np.hstack((self.eta_iter[Si][1] + S_in*2*self.dx[Si], self.eta_iter[Si], self.eta_iter[Si][-2] - self.S_t_out[Si]*2*self.dx[Si]))
    """
    d_eta_iter__dx = (eta_iter__with_ghost[2:] - eta_iter__with_ghost[:-2]) / (2*self.dx[Si])


    #A1 = (- ( (self.h/self.D) * detatmp_dx ) - 0.0816)**.5
    # HAVE TO CHECK ABS TO LET UPSTREAM QS HAPPEN
    #A1_inside_inside = -self.h[Si]*d_eta_iter__dx/self.D # - because MPM has slope down positive
    A1_inside = np.abs(self.h[Si]*d_eta_iter__dx/self.D) - 0.0816
    A1_inside[A1_inside < 0] = 0 # no transport
    #print A1_inside
    A1 = np.sign(A1_inside) * (A1_inside)**0.5

    # Minus for flipping eta(t) and eta(t+1)
    A = - self.A0[Si] * A1
    #A = 0*A+1 # Making A linear for the moment -- none of the above matters!!!
    Adt = A * self.dt
    # upstream on left -- becomes on top w/ multiplication
    l1 =  Adt/self.dx[Si]**2
    c0 =  -2*Adt/self.dx[Si]**2 + 1 # +1 for eta(t+1)
    r1 =  Adt/self.dx[Si]**2
    #r1[0]  += l1[0]
    #l1[-1] += r1[-1]

    ##############!!!!!!!!!!!! diff from newer code starts here
    #########print d_eta_iter__dx

    # Now populate tridiagonal
    l1 = np.roll(l1, -1)
    r1 = np.roll(r1, 1)

    diags = np.vstack((l1,c0,r1))
    offsets = np.array([-1,0,1])
    
    self.Adt.append(Adt)

    coeff_matrix_block = spdiags(diags, offsets, self.nx, self.nx, format='csr')
    
    return coeff_matrix_block
    
  def build_coeff_matrix(self, q_s_equilibrium = None):
    """
    This goes within of the Picard iteration
    
    Note that the first stream must always be upst%ream, and the last
    must always be downstream, due to the matrix edges
    
    For each of these that is false, use the array of adjacencies
    (self.flow_from_to) to map to get inputs from proper blocks
    This will mean using the first (or last) row of the new block, and placing
    a point in the second-to-last (or second-to-first) column of the block
    that it references
    """

    ##############################
    ### LINKS BETWEEN SEGMENTS ###
    ##############################
    
    # 1. List whether these are upstream-most or downstream-most ends
    self.at_upstream_end = []
    self.at_downstream_end = []
    # Si for segment (or stream) index
    for Si in range(len(self.eta)):
      try:
        self.at_upstream_end.append( (self.flow_from_to[:,1] != Si).all() )
      except:
        # Only one link?
        self.at_upstream_end.append(True)
      try:
        self.at_downstream_end.append( (self.flow_from_to[:,0] != Si).all() )
      except:
        # Only one link?
        self.at_downstream_end.append(True)

    ###########################
    ### BOUNDARY CONDITIONS ###
    ###########################

    # CONTRIVED TO BE AT EQUILIBRIUM ON EDGES, TEMPORARILY
    if q_s_equilibrium is not None:
      pass
    else:
      # Default updates with array through time
      q_s_equilibrium = np.array(self.sediment__discharge_per_unit_width())
    self.S_t_in = [] # flux_boundary_conditions_upstream
    self.S_t_out = [] # flux_boundary_conditions_downstream
    for Si in range(len(self.x)):
      # Upstream: so other river is the one that is above this one
      if self.at_upstream_end[Si]:
        if Si == 0:
          bc = 1.5*self.transport__slope(q_s_equilibrium[Si][0], self.h[Si][0]) #1?
        elif Si == 1:
          bc = 0.8*self.transport__slope(q_s_equilibrium[Si][0], self.h[Si][0]) #1?
        else:
          bc = 2*self.transport__slope(q_s_equilibrium[Si][0], self.h[Si][0]) #1?
      else:
        # internal boundary conditions -- what goes in, must come out
        # and at equilibrium (unforced)
        # CHANGED! USING BOUNDARY SLOPES
        _q_s = 0
        # Sf = stream from (above Si)
        for Sf in self.flow_from[Si]:
          Sf = int(Sf)
          #_z = [self.eta[Sf][-1], self.eta[Si][0]]
          #_x = [self.x[Sf][-1], self.x[Si][0]]
          #_h = np.mean([self.h[Sf][-1], self.h[Si][0]])
          # WIDTH ADJUSTMENT
          #_q_s += self.sediment__discharge_per_unit_width_at_link(_h, _z, _x) * self.b[Sf]/self.b[Si]
          _q_s += q_s_equilibrium[Sf][-1] * self.b[Sf]/self.b[Si]
        bc = self.transport__slope(_q_s, self.h[Si][0])
      self.S_t_in.append(bc)
      # So if at downstream end, then flowing to something
      # and Si is upstream
      if self.at_downstream_end[Si]:
        #bc = 1.2 * self.transport__slope(q_s_equilibrium[Si][-1], self.h[Si][-1]) #-2?
        # Just transporting out as much as it gets
        bc = self.transport__slope(self.sediment__discharge_per_unit_width()[Si][-1], self.h[Si][-1]) #-2?
      else:
        # internal boundary conditions -- what goes in, must come out
        # and at equilibrium (unforced)
        # CHANGED! USING BOUNDARY SLOPES
        _q_s = 0
        for St in self.flow_to[Si]:
          St = int(St)
          # To use the slope right at the margin
          #_z = [self.eta[Si][-1], self.eta[St][0]]
          #_x = [self.x[Si][-1], self.x[St][0]]
          #_h = np.mean([self.h[Si][-1], self.h[St][0]])
          #_q_s += self.sediment__discharge_per_unit_width_at_link(_h, _z, _x)# * self.b[Si]/self.b[St]
          _q_s += q_s_equilibrium[St][0]
          #_q_s /= 2. # test!!!
        bc = self.transport__slope(_q_s, self.h[Si][-1])
      self.S_t_out.append(bc)

    ################################
    ### BACKBONE OF COEFF_MATRIX ###
    ################################
    
    coeff_matrix_blocks = []
    self.Adt = []
    for Si in range(len(self.eta)):
      coeff_matrix_blocks.append(self.build_coeff_matrix_section_core(Si))
    self.coeff_matrix = block_diag(coeff_matrix_blocks, format='csr')
    self.Adt_stack = np.hstack(self.Adt)
    
    ##############################
    ### LINKS BETWEEN SEGMENTS ###
    ##############################
      
    for Si in range(len(self.eta)):
      n = len(self.eta[Si]) # n x n matrix -- assuming that all segments are equal in cell length!!!!!!!
      # Connections at upstream ends of rivers
      if self.at_upstream_end[Si]:
        # Assuming symmetry (shortcut for constant flux, not generalized!):
        self.coeff_matrix[Si*n, Si*n+1] *= 2
      else:
        from_segments = list(self.flow_from_to[:,0][self.flow_from_to[:,1] == Si])
        for Sf in from_segments:
          # Internal symmetry: always safe to assume
          # So look at next upstream on block it comes from
          # This may be a bit of a shortcut, though!
          # I'm dividing by 2 by assuming that width doubles -- so can transport twice as much sediment
          self.coeff_matrix[Si*n, (Sf+1)*n-1] += self.coeff_matrix[(Sf+1)*n-1, (Sf+1)*n-2] * self.b[Sf]/self.b[Si]
      # Connections at downstream ends of rivers
      if self.at_downstream_end[Si]:
        self.coeff_matrix[(Si+1)*n-1, (Si+1)*n-2] *= 2
      else:
        #to_segments = list(self.flow_from_to[Si,1])
        # FOR TRIBUTARIES, LIST ISN'T NEEDED. This will work though, if I 
        # generate anastomosing networks later
        #for St in to_segments:
        St = self.flow_from_to[Si,1] # to-segment
        #print St
        # Internal symmetry: always safe to assume
        # So look at next upstream on block it comes from
        # This may be a bit of a shortcut, though!
        # last row, ghost column, to first row, first column: pick up start 
        # of next river as end of previous one
        self.coeff_matrix[(Si+1)*n-1, St*n] += self.coeff_matrix[St*n, St*n+1]

      """
      ########################
      ### NOW IT IS FORMED ###
      ########################
      plt.imshow(self.coeff_matrix.todense(), interpolation='nearest'); plt.show()
      """

  def build_RHS(self):
    RHS = self.eta_stack.copy() # This is not the same as eta_iter!
    # Add boundary conditions at very ends
    # Assuming all inputs have same transport slope -- can index this in a 
    # future code.
    for Si in range(len(self.eta)):
      n = len(self.eta[Si]) # n x n matrix -- assuming that all segments are equal in cell length!!!!!!!
      if self.at_upstream_end[Si]:
        ghost_left = self.eta_stack[Si*n+1] + self.S_t_in[Si]*2*self.dx[Si]
        RHS[Si*n] += self.Adt_stack[Si*n] * (self.eta_stack[Si*n+1] - ghost_left) / (self.dx[Si]**2)
      if self.at_downstream_end[Si]:
        ghost_right = self.eta_stack[(Si+1)*n-2] - self.S_t_out[Si]*2*self.dx[Si]
        RHS[(Si+1)*n-1] -= self.Adt_stack[(Si+1)*n-1] * (ghost_right - self.eta_stack[(Si+1)*n-2]) / (self.dx[Si]**2)
    self.RHS = RHS

  def solve(self):
    eta_iter_tmp = spsolve(self.coeff_matrix, self.RHS, use_umfpack=True) 
    self.eta_iter = []
    for Si in range(len(self.eta)):
      n = len(self.eta[Si]) # n x n matrix -- assuming that all segments are equal in cell length!!!!!!!
      # ALSO DEPENDS ON SAME NUMBER OF CELLS IN ALL
      self.eta_iter.append(eta_iter_tmp[Si*n:(Si+1)*n])

  def update(self):
    self.eta = copy.deepcopy(self.eta_iter)

  def stack_vars(self):
    """
    Turn lists of variable lists into long single arrays, for use with solvers and plotting
    """
    self.eta_stack = np.hstack(self.eta)
    self.x_stack = np.hstack(self.x) # Take this out, eventually!

  def plot_coeff_matrix(self):
    plt.figure()
    cm_d = self.coeff_matrix.todense()
    cm_d[cm_d == 0] = np.nan
    plt.imshow(cm_d, interpolation='nearest')
    plt.show()
    
  def riverplot(self, linewidth=1, plot_start=False, _colors=True):
    if _colors:
      colors = ['red', 'orange', 'green', 'blue', 'purple'] * 5
    else:
      colors = ['k']*25
    plt.figure(figsize=(12,6))
    for i in range(len(self.x)):
      plt.plot(self.x[i]/1000., self.eta[i], 'k-', color=colors[i], linewidth=linewidth)
      plt.plot(self.x[i]/1000., self.eta0[i], '-', color='0.5', linewidth=2)
    try:
      for _from, _to in self.flow_from_to:
        plt.plot([self.x[_from][-1]/1000., self.x[_to][0]/1000.], [self.eta[_from][-1], self.eta[_to][0]], 'k-', color=colors[_from], linewidth=linewidth)
        if plot_start:
          plt.plot([self.x[_from][-1]/1000., self.x[_to][0]/1000.], [self.eta0[_from][-1], self.eta0[_to][0]], '-', linewidth=2, color='.5')
    except:
      plt.plot([self.x[-1]/1000., self.x[0]/1000.], [self.eta[_from][-1], self.eta[_to][0]], 'k-', color=colors[_from], linewidth=linewidth)
    plt.xlabel('$x$ [km]', fontsize=16)
    plt.ylabel('$x$ [km]', fontsize=16)
    plt.title('Alluvial long profiles', fontsize=20)
    plt.xlim((0, 3))
    plt.ylim((0, 35))
    plt.tight_layout()
    
  def shear_stress_depth_slope(self):
    for row in self.h:
      pass


  def channel_width_Parker(self):
    pass
    
  def channel_closure_Parker(self):
    self.Q = 100. # needs to be defined outside
    self.tau_start_crit = 0.03
    S = 1E-2 # replace with diff later
    h = 1.2 * self.tau_start_crit * 1.65 * self.D/S
    z0 = self.D/8. # Guess, based on idea of D84
    u = (9.8**.5/.407) * (1.2 * self.tau_start_crit * 1.65)**.5 * np.log(h/(np.e*self.D))
    b = self.Q/(h*u)
    
    
