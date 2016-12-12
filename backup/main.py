#from __future__ import division
import numpy as np
from scipy.sparse import spdiags, block_diag
from scipy.sparse.linalg import spsolve, isolve
from matplotlib import pyplot as plt
import copy

class rnet(object):

  def __init__(self, eta_segments=None, D=None):
    self.tau_star_crit = 0.0495 # Wong and Parker
    self.rho_s = 2650.
    self.rho = 1000.
    self.g = 9.8
    self.kappa_von_Karman = 0.407
    if D:
      self.D = D
      self.k_s = 6*self.D # rough guess --> 3.5*D_84 = 6*D_50?
    self.lambda_p = 0.35
    self.z_0 = self.k_s/30. # hydraulically rough flow
    if eta_segments:
      self.eta = eta_segments
      self.nsegments = len(self.eta)
      self.eta_iter = copy.deepcopy(self.eta)
    if D:
      pass
    else:
      print "Choose a single grain size, D!"

  def sediment__discharge_per_unit_width(self, Si_values=None, i=None):
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
    if Si_values is not None:
      if type(Si_values) is list:
        pass
      else:
        Si_values = [Si_values]
    else:
      Si_values = range(len(self.eta))
    z = self.eta_iter
    q_s = []
    for Si in Si_values:
      print "Si", Si
      if i is not None:
        _S = -(z[Si][i+1] - z[Si][i-1]) / (self.x[Si][i+1] - self.x[Si][i-1]) # = -dz_dx
        q_s_inner = np.abs(self.h[Si][i] * _S/self.D) - 0.0816
        if type(q_s_inner) is list:
          q_s_inner[q_s_inner < 0] = 0 # no transport if less than threshold
        else:
          if q_s_inner < 0:
            q_s_inner = 0
        q_s_i = np.sign(_S) * 7.55 * q_s_inner**1.5
        q_s.append(q_s_i)
      else:
        _S = -(z[Si][2:] - z[Si][:-2]) / (self.x[Si][2:] - self.x[Si][:-2]) # = -dz_dx
        q_s_inner = np.abs(self.h[Si][1:-1] * _S/self.D) - 0.0816
        if type(q_s_inner) is list or type(q_s_inner) is np.ndarray:
          q_s_inner[q_s_inner < 0] = 0 # no transport if less than threshold
        else:
          if q_s_inner < 0:
            q_s_inner = 0
        q_s_i = np.sign(_S) * 7.55 * q_s_inner**1.5
        q_s.append(q_s_i)
    if len(q_s) == 1:
      q_s = q_s[0]
    return q_s

  def sediment__discharge_per_unit_width_at_link(self, h, z, x):
    _S = -(float(z[1]) - float(z[0])) / (float(x[1]) - float(x[0])) # = -dz_dx
    q_s_inner = np.abs(h*_S/self.D) - 0.0816
    if q_s_inner < 0:
      q_s_inner = 0 # no transport if less than threshold
    q_s = np.sign(_S) * 7.55 * q_s_inner**1.5
    return q_s
    
  def channel_closure_Parker(self, bc_type):
    self.boundary_conditions__copy_arrays()
    # Get elevations outside of domain
    if bc_type == q_s:
      self.boundary_conditions__sediment_discharge_transport_slope__network()
    if bc_type == Q:
      self.boundary_conditions__constant_slope()
    # Discharge through network
    self.discharge__in_network()
    # Channel geometry through network:
    self.channel_depth_Parker__in_network()
    self.channel_velocity_Parker__in_network()
    self.channel_width_Parker__in_network()
    
  def channel_depth_Parker__in_network(self):
    self.h = []
    for Si in range(len(eta)):
      h_Si = 1.2 * self.tau_star_crit * (self.rho_s - self.rho)/self.rho \
          * self.D / self.S
      self.h.append(h_Si)

  def channel_velocity_Parker__in_network(self):
    if self.h:
      self.u = []
      for Si in range(len(eta)):
        ustar_Si = (self.g * self.h[Si] * self.S)
        u_Si = ustar_Si/self.kappa_von_Karman * np.log(h[Si]/(np.e*self.z_0))
        self.u.append(u_Si)
    else:
      print "Write way to solve for velocity with h implicit"

  def channel_width_Parker__in_network(self):
    self.b = []
    for Si in range(len(eta)):
      b_Si = self.Q[Si] / (self.h[Si] * self.u[Si])
      self.b.append(b_Si)
  
  def boundary_condition__water_discharge(self, Si, Q=None):
    pass
    #if len(self.eta_with_bc[Si]) == len(self.eta[Si]):
    #  Q[]
  
  def discharge__in_network(self, Qin):
    """
    Sum discharges through network
    """
    pass
    #for Qi in Qin
    #self.Q = 
    """
    headwaters_segments = []
    for Si in range(len(self.eta)):
      headwaters_segments_Si = []
      Si_up = Si # copies value for integers
      for Si_up in next_upstream_segments:
        next_upstream_segments = flow_from_to[:,0][flow_from_to[:,1] == Si]
        if Si
      headwaters_segments.append()
    for segment, Q in self.segment_Q_in:
      self.Qin
    """
    # Mapping
    self.headwaters_flow_to_segments = []
    for Si_in, Q in self.headwaters_segments:
      _headwaters_flow_to_segments_Si = []
      to_segment = [Si_in]
      while len(to_segment):
        _headwaters_flow_to_segments_Si += list(to_segment)
        to_segment = self.flow_from_to[:,1][self.flow_from_to[:,0] == to_segment]
      self.headwaters_flow_to_segments.append(_headwaters_flow_to_segments_Si)
    # Compute discharge
    self.Q_in_segment = np.zeros(nsegments)
    for HWi in range(len(self.headwaters_flow_to_segments)):
      for Si in self.headwaters_flow_to_segments[HWi]:
        self.Q_in_segment[Si] += self.headwaters_segments[HWi,1]
    
  def boundary_conditions__copy_arrays(self):
    self.eta_with_bc = list(self.eta_iter)
    self.h_with_bc = list(self.h)
    self.b_with_bc = list(self.b)

  def boundary_conditions__constant_slope(self):
    """
    Used with channel geometry boundary conditions: start by considering 
    constant slope at boundaries (instead of changing transport slope)
    """
    if self.at_upstream_end(Si):
      self.eta_with_bc[Si] = np.hstack((eta_with_bc[Si][0], self.eta_with_bc[Si]))
    if self.at_downstream_end(Si):
      self.eta_with_bc[Si] = np.hstack((eta_with_bc[Si][0], self.eta_with_bc[Si]))

  def slope_at_each_cell_and_padded_arrays_with_ghost_cells(self):
    """
    Uses boundary conditions and internal adjacencies to compute slopes
    """
    self.S = []
    self.eta_with_boundary_values = []
    for Si in range(len(self.eta)):
      eta_with_bc_Si = self.eta_with_bc[Si]
      # Flowing to another segment?
      if (self.flow_from_to[:,0] == Si).any():
        to_Si = self.flow_from_to[:,1][self.flow_from_to[:,0] == Si]
        boundary_value = self.eta_iter[to_Si][0]
        eta_with_bc_Si = np.hstack((eta_with_bc_Si, boundary_value))
      # Getting flow from upstream segments?
      if (self.flow_from_to[:,1] == Si).any():
        _q_s = 0
        for from_Si in self.flow_from_to[:,0][self.flow_from_to[:,1] == Si]:
          #boundary_values.append(self.eta_iter[from_Si][-1])
          # NO, TRY A VALUE THAT IS REPRESENTATIVE OF THE AMOUNT OF SEDIMENT 
          # THAT ENTERS (NO TO THE BELOW)
          #_S = (self.eta_with_bc[from_Si][-1] - eta_with_bc_Si[0])/self.dx[Si]
          _q_s += self.sediment__discharge_per_unit_width(from_Si, -1)
        print _q_s
        # STRAIGHT MEAN IS NOT THE BEST WAY TO GO, BUT IS A START
        # I AVOID THIS IN MAIN SOLVER BY HAVING UPSTREAM SEGMENTS INCLUDE
        # THE LINKS... BUT HERE, FOR DEPTH/SLOPE/VELOCITY, I CAN'T...
        boundary_value = self.transport__slope(_q_s, self.h[Si][0])
        eta_with_bc_Si = np.hstack((boundary_value, eta_with_bc_Si))
      self.eta_with_boundary_values.append(eta_with_bc_Si)
      # central differencing for S
      self.S.append( -( (self.eta_with_boundary_values[Si][2:] - \
                         self.eta_with_boundary_values[Si][:-2]) ) / (2*self.dx[Si]) )
    # not yet, not with rigid q_s boundary!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #self.eta_with_bc = list(self.eta_with_boundary_values)

  def transport_slope_at_each_cell(self):
    self.S = []
    for Si in range(len(self.eta)):
      self.S.append( -( (self.eta_with_bc[Si][2:] - \
                         self.eta_with_bc[Si][:-2]) ) / (2*self.dx[Si]) )

  def boundary_conditions__sediment_discharge_transport_slope__network(self, q_s_all=None):
    self.boundary_conditions__copy_arrays()
    if q_s_all is not None:
      pass
    else:
      q_s_all = [None] * len(self.eta)
    for Si in range(len(self.eta)):
      self.boundary_condition__sediment_discharge_transport_slope(Si, q_s_all[Si])

  def boundary_condition__sediment_discharge_transport_slope(self, Si, q_s=None):
    if len(self.eta_with_bc[Si]) == len(self.eta_iter[Si]):
      if self.at_upstream_end[Si]:
        #print "Si", Si, "has an upstream boundary"
        _b = self.b[Si][0]
        _h = self.h[Si][0]
        if q_s == None:
          q_s = self.sediment__discharge_per_unit_width(Si, 0)
          #Q_s = q_s * _b
        else:
          q_s = q_s[0]
          #q_s = Q_s/_b
        self.b_with_bc[Si] = np.hstack((_b, self.b_with_bc[Si]))
        self.h_with_bc[Si] = np.hstack((_h, self.h_with_bc[Si]))
        if Si == 0:
          S_t = 1.5*self.transport__slope(q_s, _h)
        elif Si == 1:
          S_t = 1.2*self.transport__slope(q_s, _h)
        else:
          S_t = 2.*self.transport__slope(q_s, _h)
        #print S_t
        self.eta_with_bc[Si] = np.hstack((S_t*2*self.dx[Si] + \
            self.eta_with_bc[Si][1], self.eta_with_bc[Si]))
      #"""
      else:
        # edge boundary conditions -- added in once I realized that just having
        # adjacent slopes (not related to transport efficiency) didn't work.
        _q_s = 0
        for Sf in self.flow_from[Si]:
          _q_s += self.q_s_equilibrium[Si][0]
        S_t = self.transport__slope(_q_s, self.h[Si][0])
        _eta = S_t * 2 * self.dx[Si] + self.eta_iter[Si][1]
        self.eta_with_bc[Si] = np.hstack((_eta, self.eta_with_bc[Si]))
      #"""
      if self.at_downstream_end[Si]:
        #print "Si", Si, "has a downstream boundary"
        _b = self.b[Si][-1]
        _h = self.h[Si][-1]
        if q_s == None:
          q_s = self.sediment__discharge_per_unit_width(Si, -1)
          #Q_s = q_s * _b
        else:
          q_s = q_s[-1]
          #q_s = Q_s/_b
        self.b_with_bc[Si] = np.hstack((_b, self.b_with_bc[Si]))
        self.h_with_bc[Si] = np.hstack((_h, self.b_with_bc[Si]))
        S_t = self.transport__slope(q_s, self.h[Si][-1])
        self.eta_with_bc[Si] = np.hstack((self.eta_with_bc[Si], \
            S_t*-2*self.dx[Si] + self.eta_with_bc[Si][-2]))
      # Again, realizing problem with using 
      # slope_at_each_cell_and_padded_arrays_with_ghost_cells
      # at least with this rigid sediment transport-based formulation:
      # must keep transport slope as export from each river sub-catchment
      # Maybe change this later when allowing b, u, h to vary
      #"""
      else:
        S_t = self.transport__slope(self.q_s_equilibrium[Si][-1], self.h[Si][-1])
        _eta = -S_t * 2 * self.dx[Si] + self.eta_iter[Si][-2]
        self.eta_with_bc[Si] = np.hstack((self.eta_with_bc[Si], _eta))
      #"""
    else:
      print "WARNING: boundary condition has probably already been created"
      print "Doing nothing."


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
    
    #A1 = (- ( (self.h/self.D) * detatmp_dx ) - 0.0816)**.5
    # HAVE TO CHECK ABS TO LET UPSTREAM QS HAPPEN
    A1_inside = np.abs(self.h[Si]*-self.S[Si]/self.D) - 0.0816
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
    
    # Now populate tridiagonal
    l1 = np.roll(l1, -1)
    r1 = np.roll(r1, 1)

    diags = np.vstack((l1,c0,r1))
    offsets = np.array([-1,0,1])
    
    self.Adt.append(Adt)

    coeff_matrix_block = spdiags(diags, offsets, self.nx, self.nx, format='csr')
    
    return coeff_matrix_block
  
  def upstream_downstream_end_segments(self):
    # 1. List whether these are upstream-most or downstream-most ends
    self.at_upstream_end = []
    self.at_downstream_end = []
    # Si for segment (or stream) index
    for Si in range(len(self.eta)):
      self.at_upstream_end.append( (self.flow_from_to[:,1] != Si).all() )
      self.at_downstream_end.append( (self.flow_from_to[:,0] != Si).all() )
  
  def boundary_condition_sediment_flux(self):
    #########################################################
    # BC 1: SET q_s WITH SET b (and therefore h, u, for S_t #
    #########################################################
    
    self.boundary_conditions__sediment_discharge_transport_slope__network(self.q_s_equilibrium)
    #self.slope_at_each_cell_and_padded_arrays_with_ghost_cells()
    self.transport_slope_at_each_cell()
    """
    for Si in range(len(self.x)):
      # Upstream: so other river is the one that is above this one
      if self.at_upstream_end[Si]:
        if Si == 0:
          bc = 1.5*self.transport__slope(self.q_s_equilibrium[Si][0], self.h[Si][0]) #1?
        elif Si == 1:
          bc = 0.8*self.transport__slope(self.q_s_equilibrium[Si][0], self.h[Si][0]) #1?
        else:
          bc = 2*self.transport__slope(self.q_s_equilibrium[Si][0], self.h[Si][0]) #1?
      else:
        # internal boundary conditions -- what goes in, must come out
        # and at equilibrium (unforced)
        # CHANGED! USING BOUNDARY SLOPES
        _q_s = 0
        # Sf = stream from (above Si)
        for Sf in self.flow_from[Si]:
          Sf = int(Sf)
          # WIDTH ADJUSTMENT
          _q_s += self.q_s_equilibrium[Sf][-1] * self.b[Sf][-1]/self.b[Si][0]
        bc = self.transport__slope(_q_s, self.h[Si][0])
      self.S_t_in.append(bc)
      # So if at downstream end, then flowing to something
      # and Si is upstream
      if self.at_downstream_end[Si]:
        # Just transporting out as much as it gets
        bc = self.transport__slope(self.sediment__discharge_per_unit_width()[Si][-1], self.h[Si][-1]) #-2?
      else:
        # internal boundary handling
        _q_s = 0
        for St in self.flow_to[Si]:
          St = int(St)
          # To use the slope right at the margin
          _q_s += self.q_s_equilibrium[St][0]
        bc = self.transport__slope(_q_s, self.h[Si][-1])
      self.S_t_out.append(bc)
    """


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
    
    self.upstream_downstream_end_segments()
    
    ###########################
    ### BOUNDARY CONDITIONS ###
    ###########################

    # CONTRIVED TO BE AT EQUILIBRIUM ON EDGES, TEMPORARILY
    if q_s_equilibrium is not None:
      self.q_s_equilibrium = q_s_equilibrium
    else:
      # Default updates with array through time
      self.q_s_equilibrium = np.array(self.sediment__discharge_per_unit_width())
    #self.S_t_in = [] # flux_boundary_conditions_upstream
    #self.S_t_out = [] # flux_boundary_conditions_downstream
    

    ################################################################
    # BC 2: DISCHARGE AT UPSTREAM, DOWNSTREAM: SELF-FORMED CHANNEL #
    ################################################################
    #if self.bc == 'q_s_boundary_with_set_width_and_depth':
    self.boundary_condition_sediment_flux()
    #elif self.bc == 'self_adjusting_channel':
    #  pass

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
          self.coeff_matrix[Si*n, (Sf+1)*n-1] += self.coeff_matrix[(Sf+1)*n-1, (Sf+1)*n-2] * self.b[Sf][-1]/self.b[Si][0]
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
    # Reverting to old version -- concerns about looking upstream and having 
    # two separate rivers with incoming slopes
    self.transport_slope__network(self.q_s_equilibrium)
    for Si in range(len(self.eta)):
      # n = number of elements in the row
      n = len(self.eta[Si]) # n x n matrix -- assuming that all segments are equal in cell length!!!!!!!
      if self.at_upstream_end[Si]:
        ghost_left = self.eta_stack[Si*n+1] + self.S_t_in[Si]*2*self.dx[Si]
        RHS[Si*n] += self.Adt_stack[Si*n] * (self.eta_stack[Si*n+1] - ghost_left) / (self.dx[Si]**2)
      if self.at_downstream_end[Si]:
        ghost_right = self.eta_stack[(Si+1)*n-2] - self.S_t_out[Si]*2*self.dx[Si]
        RHS[(Si+1)*n-1] -= self.Adt_stack[(Si+1)*n-1] * (ghost_right - self.eta_stack[(Si+1)*n-2]) / (self.dx[Si]**2)
      """
      if self.at_upstream_end[Si]:
        RHS[Si*n] += self.Adt_stack[Si*n] * (self.eta_with_bc[Si][2] - self.eta_with_bc[Si][0]) / (self.dx[Si]**2)
      if self.at_downstream_end[Si]:
        RHS[(Si+1)*n-1] -= self.Adt_stack[(Si+1)*n-1] * (self.eta_with_bc[Si][-1] - self.eta_stack[(Si+1)*n-2]) / (self.dx[Si]**2)
      """
    self.RHS = RHS

  def solve(self):
    eta_iter_tmp = spsolve(self.coeff_matrix, self.RHS, use_umfpack=True) 
    self.eta_iter = []
    for Si in range(len(self.eta)):
      n = len(self.eta[Si]) # n x n matrix -- assuming that all segments are equal in cell length!!!!!!!
      # ALSO DEPENDS ON SAME NUMBER OF CELLS IN ALL
      self.eta_iter.append(eta_iter_tmp[Si*n:(Si+1)*n])
      
    print self.eta_with_bc[0][-1]# - self.eta_with_bc[-1][-2]
    #del self.eta_with_bc # cleanup

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
    
  def riverplot(self, linewidth=1):
    for i in range(len(self.x)):
      plt.plot(self.x[i], self.eta[i], 'k-', linewidth=linewidth)
    for _from, _to in self.flow_from_to:
      plt.plot([self.x[_from][-1], self.x[_to][0]], [self.eta[_from][-1], self.eta[_to][0]], 'k-', linewidth=linewidth)
      
  def shear_stress_depth_slope(self):
    tau_b = []
    for row in self.h:
      #tau_b.append(1000 * 9.8 * self.h * 
      pass

  def transport_slope__network(self, q_s_equilibrium=None):
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
          bc = 1*self.transport__slope(self.q_s_equilibrium[Si][0], self.h[Si][0]) #1?
        elif Si == 1:
          bc = 1*self.transport__slope(self.q_s_equilibrium[Si][0], self.h[Si][0]) #1?
        else:
          bc = 1*self.transport__slope(self.q_s_equilibrium[Si][0], self.h[Si][0]) #1?
      else:
        # internal boundary conditions -- what goes in, must come out
        # and at equilibrium (unforced)
        # CHANGED! USING BOUNDARY SLOPES
        _q_s = 0
        # Sf = stream from (above Si)
        for Sf in self.flow_from[Si]:
          Sf = int(Sf)
          # WIDTH ADJUSTMENT
          _q_s += self.q_s_equilibrium[Sf][-1] * self.b[Sf][-1]/self.b[Si][0]
        bc = self.transport__slope(_q_s, self.h[Si][0])
      self.S_t_in.append(bc)
      # So if at downstream end, then flowing to something
      # and Si is upstream
      if self.at_downstream_end[Si]:
        # Just transporting out as much as it gets
        bc = self.transport__slope(self.sediment__discharge_per_unit_width()[Si][-1], self.h[Si][-1]) #-2?
      else:
        # internal boundary handling
        _q_s = 0
        for St in self.flow_to[Si]:
          St = int(St)
          # To use the slope right at the margin
          _q_s += self.q_s_equilibrium[St][0]
        bc = self.transport__slope(_q_s, self.h[Si][-1])
      self.S_t_out.append(bc)    
