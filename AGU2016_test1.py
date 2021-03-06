import numpy as np
from scipy.sparse import spdiags, block_diag
from scipy.sparse.linalg import spsolve, isolve
from matplotlib import pyplot as plt
import copy
import time

import ThreeChannels_generalizing

reload(ThreeChannels_generalizing)
r = ThreeChannels_generalizing.rnet()
self = r

plt.ion()



# PER RIVER #
#############
self.eta = []
self.nx = int(5E1 + 1)


#######################
### INPUT VARIABLES ###
#######################

# GLOBAL UNIFORM #
##################
self.D = 50E-3 # [m] [uniform so far]
porosity = lambda_p = 0.35 # [-]
n_time_steps = 1
self.flow_from_to = np.array([[]])
self.flow_from = [[]]
self.flow_to = [[]]
self.b = [20]
self.segment_Q_in = self.headwaters_segments = np.array([[0,40]])
self.nsegments = len(self.flow_from)


xmax = 5E3
self.B = 100 * np.ones(self.nx)
S = 1E-2
self.dt = 3.15E2

self.x = []
self.dx = []
self.h = []
self.eta = []
# Multiple rivers
for Si in range(len(self.flow_to)):
  self.x.append(np.linspace(0, xmax, self.nx))
  self.dx.append(np.mean(np.diff(self.x[-1]))) # Special case of uniform grid spacing
  self.h.append(3. * np.ones(self.nx)) # specific case of 2 m depth everywhere
for row in self.x:
  self.eta.append( -S * row + np.max(self.x)*S )
  self.eta[-1] = np.round(self.eta[-1], 6) # coarse trick to rmv floating point issues

self.eta0 = copy.deepcopy(self.eta)

#########################
### DERIVED VARIABLES ###
#########################

self.nts = np.linspace(0, n_time_steps, n_time_steps+1) # start at 1 below, t0 is initial
self.A0 = []
for Si in range(len(self.x)):
  self.A0.append( 11.325 / (1 - lambda_p) * self.h[Si]/self.D )
#q_s_in = 0.69623693 # [m^3 s^{-1}]
# q_s for equilibrium in each channel; used for transport slope upstream
# boundary conditions

q_s_equilibrium = np.array(self.sediment__discharge_per_unit_width())

#fig = plt.figure()
#ax = fig.add_subplot(111)
#for row in self.eta:
#  row += 10
z_terrace_start = []
z_terrace_mid = []
q_s_terrace_start = []
q_s_terrace_mid = []
for realization in range(1):
  for ts in range(500): # self.nts
    # 3 iterations is usually good; nothing special about it, though.
    self.eta_iter = copy.deepcopy(self.eta) # For iteration
    self.stack_vars()
    #self._q_s_in = q_s_equilibrium[0][0] * 2
    self._q_s_in_multiplier = 1 + 0.1*np.random.randn() - 0.02
    #if ts < 100:
    #  self._q_s_in_multiplier = 2
    #else:
    #  self._q_s_in_multiplier = 1
    z_terrace_start.append(self.eta[0][1])
    z_terrace_mid.append(self.eta[0][self.nx/2])
    q_s_terrace_start.append(self.sediment__discharge_per_unit_width_at_link \
        (self.h[0][1], [self.eta[0][0], self.eta[0][2]], [self.x[0][0], self.x[0][2]]))
    q_s_terrace_mid.append(self.sediment__discharge_per_unit_width_at_link \
        (self.h[0][self.nx/2], [self.eta[0][self.nx/2-1], self.eta[0][self.nx/2+1]], \
         [self.x[0][self.nx/2 -1], self.x[0][self.nx/2 + 1]]))
    for iter_i in range(1):
      self.build_coeff_matrix(q_s_equilibrium)
      self.build_RHS()
      #print np.max(np.hstack(self.eta_iter))
      self.solve()
    self.update()
    #ax.clear()
    if ts % 5 == 0:
      self.riverplot(linewidth=.5, _num=1)
      #plt.ylim((0,40))
      plt.draw()
      #plt.pause(0.1)

self.stack_vars()

#self.plot_coeff_matrix()

#plt.ylim((0,40))
#self.riverplot(linewidth=.5, plot_start=True, _num=1)
#self.riverplot(linewidth=4, plot_start=True, _num=1)
#plt.show()
plt.figure(num=1)
plt.plot(self.x[0]/1000., self.eta0[0], color='r', linewidth=6)
plt.xlim((0,3.4))
plt.ylim((14,60))
plt.show()
  

