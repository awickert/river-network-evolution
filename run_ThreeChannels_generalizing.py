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
self.nx = 1E1+1


#######################
### INPUT VARIABLES ###
#######################

# GLOBAL UNIFORM #
##################
self.D = 200E-3 # [m] [uniform so far]
porosity = lambda_p = 0.35 # [-]
n_time_steps = 10
self.flow_from_to = np.array([[0,2],[1,2]])
self.flow_from = [[], [], [0,1]]
self.flow_to = [[2], [2], []]
self.b = [20, 20, 40]
self.segment_Q_in = self.headwaters_segments = np.array([[0,40],[1,20]])
self.nsegments = len(self.flow_from)

#self.flow_from_to = np.array([[0,1]])
#self.flow_from = [[], [0]]
#self.flow_to = [[1], []]

self.flow_from_to = np.array([[0,2],[1,2],[2,4],[3,4]])
self.flow_from = [[], [], [0,1], [], [2,3]]
self.flow_to = [[2], [2], [4], [4], []]
self.b = [20, 20, 40, 20, 60]
#self.b = [20, 30, 50, 10, 60]
self.segment_Q_in = self.headwaters_segments = np.array([[0,40],[1,20],[3,50]])

"""
self.flow_from_to = np.array([[0,1]])
self.flow_from = [[], [0]]
self.flow_to = [[1], []]
#self.b = [20, 20, 40, 20, 60]
self.b = [20, 20]
"""

# 3 rivers -- would often pull them in from GIS
# Keep everything uniform for starters
xmax = 1E3
self.B = 100 * np.ones(self.nx)
S = 1E-2
self.dt = 3.15E3

self.x = []
self.dx = []
self.h = []
self.eta = []
# Multiple rivers
for Si in range(len(self.flow_to)):
  self.x.append(np.linspace(0, xmax, self.nx))
  self.dx.append(np.mean(np.diff(self.x[-1]))) # Special case of uniform grid spacing
  self.h.append(2. * np.ones(self.nx)) # specific case of 2 m depth everywhere
#self.x[-1] += self.x[-2][-1] + self.dx[-1] #Very specific to this 3-river set here
self.x[-3] += self.x[1][-1] + self.dx[-1] #Very specific to this 5-river set here
self.x[-2] += self.x[1][-1] + self.dx[-1] #Very specific to this 5-river set here
self.x[-1] += self.x[2][-1] + self.dx[-1] #Very specific to this 5-river set here
#self.x[-1] += self.x[-2][-1] + self.dx[-1] #Very specific to this 2-river set here
for row in self.x:
  self.eta.append( -S * row + np.max(self.x)*S )
  self.eta[-1] = np.round(self.eta[-1], 6) # coarse trick to rmv floating point issues

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

    
#q_s_out = whatever it has to be to transport out as much material as it receives

q_s_equilibrium = np.array(self.sediment__discharge_per_unit_width())


#print np.mean(eta)


# Ignoring for now -- for iterating

# Assuming in order: so flow_from is really irrelavant; flow_to is the important part

fig = plt.figure()
plt.ylim((0,50))
ax = plt.subplot(111)
#for row in self.eta:
#  row += 10
for ts in range(100):#self.nts:
  # 3 iterations is usually good; nothing special about it, though.
  self.eta_iter = copy.deepcopy(self.eta) # For iteration
  self.stack_vars()
  for iter_i in range(5):
    self.build_coeff_matrix(q_s_equilibrium)
    self.build_RHS()
    print np.max(np.hstack(self.eta_iter))
    self.solve()
  self.update()
  ax.clear()
  if ts % 5 == 0:
    self.riverplot(linewidth=2)
    #plt.ylim((0,40))
    #plt.draw()
    plt.pause(0.0001)
self.stack_vars()

#self.plot_coeff_matrix()

#plt.ylim((0,40))
self.riverplot(linewidth=4)
plt.show()
  
  

