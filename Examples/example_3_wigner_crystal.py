'''
Program: example_3_wigner_crystal.py
Author: Tommaso Comparin

This program produces a snapshot in the Wigner-crystal phase (for a
large value of the plasma-analogy charge m).
'''

from __future__ import print_function
import numpy
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
from lib_laughlin_metropolis import main_laughlin_mc


# Define parameters
N = 45
m = 120.0                               # Particle charge (in the plasma analogy)
Nqh = 0                                 # Number of quasiholes
xqh = numpy.array([[]])                 # Coordinates of quasiholes
delta = 0.2                             # Maximum displacement in Monte Carlo moves
ID = 'N%03i_Nqh%i_m%f' % (N, Nqh, m)    # Identifier of this run
nsteps = 500000                         # Number of Monte Carlo moves
print('Now starting the Monte Carlo run.\nWith nsteps=500000 and N=45,' +
      ' this should take about 40 seconds on a laptop.\n')

# Perform Monte Carlo run
main_laughlin_mc(N, m, Nqh, xqh, delta, nsteps,
                 skip_for_rsq=100, skip_for_xy_hist=0,
                 ContinuePreviousRun=0, xmax=150.0)

# Load and rescale data
x, y = numpy.loadtxt('data_%s_config.dat' % ID, unpack=True)
lB = 1.0 / numpy.sqrt(2.0)
x /= lB
y /= lB

# Show the final configuration
fig, ax = plt.subplots(1, 1)
ax.scatter(x, y)
ax.set_aspect(1)
ax.set_xlabel('$x / l_B$')
ax.set_ylabel('$y / l_B$')
ax.set_title('$N=%i$, $m=%i$' % (N, m))
plt.savefig('fig_example_3_wigner_crystal.pdf', bbox_inches='tight')
plt.show()
