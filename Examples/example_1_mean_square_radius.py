'''
Program: example_1_mean_square_radius.py
Author: Tommaso Comparin

This program performs several Monte Carlo runs and computes the
corresponding mean square radius.
'''

from __future__ import print_function
import numpy
import sys
sys.path.append('..')
from lib_laughlin_metropolis import main_laughlin_mc


N = 10                                        # Number of particles
m = 2.0                                       # Particle charge (in the plasma analogy)
Nqh = 2                                       # Number of quasiholes
xqh = numpy.array([[0.0, 0.0], [1.0, 0.0]])   # Coordinates of quasiholes
delta = 0.5                                   # Maximum displacement in Monte Carlo moves
n_runs = 10                                   # Number of Monte Carlo runs

for run in xrange(n_runs):
    # Thermalization
    nsteps = 10000
    main_laughlin_mc(N, m, Nqh, xqh, delta, nsteps,
                     skip_for_rsq=100, skip_for_xy_hist=0,
                     ContinuePreviousRun=0, xmax=20.0)
    # Measurements
    nsteps = 200000
    main_laughlin_mc(N, m, Nqh, xqh, delta, nsteps,
                     skip_for_rsq=100, skip_for_xy_hist=0,
                     ContinuePreviousRun=1)

    # Compute average of r^2
    ID = 'N%03i_Nqh%i_m%f' % (N, Nqh, m)
    rsq = numpy.loadtxt('data_%s_rsq.dat' % ID)
    rsq_mean = rsq.mean()
    print('Run %2i, <r^2>=%f' % (run, rsq_mean))
