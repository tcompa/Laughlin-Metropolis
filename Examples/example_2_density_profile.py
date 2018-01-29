from __future__ import print_function

import json
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
sys.path.append('..')
from lib_laughlin_metropolis import main_laughlin_mc


# Part 1/3: Monte Carlo run
## Define parameters
N = 20                                        # Number of particles
m = 2.0                                       # Particle charge (in the plasma analogy)
Nqh = 2                                       # Number of quasiholes
xqh = numpy.array([[-2.5, 0.0], [2.5, 0.0]])  # Coordinates of quasiholes
delta = 0.4                                   # Maximum displacement in Monte Carlo moves
ID = 'N%03i_Nqh%i_m%f' % (N, Nqh, m)          # Identifier of this run
nsteps = 500000                               # Number of Monte Carlo moves
print('Now starting the Monte Carlo run.\nWith nsteps=500000, this should ' + 
      'take about 10-15 seconds on a laptop.\nFor smaller values of nsteps ' +
      'the calculation is faster (but the histogram is more noisy).')
## Thermalization run
main_laughlin_mc(N, m, Nqh, xqh, delta, 10000,
                 skip_for_rsq=0, skip_for_xy_hist=0,
                 ContinuePreviousRun=0, xmax=30.0)
## Production run
main_laughlin_mc(N, m, Nqh, xqh, delta, nsteps,
                 skip_for_rsq=0, skip_for_xy_hist=5,
                 ContinuePreviousRun=1)

# Part 2/3: Load Monte Carlo data and build density-profile histogram
## Load histogram data
H = numpy.loadtxt('data_%s_xy_hist.dat' % ID).T
with open('data_%s_xy_hist_params.json' % ID, 'r') as parameter_file:
    H_params = json.load(parameter_file)
    binwidth = H_params['binwidth']
    xmax_hist = H_params['xmax_hist']
    nbins = H_params['nbins']
xedges = numpy.linspace(-xmax_hist, xmax_hist, nbins + 1)
yedges = numpy.linspace(-xmax_hist, xmax_hist, nbins + 1)
## Find empty bins and set the corresponding histogram values to NaN
indices_empty_bins = H < 1
H[indices_empty_bins] = numpy.nan
## Normalize density profile
tot_measures = numpy.nansum(H)
H /= (float(tot_measures) * binwidth ** 2)
H *= N
## Rescale all lengths by lB=1/sqrt(2)
lB = 1.0 / numpy.sqrt(2.0)
xedges /= lB
yedges /= lB
xqh /= lB
H *= (lB ** 2)
print('Note: All length scales are being rescaled by lB=1/sqrt(2).')

# Part 3/3: Show density-profile histogram
## Plot histogram
fig, ax = plt.subplots(1, 1, figsize=(5.0, 4.0))
hb = ax.imshow(H, interpolation='nearest', origin='low', vmin=0.0,
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
## Add colorbar
divider = make_axes_locatable(ax)
cbar = plt.colorbar(hb, cax=divider.append_axes('right', '5%', pad='3%'))
## Plot quasiholes
ax.scatter(xqh[:, 0], xqh[:, 1], c='r', edgecolor='k', s=30, zorder=10)
## Finalize plot
xmax_plot = 12.5
ax.set_xlim(-xmax_plot, xmax_plot)
ax.set_ylim(-xmax_plot, xmax_plot)
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel('$x / l_B$')
ax.set_ylabel('$y / l_B$')
ax.set_title('$N=%i$, $m=%i$' % (N, m))
plt.savefig('fig_example_2_density_profile.pdf', bbox_inches='tight')
plt.show()
