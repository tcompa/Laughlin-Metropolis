#cython: language_level=3

'''
Program: lib_laughlin_metropolis.pyx
Authors: Tommaso Comparin (tommaso.comparin@unitn.it),
         Elia Macaluso (elia.macaluso@unitn.it)

This program implements Monte Carlo sampling for the classical-plasma
analogy of the Laughlin wave function. The probability distribution
which is sampled is defined in the README.md file, which also includes
a note on the physical units used in this program.
'''

import os
import sys
import time
import socket
import random

import cython
import numpy


cdef extern from 'math.h':
    double c_exp 'exp' (double)
    double c_log 'log' (double)
    double c_sqrt 'sqrt' (double)
    double c_pow 'pow' (double, double)


cdef double dist(double x1, double y1, double x2, double y2):
    ''' Compute distance between (x1, y1) and (x2, y2). '''
    return c_sqrt(c_pow(x1 - x2,  2) + c_pow(y1 - y2, 2))


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
cdef double compute_beta_energy(double [:, :] x, int N,
                                double [:, :] xqh, int Nqh,
                                double m):
    '''
    Compute beta*energy for the system configuration formed by the
    N-particle configuration x and the Nqh-quasihole configuration xqh.
    '''
    cdef double sqrt_two = c_sqrt(2.0)
    cdef double Etrap = 0.0, E2body = 0.0, Eqh = 0.0
    cdef int i_part, j_part, i_qh
    for i_part in xrange(N):
        Etrap += (c_pow(x[i_part, 0], 2) + c_pow(x[i_part, 1], 2))
        for i_qh in range(Nqh):
            dist_qh = dist(x[i_part, 0], x[i_part, 1], xqh[i_qh, 0], xqh[i_qh, 1])
            Eqh += c_log(sqrt_two * dist_qh)
        for j_part in xrange(i_part):
            dist_ij = dist(x[i_part, 0], x[i_part, 1], x[j_part, 0], x[j_part, 1])
            E2body += c_log(sqrt_two * dist_ij)
    Eqh *= (-2.0)
    E2body *= (-2.0 * m)
    return (Etrap + Eqh + E2body)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.initializedcheck(False)
def main_laughlin_mc(int N, double m, int Nqh, double [:, :] xqh, double delta,
                     int nsteps,
                     int skip_for_rsq=1,
                     int skip_for_xy_hist=1, int nbins=200, double xmax_hist=15.0,
                     double xmax=100.0,
                     int ContinuePreviousRun=1):
    '''
    Perform a Monte Carlo run.

    Parameters
    ----------
    N : int
        Number of particles
    m : double
        Particle charge in the plasma analogy (such that beta=2/m)
    Nqh : int
        Number of quasiholes
    xqh : two-dimensional array, shape (Nqh, 2)
        Coordinates of the Nqh quasiholes
    delta : double
        Maximum displacement per particle in a Monte Carlo move.
    skip_for_rsq : int 
        Number of Monte Carlo steps to skip before measuring r^2;
        if skip_for_rsq=0, r^2 is never measured.
    skip_for_xy_hist : int 
        Number of Monte Carlo steps to skip before adding a sample to
        the (x, y) histogram; if skip_for_xy_hist=0, the (x, y)
        histogram is not built.
    nbins : int
        Number of bins for the (x, y) histogram.
    xmax_hist : double
        Maximum value of (x, y) particle coordinates for the (x, y)
        histogram.
    ContinuePreviousRun : int
        If ContinuePreviousRun>0, try starting from an existing configuration.
    xmax : double
        Maximum value of (x, y) particle coordinates in the initial
        state (if not starting from an existing configuration).
    '''

    # Pick a seed for the random-number generator
    seed = random.randint(0, 131231323)
    random.seed(seed)

    # Set/initialize useful variables
    cdef int n_acc = 0
    cdef int step, i_part
    cdef double rsq, beta_energy

    # Variables for (x, y) histogram
    nbins += (nbins % 2)      # make nbins an even number
    cdef long int [:, :] hist_xy = numpy.zeros((nbins, nbins), dtype=numpy.int64)
    cdef double binwidth = 2.0 * xmax_hist / float(nbins)
    cdef int xbin, ybin

    # Prepare some output files
    ID = 'N%03i_Nqh%i_m%f' % (N, Nqh, m)
    if skip_for_rsq > 0:
        out_rsq = open('data_%s_rsq.dat' % ID, 'w')
    if skip_for_xy_hist > 0:
        with open('data_%s_xy_hist_params.json' % ID, 'w') as out_tmp:
            out_tmp.write('{\n"nbins": %i,\n "xmax_hist": %.8f,\n "binwidth": %.8f\n}' % (nbins, xmax_hist, binwidth))
    fparams = 'data_%s_params.dat' % ID
    time_start = time.process_time()
    with open(fparams, 'w') as out_params:
        out_params.write('# ' + time.strftime("%a, %d %b %Y %H:%M:%S %Z") + '\n')
        out_params.write('# Running on %s\n' % socket.gethostname())
        out_params.write('# Seed for random generator: %i\n' % seed)
        out_params.write('# nbins: %i\n' % nbins)
        out_params.write('# xmax_hist: %.8f\n' % xmax_hist)
        out_params.write('# binwidth: %.8f\n' % binwidth)

    # Declare and initialize particle configuration and its energy
    file_config = 'data_%s_config.dat' % ID
    cdef double [:, :] x = numpy.zeros((N, 2))
    cdef double[:, :] xnew = numpy.zeros((N, 2))
    if ContinuePreviousRun and os.path.isfile(file_config):
        with open(fparams, 'a') as out_params:
            out_params.write('# Reading initial configuration from %s.\n' % file_config)
        x = numpy.loadtxt(file_config)
    else:
        with open(fparams, 'a') as out_params:
            out_params.write('# Generating random initial configuration with xmax=%f.' % xmax)
        for i_part in xrange(N):
            x[i_part, 0] = random.uniform(-xmax, xmax)
            x[i_part, 1] = random.uniform(-xmax, xmax)
    cdef double betaE = compute_beta_energy(x, N, xqh, Nqh, m)
    cdef double deltabetaE

    # Monte Carlo loop
    for step in xrange(nsteps):

        # Metropolis move
        for i_part in xrange(N):
            xnew[i_part, 0] = x[i_part, 0] + random.uniform(-delta, delta)
            xnew[i_part, 1] = x[i_part, 1] + random.uniform(-delta, delta)
        deltabetaE = compute_beta_energy(xnew, N, xqh, Nqh, m) - betaE
        if deltabetaE < 0.0 or random.random() < c_exp(-deltabetaE):
            n_acc += 1
            betaE += deltabetaE
            x[:, :] = xnew

        # Take sample of r^2
        if skip_for_rsq > 0 and step % skip_for_rsq == 0:
            rsq = 0.0
            for i_part in xrange(N):
                rsq += x[i_part, 0] * x[i_part, 0] + x[i_part, 1] * x[i_part, 1]
            rsq /= float(N)
            out_rsq.write('%.10f\n' % rsq)

        # Take samples for (x, y) histogram
        if skip_for_xy_hist > 0 and step % skip_for_xy_hist == 0:
            for i_part in xrange(N):
                if abs(x[i_part, 0]) < xmax_hist and abs(x[i_part, 1]) < xmax_hist:
                    xbin = int((x[i_part, 0] + xmax_hist) / binwidth)
                    ybin = int((x[i_part, 1] + xmax_hist) / binwidth)
                    if xbin < 0 or xbin >= nbins or ybin < 0 or ybin >= nbins:
                        msg = 'ERROR in (x, y) histogram\n'
                        msg += 'x=%f\ny=%f\n' % (x[i_part, 0], x[i_part, 1])
                        msg += 'xb=%i\nyb=%i\n' % (xbin, ybin)
                        msg += 'nb=%i' % nbins
                        sys.exit(msg)
                    hist_xy[xbin, ybin] += 1

        # Update output every 1e5 steps, or when the run reached the last step
        if (step % 100000 * N == 0) or (step == nsteps - 1):
            if skip_for_rsq > 0:
                out_rsq.flush()
            if skip_for_xy_hist > 0:
                out_xy_hist = open('data_%s_xy_hist.dat' % ID, 'w')
                for xbin in xrange(nbins):
                    for ybin in xrange(nbins):
                        out_xy_hist.write('%i ' % hist_xy[xbin, ybin])
                    out_xy_hist.write('\n')
                out_xy_hist.close()

    # Save some more output
    numpy.savetxt(file_config, x)
    with open(fparams, 'a') as out_params:
        out_params.write('# Output\n')
        out_params.write('  Acceptance_ratio: %f [%i/%i]\n' % (n_acc / float(nsteps), n_acc, nsteps))
        elapsed = time.process_time() - time_start
        out_params.write('  Elapsed_time: %.2f s %.3e s/step\n#\n' % (elapsed, elapsed / nsteps))
