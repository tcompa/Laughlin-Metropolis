import numpy
import sys
sys.path.append('..')
from lib_laughlin_metropolis import main_laughlin_mc


def test_start_from_scratch():
    N = 11
    m = 2.0
    Nqh = 3
    xqh = numpy.random.uniform(-1.0, 1.0, (Nqh, 2))
    delta = 0.1
    nsteps = 1000
    main_laughlin_mc(N, m, Nqh, xqh, delta, nsteps,
                     skip_for_rsq=0, skip_for_xy_hist=0,
                     ContinuePreviousRun=0, xmax=100.0)


def test_continue_previous_run():
    N = 11
    m = 2.0
    Nqh = 3
    xqh = numpy.random.uniform(-1.0, 1.0, (Nqh, 2))
    delta = 0.1
    nsteps = 1000
    main_laughlin_mc(N, m, Nqh, xqh, delta, nsteps,
                     skip_for_rsq=0, skip_for_xy_hist=0,
                     ContinuePreviousRun=0, xmax=100.0)
    main_laughlin_mc(N, m, Nqh, xqh, delta, nsteps,
                     skip_for_rsq=0, skip_for_xy_hist=0,
                     ContinuePreviousRun=1)
