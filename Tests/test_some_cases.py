import numpy
import sys
sys.path.append('..')
from lib_laughlin_metropolis import main_laughlin_mc


def run_one_specific_case(N=None, m=None, Nqh=None, xqh=None, delta=None,
                          known_rsq_mean=None, err_known_rsq_mean=None):
    '''
    Run a full Monte Carlo calculation (thermalization + measurements)
    and compare its result for <r^2> with known values. The test fails
    if the difference is larger than 5 standard errors.
    Note: Monte Carlo parameters (delta and skip_for_rsq) need to be
    chosen carefully, to avoid an underestimate of the error.
    '''

    # Thermalization run
    nsteps = 10000
    main_laughlin_mc(N, m, Nqh, xqh, delta, nsteps,
                     skip_for_rsq=100, skip_for_xy_hist=0,
                     ContinuePreviousRun=0, xmax=20.0)
    # Measurement run
    nsteps = 1000000
    main_laughlin_mc(N, m, Nqh, xqh, delta, nsteps,
                     skip_for_rsq=100, skip_for_xy_hist=0,
                     ContinuePreviousRun=1)

    # Analysis and comparison with known result
    ID = 'N%03i_Nqh%i_m%f' % (N, Nqh, m)
    rsq = numpy.loadtxt('data_%s_rsq.dat' % ID)
    rsq_mean = rsq.mean()
    variance_rsq = (rsq ** 2).mean() - rsq_mean ** 2
    err_rsq_mean = (variance_rsq / float(len(rsq - 1))) ** 0.5
    sigma_diff = abs(rsq_mean - known_rsq_mean)
    sigma_diff /= (err_rsq_mean ** 2 + err_known_rsq_mean ** 2) ** 0.5
    assert sigma_diff < 5.0


def test_case_1():
    xqh = numpy.array([[0.0, 0.0], [0.0, 0.3]])
    known_rsq_mean = 6.9951
    err_known_rsq_mean = 0.0037
    run_one_specific_case(N=5, m=2.0, Nqh=2, xqh=xqh, delta=0.6,
                          known_rsq_mean=known_rsq_mean,
                          err_known_rsq_mean=err_known_rsq_mean)


def test_case_2():
    xqh = numpy.array([[1.0, 1.0], [20.0, -20.0], [2.0, 0.0]])
    known_rsq_mean = 13.2948
    err_known_rsq_mean = 0.0041
    run_one_specific_case(N=8, m=3.0, Nqh=3, xqh=xqh, delta=0.5,
                          known_rsq_mean=known_rsq_mean,
                          err_known_rsq_mean=err_known_rsq_mean)


def test_case_3():
    xqh = numpy.array([[]])
    known_rsq_mean = 300.994
    err_known_rsq_mean = 0.021
    run_one_specific_case(N=7, m=100.0, Nqh=0, xqh=xqh, delta=0.5,
                          known_rsq_mean=known_rsq_mean,
                          err_known_rsq_mean=err_known_rsq_mean)
