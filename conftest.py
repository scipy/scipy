from scpdt.conftest import dt_config

from contextlib import contextmanager
import warnings


@contextmanager
def warnings_errors_and_rng(test):
    """Temporarily turn (almost) all warnings to errors, also control the rng seed
    Some functions are allowed to emit specific warnings though;
    these are listed explicitly.
    """
    from scipy import integrate

    known_warnings = dict()

    # these functions are known to emit "divide by zero" RuntimeWarnings
    divide_by_zero = ['scipy.linalg.norm', 'scipy.ndimage.center_of_mass',
        'scipy.signal.bohman', 'scipy.signal.windows.bohman',
        'scipy.signal.cosine', 'scipy.signal.windows.cosine',
    ]
    for name in divide_by_zero:
        known_warnings[name] = dict(category=RuntimeWarning,
                                    message='divide by zero')

    # the funcions are known to emit IntergrationWarnings
    integration_w = ['scipy.special.ellip_normal',
                     'scipy.special.ellip_harm_2',
    ]
    for name in integration_w:
        known_warnings[name] = dict(category=integrate.IntegrationWarning,
                                    message='The occurrence of roundoff')


    # additional one-off warnings to filter
    dct = {
        'scipy.stats.anderson_ksamp':
            dict(category=UserWarning, message='p-value capped:'),
        # tutorials
        'linalg.rst':
            dict(message='the matrix subclass is not',
                 category=PendingDeprecationWarning),
        'stats.rst':
            dict(message='The maximum number of subdivisions',
                 category=integrate.IntegrationWarning),
    }
    known_warnings.update(dct)

    # Now, the mean of the matter: filter warnings according to the dict above,
    # also control the random seed for each doctest.

    # XXX: this matches the refguide-check behavior, but is a tad strange:
    # makes sure that the seed the old-fashioned np.random* methods is *NOT*
    # reproducible but the new-style `default_rng()` *IS* repoducible.
    # Should these two be either both repro or both not repro?

    from scipy._lib._util import _fixed_default_rng
    import numpy as np
    with _fixed_default_rng():
        np.random.seed(None)
        with warnings.catch_warnings(record=True) as w:
            if test.name in known_warnings:
                warnings.filterwarnings('ignore', **known_warnings[test.name])      
                yield
            else:
                warnings.simplefilter('error', Warning)
                yield


dt_config.user_context_mgr = warnings_errors_and_rng

# printing the Romberg extrap table is flaky, adds blanklines on some platforms
dt_config.stopwords.add('integrate.romb(y, show=True)')

# these names are known to fail doctesting and we like to keep it that way
# e.g. sometimes pseudocode is acceptable etc
# TODO: we need to remove duplicate functions from the list once pytest collection is tightened: 
# https://github.com/ev-br/scpdt/issues/102
dt_config.skiplist = [
    'scipy.stats.kstwobign',  # inaccurate cdf or ppf
    'scipy.stats._continuous_distns.kstwobign',
    'scipy.stats.levy_stable',
    'scipy.special.sinc',
    'scipy.special._basic.sinc',  # comes from numpy
    'scipy.fft.fftfreq',
    'scipy.fft.rfftfreq',
    'scipy.fft.fftshift',
    'scipy.fft.ifftshift',
    'scipy.fftpack.fftfreq',
    'scipy.fftpack._helper.fftfreq',
    'scipy.fftpack.fftshift',
    'scipy.fftpack._helper.fftshift',
    'scipy.fftpack.ifftshift',
    'scipy.fftpack._helper.ifftshift',
    'scipy.integrate.trapezoid',
    'scipy.integrate._quadrature.trapezoid',
    'scipy.linalg.LinAlgError',
    'scipy.linalg._misc.LinAlgError',
    'scipy.signal.bspline',
    'scipy.signal._bsplines.bspline',
    'scipy.signal.cubic',
    'scipy.signal._bsplines.cubic',
    'scipy.signal.quadratic',
    'scipy.signal._bsplines.quadratic',
    'scipy.optimize.show_options',
    'scipy.optimize._optimize.show_options',
]