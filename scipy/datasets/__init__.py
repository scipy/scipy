"""
================================
Datasets (:mod:`scipy.datasets`)
================================

.. currentmodule:: scipy.datasets

.. autosummary::
   :toctree: generated/

   ascent - Get example image for processing
   face - Get example image for processing
   electrocardiogram - Load an example of a one-dimensional signal.


Usage of Datasets
^^^^^^^^^^^^^^^^^

SciPy dataset methods can be simply called as follows:
``'<dataset-name>()'`` This downloads the dataset files over
the network once, and saves the cache, before returning
a `numpy.ndarray` object representing the dataset.

Note that the return data structure and data type might be different
for different dataset methods. For a more detailed
example on usage, please look into the
particular dataset method documentation above.


How dataset retrieval and storage works?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SciPy dataset files are stored within individual github repositories under
the SciPy GitHub organization, following a naming convention as
``'dataset-<name>'``, for example `scipy.datasets.face` files live at
https://github.com/scipy/dataset-face.
The `scipy.datasets` submodule utilizes and depends on
`Pooch <https://www.fatiando.org/pooch/latest/>`_, a Python package built to
simplify fetching data files. Pooch uses these repos to retrieve
the respective dataset files when calling the dataset function.

A registry of all the datasets, essentially a mapping of filenames
with their SHA256 hash and repo urls are maintained,
which Pooch uses to handle and verify the downloads
on function call. After downloading the dataset once, the files
are saved in the system cache directory under ``'scipy-data'``.

Dataset cache locations may vary on different platforms.

For Mac::

    '~/Library/Caches/scipy-data'

For Unix::

    '~/.cache/scipy-data' or the value of the XDG_CACHE_HOME environment
    variable, if defined.

For Windows::

    'C:\\Users\\<user>\\AppData\\Local\\<AppAuthor>\\scipy-data\\Cache'


In environments with constrained network connectivity for various security
reasons or on systems without continuous internet connections, one may
manually load the cache of the datasets by placing the contents of the dataset
repo in the above mentioned cache directory to avoid fetching dataset errors
without the internet connectivity.

"""
import warnings

try:
    # https://github.com/scipy/scipy/pull/15607#issuecomment-1176457275
    # TODO: Remove warning filter after next certifi release
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=DeprecationWarning)
        import pooch
except ImportError:
    pooch = None

msg = (
    "Missing optional dependency 'pooch' required for scipy.datasets module. "
    "Please use pip or conda to install 'pooch'."
)
if pooch is None:
    raise ImportError(msg)


from ._fetchers import face, ascent, electrocardiogram  # noqa: E402
__all__ = ['ascent', 'electrocardiogram', 'face']


from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
