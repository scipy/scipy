# -*- coding: utf-8 -*-
import sys, os, re
import glob
from datetime import date
import warnings

import numpy as np

# Currently required to build scipy.fft docs
os.environ['_SCIPY_BUILDING_DOC'] = 'True'

# Check Sphinx version
import sphinx
if sphinx.__version__ < "2.0":
    raise RuntimeError("Sphinx 2.0 or newer required")

needs_sphinx = '2.0'

# Workaround for sphinx-doc/sphinx#6573
# ua._Function should not be treated as an attribute
from sphinx.util import inspect
import scipy._lib.uarray as ua
from scipy.stats._distn_infrastructure import rv_generic  # noqa: E402
from scipy.stats._multivariate import multi_rv_generic  # noqa: E402
old_isdesc = inspect.isdescriptor
inspect.isdescriptor = (lambda obj: old_isdesc(obj)
                        and not isinstance(obj, ua._Function))

# -----------------------------------------------------------------------------
# General configuration
# -----------------------------------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

import numpydoc.docscrape as np_docscrape  # noqa:E402

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.intersphinx',
    'numpydoc',
    'sphinx_panels',
    'scipyoptdoc',
    'doi_role',
    'matplotlib.sphinxext.plot_directive',
    'sphinx_tabs.tabs',
]

# Determine if the matplotlib has a recent enough version of the
# plot_directive.
from matplotlib.sphinxext import plot_directive
if plot_directive.__version__ < 2:
    raise RuntimeError("You need a recent enough version of matplotlib")
# Do some matplotlib config in case users have a matplotlibrc that will break
# things
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.ioff()

# sphinx-panels shouldn't add bootstrap css since the pydata-sphinx-theme
# already loads it
panels_add_bootstrap_css = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The main toctree document.
master_doc = 'index'

# General substitutions.
project = 'SciPy'
copyright = '2008-%s, The SciPy community' % date.today().year

# The default replacements for |version| and |release|, also used in various
# other places throughout the built documents.
import scipy
version = re.sub(r'\.dev-.*$', r'.dev', scipy.__version__)
release = scipy.__version__

if os.environ.get('CIRCLE_JOB', False) and \
        os.environ.get('CIRCLE_BRANCH', '') != 'main':
    version = os.environ['CIRCLE_BRANCH']
    release = version

print("%s (VERSION %s)" % (project, version))

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
today_fmt = '%B %d, %Y'

# List of documents that shouldn't be included in the build.
#unused_docs = []

# The reST default role (used for this markup: `text`) to use for all documents.
default_role = "autolink"

# List of directories, relative to source directories, that shouldn't be searched
# for source files.
exclude_dirs = []

# If true, '()' will be appended to :func: etc. cross-reference text.
add_function_parentheses = False

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
# show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
# pygments_style = 'sphinx'

# Ensure all our internal links work
nitpicky = True
nitpick_ignore = [
    # This ignores errors for classes (OptimizeResults, sparse.dok_matrix)
    # which inherit methods from `dict`. missing references to builtins get
    # ignored by default (see https://github.com/sphinx-doc/sphinx/pull/7254),
    # but that fix doesn't work for inherited methods.
    ("py:class", "a shallow copy of D"),
    ("py:class", "a set-like object providing a view on D's keys"),
    ("py:class", "a set-like object providing a view on D's items"),
    ("py:class", "an object providing a view on D's values"),
    ("py:class", "None.  Remove all items from D."),
    ("py:class", "(k, v), remove and return some (key, value) pair as a"),
    ("py:class", "None.  Update D from dict/iterable E and F."),
    ("py:class", "v, remove specified key and return the corresponding value."),
]

exclude_patterns = [  # glob-style

]

# be strict about warnings in our examples, we should write clean code
# (exceptions permitted for pedagogical purposes below)
warnings.resetwarnings()
warnings.filterwarnings('error')
# allow these and show them
warnings.filterwarnings('default', module='sphinx')  # internal warnings
# global weird ones that can be safely ignored
for key in (
        r"OpenSSL\.rand is deprecated",  # OpenSSL package in linkcheck
        r"distutils Version",  # distutils
        ):
    warnings.filterwarnings(  # deal with other modules having bad imports
        'ignore', message=".*" + key, category=DeprecationWarning)
warnings.filterwarnings(  # matplotlib<->pyparsing issue
    'ignore', message="Exception creating Regex for oneOf.*",
    category=SyntaxWarning)
# warnings in examples (mostly) that we allow
# TODO: eventually these should be eliminated!
for key in (
        'invalid escape sequence',  # numpydoc 0.8 has some bad escape chars
        'The integral is probably divergent',  # stats.mielke example
        'underflow encountered in square',  # signal.filtfilt underflow
        'underflow encountered in multiply',  # scipy.spatial.HalfspaceIntersection
        'underflow encountered in nextafter',  # tuterial/interpolate.rst
        # stats.skewnorm, stats.norminvgauss, stats.gaussian_kde,
        # tutorial/stats.rst (twice):
        'underflow encountered in exp',
        ):
    warnings.filterwarnings(
        'once', message='.*' + key)

# -----------------------------------------------------------------------------
# HTML output
# -----------------------------------------------------------------------------

html_theme = 'pydata_sphinx_theme'

html_logo = '_static/logo.svg'
html_favicon = '_static/favicon.ico'

html_theme_options = {
  "logo_link": "index",
  "github_url": "https://github.com/scipy/scipy",
  "navbar_end": ["version-switcher", "navbar-icon-links"],
  "switcher": {
      "json_url": "https://scipy.github.io/devdocs/_static/version_switcher.json",
      "version_match": version,
  }
}

if 'dev' in version:
    html_theme_options["switcher"]["version_match"] = "dev"

if 'versionwarning' in tags:
    # Specific to docs.scipy.org deployment.
    # See https://github.com/scipy/docs.scipy.org/blob/main/_static/versionwarning.js_t
    src = ('var script = document.createElement("script");\n'
           'script.type = "text/javascript";\n'
           'script.src = "/doc/_static/versionwarning.js";\n'
           'document.head.appendChild(script);');
    html_context = {
        'VERSIONCHECK_JS': src
    }
    html_js_files = ['versioncheck.js']

html_title = "%s v%s Manual" % (project, version)
html_static_path = ['_static']
html_last_updated_fmt = '%b %d, %Y'

html_css_files = [
    "scipy.css",
]

# html_additional_pages = {
#     'index': 'indexcontent.html',
# }
html_additional_pages = {}
html_use_modindex = True
html_domain_indices = False
html_copy_source = False
html_file_suffix = '.html'

htmlhelp_basename = 'scipy'

mathjax_path = "scipy-mathjax/MathJax.js?config=scipy-mathjax"

# -----------------------------------------------------------------------------
# Intersphinx configuration
# -----------------------------------------------------------------------------
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/devdocs', None),
    'neps': ('https://numpy.org/neps', None),
    'matplotlib': ('https://matplotlib.org/stable', None),
    'asv': ('https://asv.readthedocs.io/en/stable/', None),
}


# -----------------------------------------------------------------------------
# Numpy extensions
# -----------------------------------------------------------------------------

# If we want to do a phantom import from an XML file for all autodocs
phantom_import_file = 'dump.xml'

# Generate plots for example sections
numpydoc_use_plots = True
np_docscrape.ClassDoc.extra_public_methods = [  # should match class.rst
    '__call__', '__mul__', '__getitem__', '__len__',
]

# -----------------------------------------------------------------------------
# Autosummary
# -----------------------------------------------------------------------------

autosummary_generate = True


# -----------------------------------------------------------------------------
# Autodoc
# -----------------------------------------------------------------------------

autodoc_default_options = {
    'inherited-members': None,
}
autodoc_typehints = 'none'


# -----------------------------------------------------------------------------
# Coverage checker
# -----------------------------------------------------------------------------
coverage_ignore_modules = r"""
    """.split()
coverage_ignore_functions = r"""
    test($|_) (some|all)true bitwise_not cumproduct pkgload
    generic\.
    """.split()
coverage_ignore_classes = r"""
    """.split()

coverage_c_path = []
coverage_c_regexes = {}
coverage_ignore_c_items = {}


#------------------------------------------------------------------------------
# Matplotlib plot_directive options
#------------------------------------------------------------------------------

plot_pre_code = """
import numpy as np
np.random.seed(123)
"""
plot_include_source = True
plot_formats = [('png', 96)]
plot_html_show_formats = False
plot_html_show_source_link = False

import math
phi = (math.sqrt(5) + 1)/2

font_size = 13*72/96.0  # 13 px

plot_rcparams = {
    'font.size': font_size,
    'axes.titlesize': font_size,
    'axes.labelsize': font_size,
    'xtick.labelsize': font_size,
    'ytick.labelsize': font_size,
    'legend.fontsize': font_size,
    'figure.figsize': (3*phi, 3),
    'figure.subplot.bottom': 0.2,
    'figure.subplot.left': 0.2,
    'figure.subplot.right': 0.9,
    'figure.subplot.top': 0.85,
    'figure.subplot.wspace': 0.4,
    'text.usetex': False,
}

# -----------------------------------------------------------------------------
# Source code links
# -----------------------------------------------------------------------------

import re
import inspect
from os.path import relpath, dirname

for name in ['sphinx.ext.linkcode', 'linkcode', 'numpydoc.linkcode']:
    try:
        __import__(name)
        extensions.append(name)
        break
    except ImportError:
        pass
else:
    print("NOTE: linkcode extension not found -- no links to source generated")

def linkcode_resolve(domain, info):
    """
    Determine the URL corresponding to Python object
    """
    if domain != 'py':
        return None

    modname = info['module']
    fullname = info['fullname']

    submod = sys.modules.get(modname)
    if submod is None:
        return None

    obj = submod
    for part in fullname.split('.'):
        try:
            obj = getattr(obj, part)
        except Exception:
            return None

    # Use the original function object if it is wrapped.
    obj = getattr(obj, "__wrapped__", obj)
    # SciPy's distributions are instances of *_gen. Point to this
    # class since it contains the implementation of all the methods.
    if isinstance(obj, (rv_generic, multi_rv_generic)):
        obj = obj.__class__
    try:
        fn = inspect.getsourcefile(obj)
    except Exception:
        fn = None
    if not fn:
        try:
            fn = inspect.getsourcefile(sys.modules[obj.__module__])
        except Exception:
            fn = None
    if not fn:
        return None

    try:
        source, lineno = inspect.getsourcelines(obj)
    except Exception:
        lineno = None

    if lineno:
        linespec = "#L%d-L%d" % (lineno, lineno + len(source) - 1)
    else:
        linespec = ""

    startdir = os.path.abspath(os.path.join(dirname(scipy.__file__), '..'))
    fn = relpath(fn, start=startdir).replace(os.path.sep, '/')

    if fn.startswith('scipy/'):
        m = re.match(r'^.*dev0\+([a-f0-9]+)$', scipy.__version__)
        if m:
            return "https://github.com/scipy/scipy/blob/%s/%s%s" % (
                m.group(1), fn, linespec)
        elif 'dev' in scipy.__version__:
            return "https://github.com/scipy/scipy/blob/main/%s%s" % (
                fn, linespec)
        else:
            return "https://github.com/scipy/scipy/blob/v%s/%s%s" % (
                scipy.__version__, fn, linespec)
    else:
        return None
