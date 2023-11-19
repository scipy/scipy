import math
import os
from os.path import relpath, dirname
import re
import sys
import warnings
from datetime import date
from docutils import nodes
from docutils.parsers.rst import Directive

import matplotlib
import matplotlib.pyplot as plt
from numpydoc.docscrape_sphinx import SphinxDocString
from sphinx.util import inspect

import scipy
from scipy._lib._util import _rng_html_rewrite
# Workaround for sphinx-doc/sphinx#6573
# ua._Function should not be treated as an attribute
import scipy._lib.uarray as ua
from scipy.stats._distn_infrastructure import rv_generic
from scipy.stats._multivariate import multi_rv_generic


old_isdesc = inspect.isdescriptor
inspect.isdescriptor = (lambda obj: old_isdesc(obj)
                        and not isinstance(obj, ua._Function))

# Currently required to build scipy.fft docs
os.environ['_SCIPY_BUILDING_DOC'] = 'True'

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
    'sphinx_design',
    'scipyoptdoc',
    'doi_role',
    'matplotlib.sphinxext.plot_directive',
    'myst_nb',
]


# Do some matplotlib config in case users have a matplotlibrc that will break
# things
matplotlib.use('agg')
plt.ioff()

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
version = re.sub(r'\.dev.*$', r'.dev', scipy.__version__)
release = version

if os.environ.get('CIRCLE_JOB', False) and \
        os.environ.get('CIRCLE_BRANCH', '') != 'main':
    version = os.environ['CIRCLE_BRANCH']
    release = version

print(f"{project} (VERSION {version})")

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
# docutils warnings when using notebooks (see gh-17322)
# these will hopefully be removed in the near future
for key in (
    r"The frontend.OptionParser class will be replaced",
    r"The frontend.Option class will be removed",
    ):
    warnings.filterwarnings('ignore', message=key, category=DeprecationWarning)
warnings.filterwarnings(
    'ignore',
    message=r'.*is obsoleted by Node.findall()',
    category=PendingDeprecationWarning,
)
warnings.filterwarnings(
    'ignore',
    message=r'There is no current event loop',
    category=DeprecationWarning,
)
# TODO: remove after gh-19228 resolved:
warnings.filterwarnings(
    'ignore',
    message=r'.*path is deprecated.*',
    category=DeprecationWarning,
)

# -----------------------------------------------------------------------------
# HTML output
# -----------------------------------------------------------------------------

html_theme = 'pydata_sphinx_theme'

html_logo = '_static/logo.svg'
html_favicon = '_static/favicon.ico'

html_theme_options = {
  "github_url": "https://github.com/scipy/scipy",
  "twitter_url": "https://twitter.com/SciPy_team",
  "navbar_end": ["theme-switcher", "version-switcher", "navbar-icon-links"],
  "switcher": {
      "json_url": "https://scipy.github.io/devdocs/_static/version_switcher.json",
      "version_match": version,
  }
}

if 'dev' in version:
    html_theme_options["switcher"]["version_match"] = "development"

if 'versionwarning' in tags:  # noqa: F821
    # Specific to docs.scipy.org deployment.
    # See https://github.com/scipy/docs.scipy.org/blob/main/_static/versionwarning.js_t
    src = ('var script = document.createElement("script");\n'
           'script.type = "text/javascript";\n'
           'script.src = "/doc/_static/versionwarning.js";\n'
           'document.head.appendChild(script);')
    html_context = {
        'VERSIONCHECK_JS': src
    }
    html_js_files = ['versioncheck.js']

html_title = f"{project} v{version} Manual"
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
    'statsmodels': ('https://www.statsmodels.org/stable', None),
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

# maps functions with a name same as a class name that is indistinguishable
# Ex: scipy.signal.czt and scipy.signal.CZT or scipy.odr.odr and scipy.odr.ODR
# Otherwise, the stubs are overwritten when the name is same for
# OS (like MacOS) which has a filesystem that ignores the case
# See https://github.com/sphinx-doc/sphinx/pull/7927
autosummary_filename_map = {
    "scipy.odr.odr": "odr-function",
    "scipy.signal.czt": "czt-function",
}


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
import warnings
for key in (
        'lsim2 is deprecated',  # Deprecation of scipy.signal.lsim2
        'impulse2 is deprecated',  # Deprecation of scipy.signal.impulse2
        'step2 is deprecated',  # Deprecation of scipy.signal.step2
        'interp2d` is deprecated',  # Deprecation of scipy.interpolate.interp2d
        'scipy.misc',  # scipy.misc deprecated in v1.10.0; use scipy.datasets
        'kurtosistest only valid',  # intentionally "bad" excample in docstring
        'scipy.signal.daub is deprecated',
        'scipy.signal.qmf is deprecated',
        'scipy.signal.cascade is deprecated',
        'scipy.signal.morlet is deprecated',
        'scipy.signal.morlet2 is deprecated',
        'scipy.signal.ricker is deprecated',
        'scipy.signal.cwt is deprecated',
        ):
    warnings.filterwarnings(action='ignore', message='.*' + key + '.*')

import numpy as np
np.random.seed(123)
"""

plot_include_source = True
plot_formats = [('png', 96)]
plot_html_show_formats = False
plot_html_show_source_link = False

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
# Notebook tutorials with MyST-NB
# -----------------------------------------------------------------------------

nb_execution_mode = "auto"

# -----------------------------------------------------------------------------
# Source code links
# -----------------------------------------------------------------------------

# Not the same as from sphinx.util import inspect and needed here
import inspect  # noqa: E402

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
    while hasattr(obj, "__wrapped__"):
        obj = obj.__wrapped__
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
        base_url = "https://github.com/scipy/scipy/blob"
        if m:
            return f"{base_url}/{m.group(1)}/{fn}{linespec}"
        elif 'dev' in scipy.__version__:
            return f"{base_url}/main/{fn}{linespec}"
        else:
            return f"{base_url}/v{scipy.__version__}/{fn}{linespec}"
    else:
        return None


# Tell overwrite numpydoc's logic to render examples containing rng.
SphinxDocString._str_examples = _rng_html_rewrite(
    SphinxDocString._str_examples
)


class LegacyDirective(Directive):
    """
    Adapted from docutils/parsers/rst/directives/admonitions.py

    Uses a default text if the directive does not have contents. If it does,
    the default text is concatenated to the contents.

    """
    has_content = True
    node_class = nodes.admonition
    optional_arguments = 1

    def run(self):
        try:
            obj = self.arguments[0]
        except IndexError:
            # Argument is empty; use default text
            obj = "submodule"
        text = (f"This {obj} is considered legacy and will no longer receive "
                "updates. This could also mean it will be removed in future "
                "SciPy versions.")

        try:
            self.content[0] = text+" "+self.content[0]
        except IndexError:
            # Content is empty; use the default text
            source, lineno = self.state_machine.get_source_and_line(
                self.lineno
            )
            self.content.append(
                text,
                source=source,
                offset=lineno
            )
        text = '\n'.join(self.content)
        # Create the admonition node, to be populated by `nested_parse`
        admonition_node = self.node_class(rawsource=text)
        # Set custom title
        title_text = "Legacy"
        textnodes, _ = self.state.inline_text(title_text, self.lineno)
        title = nodes.title(title_text, '', *textnodes)
        # Set up admonition node
        admonition_node += title
        # Select custom class for CSS styling
        admonition_node['classes'] = ['admonition-legacy']
        # Parse the directive contents
        self.state.nested_parse(self.content, self.content_offset,
                                admonition_node)
        return [admonition_node]


def setup(app):
    app.add_directive("legacy", LegacyDirective)
