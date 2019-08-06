# -*- coding: utf-8 -*-
from __future__ import print_function
import sys, os, re
import glob
from datetime import date
import warnings

import numpy as np

# Currently required to build scipy.fft docs
os.environ['_SCIPY_BUILDING_DOC'] = 'True'

# Check Sphinx version
import sphinx
if sphinx.__version__ < "1.6":
    raise RuntimeError("Sphinx 1.6 or newer required")

needs_sphinx = '1.6'

# Workaround for sphinx-doc/sphinx#6573
# ua._Function should not be treated as an attribute
from sphinx.util import inspect
import scipy._lib.uarray as ua
old_isattrdesc = inspect.isattributedescriptor
inspect.isattributedescriptor = (lambda obj: old_isattrdesc(obj)
                                 and not isinstance(obj, ua._Function))

# -----------------------------------------------------------------------------
# General configuration
# -----------------------------------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.

sys.path.insert(0, os.path.abspath('../sphinxext'))
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

import numpydoc.docscrape as np_docscrape  # noqa:E402

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.intersphinx',
    'numpydoc',
    'scipyoptdoc',
    'doi_role',
    'matplotlib.sphinxext.plot_directive',
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

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General substitutions.
project = 'SciPy'
copyright = '2008-%s, The SciPy community' % date.today().year

# The default replacements for |version| and |release|, also used in various
# other places throughout the built documents.
import scipy
version = re.sub(r'\.dev-.*$', r'.dev', scipy.__version__)
release = scipy.__version__

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
show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# Ensure all our internal links work
nitpicky = True
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
        "'U' mode is deprecated",  # sphinx io
        "OpenSSL\.rand is deprecated",  # OpenSSL package in linkcheck
        "Using or importing the ABCs from",  # 3.5 importlib._bootstrap
        ):
    warnings.filterwarnings(  # deal with other modules having bad imports
        'ignore', message=".*" + key, category=DeprecationWarning)
# warnings in examples (mostly) that we allow
# TODO: eventually these should be eliminated!
for key in (
        'invalid escape sequence',  # numpydoc 0.8 has some bad escape chars
        '\nWARNING. The coefficients',  # interpolate.LSQSphereBivariateSpline
        'The integral is probably divergent',  # stats.mielke example
        'underflow encountered in square',  # signal.filtfilt underflow
        'slepian is deprecated',  # signal.windows.slepian deprecation
        'underflow encountered in multiply',  # scipy.spatial.HalfspaceIntersection
        '`frechet_l` is deprecated',  # stats.frechet_l
        '`frechet_r` is deprecated',  # stats.frechet_r
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

themedir = os.path.join(os.pardir, 'scipy-sphinx-theme', '_theme')
if os.path.isdir(themedir):
    html_theme = 'scipy'
    html_theme_path = [themedir]

    if 'scipyorg' in tags:
        # Build for the scipy.org website
        html_theme_options = {
            "edit_link": True,
            "sidebar": "right",
            "scipy_org_logo": True,
            "rootlinks": [("https://scipy.org/", "SciPy.org"),
                          ("https://docs.scipy.org/", "Docs")]
        }
    else:
        # Default build
        html_theme_options = {
            "edit_link": False,
            "sidebar": "left",
            "scipy_org_logo": False,
            "rootlinks": []
        }
        html_logo = '_static/scipyshiny_small.png'
        html_sidebars = {'index': ['indexsidebar.html', 'searchbox.html']}
else:
    # Build without scipy.org sphinx theme present
    if 'scipyorg' in tags:
        raise RuntimeError("Get the scipy-sphinx-theme first, "
                           "via git submodule init & update")
    else:
        html_style = 'scipy_fallback.css'
        html_logo = '_static/scipyshiny_small.png'
        html_sidebars = {'index': ['indexsidebar.html', 'searchbox.html']}

html_title = "%s v%s Reference Guide" % (project, version)
html_static_path = ['_static']
html_last_updated_fmt = '%b %d, %Y'

html_additional_pages = {}
html_domain_indices = True
html_copy_source = False
html_file_suffix = '.html'

htmlhelp_basename = 'scipy'

mathjax_path = "scipy-mathjax/MathJax.js?config=scipy-mathjax"


# -----------------------------------------------------------------------------
# LaTeX output
# -----------------------------------------------------------------------------

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
_stdauthor = 'Written by the SciPy community'
latex_documents = [
  ('index', 'scipy-ref.tex', 'SciPy Reference Guide', _stdauthor, 'manual'),
#  ('user/index', 'scipy-user.tex', 'SciPy User Guide',
#   _stdauthor, 'manual'),
]

# Not available on many systems:
latex_use_xindy = False

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
latex_domain_indices = False

# fix issues with Unicode characters
latex_engine = 'xelatex'

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    'preamble': r'''
% In the parameters etc. sections, align uniformly, and adjust label emphasis
\usepackage{expdlist}
\let\latexdescription=\description
\let\endlatexdescription=\enddescription
\renewenvironment{description}
{\renewenvironment{description}
   {\begin{latexdescription}%
    [\setleftmargin{50pt}\breaklabel\setlabelstyle{\bfseries}]%
   }%
   {\end{latexdescription}}%
 \begin{latexdescription}%
    [\setleftmargin{15pt}\breaklabel\setlabelstyle{\bfseries\itshape}]%
}%
{\end{latexdescription}}
% Fix bug in expdlist's modified \@item
\usepackage{etoolbox}
\makeatletter
\patchcmd\@item{{\@breaklabel} }{{\@breaklabel}}{}{}
% Fix bug in expdlist's way of breaking the line after long item label
\def\breaklabel{%
    \def\@breaklabel{%
        \leavevmode\par
        % now a hack because Sphinx inserts \leavevmode after term node
        \def\leavevmode{\def\leavevmode{\unhbox\voidb@x}}%
    }%
}
\makeatother

% Make Examples/etc section headers smaller and more compact
\titlespacing*{\paragraph}{0pt}{1ex}{0pt}

% Save vertical space in parameter lists and elsewhere
\makeatletter
\renewenvironment{quote}%
               {\list{}{\topsep=0pt\relax
                        \parsep \z@ \@plus\p@}%
                \item\relax}%
               {\endlist}
\makeatother
% Avoid small font size in code-blocks
\fvset{fontsize=auto}
% Use left-alignment per default in tabulary rendered tables
\newcolumntype{T}{L}
% Get some useful deeper bookmarks and table of contents in PDF
\setcounter{tocdepth}{1}
% Fix: ≠ is unknown to XeLaTeX's default font Latin Modern
\usepackage{newunicodechar}
\newunicodechar{≠}{\ensuremath{\neq}}
% Get PDF to use maximal depth bookmarks
\hypersetup{bookmarksdepth=subparagraph}
% reduce hyperref warnings
\pdfstringdefDisableCommands{%
  \let\sphinxupquote\empty
  \let\sphinxstyleliteralintitle\empty
  \let\sphinxstyleemphasis\empty
}
''',
    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',

    # benefit from  Sphinx built-in workaround of LaTeX's list limitations
    'maxlistdepth': '12',

    # reduce TeX warnings about underfull boxes in the index
    'printindex': r'\raggedright\printindex',

    # avoid potential problems arising from erroneous mark-up of the
    # \mathbf{\Gamma} type
    'passoptionstopackages': r'\PassOptionsToPackage{no-math}{fontspec}',
}


# -----------------------------------------------------------------------------
# Intersphinx configuration
# -----------------------------------------------------------------------------
intersphinx_mapping = {
    'python': ('https://docs.python.org/dev', None),
    'numpy': ('https://docs.scipy.org/doc/numpy', None),
    'matplotlib': ('https://matplotlib.org', None),
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

autosummary_generate = glob.glob("*.rst")

# -----------------------------------------------------------------------------
# Autodoc
# -----------------------------------------------------------------------------

autodoc_default_options = {
    'inherited-members': None,
}

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
plot_formats = [('png', 96), 'pdf']
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
            return "https://github.com/scipy/scipy/blob/master/%s%s" % (
                fn, linespec)
        else:
            return "https://github.com/scipy/scipy/blob/v%s/%s%s" % (
                scipy.__version__, fn, linespec)
    else:
        return None
