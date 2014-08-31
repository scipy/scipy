# -*- coding: utf-8 -*-

import sys, os, re

# Check Sphinx version
import sphinx
if sphinx.__version__ < "1.1":
    raise RuntimeError("Sphinx 1.1 or newer required")

needs_sphinx = '1.1'

# -----------------------------------------------------------------------------
# General configuration
# -----------------------------------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.

sys.path.insert(0, os.path.abspath('../sphinxext'))

extensions = ['sphinx.ext.autodoc', 'sphinx.ext.pngmath', 'numpydoc',
              'sphinx.ext.intersphinx', 'sphinx.ext.coverage',
              'sphinx.ext.autosummary']

# Determine if the matplotlib has a recent enough version of the
# plot_directive.
try:
    from matplotlib.sphinxext import plot_directive
except ImportError:
    use_matplotlib_plot_directive = False
else:
    try:
        use_matplotlib_plot_directive = (plot_directive.__version__ >= 2)
    except AttributeError:
        use_matplotlib_plot_directive = False

if use_matplotlib_plot_directive:
    extensions.append('matplotlib.sphinxext.plot_directive')
else:
    raise RuntimeError("You need a recent enough version of matplotlib")

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General substitutions.
project = 'SciPy'
copyright = '2008-2014, The Scipy community'

# The default replacements for |version| and |release|, also used in various
# other places throughout the built documents.
import scipy
version = re.sub(r'\.dev-.*$', r'.dev', scipy.__version__)
release = scipy.__version__

print "Scipy (VERSION %s)" % (version,)

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
            "rootlinks": [("http://scipy.org/", "Scipy.org"),
                          ("http://docs.scipy.org/", "Docs")]
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
        html_sidebars = {'index': 'indexsidebar.html'}
else:
    # Build without scipy.org sphinx theme present
    if 'scipyorg' in tags:
        raise RuntimeError("Get the scipy-sphinx-theme first, "
                           "via git submodule init & update")
    else:
        html_style = 'scipy_fallback.css'
        html_logo = '_static/scipyshiny_small.png'
        html_sidebars = {'index': 'indexsidebar.html'}

html_title = "%s v%s Reference Guide" % (project, version)
html_static_path = ['_static']
html_last_updated_fmt = '%b %d, %Y'

html_additional_pages = {}
html_use_modindex = True
html_copy_source = False
html_file_suffix = '.html'

htmlhelp_basename = 'scipy'

pngmath_use_preview = True
pngmath_dvipng_args = ['-gamma', '1.5', '-D', '96', '-bg', 'Transparent']


# -----------------------------------------------------------------------------
# LaTeX output
# -----------------------------------------------------------------------------

# The paper size ('letter' or 'a4').
#latex_paper_size = 'letter'

# The font size ('10pt', '11pt' or '12pt').
#latex_font_size = '10pt'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
_stdauthor = 'Written by the SciPy community'
latex_documents = [
  ('index', 'scipy-ref.tex', 'SciPy Reference Guide', _stdauthor, 'manual'),
#  ('user/index', 'scipy-user.tex', 'SciPy User Guide',
#   _stdauthor, 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# Additional stuff for the LaTeX preamble.
latex_preamble = r'''
\usepackage{amsmath}

\DeclareUnicodeCharacter{00A0}{\nobreakspace}

% In the parameters etc. sections, align uniformly, and adjust label emphasis
\usepackage{expdlist}
\let\latexdescription=\description
\let\endlatexdescription=\enddescription
\renewenvironment{description}%
{\begin{latexdescription}[\setleftmargin{60pt}\breaklabel\setlabelstyle{\bfseries\itshape}]}%
{\end{latexdescription}}

% Make Examples/etc section headers smaller and more compact
\makeatletter
\titleformat{\paragraph}{\normalsize\normalfont\bfseries\itshape}%
            {\py@NormalColor}{0em}{\py@NormalColor}{\py@NormalColor}
\titlespacing*{\paragraph}{0pt}{1ex}{0pt}
\makeatother

% Save vertical space in parameter lists and elsewhere
\makeatletter
\renewenvironment{quote}%
               {\list{}{\topsep=0pt%
                        \parsep \z@ \@plus\p@}%
                \item\relax}%
               {\endlist}
\makeatother

% Fix footer/header
\renewcommand{\chaptermark}[1]{\markboth{\MakeUppercase{\thechapter.\ #1}}{}}
\renewcommand{\sectionmark}[1]{\markright{\MakeUppercase{\thesection.\ #1}}}
'''

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
latex_use_modindex = False


# -----------------------------------------------------------------------------
# Intersphinx configuration
# -----------------------------------------------------------------------------
intersphinx_mapping = {
        'http://docs.python.org/dev': None,
        'http://docs.scipy.org/doc/numpy': None,
}


# -----------------------------------------------------------------------------
# Numpy extensions
# -----------------------------------------------------------------------------

# If we want to do a phantom import from an XML file for all autodocs
phantom_import_file = 'dump.xml'

# Generate plots for example sections
numpydoc_use_plots = True

# -----------------------------------------------------------------------------
# Autosummary
# -----------------------------------------------------------------------------

if sphinx.__version__ >= "0.7":
    import glob
    autosummary_generate = glob.glob("*.rst")

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
# Plot
#------------------------------------------------------------------------------
plot_pre_code = """
import numpy as np
np.random.seed(123)
"""
plot_include_source = True
plot_formats = [('png', 96), 'pdf']
plot_html_show_formats = False

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

if not use_matplotlib_plot_directive:
    import matplotlib
    matplotlib.rcParams.update(plot_rcparams)

# -----------------------------------------------------------------------------
# Source code links
# -----------------------------------------------------------------------------

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
    print "NOTE: linkcode extension not found -- no links to source generated"

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
        except:
            return None

    try:
        fn = inspect.getsourcefile(obj)
    except:
        fn = None
    if not fn:
        try:
            fn = inspect.getsourcefile(sys.modules[obj.__module__])
        except:
            fn = None
    if not fn:
        return None

    try:
        source, lineno = inspect.findsource(obj)
    except:
        lineno = None

    if lineno:
        linespec = "#L%d" % (lineno + 1)
    else:
        linespec = ""

    fn = relpath(fn, start=dirname(scipy.__file__))

    if 'dev' in scipy.__version__:
        return "http://github.com/scipy/scipy/blob/master/scipy/%s%s" % (
           fn, linespec)
    else:
        return "http://github.com/scipy/scipy/blob/v%s/scipy/%s%s" % (
           scipy.__version__, fn, linespec)
