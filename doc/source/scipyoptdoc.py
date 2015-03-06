"""
===========
scipyoptdoc
===========

Proper docstrings for scipy.optimize.minimize et al.

Usage::

    .. scipy-optimize:function:: scipy.optimize.minimize
       :impl: scipy.optimize.optimize._minimize_nelder_mead
       :method: Nelder-Mead

Produces output similar to autodoc, except

- The docstring is obtained from the 'impl' function
- The call signature is mangled so that the default values for method keyword
  and options dict are substituted
- 'Parameters' section is replaced by 'Options' section
- See Also link to the actual function documentation is inserted

"""
from __future__ import division, absolute_import, print_function

import os, sys, re, pydoc
import sphinx
import inspect
import collections
import textwrap
import pydoc
import inspect

if sphinx.__version__ < '1.0.1':
    raise RuntimeError("Sphinx 1.0.1 or newer is required")

from numpydoc.numpydoc import mangle_docstrings
from docutils.parsers.rst import Directive
from docutils.statemachine import ViewList
from sphinx.domains.python import PythonDomain


if sys.version_info[0] >= 3:
    sixu = lambda s: s
else:
    sixu = lambda s: unicode(s, 'unicode_escape')


def setup(app):
    app.add_domain(ScipyOptimizeInterfaceDomain)


def _option_required_str(x):
    if not x:
        raise ValueError("value is required")
    return str(x)


def _import_object(name):
    parts = name.split('.')
    module_name = '.'.join(parts[:-1])
    __import__(module_name)
    obj = getattr(sys.modules[module_name], parts[-1])
    return obj


class ScipyOptimizeInterfaceDomain(PythonDomain):
    name = 'scipy-optimize'

    def __init__(self, *a, **kw):
        super(ScipyOptimizeInterfaceDomain, self).__init__(*a, **kw)
        self.directives = dict(self.directives)
        self.directives['function'] = wrap_mangling_directive(self.directives['function'])


BLURB = """
.. seealso:: For documentation for the rest of the parameters, see `%s`
"""

def wrap_mangling_directive(base_directive):
    class directive(base_directive):
        def run(self):
            env = self.state.document.settings.env

            # Interface function
            name = self.arguments[0].strip()
            obj = _import_object(name)
            args, varargs, keywords, defaults = inspect.getargspec(obj)

            # Implementation function
            impl_name = self.options['impl']
            impl_obj = _import_object(impl_name)
            impl_args, impl_varargs, impl_keywords, impl_defaults = inspect.getargspec(impl_obj)

            # Format signature taking implementation into account
            args = list(args)
            defaults = list(defaults)

            def set_default(arg, value):
                j = args.index(arg)
                defaults[len(defaults) - (len(args) - j)] = value

            def remove_arg(arg):
                if arg not in args:
                    return
                j = args.index(arg)
                if j < len(args) - len(defaults):
                    del args[j]
                else:
                    del defaults[len(defaults) - (len(args) - j)]
                    del args[j]

            options = []
            for j, opt_name in enumerate(impl_args):
                if opt_name in args:
                    continue
                if j >= len(impl_args) - len(impl_defaults):
                    options.append((opt_name, impl_defaults[len(impl_defaults) - (len(impl_args) - j)]))
                else:
                    options.append((opt_name, None))
            set_default('options', dict(options))
            set_default('method', self.options['method'].strip())

            for arg in list(args):
                if arg not in impl_args and arg not in ('fun', 'x0', 'args', 'tol',
                                                        'callback', 'method', 'options'):
                    remove_arg(arg)

            signature = inspect.formatargspec(args, varargs, keywords, defaults)

            # Produce output
            self.options['noindex'] = True
            self.arguments[0] = name + signature
            lines = textwrap.dedent(pydoc.getdoc(impl_obj)).splitlines()
            new_lines = []
            for line in lines:
                if line.strip() == 'Options':
                    new_lines.append("Other Parameters")
                elif line.strip() == "-"*len('Options'):
                    new_lines.append("-"*len("Other Parameters"))
                else:
                    new_lines.append(line)
            mangle_docstrings(env.app, 'function', name, None, None, new_lines)
            lines = new_lines
            new_lines = []
            for line in lines:
                if line.strip() == ':Other Parameters:':
                    new_lines.extend((BLURB % (name,)).splitlines())
                    new_lines.append(':Options:')
                else:
                    new_lines.append(line)
            self.content = ViewList(new_lines, self.content.parent)
            return base_directive.run(self)

        option_spec = dict(base_directive.option_spec)
        option_spec['impl'] = _option_required_str
        option_spec['method'] = _option_required_str

    return directive
