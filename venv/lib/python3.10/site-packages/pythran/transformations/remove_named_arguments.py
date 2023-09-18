""" RemoveNamedArguments turns named arguments into regular ones.  """

from pythran.analyses import Aliases
from pythran.passmanager import Transformation
from pythran.syntax import PythranSyntaxError
from pythran.tables import MODULES

import gast as ast
from copy import deepcopy


def handle_special_calls(func_alias, node):
    if func_alias is MODULES['numpy']['arange']:
        if len(node.args) == 1:
            node.args.insert(0, ast.Constant(0, None))


class RemoveNamedArguments(Transformation):
    '''
    Replace call with named arguments to regular calls

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> code = 'def foo(x, y): return x + y\\ndef bar(z): return foo(y=z, x=0)'
    >>> node = ast.parse(code)
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(RemoveNamedArguments, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(x, y):
        return (x + y)
    def bar(z):
        return foo(0, z)
    '''

    def __init__(self):
        super(RemoveNamedArguments, self).__init__(Aliases)

    def handle_keywords(self, func, node, offset=0):
        '''
        Gather keywords to positional argument information

        Assumes the named parameter exist, raises a KeyError otherwise
        '''
        func_argument_names = {}
        for i, arg in enumerate(func.args.args[offset:]):
            assert isinstance(arg, ast.Name)
            func_argument_names[arg.id] = i

        func_argument_kwonly_names = {}
        for i, arg in enumerate(func.args.kwonlyargs):
            assert isinstance(arg, ast.Name)
            func_argument_kwonly_names[arg.id] = i

        nargs = len(func.args.args) - offset
        defaults = func.args.defaults
        keywords = {func_argument_names[kw.arg]: kw.value
                    for kw in node.keywords
                    if kw.arg not in func_argument_kwonly_names}

        keywords_only = []
        nb_kw = len(node.keywords)
        for i, kw in enumerate(list(reversed(node.keywords))):
            if kw.arg in func_argument_kwonly_names:
                keywords_only.append((func_argument_kwonly_names[kw.arg],
                                      kw.value))
                node.keywords.pop(nb_kw - i - 1)

        keywords_only = [v for _, v in sorted(keywords_only)]

        extra_keyword_offset = max(keywords.keys()) if keywords else 0
        node.args.extend([None] * (1 + extra_keyword_offset - len(node.args)))

        replacements = {}
        for index, arg in enumerate(node.args):
            if arg is None:
                if index in keywords:
                    replacements[index] = deepcopy(keywords[index])
                else:  # must be a default value
                    replacements[index] = deepcopy(defaults[index - nargs])

        if not keywords_only:
            return replacements

        node.args.append(ast.Call(
            ast.Attribute(
                ast.Attribute(
                    ast.Name("builtins", ast.Load(), None, None),
                    "pythran",
                    ast.Load()),
                "kwonly",
                ast.Load()),
            [], [])
        )
        node.args.extend(keywords_only)

        return replacements

    def visit_Call(self, node):
        if node.keywords:
            self.update = True

            aliases = self.aliases[node.func]
            assert aliases, "at least one alias"

            # all aliases should have the same structural type...
            # call to self.handle_keywords raises an exception otherwise
            try:
                replacements = {}
                for func_alias in aliases:
                    handle_special_calls(func_alias, node)

                    if func_alias is None:  # aliasing computation failed
                        pass
                    elif isinstance(func_alias, ast.Call):  # nested function
                        # func_alias looks like functools.partial(foo, a)
                        # so we reorder using alias for 'foo'
                        offset = len(func_alias.args) - 1
                        call = func_alias.args[0]
                        for func_alias in self.aliases[call]:
                            replacements = self.handle_keywords(func_alias,
                                                                node, offset)
                    else:
                        replacements = self.handle_keywords(func_alias, node)

                # if we reach this point, we should have a replacement
                # candidate, or nothing structural typing issues would have
                # raised an exception in handle_keywords
                if replacements:
                    for index, value in replacements.items():
                        node.args[index] = value
                    node.keywords = []

            except KeyError as ve:
                err = ("function uses an unknown (or unsupported) keyword "
                       "argument `{}`".format(ve.args[0]))
                raise PythranSyntaxError(err, node)
        return self.generic_visit(node)
