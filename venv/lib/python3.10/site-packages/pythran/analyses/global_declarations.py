""" GlobalDeclarations gathers top-level declarations. """

import gast as ast
from pythran.passmanager import ModuleAnalysis

class GlobalDeclarations(ModuleAnalysis):

    """ Gather all kind of identifier defined at global scope.

    >>> import gast as ast
    >>> from pythran import passmanager
    >>> from pythran.analyses import GlobalDeclarations
    >>> node = ast.parse('''
    ... import math
    ... import math as maths
    ... from math import cos
    ... c = 12
    ... def foo(a):
    ...     b = a + 1''')
    >>> pm = passmanager.PassManager("test")
    >>> sorted(pm.gather(GlobalDeclarations, node).keys())
    ['c', 'cos', 'foo', 'math', 'maths']

    """

    def __init__(self):
        """ Result is an identifier with matching definition. """
        self.result = dict()
        super(GlobalDeclarations, self).__init__()

    def visit_FunctionDef(self, node):
        """ Import module define a new variable name. """
        self.result[node.name] = node

    def visit_Import(self, node):
        for alias in node.names:
            self.result[alias.asname or alias.name] = alias

    def visit_ImportFrom(self, node):
        for alias in node.names:
            self.result[alias.asname or alias.name] = alias

    def visit_Name(self, node):
        if isinstance(node.ctx, ast.Store):
            self.result[node.id] = node
