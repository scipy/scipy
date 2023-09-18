""" Globals computes the value of globals(). """

from pythran.analyses.global_declarations import GlobalDeclarations
from pythran.passmanager import ModuleAnalysis


class Globals(ModuleAnalysis):
    def __init__(self):
        self.result = set()
        super(Globals, self).__init__(GlobalDeclarations)

    def visit_Module(self, node):
        self.result = {'builtins',
                       '__dispatch__'}.union(self.global_declarations.keys())
