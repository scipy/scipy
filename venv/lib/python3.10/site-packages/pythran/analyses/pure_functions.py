"""
PureFunctions lists functions without side-effects.
"""

from pythran.analyses.argument_effects import ArgumentEffects
from pythran.analyses.global_effects import GlobalEffects
from pythran.passmanager import ModuleAnalysis


class PureFunctions(ModuleAnalysis):
    '''Yields the set of pure functions'''

    def __init__(self):
        self.result = set()
        super(PureFunctions, self).__init__(ArgumentEffects, GlobalEffects)


    def prepare(self, node):
        super(PureFunctions, self).prepare(node)
        self.result = {func for func, ae in self.argument_effects.items()
                       if func not in self.global_effects and not any(ae)}
