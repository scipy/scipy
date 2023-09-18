""" UsedDefChain build used-define chains analysis for each variable. """

from pythran.passmanager import ModuleAnalysis
import pythran.metadata as md

import beniget


class ExtendedDefUseChains(beniget.DefUseChains):

    def unbound_identifier(self, name, node):
        # don't warn on unbound identifier
        pass

    def visit(self, node):
        # be aware of metadata
        md.visit(self, node)
        return super(ExtendedDefUseChains, self).visit(node)


class UseDefChains(ModuleAnalysis):

    """
    Build use-define chains analysis for each variable.
    """

    def __init__(self):
        self.result = None
        super(UseDefChains, self).__init__(DefUseChains)

    def visit_Module(self, node):
        udc = beniget.UseDefChains(self.def_use_chains)
        self.result = udc.chains


class DefUseChains(ModuleAnalysis):

    """
    Build define-use-define chains analysis for each variable.
    """

    def __init__(self):
        self.result = None
        super(DefUseChains, self).__init__()

    def visit_Module(self, node):
        duc = ExtendedDefUseChains()
        duc.visit(node)
        self.result = duc
