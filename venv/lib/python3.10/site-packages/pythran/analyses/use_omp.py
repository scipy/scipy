"""
UseOMP detects if a function use OpenMP
"""

from pythran.passmanager import FunctionAnalysis


class UseOMP(FunctionAnalysis):
    """Detects if a function use openMP"""
    def __init__(self):
        self.result = False
        super(UseOMP, self).__init__()

    def visit_OMPDirective(self, _):
        self.result = True
