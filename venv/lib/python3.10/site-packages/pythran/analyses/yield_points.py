"""
YieldPoints gathers all yield points from a node
"""

from pythran.passmanager import FunctionAnalysis


class YieldPoints(FunctionAnalysis):
    '''Gathers all yield points of a generator, if any.'''
    def __init__(self):
        self.result = list()
        super(YieldPoints, self).__init__()

    def visit_Yield(self, node):
        self.result.append(node)
