"""
Identifiers gathers all identifiers used in a node
"""

from pythran.passmanager import NodeAnalysis


class Identifiers(NodeAnalysis):
    """Gather all identifiers used throughout a node."""
    def __init__(self):
        self.result = set()
        super(Identifiers, self).__init__()

    def visit_Name(self, node):
        self.result.add(node.id)

    def visit_FunctionDef(self, node):
        self.result.add(node.name)
        self.generic_visit(node)

    def visit_ImportFrom(self, node):
        self.generic_visit(node)
        self.result.add(node.module)

    def visit_alias(self, node):
        self.result.add(node.name)
        if node.asname:
            self.result.add(node.asname)
