""" NormalizeReturn adds return statement where relevant. """

from pythran.analyses import CFG, YieldPoints
from pythran.passmanager import Transformation

import gast as ast


class NormalizeReturn(Transformation):
    '''
    Adds Return statement when they are implicit,
    and adds the None return value when not set

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("def foo(y): print(y)")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(NormalizeReturn, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(y):
        print(y)
        return builtins.None
    '''

    def __init__(self):
        super(NormalizeReturn, self).__init__(CFG)

    def visit_FunctionDef(self, node):
        self.yield_points = self.gather(YieldPoints, node)
        for stmt in node.body:
            self.visit(stmt)
        # Look for nodes that have no successors; the predecessors of
        # the special NIL node are those AST nodes that end control flow
        # without a return statement.
        for n in self.cfg.predecessors(CFG.NIL):
            if not isinstance(n, (ast.Return, ast.Raise)):
                self.update = True
                if self.yield_points:
                    node.body.append(ast.Return(None))
                else:
                    none = ast.Attribute(
                        ast.Name("builtins", ast.Load(), None, None),
                        'None',
                        ast.Load())
                    node.body.append(ast.Return(none))
                break

        return node

    def visit_Return(self, node):
        if not node.value and not self.yield_points:
            none = ast.Attribute(ast.Name("builtins", ast.Load(), None, None),
                                 'None', ast.Load())
            node.value = none
            self.update = True
        return node
