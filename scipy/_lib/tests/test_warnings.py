"""
Tests which scan for certain occurrences in the code, they may not find
all of these occurrences but should catch almost all. This file was adapted
from numpy.
"""


from __future__ import division, absolute_import, print_function


import sys
if sys.version_info >= (3, 4):
    from pathlib import Path
    import ast
    import tokenize
    import scipy
    from numpy.testing import run_module_suite
    from numpy.testing.decorators import slow


    class ParseCall(ast.NodeVisitor):
        def __init__(self):
            self.ls = []

        def visit_Attribute(self, node):
            ast.NodeVisitor.generic_visit(self, node)
            self.ls.append(node.attr)

        def visit_Name(self, node):
            self.ls.append(node.id)


    class FindFuncs(ast.NodeVisitor):
        def __init__(self, filename):
            super().__init__()
            self.__filename = filename

        def visit_Call(self, node):
            p = ParseCall()
            p.visit(node.func)
            ast.NodeVisitor.generic_visit(self, node)

            if p.ls[-1] == 'simplefilter' or p.ls[-1] == 'filterwarnings':
                if node.args[0].s == "ignore":
                    raise AssertionError(
                        "ignore filter should not be used; found in "
                        "{} on line {}".format(self.__filename, node.lineno))

            if p.ls[-1] == 'warn' and (
                    len(p.ls) == 1 or p.ls[-2] == 'warnings'):

                if "_lib/tests/test_warnings.py" is self.__filename:
                    # This file
                    return

                # See if stacklevel exists:
                # if len(node.args) == 3:
                #     return
                # args = {kw.arg for kw in node.keywords}
                # if "stacklevel" in args:
                #     return
                # raise AssertionError(
                #     "warnings should have an appropriate stacklevel; found in "
                #     "{} on line {}".format(self.__filename, node.lineno))


    @slow
    def test_warning_calls():
        # combined "ignore" and stacklevel error
        base = Path(scipy.__file__).parent

        for path in base.rglob("*.py"):
            # There is still one missing occurance in optimize.py,
            # this is one that should be fixed and this removed then.
            if path == base / "optimize" / "optimize.py":
                continue
            # use tokenize to auto-detect encoding on systems where no
            # default encoding is defined (e.g. LANG='C')
            with tokenize.open(str(path)) as file:
                tree = ast.parse(file.read())
                FindFuncs(path).visit(tree)


    if __name__ == "__main__":
        run_module_suite()
