""" RemoveNestedFunctions turns nested function into top-level functions. """

from pythran.analyses import GlobalDeclarations, ImportedIds
from pythran.passmanager import Transformation
from pythran.tables import MODULES
from pythran.conversion import mangle

import pythran.metadata as metadata

import gast as ast


class _NestedFunctionRemover(ast.NodeTransformer):
    def __init__(self, parent):
        ast.NodeTransformer.__init__(self)
        self.parent = parent
        self.identifiers = set(self.global_declarations.keys())

    def __getattr__(self, attr):
        return getattr(self.parent, attr)

    def visit_FunctionDef(self, node):
        self.update = True
        if MODULES['functools'] not in self.global_declarations.values():
            import_ = ast.Import([ast.alias('functools', mangle('functools'))])
            self.ctx.module.body.insert(0, import_)
            functools_module = MODULES['functools']
            self.global_declarations[mangle('functools')] = functools_module

        self.ctx.module.body.append(node)

        former_name = node.name
        seed = 0
        new_name = "pythran_{}{}"

        while new_name.format(former_name, seed) in self.identifiers:
            seed += 1

        new_name = new_name.format(former_name, seed)
        self.identifiers.add(new_name)

        ii = self.gather(ImportedIds, node)
        binded_args = [ast.Name(iin, ast.Load(), None, None)
                       for iin in sorted(ii)]
        node.args.args = ([ast.Name(iin, ast.Param(), None, None)
                           for iin in sorted(ii)] +
                          node.args.args)

        metadata.add(node, metadata.Local())

        class Renamer(ast.NodeTransformer):
            def visit_Call(self, node):
                self.generic_visit(node)
                if (isinstance(node.func, ast.Name) and
                        node.func.id == former_name):
                    node.func.id = new_name
                    node.args = (
                        [ast.Name(iin, ast.Load(), None, None)
                         for iin in sorted(ii)] +
                        node.args
                        )
                return node
        Renamer().visit(node)

        node.name = new_name
        self.global_declarations[node.name] = node
        proxy_call = ast.Name(new_name, ast.Load(), None, None)

        new_node = ast.Assign(
            [ast.Name(former_name, ast.Store(), None, None)],
            ast.Call(
                ast.Attribute(
                    ast.Name(mangle('functools'), ast.Load(), None, None),
                    "partial",
                    ast.Load()
                    ),
                [proxy_call] + binded_args,
                [],
                ),
            None)

        self.generic_visit(node)
        return new_node


class RemoveNestedFunctions(Transformation):

    """
    Replace nested function by top-level functions.

    Also add a call to a bind intrinsic that
    generates a local function with some arguments binded.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("def foo(x):\\n def bar(y): return x+y\\n bar(12)")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(RemoveNestedFunctions, node)
    >>> print(pm.dump(backend.Python, node))
    import functools as __pythran_import_functools
    def foo(x):
        bar = __pythran_import_functools.partial(pythran_bar0, x)
        bar(12)
    def pythran_bar0(x, y):
        return (x + y)
    """
    def __init__(self):
        super(RemoveNestedFunctions, self).__init__(GlobalDeclarations)

    def visit_Module(self, node):
        # keep original node as it's updated by _NestedFunctionRemover
        for stmt in node.body:
            self.visit(stmt)
        return node

    def visit_FunctionDef(self, node):
        nfr = _NestedFunctionRemover(self)
        node.body = [nfr.visit(stmt) for stmt in node.body]
        self.update |= nfr.update
        return node
