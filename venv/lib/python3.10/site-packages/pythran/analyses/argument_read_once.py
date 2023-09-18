""" ArgumentReadOnce counts the usages of each argument of each function. """

from pythran.analyses.aliases import Aliases
from pythran.analyses.global_declarations import GlobalDeclarations
from pythran.passmanager import ModuleAnalysis
from pythran.tables import MODULES
import pythran.intrinsic as intrinsic

import gast as ast
from functools import reduce


class ArgumentReadOnce(ModuleAnalysis):

    """
    Counts the usages of each argument of each function.

    Attributes
    ----------
    result : {FunctionEffects}
        Number of use for each argument of each function.
    node_to_functioneffect : {???: ???}
        FunctionDef ast node to function effect binding.
    """

    class FunctionEffects(object):
        def __init__(self, node):
            self.func = node
            self.dependencies = lambda ctx: 0
            if isinstance(node, ast.FunctionDef):
                self.read_effects = [-1] * len(node.args.args)
            elif isinstance(node, intrinsic.Intrinsic):
                self.read_effects = [
                    1 if isinstance(x, intrinsic.ReadOnceEffect)
                    else 2 for x in node.argument_effects]
            elif isinstance(node, ast.alias):
                self.read_effects = []
            else:
                raise NotImplementedError

    class ConstructorEffects(object):
        def __init__(self, node):
            self.func = node
            self.dependencies = lambda ctx: 0
            self.read_effects = [0]

    class Context(object):
        def __init__(self, function, index, path, global_dependencies):
            self.function = function
            self.index = index
            self.path = path
            self.global_dependencies = global_dependencies

    def __init__(self):
        """ Basic initialiser for class attributes. """
        self.result = set()
        self.node_to_functioneffect = dict()
        super(ArgumentReadOnce, self).__init__(Aliases, GlobalDeclarations)

    def prepare(self, node):
        """
        Initialise arguments effects as this analysis in inter-procedural.

        Initialisation done for Pythonic functions and default values set for
        user defined functions.
        """
        super(ArgumentReadOnce, self).prepare(node)
        # global functions init
        for n in self.global_declarations.values():
            fe = ArgumentReadOnce.FunctionEffects(n)
            self.node_to_functioneffect[n] = fe
            self.result.add(fe)

        # Pythonic functions init
        def save_effect(module):
            """ Recursively save read once effect for Pythonic functions. """
            for intr in module.values():
                if isinstance(intr, dict):  # Submodule case
                    save_effect(intr)
                else:
                    fe = ArgumentReadOnce.FunctionEffects(intr)
                    self.node_to_functioneffect[intr] = fe
                    self.result.add(fe)
                    if isinstance(intr, intrinsic.Class):  # Class case
                        save_effect(intr.fields)

        for module in MODULES.values():
            save_effect(module)

    def run(self, node):
        result = super(ArgumentReadOnce, self).run(node)
        for fun in result:
            for i in range(len(fun.read_effects)):
                self.recursive_weight(fun, i, set())
        self.result = {f.func: f.read_effects for f in result}
        return self.result

    def recursive_weight(self, function, index, predecessors):
        # TODO : Find out why it happens in some cases
        if len(function.read_effects) <= index:
            return 0
        if function.read_effects[index] == -1:
            # In case of recursive/cyclic calls
            cycle = function in predecessors
            predecessors.add(function)
            if cycle:
                function.read_effects[index] = 2 * function.dependencies(
                    ArgumentReadOnce.Context(function, index,
                                             predecessors, False))
            else:
                function.read_effects[index] = function.dependencies(
                    ArgumentReadOnce.Context(function, index,
                                             predecessors, True))
        return function.read_effects[index]

    def argument_index(self, node):
        while isinstance(node, ast.Subscript):
            node = node.value
        if node in self.aliases:
            for n_alias in self.aliases[node]:
                try:
                    return self.current_function.func.args.args.index(n_alias)
                except ValueError:
                    pass
        return -1

    def local_effect(self, node, effect):
        index = self.argument_index(node)
        return lambda ctx: effect if index == ctx.index else 0

    def generic_visit(self, node):
        lambdas = [self.visit(child) for child in ast.iter_child_nodes(node)]
        return lambda ctx: sum(l(ctx) for l in lambdas)

    def visit_FunctionDef(self, node):
        self.current_function = self.node_to_functioneffect[node]
        assert self.current_function in self.result
        self.current_function.dependencies = self.generic_visit(node)

    def visit_Return(self, node):
        dep = self.generic_visit(node)
        if isinstance(node.value, ast.Name):
            local = self.local_effect(node.value, 2)
            return lambda ctx: dep(ctx) + local(ctx)
        else:
            return dep

    def visit_Assign(self, node):
        dep = self.generic_visit(node)
        local = [self.local_effect(t, 2) for t in node.targets
                 if isinstance(t, ast.Subscript)]
        return lambda ctx: dep(ctx) + sum(l(ctx) for l in local)

    def visit_AugAssign(self, node):
        dep = self.generic_visit(node)
        local = self.local_effect(node.target, 2)
        return lambda ctx: dep(ctx) + local(ctx)

    def visit_For(self, node):
        iter_local = self.local_effect(node.iter, 1)
        iter_deps = self.visit(node.iter)
        body_deps = [self.visit(stmt) for stmt in node.body]
        else_deps = [self.visit(stmt) for stmt in node.orelse]
        return lambda ctx: iter_local(ctx) + iter_deps(ctx) + 2 * sum(
            l(ctx) for l in body_deps) + sum(l(ctx) for l in else_deps)

    def visit_While(self, node):
        test_deps = self.visit(node.test)
        body_deps = [self.visit(stmt) for stmt in node.body]
        else_deps = [self.visit(stmt) for stmt in node.orelse]
        return lambda ctx: test_deps(ctx) + 2 * sum(
            l(ctx) for l in body_deps) + sum(l(ctx) for l in else_deps)

    def visit_If(self, node):
        test_deps = self.visit(node.test)
        body_deps = [self.visit(stmt) for stmt in node.body]
        else_deps = [self.visit(stmt) for stmt in node.orelse]
        return lambda ctx: test_deps(ctx) + max(sum(
            l(ctx) for l in body_deps), sum(l(ctx) for l in else_deps))

    def visit_Call(self, node):
        l0 = self.generic_visit(node)
        index_corres = dict()
        func = None
        for i, arg in enumerate(node.args):
            n = self.argument_index(arg)
            if n >= 0:
                func_aliases = self.aliases[node.func]

                # expand argument if any
                func_aliases = reduce(
                    lambda x, y: x + (
                        # all functions
                        list(self.node_to_functioneffect.keys())
                        if (isinstance(y, ast.Name) and
                            self.argument_index(y) >= 0)
                        else [y]),
                    func_aliases,
                    list())

                for func_alias in func_aliases:
                    # special hook for binded functions
                    if isinstance(func_alias, ast.Call):
                        bound_name = func_alias.args[0].id
                        func_alias = self.global_declarations[bound_name]

                    if func_alias is intrinsic.UnboundValue:
                        continue

                    if func_alias not in self.node_to_functioneffect:
                        continue

                    func = self.node_to_functioneffect[func_alias]
                    index_corres[n] = i

        def merger(ctx):
            base = l0(ctx)
            if (ctx.index in index_corres) and ctx.global_dependencies:
                rec = self.recursive_weight(func, index_corres[ctx.index],
                                            ctx.path)
            else:
                rec = 0
            return base + rec

        return merger

    def visit_Subscript(self, node):
        dep = self.generic_visit(node)
        local = self.local_effect(node.value, 2)
        return lambda ctx: dep(ctx) + local(ctx)

    def visit_comprehension(self, node):
        dep = self.generic_visit(node)
        local = self.local_effect(node.iter, 1)
        return lambda ctx: dep(ctx) + local(ctx)
