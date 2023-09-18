'''
This module performs the return type inference, according to symbolic types,
   It then reorders function declarations according to the return type deps.
    * type_all generates a node -> type binding
'''

from pythran.analyses import LazynessAnalysis, StrictAliases, YieldPoints
from pythran.analyses import LocalNodeDeclarations, Immediates, RangeValues
from pythran.analyses import Ancestors
from pythran.config import cfg
from pythran.cxxtypes import TypeBuilder, ordered_set
from pythran.intrinsic import UserFunction, Class
from pythran.passmanager import ModuleAnalysis
from pythran.tables import operator_to_lambda, MODULES
from pythran.types.conversion import pytype_to_ctype
from pythran.types.reorder import Reorder
from pythran.utils import attr_to_path, cxxid, isnum, isextslice

from collections import defaultdict
from functools import reduce
import gast as ast
from itertools import islice
from copy import deepcopy


class UnboundableRValue(Exception):
    pass


class Types(ModuleAnalysis):

    """ Infer symbolic type for all AST node. """

    def __init__(self):

        self.max_seq_size = cfg.getint('typing',
                                       'max_heterogeneous_sequence_size')

        class TypeResult(dict):
            def __init__(self):
                self.builder = TypeBuilder()

            def copy(self):
                other = TypeResult()
                other.update(self.items())
                other.builder = self.builder
                return other

        self.result = TypeResult()
        self.builder = self.result.builder
        self.result["bool"] = self.builder.NamedType("bool")
        self.combiners = defaultdict(UserFunction)
        self.current_global_declarations = dict()
        self.max_recompute = 1  # max number of use to be lazy
        ModuleAnalysis.__init__(self, Reorder, StrictAliases, LazynessAnalysis,
                                Immediates, RangeValues, Ancestors)
        self.curr_locals_declaration = None

    def combined(self, *types):
        all_types = ordered_set()
        for ty in types:
            if isinstance(ty, self.builder.CombinedTypes):
                all_types.extend(ty.types)
            elif ty is not self.builder.UnknownType:
                all_types.append(ty)

        if len(all_types) == 0:
            return self.builder.UnknownType
        elif len(all_types) == 1:
            return next(iter(all_types))
        else:
            return self.builder.CombinedTypes(*all_types)



    def prepare(self, node):
        """
        Initialise values to prepare typing computation.

        Reorder functions to avoid dependencies issues and prepare typing
        computation setting typing values for Pythonic functions.
        """

        def register(name, module):
            """ Recursively save function typing and combiners for Pythonic."""
            for fname, function in module.items():
                if isinstance(function, dict):
                    register(name + "::" + fname, function)
                else:
                    tname = 'pythonic::{0}::functor::{1}'.format(name, fname)
                    self.result[function] = self.builder.NamedType(tname)
                    self.combiners[function] = function
                    if isinstance(function, Class):
                        register(name + "::" + fname, function.fields)

        for mname, module in MODULES.items():
            register(mname, module)
        super(Types, self).prepare(node)

    def run(self, node):
        super(Types, self).run(node)
        for head in self.current_global_declarations.values():
            if head not in self.result:
                self.result[head] = "pythonic::types::none_type"
        return self.result

    def register(self, fname, nid, ptype):
        """register ptype as a local typedef"""
        # Too many of them leads to memory burst
        if len(self.typedefs[fname, nid]) < cfg.getint('typing', 'max_combiner'):
            self.typedefs[fname, nid].append(ptype)
            return True
        return False

    def node_to_id(self, n, depth=()):
        if isinstance(n, ast.Name):
            return (n.id, depth)
        elif isinstance(n, ast.Subscript):
            if isinstance(n.slice, ast.Slice):
                return self.node_to_id(n.value, depth)
            else:
                index = n.slice.value if isnum(n.slice) else None
                return self.node_to_id(n.value, depth + (index,))
        # use alias information if any
        elif isinstance(n, ast.Call):
            for alias in self.strict_aliases[n]:
                if alias is n:  # no specific alias info
                    continue
                try:
                    return self.node_to_id(alias, depth)
                except UnboundableRValue:
                    continue
        raise UnboundableRValue()

    def isargument(self, node):
        """ checks whether node aliases to a parameter."""
        try:
            node_id, _ = self.node_to_id(node)
            return (node_id in self.name_to_nodes and
                    any(isinstance(n, ast.Name) and
                        isinstance(n.ctx, ast.Param)
                        for n in self.name_to_nodes[node_id]))
        except UnboundableRValue:
            return False

    def combine(self, node, op, othernode):
        """
        Change `node` typing with combination of `node` and `othernode`.
        """
        if self.result[othernode] is self.builder.UnknownType:
            if node not in self.result:
                self.result[node] = self.builder.UnknownType
            return

        if op is None:
            op = lambda x: x

        node_aliases = ordered_set([node])
        node_aliases.extend(self.strict_aliases.get(node, ()))
        for a in node_aliases:
            self.combine_(a, op, othernode)

    def combine_(self, node, op, othernode):
        try:
            # This comes from an assignment,so we must check where the value is
            # assigned
            try:
                node_id, depth = self.node_to_id(node)
                if depth:
                    node = ast.Name(node_id, ast.Load(), None, None)
                    former_op = op

                    # update the type to reflect container nesting
                    def merge_container_type(ty, index):
                        # integral index make it possible to correctly
                        # update tuple type
                        if isinstance(index, int):
                            kty = self.builder.NamedType(
                                    'std::integral_constant<long,{}>'
                                    .format(index))
                            return self.builder.IndexableContainerType(kty,
                                                                       ty)
                        elif isinstance(index, float):
                            kty = self.builder.NamedType('double')
                            return self.builder.IndexableContainerType(kty,
                                                                       ty)
                        else:
                            # FIXME: what about other key types?
                            return self.builder.ContainerType(ty)

                    def op(*args):
                        return reduce(merge_container_type, depth,
                                      former_op(*args))

                self.name_to_nodes[node_id].append(node)
            except UnboundableRValue:
                pass

            # perform inter procedural combination
            if self.isargument(node):
                node_id, _ = self.node_to_id(node)
                if node not in self.result:
                    self.result[node] = op(self.result[othernode])
                assert self.result[node], "found an alias with a type"

                parametric_type = self.builder.PType(self.current,
                                                     self.result[othernode])

                if self.register(self.current, node_id, parametric_type):

                    current_function = self.combiners[self.current]

                    def translator_generator(args, op):
                        ''' capture args for translator generation'''
                        def interprocedural_type_translator(s, n):
                            translated_othernode = ast.Name(
                                '__fake__', ast.Load(), None, None)
                            s.result[translated_othernode] = (
                                parametric_type.instanciate(
                                    s.current,
                                    [s.result[arg] for arg in n.args]))

                            # look for modified argument
                            for p, effective_arg in enumerate(n.args):
                                formal_arg = args[p]
                                if formal_arg.id == node_id:
                                    translated_node = effective_arg
                                    break
                            try:
                                s.combine(translated_node,
                                          op,
                                          translated_othernode)
                            except NotImplementedError:
                                pass
                                # this may fail when the effective
                                # parameter is an expression
                            except UnboundLocalError:
                                pass
                                # this may fail when translated_node
                                # is a default parameter
                        return interprocedural_type_translator

                    translator = translator_generator(self.current.args.args, op)
                    current_function.add_combiner(translator)
            else:
                self.update_type(node, op, self.result[othernode])

        except UnboundableRValue:
            pass

    def update_type(self, node, ty_builder, *args):
        if ty_builder is None:
            ty, = args
        # propagate UnknowType status if one of *args is itself unknown.
        elif any(arg is self.builder.UnknownType for arg in args):
            ty = self.builder.UnknownType
        else:
            ty = ty_builder(*args)

        curr_ty = self.result.get(node, self.builder.UnknownType)
        if isinstance(curr_ty, tuple):
            return

        self.result[node] = self.combined(curr_ty, ty)

    def visit_FunctionDef(self, node):
        self.delayed_types = set()
        self.curr_locals_declaration = self.gather(
            LocalNodeDeclarations,
            node)
        self.current = node
        self.typedefs = defaultdict(list)
        self.name_to_nodes = defaultdict(ordered_set)
        for arg in node.args.args:
            self.name_to_nodes[arg.id].append(arg)

        self.yield_points = self.gather(YieldPoints, node)

        self.generic_visit(node)

        for delayed_node in self.delayed_types:
            delayed_type = self.result[delayed_node]
            all_types = ordered_set(self.result[ty] for ty in
                                    self.name_to_nodes[delayed_node.id])
            final_type = self.combined(*all_types)
            delayed_type.final_type = final_type

        # propagate type information through all aliases
        for name, nodes in self.name_to_nodes.items():
            all_types = ordered_set(self.result[ty] for ty in nodes)
            final_type = self.combined(*all_types)
            for n in nodes:
                if n not in self.delayed_types:
                    self.result[n] = final_type
        self.current_global_declarations[node.name] = node

        # return type may be unset if the function always raises
        return_type = self.result.get(
            node,
            self.builder.NamedType("pythonic::types::none_type"))

        self.result[node] = (
            self.builder.Returnable(return_type),
            reduce(list.__add__, self.typedefs.values(), []))
        for k in self.curr_locals_declaration:
            self.result[k] = self.get_qualifier(k)(self.result[k])

    def get_qualifier(self, node):
        lazy_res = self.lazyness_analysis[node.id]
        return (self.builder.Lazy
                if lazy_res <= self.max_recompute
                else self.builder.Assignable)

    def visit_Return(self, node):
        """ Compute return type and merges with others possible return type."""
        self.generic_visit(node)
        # No merge are done if the function is a generator.
        if not self.yield_points:
            assert node.value, "Values were added in each return statement."
            self.update_type(self.current, None, self.result[node.value])

    def visit_Yield(self, node):
        """ Compute yield type and merges it with others yield type. """
        self.generic_visit(node)
        self.update_type(self.current, None, self.result[node.value])

    def visit_Assign(self, node):
        self.visit(node.value)
        for t in node.targets:
            self.combine(t, None, node.value)
            if t in self.curr_locals_declaration:
                self.result[t] = self.get_qualifier(t)(self.result[t])
            if isinstance(t, ast.Subscript):
                if self.visit_AssignedSubscript(t):
                    for alias in self.strict_aliases[t.value]:
                        fake = ast.Subscript(alias, t.slice, ast.Store())
                        self.combine(fake, None, node.value)

    def visit_AugAssign(self, node):
        self.visit(node.value)

        if isinstance(node.target, ast.Subscript):
            if self.visit_AssignedSubscript(node.target):
                for alias in self.strict_aliases[node.target.value]:
                    fake = ast.Subscript(alias, node.target.slice, ast.Store())
                    self.combine(fake, None, node.value)
        else:
            self.combine(node.target, None, node.value)


    def visit_For(self, node):
        self.visit(node.iter)
        self.combine(node.target, self.builder.IteratorContentType, node.iter)
        for n in node.body + node.orelse:
            self.visit(n)

    def visit_BoolOp(self, node):
        """
        Merge BoolOp operand type.

        BoolOp are "and" and "or" and may return any of these results so all
        operands should have the combinable type.
        """
        # Visit subnodes
        self.generic_visit(node)
        # Merge all operands types.
        for value in node.values:
            self.update_type(node, None, self.result[value])

    def visit_BinOp(self, node):
        self.generic_visit(node)
        self.update_type(node,
                         self.builder.ExpressionType,
                         operator_to_lambda[type(node.op)],
                         self.result[node.left],
                         self.result[node.right])

    def visit_UnaryOp(self, node):
        self.generic_visit(node)
        self.update_type(node,
                         self.builder.ExpressionType,
                         operator_to_lambda[type(node.op)],
                         self.result[node.operand])

    def visit_IfExp(self, node):
        self.generic_visit(node)
        self.update_type(node, None, self.result[node.body])
        self.update_type(node, None, self.result[node.orelse])

    def visit_Compare(self, node):
        self.generic_visit(node)
        all_compare = zip(node.ops,
                          [node.left] + node.comparators[:-1],
                          node.comparators)

        for op, left, right in all_compare:
            self.update_type(node,
                             self.builder.ExpressionType,
                             operator_to_lambda[type(op)],
                             self.result[left],
                             self.result[right])

    def visit_Call(self, node):
        self.generic_visit(node)

        func = node.func

        for alias in self.strict_aliases[func]:
            # this comes from a bind
            if isinstance(alias, ast.Call):
                a0 = alias.args[0]
                # by construction of the bind construct
                assert len(self.strict_aliases[a0]) == 1
                bounded_function = next(iter(self.strict_aliases[a0]))
                fake_name = deepcopy(a0)
                fake_node = ast.Call(fake_name, alias.args[1:] + node.args,
                                     [])
                self.combiners[bounded_function].combiner(self, fake_node)

            # handle backward type dependencies from function calls
            else:
                self.combiners[alias].combiner(self, node)

        UnknownType = self.builder.UnknownType

        # recurring nightmare
        def last_chance():
            # maybe we can get saved if we have a hint about
            # the called function return type
            for alias in self.strict_aliases[func]:
                if alias is self.current and alias in self.result:
                    # great we have a (partial) type information
                    self.result[node] = self.result[alias]
                    return
            self.result[node] = UnknownType

        if self.result[node.func] is UnknownType:
            return last_chance()

        if any(self.result[arg] is UnknownType for arg in node.args):
            return last_chance()

        # special handler for getattr: use the attr name as an enum member
        if (isinstance(func, ast.Attribute) and func.attr == 'getattr'):
            self.update_type(node,
                             self.builder.GetAttr,
                             self.result[node.args[0]],
                             node.args[1].value)
        # default behavior
        else:
            self.update_type(node,
                             self.builder.ReturnType,
                             self.result[func],
                             *[self.result[arg] for arg in node.args])

    def visit_Constant(self, node):
        """ Set the pythonic constant type. """
        ty = type(node.value)
        if ty is str and len(node.value) == 1:
            sty = 'pythonic::types::chr'
        else:
            sty = pytype_to_ctype(ty)
        if node in self.immediates:
            sty = "std::integral_constant<%s, %s>" % (sty,
                                                      str(node.value).lower())
        self.result[node] = self.builder.NamedType(sty)

    def visit_Attribute(self, node):
        """ Compute typing for an attribute node. """
        obj, path = attr_to_path(node)
        # If no type is given, use a decltype
        if obj.isliteral():
            typename = pytype_to_ctype(obj.signature)
            self.result[node] = self.builder.NamedType(typename)
        else:
            path = '::'.join(map(cxxid, path)) + '{}'
            self.result[node] = self.builder.DeclType(path)

    def visit_Slice(self, node):
        """
        Set slicing type using continuous information if provided.

        Also visit subnodes as they may contains relevant typing information.
        """
        self.generic_visit(node)
        nstep = node.step
        if nstep is None or (isnum(nstep) and nstep.value > 0):
            if nstep is None or nstep.value == 1:
                if all(self.range_values[p].low >= 0
                       for p in (node.lower, node.upper)):
                    ntype = "pythonic::types::fast_contiguous_slice"
                else:
                    ntype = "pythonic::types::contiguous_slice"
            else:
                ntype = "pythonic::types::cstride_slice<{}>".format(nstep.value)
            self.result[node] = self.builder.NamedType(ntype)
        else:
            self.result[node] = self.builder.NamedType(
                'pythonic::types::slice')

    def stores_to(self, node):
        ancestors = self.ancestors[node] + (node,)
        stmt_indices = [i for i, n in enumerate(ancestors)
                        if isinstance(n, (ast.Assign, ast.For))]
        if not stmt_indices:
            return True

        stmt_index = stmt_indices[-1]

        if isinstance(ancestors[stmt_index], ast.Assign):
            return ancestors[stmt_index + 1] is ancestors[stmt_index].value
        else:
            return ancestors[stmt_index + 1] is not ancestors[stmt_index].target

    def visit_Subscript(self, node):
        self.visit(node.value)
        if self.stores_to(node):
            self.result[node.value] = self.builder.AddConst(self.result[node.value])
        # type of a[1:2, 3, 4:1] is the type of: declval(a)(slice, long, slice)
        if isextslice(node.slice):
            self.visit(node.slice)
            self.update_type(node,
                             self.builder.ExpressionType,
                             lambda a, *b: "{0}({1})".format(a, ", ".join(b)),
                             self.result[node.value],
                             *[self.result[d] for d in node.slice.elts])
        elif isnum(node.slice) and node.slice.value >= 0:
            # type of a[2] is the type of an elements of a
            # this special case is to make type inference easier
            # for the back end compiler
            self.update_type(node,
                             self.builder.ElementType,
                             node.slice.value,
                             self.result[node.value])
        else:
            # type of a[i] is the return type of the matching function
            self.visit(node.slice)
            self.update_type(node,
                             self.builder.ExpressionType,
                             "{0}[{1}]".format,
                             self.result[node.value],
                             self.result[node.slice])

    def visit_AssignedSubscript(self, node):
        if isinstance(node.slice, ast.Slice):
            return False
        elif isextslice(node.slice):
            return False
        else:
            self.visit(node.slice)
            self.combine(node.value, self.builder.IndexableType, node.slice)
            return True

    def delayed(self, node):
        fallback_type = self.combined(*[self.result[n] for n in
                                        self.name_to_nodes[node.id]])
        self.delayed_types.add(node)
        return self.builder.LType(fallback_type, node)

    def visit_Name(self, node):
        if node.id in self.name_to_nodes:
            self.result[node] = self.delayed(node)
        elif node.id in self.current_global_declarations:
            newtype = self.builder.NamedType(
                self.current_global_declarations[node.id].name)
            if node not in self.result:
                self.result[node] = newtype
        else:
            self.result[node] = self.builder.UnknownType

    def visit_List(self, node):
        """ Define list type from all elements type (or empty_list type). """
        self.generic_visit(node)
        if node.elts:
            for elt in node.elts[:self.max_seq_size]:
                self.update_type(node, self.builder.ListType, self.result[elt])
        else:
            self.update_type(node,
                             self.builder.NamedType,
                             "pythonic::types::empty_list")

    def visit_Set(self, node):
        """ Define set type from all elements type (or empty_set type). """
        self.generic_visit(node)
        if node.elts:
            for elt in node.elts[:self.max_seq_size]:
                self.update_type(node, self.builder.SetType, self.result[elt])
        else:
            self.update_type(node, self.builder.NamedType,
                "pythonic::types::empty_set")

    def visit_Dict(self, node):
        """ Define set type from all elements type (or empty_dict type). """
        self.generic_visit(node)
        if node.keys:
            for key, value in islice(zip(node.keys, node.values),
                                     self.max_seq_size):
                self.update_type(node,
                                 self.builder.DictType,
                                 self.result[key],
                                 self.result[value])
        else:
            self.update_type(node,
                             self.builder.NamedType,
                             "pythonic::types::empty_dict")

    def visit_ExceptHandler(self, node):
        if node.type and node.name:
            if not isinstance(node.type, ast.Tuple):
                tname = self.builder.NamedType(
                    'pythonic::types::{0}'.format(node.type.attr))
                self.result[node.type] = tname
                self.combine(node.name, None, node.type)
        for n in node.body:
            self.visit(n)

    def visit_Tuple(self, node):
        self.generic_visit(node)
        types = [self.result[elt] for elt in node.elts]
        self.update_type(node,
                         self.builder.TupleType,
                         *types)

    def visit_arguments(self, node):
        for i, arg in enumerate(node.args):
            self.update_type(arg, self.builder.ArgumentType, i)
        for n in node.defaults:
            self.visit(n)
