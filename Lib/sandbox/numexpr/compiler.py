
import sys
import numpy

import interpreter, expressions

typecode_to_kind = {'b': 'bool', 'i': 'int', 'f': 'float', 'c': 'complex', 'n' : 'none'}
kind_to_typecode = {'bool': 'b', 'int': 'i', 'float': 'f', 'complex': 'c', 'none' : 'n'}
type_to_kind = expressions.type_to_kind
kind_to_type = expressions.kind_to_type

class ASTNode(object):
    """Abstract Syntax Tree node.

    Members:

    astType      -- type of node (op, constant, variable, raw, or alias)
    astKind      -- the type of the result (bool, float, etc.)
    value        -- value associated with this node.
                    An opcode, numerical value, a variable name, etc.
    children     -- the children below this node
    reg          -- the register assigned to the result for this node.
    """
    cmpnames = ['astType', 'astKind', 'value', 'children']
    def __init__(self, astType='generic', astKind='unknown',
                 value=None, children=()):
        object.__init__(self)
        self.astType = astType
        self.astKind = astKind
        self.value = value
        self.children = tuple(children)
        self.reg = None

    def __eq__(self, other):
        if self.astType == 'alias':
            self = self.value
        if other.astType == 'alias':
            other = other.value
        if not isinstance(other, ASTNode):
            return False
        for name in self.cmpnames:
            if getattr(self, name) != getattr(other, name):
                return False
        return True

    def __hash__(self):
        if self.astType == 'alias':
            self = self.value
        hashval = 0
        for name in self.cmpnames:
            hashval ^= hash(getattr(self, name))
        return hashval

    def __str__(self):
        return 'AST(%s, %s, %s, %s, %s)' % (self.astType, self.astKind,
                                            self.value, self.children, self.reg)
    def __repr__(self): return '<AST object at %s>' % id(self)

    def key(self):
        return (self.astType, self.astKind, self.value, self.children)

    def typecode(self):
        return kind_to_typecode[self.astKind]

    def postorderWalk(self):
        for c in self.children:
            for w in c.postorderWalk():
                yield w
        yield self

    def allOf(self, *astTypes):
        astTypes = set(astTypes)
        for w in self.postorderWalk():
            if w.astType in astTypes:
                yield w

def expressionToAST(ex):
    """Take an expression tree made out of expressions.ExpressionNode,
    and convert to an AST tree.

    This is necessary as ExpressionNode overrides many methods to act
    like a number.
    """
    this_ast = ASTNode(ex.astType, ex.astKind, ex.value,
                           [expressionToAST(c) for c in ex.children])
    return this_ast


def sigPerms(s):
    """Generate all possible signatures derived by upcasting the given
    signature.
    """
    codes = 'bifc'
    if not s:
        yield ''
    elif s[0] in codes:
        start = codes.index(s[0])
        for x in codes[start:]:
            for y in sigPerms(s[1:]):
                yield x + y
    else:
        yield s

def typeCompileAst(ast):
    """Assign appropiate types to each node in the AST.

    Will convert opcodes and functions to appropiate upcast version,
    and add "cast" ops if needed.
    """
    children = list(ast.children)
    if ast.astType  == 'op':
        retsig = ast.typecode()
        basesig = ''.join(x.typecode() for x in list(ast.children))
        # Find some operation that will work on an acceptable casting of args.
        for sig in sigPerms(basesig):
            value = ast.value + '_' + retsig + sig
            if value in interpreter.opcodes:
                break
        else:
            for sig in sigPerms(basesig):
                funcname = ast.value + '_' + retsig + sig
                if funcname in interpreter.funccodes:
                    value = 'func_%s' % (retsig+sig)
                    children += [ASTNode('raw', 'none',
                                         interpreter.funccodes[funcname])]
                    break
            else:
                raise NotImplementedError(
                    "couldn't find matching opcode for '%s'"
                    % (ast.value + '_' + retsig+basesig))
        # First just cast constants, then cast variables if necessary:
        for i, (have, want)  in enumerate(zip(basesig, sig)):
            if have != want:
                kind = typecode_to_kind[want]
                if children[i].astType == 'constant':
                    children[i] = ASTNode('constant', kind, children[i].value)
                else:
                    opname = "cast"
                    children[i] = ASTNode('op', kind, opname, [children[i]])
    else:
        value = ast.value
        children = ast.children
    new_ast = ASTNode(ast.astType, ast.astKind, value,
                           [typeCompileAst(c) for c in children])
    return new_ast


class Register(object):
    """Abstraction for a register in the VM.

    Members:
    node          -- the AST node this corresponds to
    temporary     -- True if this isn't an input or output
    immediate     -- not a register, but an immediate value
    n             -- the physical register number.
                     None if no number assigned yet.
    """
    def __init__(self, astnode, temporary=False):
        self.node = astnode
        self.temporary = temporary
        self.immediate = False
        self.n = None

    def __str__(self):
        if self.temporary:
            name = 'Temporary'
        else:
            name = 'Register'
        return '%s(%s, %s, %s)' % (name, self.node.astType,
                                   self.node.astKind, self.n,)

    def __repr__(self):
        return self.__str__()


class Immediate(Register):
    """Representation of an immediate (integer) operand, instead of
    a register.
    """
    def __init__(self, astnode):
        Register.__init__(self, astnode)
        self.immediate = True

    def __str__(self):
        return 'Immediate(%d)' % (self.node.value,)

def makeExpressions(context):
    """Make private copy of the expressions module with a custom get_context().

    An attempt was made to make this threadsafe, but I can't guarantee it's
    bulletproof.
    """
    import sys, imp
    modname = __name__[:__name__.rfind('.')] + '.expressions'
    # get our own, private copy of expressions
    imp.acquire_lock()
    try:
        old = sys.modules.pop(modname)
        import expressions
        private = sys.modules.pop(modname)
        sys.modules[modname] = old
    finally:
        imp.release_lock()
    def get_context():
        return context
    private.get_context = get_context
    return private

def stringToExpression(s, types, context):
    """Given a string, convert it to a tree of ExpressionNode's.
    """
    expr = makeExpressions(context)
    # first compile to a code object to determine the names
    c = compile(s, '<expr>', 'eval')
    # make VariableNode's for the names
    names = {}
    for name in c.co_names:
        if name == "None":
            names[name] = None
        else:
            t = types.get(name, float)
            names[name] = expr.VariableNode(name, type_to_kind[t])
    names.update(expr.functions)
    # now build the expression
    ex = eval(c, names)
    if expressions.isConstant(ex):
        ex = expr.ConstantNode(ex, expressions.getKind(ex))
    return ex


def isReduction(ast):
    return ast.value.startswith('sum_') or ast.value.startswith('prod_')

def getInputOrder(ast, input_order=None):
    """Derive the input order of the variables in an expression.
    """
    variables = {}
    for a in ast.allOf('variable'):
        variables[a.value] = a
    variable_names = set(variables.keys())

    if input_order:
        if variable_names != set(input_order):
            raise ValueError("input names don't match those found in expression")
        ordered_names = input_order
    else:
        ordered_names = list(variable_names)
        ordered_names.sort()
    ordered_variables = [variables[v] for v in ordered_names]
    return ordered_variables

def convertConstantToKind(x, kind):
    return kind_to_type[kind](x)

def getConstants(ast):
    const_map = {}
    for a in ast.allOf('constant'):
        const_map[(a.astKind, a.value)] = a
    ordered_constants = const_map.keys()
    ordered_constants.sort()
    constants_order = [const_map[v] for v in ordered_constants]
    constants = [convertConstantToKind(a.value, a.astKind)
                 for a in constants_order]
    return constants_order, constants

def sortNodesByOrder(nodes, order):
    order_map = {}
    for i, (_, v, _) in enumerate(order):
        order_map[v] = i
    dec_nodes = [(order_map[n.value], n) for n in nodes]
    dec_nodes.sort()
    return [a[1] for a in dec_nodes]

def assignLeafRegisters(inodes, registerMaker):
    """Assign new registers to each of the leaf nodes.
    """
    leafRegisters = {}
    for node in inodes:
        key = node.key()
        if key in leafRegisters:
            node.reg = leafRegisters[key]
        else:
            node.reg = leafRegisters[key] = registerMaker(node)

def assignBranchRegisters(inodes, registerMaker):
    """Assign temporary registers to each of the branch nodes.
    """
    for node in inodes:
        node.reg = registerMaker(node, temporary=True)


def collapseDuplicateSubtrees(ast):
    """common subexpression elimination.
    """
    seen = {}
    aliases = []
    for a in ast.allOf('op'):
        if a in seen:
            target = seen[a]
            a.astType = 'alias'
            a.value = target
            a.children = ()
            aliases.append(a)
        else:
            seen[a] = a
    # Set values and registers so optimizeTemporariesAllocation
    # doesn't get confused
    for a in aliases:
        while a.value.astType == 'alias':
            a.value = a.value.value
        a.reg = a.value.reg


def optimizeTemporariesAllocation(ast):
    """Attempt to minimize the number of temporaries needed, by
    reusing old ones.
    """
    nodes = list(x for x in ast.postorderWalk() if x.reg.temporary)
    users_of = dict((n.reg, set()) for n in nodes)
    if nodes and nodes[-1] is not ast:
        for c in ast.children:
            if c.reg.temporary:
                users_of[c.reg].add(ast)
    for n in reversed(nodes):
        for c in n.children:
            if c.reg.temporary:
                users_of[c.reg].add(n)
    unused = {'bool' : set(), 'int' : set(),
              'float' : set(), 'complex' : set()}
    for n in nodes:
        for reg, users in users_of.iteritems():
            if n in users:
                users.remove(n)
                if not users:
                    unused[reg.node.astKind].add(reg)
        if unused[n.astKind]:
            reg = unused[n.astKind].pop()
            users_of[reg] = users_of[n.reg]
            n.reg = reg


def setOrderedRegisterNumbers(order, start):
    """Given an order of nodes, assign register numbers.
    """
    for i, node in enumerate(order):
        node.reg.n = start + i
    return start + len(order)

def setRegisterNumbersForTemporaries(ast, start):
    """Assign register numbers for temporary registers, keeping track of
    aliases and handling immediate operands.
    """
    seen = 0
    signature = ''
    aliases = []
    for node in ast.postorderWalk():
        if node.astType == 'alias':
            aliases.append(node)
            node = node.value
        if node.reg.immediate:
            node.reg.n = node.value
            continue
        reg = node.reg
        if reg.n < 0:
            reg.n = start + seen
            seen += 1
            signature += reg.node.typecode()
    for node in aliases:
        node.reg = node.value.reg
    return start + seen, signature

def convertASTtoThreeAddrForm(ast):
    """Convert an AST to a three address form.

    Three address form is (op, reg1, reg2, reg3), where reg1 is the
    destination of the result of the instruction.

    I suppose this should be called three register form, but three
    address form is found in compiler theory.
    """
    program = []
    for node in ast.allOf('op'):
        children = node.children
        instr = (node.value, node.reg) \
                + tuple([c.reg for c in children])
        program.append(instr)
    return program

def compileThreeAddrForm(program):
    """Given a three address form of the program, compile it a string that
    the VM understands.
    """
    def nToChr(reg):
        if reg is None:
            return '\xff'
        elif reg.n < 0:
            raise ValueError("negative value for register number %s" % (reg.n,))
        else:
            return chr(reg.n)

    def quadrupleToString(opcode, store, a1=None, a2=None):
        cop = chr(interpreter.opcodes[opcode])
        cs = nToChr(store)
        ca1 = nToChr(a1)
        ca2 = nToChr(a2)
        return cop + cs + ca1 + ca2

    def toString(*args):
        while len(args) < 4:
            args += (None,)
        opcode, store, a1, a2 = args[0:4]
        s = quadrupleToString(opcode, store, a1, a2)
        l = [s]
        args = args[4:]
        while args:
            s = quadrupleToString('noop', *args[:3])
            l.append(s)
            args = args[3:]
        return ''.join(l)

    prog_str = ''.join([toString(*t) for t in program])
    return prog_str

context_info = [
    ('optimization', ('none', 'moderate', 'aggressive'), 'aggressive'),
               ]
def getContext(map):
    context = {}
    for name, allowed, default in context_info:
        value = map.pop(name, default)
        if value in allowed:
            context[name] = value
        else:
            raise ValueError("'%s' must be one of %s" % (name, allowed))
    if map:
        raise ValueError("Unknown keyword argument '%s'" % map.popitem()[0])
    return context


def precompile(ex, signature=(), copy_args=(), **kwargs):
    """Compile the expression to an intermediate form.
    """
    types = dict(signature)
    input_order = [name for (name, type) in signature]
    context = getContext(kwargs)

    if isinstance(ex, str):
        ex = stringToExpression(ex, types, context)

    # the AST is like the expression, but the node objects don't have
    # any odd interpretations

    ast = expressionToAST(ex)

    # Add a copy for strided or unaligned arrays
    for a in ast.postorderWalk():
        if a.astType == "variable" and a.value in copy_args:
            newVar = ASTNode(*a.key())
            a.astType, a.value, a.children = ('op', 'copy', (newVar,))

    if ex.astType not in ('op'):
        ast = ASTNode('op', value='copy', astKind=ex.astKind, children=(ast,))

    ast = typeCompileAst(ast)

    reg_num = [-1]
    def registerMaker(node, temporary=False):
        reg = Register(node, temporary=temporary)
        reg.n = reg_num[0]
        reg_num[0] -= 1
        return reg

    assignLeafRegisters(ast.allOf('raw'), Immediate)
    assignLeafRegisters(ast.allOf('variable', 'constant'), registerMaker)
    assignBranchRegisters(ast.allOf('op'), registerMaker)

    collapseDuplicateSubtrees(ast)

    input_order = getInputOrder(ast, input_order)
    constants_order, constants = getConstants(ast)
    
    if isReduction(ast):
        ast.reg.temporary = False
    
    optimizeTemporariesAllocation(ast)
    
    ast.reg.temporary = False
    r_output = 0
    ast.reg.n = 0
    
    r_inputs = r_output + 1
    r_constants = setOrderedRegisterNumbers(input_order, r_inputs)
    r_temps = setOrderedRegisterNumbers(constants_order, r_constants)
    r_end, tempsig = setRegisterNumbersForTemporaries(ast, r_temps)

    threeAddrProgram = convertASTtoThreeAddrForm(ast)
    input_names = tuple([a.value for a in input_order])
    signature = ''.join(types.get(x, float).__name__[0] for x in input_names)

    return threeAddrProgram, signature, tempsig, constants, input_names


def numexpr(ex, signature=(), copy_args=(), **kwargs):
    """Compile an expression built using E.<variable> variables to a function.

    ex can also be specified as a string "2*a+3*b".

    The order of the input variables and their types can be specified using the
    signature parameter, which is a list of (name, type) pairs.

    """
    threeAddrProgram, inputsig, tempsig, constants, input_names = \
                      precompile(ex, signature, copy_args, **kwargs)
    program = compileThreeAddrForm(threeAddrProgram)
    return interpreter.NumExpr(inputsig, tempsig, program, constants,
                               input_names)


def disassemble(nex):
    """Given a NumExpr object, return a list which is the program
    disassembled.
    """
    rev_opcodes = {}
    for op in interpreter.opcodes:
        rev_opcodes[interpreter.opcodes[op]] = op
    r_constants = 1 + len(nex.signature)
    r_temps = r_constants + len(nex.constants)
    def getArg(pc, offset):
        arg = ord(nex.program[pc+offset])
        op = rev_opcodes.get(ord(nex.program[pc]))
        code = op.split('_')[1][offset-1]
        if arg == 255:
            return None
        if code != 'n':
            if arg == 0:
                return 'r0'
            elif arg < r_constants:
                return 'r%d[%s]' % (arg, nex.input_names[arg-1])
            elif arg < r_temps:
                return 'c%d[%s]' % (arg, nex.constants[arg - r_constants])
            else:
                return 't%d' % (arg,)
        else:
            return arg
    source = []
    for pc in range(0, len(nex.program), 4):
        op = rev_opcodes.get(ord(nex.program[pc]))
        dest = getArg(pc, 1)
        arg1 = getArg(pc, 2)
        arg2 = getArg(pc, 3)
        source.append( (op, dest, arg1, arg2) )
    return source


def getType(a):
    t = a.dtype.type
    if issubclass(t, numpy.bool_):
        return bool
    if issubclass(t, numpy.integer):
        return int
    if issubclass(t, numpy.floating):
        return float
    if issubclass(t, numpy.complexfloating):
        return complex
    raise ValueError("unkown type %s" % a.dtype.name)


def getExprNames(text, context):
    ex = stringToExpression(text, {}, context)
    ast = expressionToAST(ex)
    input_order = getInputOrder(ast, None)
    return [a.value for a in input_order]


_names_cache = {}
_numexpr_cache = {}
def evaluate(ex, local_dict=None, global_dict=None, **kwargs):
    """Evaluate a simple array expression elementwise.

    ex is a string forming an expression, like "2*a+3*b". The values for "a"
    and "b" will by default be taken from the calling function's frame
    (through use of sys._getframe()). Alternatively, they can be specifed
    using the 'local_dict' or 'global_dict' arguments.

    Not all operations are supported, and only real
    constants and arrays of floats currently work.
    """
    if not isinstance(ex, str):
        raise ValueError("must specify expression as a string")
    # Get the names for this expression
    expr_key = (ex, tuple(sorted(kwargs.items())))
    if expr_key not in _names_cache:
        context = getContext(kwargs)
        _names_cache[expr_key] = getExprNames(ex, context)
    names = _names_cache[expr_key]
    # Get the arguments based on the names.
    call_frame = sys._getframe(1)
    if local_dict is None:
        local_dict = call_frame.f_locals
    if global_dict is None:
        global_dict = call_frame.f_globals
    arguments = []
    copy_args = []
    for name in names:
        try:
            a = local_dict[name]
        except KeyError:
            a = global_dict[name]
        # byteswapped arrays are taken care of in the extension.
        arguments.append(numpy.asarray(a)) # don't make a data copy, if possible
        if (hasattr(a, "flags") and            # numpy object
            (not a.flags.contiguous or
             not a.flags.aligned)):
            copy_args.append(name)    # do a copy to temporary

    # Create a signature
    signature = [(name, getType(arg)) for (name, arg) in zip(names, arguments)]
    # Look up numexpr if possible. copy_args *must* be added to the key,
    # just in case a non-copy expression is already in cache.
    numexpr_key = expr_key + (tuple(signature),) + tuple(copy_args)
    try:
        compiled_ex = _numexpr_cache[numexpr_key]
    except KeyError:
        compiled_ex = _numexpr_cache[numexpr_key] = \
                      numexpr(ex, signature, copy_args, **kwargs)
    return compiled_ex(*arguments)


