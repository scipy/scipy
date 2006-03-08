
import sys
import numpy
from numexpr import interpreter

class ASTNode(object):
    def __init__(self, astType='generic', value=None, children=()):
        object.__init__(self)
        self.astType = astType
        self.value = value
        self.children = tuple(children)
        self.reg = None

    def __str__(self):
        return 'AST(%s, %s, %s, %s)' % (self.astType, self.value,
                                        self.children, self.reg)
    __repr__ = __str__

    def key(self):
        return (self.astType, self.value, self.children)

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
    this_ast = ASTNode(ex.astType, ex.value,
                           [expressionToAST(c) for c in ex.children])
    return this_ast

class ConstantNumexpr(object):
    """Used to wrap a constant when numexpr returns a constant value.
    """
    def __init__(self, value):
        self._value = value
        self.input_names = []
    def __call__(self, *args, **kargs):
        return self._value
    def __str__(self):
        return "%s(%s)" % (self.__class__.__name__, self._value)
    __repr__ = __str__


class Register(object):
    def __init__(self, astnode, temporary=False):
        self.node = astnode
        self.temporary = temporary
        self.n = None

    def __str__(self):
        if self.temporary:
            name = 'Temporary'
        else:
            name = 'Register'
        return '%s(%s, %s)' % (name, self.node.astType, self.n,)

    def __repr__(self):
        return self.__str__()


def stringToExpression(s):
    from numexpr.expressions import E, ConstantNode, functions
    # first compile to a code object to determine the names
    c = compile(s, '<expr>', 'eval')
    # make VariableNode's for the names
    names = {}
    for name in c.co_names:
        names[name] = getattr(E, name)
    names.update(functions)
    # now build the expression
    ex = eval(c, names)
    if isinstance(ex, (int, float)):
        ex = ConstantNode(ex)
    return ex


def getInputOrder(ast, input_order=None):
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

def getConstants(ast):
    const_map = {}
    for a in ast.allOf('constant'):
        const_map[a.value] = a
    ordered_constants = const_map.keys()
    ordered_constants.sort()
    constants_order = [const_map[v] for v in ordered_constants]
    constants = numpy.array([a.value for a in constants_order], dtype=float)
    return constants_order, constants

def sortNodesByOrder(nodes, order):
    order_map = {}
    for i, (_, v, _) in enumerate(order):
        order_map[v] = i
    dec_nodes = [(order_map[n.value], n) for n in nodes]
    dec_nodes.sort()
    return [a[1] for a in dec_nodes]

def assignLeafRegisters(inodes, registerMaker):
    leafRegisters = {}
    for node in inodes:
        key = node.key()
        if key in leafRegisters:
            node.reg = leafRegisters[key]
        else:
            node.reg = leafRegisters[key] = registerMaker(node)

def assignBranchRegisters(inodes, registerMaker):
    for node in inodes:
        node.reg = registerMaker(node, temporary=True)

def optimizeTemporariesAllocation(ast):
    for a in ast.allOf('op', 'function'):
        # put result in one of the operand temporaries if there is one
        for c in a.children:
            if c.reg.temporary:
                a.reg = c.reg
                break

def setOrderedRegisterNumbers(order, start):
    for i, node in enumerate(order):
        node.reg.n = start + i
    return start + len(order)

def setRegisterNumbersForTemporaries(ast, start):
    seen = 0
    for node in ast.postorderWalk():
        if node.astType == 'raw':
            node.reg.n = node.value
        reg = node.reg
        if reg.n < 0:
            reg.n = start + seen
            seen += 1
    return start + seen

def convertASTtoThreeAddrForm(ast):
    program = []
    for node in ast.allOf('op', 'function'):
        children = node.children
        instr = (node.value, node.reg.n) \
                + tuple([c.reg.n for c in children])
        program.append(instr)
    return program

def compileThreeAddrForm(program):
    def nToChr(n):
        if n is None:
            return '\xff'
        elif n < 0:
            raise ValueError("negative value for register number %s" % (n,))
        else:
            return chr(n)

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

def numexpr(ex, input_order=None, precompiled=False, debug=False):
    """Compile an expression built using E.<variable> variables to a function.

    ex can also be specified as a string "2*a+3*b".

    The order of the input variables can be specified using the input_order
    parameter.

    If precompiled is set to True, a nonusable, easy-to-read representation of
    the bytecode program used is returned instead.
    """
    if isinstance(ex, str):
        ex = stringToExpression(ex)
    if ex.astType == 'constant':
        return ConstantNumexpr(ex._value)

    # the AST is like the expression, but the node objects don't have
    # any odd interpretations
    ast = expressionToAST(ex)
    if ex.astType not in ('op', 'function'):
        ast = ASTNode('op', value='copy', children=(ast,))

    input_order = getInputOrder(ast, input_order)
    constants_order, constants = getConstants(ast)

    reg_num = [-1]
    def registerMaker(node, temporary=False):
        reg = Register(node, temporary=temporary)
        reg.n = reg_num[0]
        reg_num[0] -= 1
        return reg

    assignLeafRegisters(ast.allOf('variable', 'constant', 'raw'), registerMaker)
    assignBranchRegisters(ast.allOf('op', 'function'), registerMaker)

    if debug:
        print 'Before temporaries optimization'
        print ast

    optimizeTemporariesAllocation(ast)

    r_output = 0
    ast.reg.n = 0
    ast.reg.temporary = False
    r_inputs = r_output + 1
    r_constants = setOrderedRegisterNumbers(input_order, r_inputs)
    r_temps = setOrderedRegisterNumbers(constants_order, r_constants)
    r_end = setRegisterNumbersForTemporaries(ast, r_temps)
    n_temps = r_end - r_temps

    if debug:
        print 'After register number assignments'
        print ast

    threeAddrProgram = convertASTtoThreeAddrForm(ast)

    if precompiled:
        return threeAddrProgram

    program = compileThreeAddrForm(threeAddrProgram)
    input_names = tuple([a.value for a in input_order])
    nex = interpreter.NumExpr(n_inputs=len(input_order),
                              n_temps=n_temps,
                              program=program,
                              constants=constants,
                              input_names=input_names)
    return nex

def disassemble(nex):
    rev_opcodes = {}
    for op in interpreter.opcodes:
        rev_opcodes[interpreter.opcodes[op]] = op
    r_constants = 1 + nex.n_inputs
    r_temps = r_constants + len(nex.constants)
    def get_arg(pc):
        arg = ord(nex.program[pc])
        if arg == 0:
            return 'r0'
        elif arg == 255:
            return None
        elif arg < r_constants:
            return 'r%d[%s]' % (arg, nex.input_names[arg-1])
        elif arg < r_temps:
            return 'c%d[%s]' % (arg, nex.constants[arg - r_constants])
        else:
            return 't%d' % (arg,)
    source = []
    for pc in range(0, len(nex.program), 4):
        op = rev_opcodes.get(ord(nex.program[pc]))
        dest = get_arg(pc+1)
        arg1 = get_arg(pc+2)
        arg2 = get_arg(pc+3)
        source.append( (op, dest, arg1, arg2) )
    return source

_numexpr_cache = {}
def evaluate(ex, local_dict=None, global_dict=None):
    """Evaluate a simple array expression elementwise.

    ex is a string forming an expression, like "2*a+3*b". The values for "a"
    and "b" will by default be taken from the calling function's frame
    (through use of sys._getframe()). Alternatively, they can be specifed
    using the 'local_dict' or 'global_dict' arguments.

    Not all operations are supported, and only on real
    constants and arrays of floats currently work..
    """
    if not isinstance(ex, str):
        raise ValueError("must specify expression as a string")
    try:
        compiled_ex = _numexpr_cache[ex]
    except KeyError:
        compiled_ex = _numexpr_cache[ex] = numexpr(ex)
    call_frame = sys._getframe().f_back
    if local_dict is None:
        local_dict = call_frame.f_locals
    if global_dict is None:
        global_dict = call_frame.f_globals
    arguments = []
    for name in compiled_ex.input_names:
        try:
            a = local_dict[name]
        except KeyError:
            a = global_dict[name]
        arguments.append(a)
    return compiled_ex(*arguments)
