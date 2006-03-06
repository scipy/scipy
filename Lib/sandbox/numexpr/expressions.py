__all__ = ['E', 'numexpr', 'evaluate']

import sys
import operator
import numpy
from numexpr import interpreter

class Expression(object):
    def __init__(self):
        object.__init__(self)

    def __getattr__(self, name):
        if name.startswith('_'):
            return self.__dict__[name]
        else:
            return VariableNode(name)

E = Expression()

class Register(object):
    def __init__(self, n, name=None, temporary=False, constant=False,
                 value=None):
        self.n = n
        self.name = name
        self.temporary = temporary
        self.constant = constant
        self.value = value

    def __str__(self):
        if self.temporary:
            return 'Temporary(%d)' % (self.n,)
        elif self.constant:
            return 'Constant(%d, %s)' % (self.n, self.value)
        elif self.name:
            return 'Register(%d, %s)' % (self.n, self.name)
        else:
            return 'Register(%d)' % (self.n,)
    __repr__ = __str__




def string2expression(s):
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

def numexpr(ex, input_order=None, precompiled=False):
    """Compile an expression built using E.<variable> variables to a function.

    ex can also be specified as a string "2*a+3*b".

    The order of the input variables can be specified using the input_order
    parameter.

    If precompiled is set to True, a nonusable, easy-to-read representation of
    the bytecode program used is returned instead.
    """
    if isinstance(ex, str):
        ex = string2expression(ex)
    if isinstance(ex, ConstantNode):
        return ConstantNumexpr(ex._value)
    if not isinstance(ex, OpNode):
        ex = OpNode('copy', (ex,))

    # canonicalize the variables
    node_map = {}
    name_map = {}
    for a in ex.walk(instances_of=VariableNode):
        if a._name in name_map:
            node_map[a] = name_map[a._name]
        else:
            name_map[a._name] = a
            node_map[a] = a

    if input_order:
        for name in input_order:
            if name not in name_map:
                raise ValueError("input name %r not found in expression" %
                                 (name,))
        for name in name_map:
            if name not in input_order:
                raise ValueError("input name %r not specified in order" %
                                 (name,))
    else:
        input_order = name_map.keys()
        input_order.sort()
    n_inputs = len(input_order)

    # canonicalize the constants
    const_map = {}
    for a in ex.walk(instances_of=ConstantNode):
        if a._value in const_map:
            node_map[a] = const_map[a._value]
        else:
            const_map[a._value] = a
            node_map[a] = a

    # coerce equal nodes to one node when walking the tree
    def walk(instances_of=None):
        for a in ex.walk(instances_of=instances_of):
            yield node_map.get(a, a)

    # assign registers
    registers = {}
    # variables are 1 .. #of variables
    for i, name in enumerate(input_order):
        v = name_map[name]
        registers[v] = Register(i+1, name=name)
    # assign temporaries, without registers set yet
    for i, a in enumerate(walk(instances_of=OpNode)):
        registers[a] = Register(None, temporary=True)
    # assign constant locations
    constants = numpy.empty((len(const_map),), dtype=float)
    for i, value in enumerate(const_map):
        c = const_map[value]
        registers[c] = Register(i, constant=True, value=value)
        constants[i] = value

    for a in walk(instances_of=OpNode):
        # put result in one of the operand temporaries if there is one
        for arg in a.walk(instances_of=ExpressionNode):
            arg = node_map.get(arg, arg)
            if registers[arg].temporary:
                registers[a] = registers[arg]
                break

    # output is register 0
    registers[ex].n = 0
    registers[ex].temporary = False

    seen_temps = set()
    def reg(a):
        if isinstance(a, RawNode):
            return Register(a._value)
        a = node_map.get(a, a)
        m = registers[a]
        if m not in seen_temps and m.temporary:
            m.n = 1 + n_inputs + len(seen_temps)
            seen_temps.add(m)
        return m
    def to_string(opcode, store, a1, a2):
        cop = chr(interpreter.opcodes[opcode])
        cs = chr(store.n)
        if a1 is None:
            ca1 = chr(0)
        else:
            ca1 = chr(a1.n)
        if a2 is None:
            ca2 = chr(0)
        else:
            ca2 = chr(a2.n)
        return cop + cs + ca1 + ca2
    program = []
    for a in ex.walk(instances_of=OpNode):
        if len(a._args) == 1:
            program.append( (a._opcode, reg(a), reg(a._args[0]), None) )
        elif len(a._args) == 2:
            program.append( (a._opcode, reg(a),
                             reg(a._args[0]), reg(a._args[1])))
        elif len(a._args) == 3:
            program.append( (a._opcode, reg(a),
                             reg(a._args[0]), reg(a._args[1])))
            program.append( ('noop', reg(a._args[2]), None, None) )
        else:
            raise ValueError('too many arguments')

    if precompiled:
        return program

    prog_str = ''.join([to_string(*t) for t in program])
    n_temps = len(seen_temps)
    nex = interpreter.NumExpr(n_inputs=n_inputs, n_temps=n_temps,
                              program=prog_str, constants=constants,
                              input_names=tuple(input_order))
    return nex

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

def binop(opname, reversed=False):
    def operation(self, other):
        if isinstance(other, (int, float)):
            other = ConstantNode(other)
        elif not isinstance(other, ExpressionNode):
            return NotImplemented
        if reversed:
            return OpNode(opname, (other, self))
        else:
            return OpNode(opname, (self, other))
    return operation


def compareop(opname):
    reverse = {'lt' : 'gt',
               'le' : 'ge',
               'eq' : 'eq',
               'ne' : 'ne'}
    def operation(self, other):
        if isinstance(other, (int, float)):
            other = ConstantNode(other)
        elif not isinstance(other, ExpressionNode):
            return NotImplemented
        if isinstance(self, ConstantNode) and isinstance(other, ConstantNode):
            return ConstantNode(getattr(self.__val__, "__%s__" % opname)(other))
        if isinstance(self, ConstantNode):
            return OpNode(reverse[opname] + '_c', (other, self))
        if isinstance(other, ConstantNode):
            return OpNode(opname + '_c', (other, self))
        if opname in reverse:
            return OpNode(reverse[opname], (other, self))
        return OpNode(opname, (self, other))
    return operation

def func(func):
    def function(*args):
        args = list(args)
        constant = True
        expression = True
        for i, x in enumerate(args):
            if isinstance(x, (int, float)):
                args[i] = x = ConstantNode(x)
            if not isinstance(x, ConstantNode):
                constant = False
            if not isinstance(x, ExpressionNode):
                expression = False
        if constant:
            return ConstantNode(func(*[x._value for x in args]))
        elif not expression:
            return NotImplemented
        elif (func.__name__ == 'where') and isinstance(args[0], ConstantNode):
            raise ValueError("too many dimensions")
        return FuncNode(func.__name__, args)
    return function

functions = {
    'sin' : func(numpy.sin),
    'cos' : func(numpy.cos),
    'tan' : func(numpy.tan),
    'sinh' : func(numpy.sinh),
    'cosh' : func(numpy.cosh),
    'tanh' : func(numpy.tanh),
    
    'arctan2' : func(numpy.arctan2),
    'fmod' : func(numpy.fmod),
    
    'where' : func(numpy.where)
            }


class ExpressionNode(object):
    def __init__(self):
        object.__init__(self)

    def __repr__(self):
        return self.__str__()

    __add__ = __radd__ = binop('add')

    def __neg__(self):
        return OpNode('neg', (self,))
    def __pos__(self):
        return self

    __sub__ = binop('sub')
    __rsub__ = binop('sub', reversed=True)
    __mul__ = __rmul__ = binop('mul')
    __div__ =  binop('div')
    __rdiv__ = binop('div', reversed=True)
    __pow__ = binop('pow')
    __rpow__ = binop('pow', reversed=True)
    __mod__ = binop('mod')
    __rmod__ = binop('mod', reversed=True)

    __gt__ = compareop('gt')
    __ge__ = compareop('ge')
    __eq__ = compareop('eq')
    __ne__ = compareop('ne')
    __lt__ = compareop('lt')
    __le__ = compareop('le')

    def _walk(self):
        yield self

    def walk(self, instances_of=None):
        if instances_of is not None:
            for a in self._walk():
                if isinstance(a, instances_of):
                    yield a
        else:
            for a in self._walk():
                yield a

class VariableNode(ExpressionNode):
    def __init__(self, name):
        ExpressionNode.__init__(self)
        self._name = name

    def __str__(self):
        return 'VariableNode(%s)' % (self._name,)

def optimize_constants(name, op):
    def operation(self, other):
        if isinstance(other, (int, float)):
            return ConstantNode(op(self._value, other))
        elif isinstance(other, ConstantNode):
            return ConstantNode(op(self._value, other._value))
        else:
            a = getattr(ExpressionNode, name)
            return a(self, other)
    return operation

class RawNode(object):
    """Used to pass raw integers to interpreter. 
    For instance, for selecting what function to use in func1.
    Purposely don't inherit from ExpressionNode, since we don't wan't 
    this to be used for anything but being walked.
    """
    def __init__(self, value):
        self._value = value
    def __str__(self):
        return 'RawNode(%s)' % (self._value,)
    __repr__ = __str__
    def walk(self, instances_of=None):
        if instances_of is not None:
            yield self
        
class ConstantNode(ExpressionNode):
    def __init__(self, value):
        ExpressionNode.__init__(self)
        self._value = value

    def __str__(self):
        return 'ConstantNode(%s)' % (self._value,)

    __add__ = optimize_constants('__add__', operator.add)
    __sub__ = optimize_constants('__sub__', operator.sub)
    __mul__ = optimize_constants('__mul__', operator.mul)
    __div__ = optimize_constants('__div__', operator.div)
    __pow__ = optimize_constants('__pow__', operator.pow)
    __mod__ = optimize_constants('__mod__', operator.mod)

    def __neg__(self):
        return ConstantNode(-self._value)
    

class OpNode(ExpressionNode):
    def __init__(self, opcode, args):
        ExpressionNode.__init__(self)
        if len(args) == 2:
            if isinstance(args[0], ConstantNode):
                if opcode in ('pow', 'mod'):
                    cp = OpNode('copy', (RawNode(0), args[0]))
                    args = (cp, args[1])
                else:
                    opcode += '_c'
                    args = (args[1], args[0])
            elif isinstance(args[1], ConstantNode):
                # constant goes last, although the op is constant - value
                if opcode == 'sub':
                    opcode = 'add_c'
                    args = args[0], -args[1]
                elif opcode == 'div':
                    opcode = 'mul_c'
                    args = args[0], ConstantNode(1./args[1]._value)
                else:
                    opcode += '_c'
        self._opcode = opcode
        self._args = args

    def __str__(self):
        return 'OpNode(%r, %s)' % (self._opcode, self._args)

    def _walk(self):
        for a in self._args:
            for w in a.walk():
                yield w
        yield self



class FuncNode(OpNode):
    
    def __init__(self, opcode, args):
        ExpressionNode.__init__(self)
        # There are three cases:
        #    1. Function is inline and can be called without copying constants
        #    2. Function is inline but constants must be copied
        #    3. Function is called from function table. Constants always copied.
        inline = (opcode in interpreter.opcodes)
        has_constants = bool([x for x in args if isinstance(x, ConstantNode)])
        copy_constants = True
        if inline:
            if has_constants:
                sig = '_'+''.join('xc'[isinstance(x, ConstantNode)] for x in args)
                if (opcode+sig) in interpreter.opcodes: 
                    opcode += sig
                    copy_constants = False
        else:
            args = args + [RawNode(interpreter.funccodes[opcode])]
            opcode = 'func_%s' % (len(args)-1)
        self._opcode = opcode
        if copy_constants:
            for i, x in enumerate(args):
                if isinstance(x, ConstantNode):
                    args[i] = OpNode('copy', (RawNode(0), x))            
        # Put all constants and RawNodes last.
        CNode = (ConstantNode, RawNode)
        vars = [x for x in args if not isinstance(x, CNode)]
        consts = [x for x in args if isinstance(x, CNode)]
        self._args = vars + consts
   
    def __str__(self):
        return 'FuncNode(%r, %s)' % (self._opcode, self._args)

    def _walk(self):
        for a in self._args:
            for w in a.walk():
                yield w
        yield self
