__all__ = ['E']

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

# helper functions for creating __magic__ methods

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
            return OpNode(reverse[opname], (other, self))
        if isinstance(other, ConstantNode):
            return OpNode(opname, (other, self))
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
            return ConstantNode(func(*[x.value for x in args]))
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
    astType = 'generic'
    def __init__(self, value=None, children=None):
        object.__init__(self)
        self.value = value
        if children is None:
            self.children = ()
        else:
            self.children = tuple(children)

    def __str__(self):
        return '%s(%s, %s)' % (self.__class__.__name__, self.value,
                               self.children)

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

class LeafNode(ExpressionNode):
    leafNode = True

class VariableNode(LeafNode):
    astType = 'variable'
    def __init__(self, value=None, children=None):
        ExpressionNode.__init__(self, value=value)

class RawNode(object):
    """Used to pass raw integers to interpreter.
    For instance, for selecting what function to use in func1.
    Purposely don't inherit from ExpressionNode, since we don't wan't
    this to be used for anything but being walked.
    """
    astType = 'raw'
    def __init__(self, value):
        self.value = value
        self.children = ()
    def __str__(self):
        return 'RawNode(%s)' % (self.value,)
    __repr__ = __str__

def optimize_constants(name, op):
    def operation(self, other):
        if isinstance(other, (int, float)):
            return ConstantNode(op(self.value, other))
        elif isinstance(other, ConstantNode):
            return ConstantNode(op(self.value, other.value))
        else:
            a = getattr(ExpressionNode, name)
            return a(self, other)
    return operation

class ConstantNode(LeafNode):
    astType = 'constant'
    def __init__(self, value=None, children=None):
        LeafNode.__init__(self, value=float(value))

    __add__ = optimize_constants('__add__', operator.add)
    __sub__ = optimize_constants('__sub__', operator.sub)
    __mul__ = optimize_constants('__mul__', operator.mul)
    __div__ = optimize_constants('__div__', operator.div)
    __pow__ = optimize_constants('__pow__', operator.pow)
    __mod__ = optimize_constants('__mod__', operator.mod)

    def __neg__(self):
        return ConstantNode(-self.value)

class OpNode(ExpressionNode):
    astType = 'op'
    def __init__(self, opcode=None, args=None):
        if len(args) == 2:
            if isinstance(args[0], ConstantNode):
                if opcode == 'div':
                    opcode += '_c'
                    args = (args[1], args[0])
            elif isinstance(args[1], ConstantNode):
                # constant goes last, although the op is constant - value
                if opcode == 'sub':
                    opcode = 'add'
                    args = args[0], -args[1]
                elif opcode == 'div':
                    opcode = 'mul'
                    args = args[0], ConstantNode(1./args[1].value)
        ExpressionNode.__init__(self, value=opcode, children=args)

class FuncNode(OpNode):
    astType = 'function'
    def __init__(self, opcode=None, args=None):
        # There are three cases:
        #    1. Function is inline and can be called without copying constants
        #    2. Function is inline but constants must be copied
        #    3. Function is called from function table. Constants always copied.
        if opcode not in interpreter.opcodes:
            args = args + [RawNode(interpreter.funccodes[opcode])]
            opcode = 'func_%s' % (len(args)-1)
        OpNode.__init__(self, opcode, args)
