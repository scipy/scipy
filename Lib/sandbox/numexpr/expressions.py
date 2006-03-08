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
            self, other = other, self
        if isinstance(self, ConstantNode) and isinstance(other, ConstantNode):
            return ConstantNode(getattr(self.value, "__%s__" % opname)(other.value))
        else:
            return BinopNode(opname, (self, other))
    return operation

def func(func):
    def function(*args):
        args = list(args)
        constant = True
        for i, x in enumerate(args):
            if isinstance(x, (int, float)):
                args[i] = x = ConstantNode(x)
            if not isinstance(x, ConstantNode):
                constant = False
            if not isinstance(x, ExpressionNode):
                return NotImplemented
        if constant:
            return ConstantNode(func(*[x.value for x in args])) 
        elif (func.__name__ == 'where') and isinstance(args[0], ConstantNode):
            # TAH: Would be nice to put this kind of special casing into the 
            # functions themselves at some point. Particularly if we get a bunch
            # more functions that require some checking.
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

    __gt__ = binop('gt')
    __ge__ = binop('ge')
    __eq__ = binop('eq')
    __ne__ = binop('ne')
    __lt__ = binop('gt', reversed=True)
    __le__ = binop('ge', reversed=True)

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


class ConstantNode(LeafNode):
    astType = 'constant'
    def __init__(self, value=None, children=None):
        LeafNode.__init__(self, value=float(value))
    def __neg__(self):
        return ConstantNode(-self.value)

class OpNode(ExpressionNode):
    astType = 'op'
    def __init__(self, opcode=None, args=None):
        ExpressionNode.__init__(self, value=opcode, children=args)    

class BinopNode(OpNode):
    def __init__(self, opcode=None, args=None):
        if (opcode == 'div') and isinstance(args[1], ConstantNode):
            opcode = 'mul'
            args = args[0], ConstantNode(1./args[1].value)
        OpNode.__init__(self, opcode, args)

class FuncNode(OpNode):
    def __init__(self, opcode=None, args=None):
        if opcode not in interpreter.opcodes:
            args = list(args) + [RawNode(interpreter.funccodes[opcode])]
            opcode = 'func_%s' % (len(args)-1)
        ExpressionNode.__init__(self, value=opcode, children=args)

