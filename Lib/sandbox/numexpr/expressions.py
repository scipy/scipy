__all__ = ['E']

import operator
import numpy
import sys

import interpreter

class Expression(object):
    def __init__(self):
        object.__init__(self)

    def __getattr__(self, name):
        if name.startswith('_'):
            return self.__dict__[name]
        else:
            return VariableNode(name, 'float')

E = Expression()


def get_context():
    """Context used to evaluate expression. Typically overridden in compiler."""
    return {}
def get_optimization():
    return get_context().get('optimization', 'none')

# helper functions for creating __magic__ methods
def ophelper(f):
    def func(*args):
        args = list(args)
        for i, x in enumerate(args):
            if isConstant(x):
                args[i] = x = ConstantNode(x)
            if not isinstance(x, ExpressionNode):
                return NotImplemented
        return f(*args)
    func.__name__ = f.__name__
    func.__doc__ = f.__doc__
    func.__dict__.update(f.__dict__)
    return func

def allConstantNodes(args):
    "returns True if args are all ConstantNodes."
    for x in args:
        if not isinstance(x, ConstantNode):
            return False
    return True

def isConstant(ex):
    "Returns True if ex is a constant scalar of an allowed type."
    return isinstance(ex, (bool, int, float, complex))

type_to_kind = {bool: 'bool', int: 'int', float: 'float', complex: 'complex'}
kind_to_type = {'bool': bool, 'int': int, 'float': float, 'complex': complex}
kind_rank = ['bool', 'int', 'float', 'complex', 'none']

def commonKind(nodes):
    n = -1
    for x in nodes:
        n = max(n, kind_rank.index(x.astKind))
    return kind_rank[n]

def bestConstantType(x):
    for converter in bool, int, float, complex:
        try:
            y = converter(x)
        except StandardError, err:
            continue
        if x == y:
            return converter

def getKind(x):
    converter = bestConstantType(x)
    return type_to_kind[converter]

def binop(opname, reversed=False, kind=None):
    @ophelper
    def operation(self, other):
        if reversed:
            self, other = other, self
        if allConstantNodes([self, other]):
            return ConstantNode(getattr(self.value, "__%s__" % opname)(other.value))
        else:
            return OpNode(opname, (self, other), kind=kind)
    return operation

def func(func, minkind=None):
    @ophelper
    def function(*args):
        if allConstantNodes(args):
            return ConstantNode(func(*[x.value for x in args]))
        kind = commonKind(args)
        if minkind and kind_rank.index(minkind) > kind_rank.index(kind):
            kind = minkind
        return FuncNode(func.__name__, args, kind)
    return function

@ophelper
def where_func(a, b, c):
    if isinstance(a, ConstantNode):
        raise ValueError("too many dimensions")
    if allConstantNodes([a,b,c]):
        return ConstantNode(numpy.where(a, b, c))
    return FuncNode('where', [a,b,c])

def encode_axis(axis):
    if isinstance(axis, ConstantNode):
        axis = axis.value
    if axis is None:
        axis = interpreter.allaxes
    else:
        if axis < 0:
            axis = interpreter.maxdims - axis
        if axis > 254:
            raise ValueError("cannot encode axis")
    return RawNode(axis)

def sum_func(a, axis=-1):
    axis = encode_axis(axis)
    if isinstance(a, ConstantNode):
        return a
    if isinstance(a, (bool, int, float, complex)):
        a = ConstantNode(a)
    kind = a.astKind
    if kind == 'bool':
        kind = 'int'
    return FuncNode('sum', [a, axis], kind=kind)

def prod_func(a, axis=-1):
    axis = encode_axis(axis)
    if isinstance(a, (bool, int, float, complex)):
        a = ConstantNode(a)
    if isinstance(a, ConstantNode):
        return a
    kind = a.astKind
    if kind == 'bool':
        kind = 'int'
    return FuncNode('prod', [a, axis], kind=kind)

@ophelper
def div_op(a, b):
    if get_optimization() in ('moderate', 'aggressive'):
        if isinstance(b, ConstantNode) and (a.astKind == b.astKind) and a.astKind in ('float', 'complex'):
            return OpNode('mul', [a, ConstantNode(1./b.value)])
    return OpNode('div', [a,b])

@ophelper
def pow_op(a, b):
    if allConstantNodes([a,b]):
        return ConstantNode(a**b)
    if isinstance(b, ConstantNode):
        x = b.value
        if get_optimization() == 'aggressive':
            RANGE = 50 # Approximate break even point with pow(x,y)
            # Optimize all integral and half integral powers in [-RANGE, RANGE]
            # Note: for complex numbers RANGE could be larger.
            if (int(2*x) == 2*x) and (-RANGE <= abs(x) <= RANGE):
                n = int(abs(x))
                ishalfpower = int(abs(2*x)) % 2
                def multiply(x, y):
                    if x is None: return y
                    return OpNode('mul', [x, y])
                r = None
                p = a
                mask = 1
                while True:
                    if (n & mask):
                        r = multiply(r, p)
                    mask <<= 1
                    if mask > n:
                        break
                    p = OpNode('mul', [p,p])
                if ishalfpower:
                    kind = commonKind([a])
                    if kind == 'int': kind = 'float'
                    r = multiply(r, OpNode('sqrt', [a], kind))
                if r is None:
                    r = OpNode('ones_like', [a])
                if x < 0:
                    r = OpNode('div', [ConstantNode(1), r])
                return r
        if get_optimization() in ('moderate', 'aggressive'):
            if x == -1:
                return OpNode('div', [ConstantNode(1),a])
            if x == 0:
                return FuncNode('ones_like', [a])
            if x == 0.5:
                kind = a.astKind
                if kind == 'int': kind = 'float'
                return FuncNode('sqrt', [a], kind=kind)
            if x == 1:
                return a
            if x == 2:
                return OpNode('mul', [a,a])
    return OpNode('pow', [a,b])

functions = {
    'copy' : func(numpy.copy),
    'ones_like' : func(numpy.ones_like),
    'sin' : func(numpy.sin, 'float'),
    'cos' : func(numpy.cos, 'float'),
    'tan' : func(numpy.tan, 'float'),
    'sqrt' : func(numpy.sqrt, 'float'),

    'sinh' : func(numpy.sinh, 'float'),
    'cosh' : func(numpy.cosh, 'float'),
    'tanh' : func(numpy.tanh, 'float'),

    'arctan2' : func(numpy.arctan2, 'float'),
    'fmod' : func(numpy.fmod, 'float'),

    'where' : where_func,

    'complex' : func(complex, 'complex'),
        
    'sum' : sum_func,
    'prod' : prod_func,
    }


class ExpressionNode(object):
    """An object that represents a generic number object.

    This implements the number special methods so that we can keep
    track of how this object has been used.
    """
    astType = 'generic'

    def __init__(self, value=None, kind=None, children=None):
        object.__init__(self)
        self.value = value
        if kind is None:
            kind = 'none'
        self.astKind = kind
        if children is None:
            self.children = ()
        else:
            self.children = tuple(children)

    def get_real(self):
        if self.astType == 'constant':
            return ConstantNode(complex(self.value).real)
        return OpNode('real', (self,), 'float')
    real = property(get_real)

    def get_imag(self):
        if self.astType == 'constant':
            return ConstantNode(complex(self.value).imag)
        return OpNode('imag', (self,), 'float')
    imag = property(get_imag)

    def __str__(self):
        return '%s(%s, %s, %s)' % (self.__class__.__name__, self.value,
                                   self.astKind, self.children)
    def __repr__(self):
        return self.__str__()

    def __neg__(self):
        return OpNode('neg', (self,))
    def __invert__(self):
        return OpNode('invert', (self,))
    def __pos__(self):
        return self

    __add__ = __radd__ = binop('add')
    __sub__ = binop('sub')
    __rsub__ = binop('sub', reversed=True)
    __mul__ = __rmul__ = binop('mul')
    __div__ =  div_op
    __rdiv__ = binop('div', reversed=True)
    __pow__ = pow_op
    __rpow__ = binop('pow', reversed=True)
    __mod__ = binop('mod')
    __rmod__ = binop('mod', reversed=True)

    # boolean operations

    __and__ = binop('and', kind='bool')
    __or__ = binop('or', kind='bool')

    __gt__ = binop('gt', kind='bool')
    __ge__ = binop('ge', kind='bool')
    __eq__ = binop('eq', kind='bool')
    __ne__ = binop('ne', kind='bool')
    __lt__ = binop('gt', reversed=True, kind='bool')
    __le__ = binop('ge', reversed=True, kind='bool')



class LeafNode(ExpressionNode):
    leafNode = True

class VariableNode(LeafNode):
    astType = 'variable'
    def __init__(self, value=None, kind=None, children=None):
        ExpressionNode.__init__(self, value=value, kind=kind)

class RawNode(object):
    """Used to pass raw integers to interpreter.
    For instance, for selecting what function to use in func1.
    Purposely don't inherit from ExpressionNode, since we don't wan't
    this to be used for anything but being walked.
    """
    astType = 'raw'
    astKind = 'none'
    def __init__(self, value):
        self.value = value
        self.children = ()
    def __str__(self):
        return 'RawNode(%s)' % (self.value,)
    __repr__ = __str__


class ConstantNode(LeafNode):
    astType = 'constant'
    def __init__(self, value=None, children=None):
        kind = getKind(value)
        LeafNode.__init__(self, value=value, kind=kind)
    def __neg__(self):
        return ConstantNode(-self.value)
    def __invert__(self):
        return ConstantNode(~self.value)

class OpNode(ExpressionNode):
    astType = 'op'
    def __init__(self, opcode=None, args=None, kind=None):
        if (kind is None) and (args is not None):
            kind = commonKind(args)
        ExpressionNode.__init__(self, value=opcode, kind=kind, children=args)

class FuncNode(OpNode):
    def __init__(self, opcode=None, args=None, kind=None):
        if (kind is None) and (args is not None):
            kind = commonKind(args)
        OpNode.__init__(self, opcode, args, kind)
