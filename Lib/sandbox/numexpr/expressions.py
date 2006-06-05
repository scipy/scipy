__all__ = ['E']

import operator
import numpy
from numexpr import interpreter
import sys

# XXX Is there any reason to keep Expression around?
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
            if isinstance(x, (int, float, complex)):
                args[i] = x = ConstantNode(x)
            if not isinstance(x, ExpressionNode):
                return NotImplemented
        return f(*args)
    return func

def all_constant(args):
    "returns true if args are all constant. Converts scalars to ConstantNodes"
    for x in args:
        if not isinstance(x, ConstantNode):
            return False
    return True

kind_rank = ['int', 'float', 'complex', 'none']
def common_kind(nodes):
    n = -1
    for x in nodes:
        n = max(n, kind_rank.index(x.astKind))
    return kind_rank[n]

def binop(opname, reversed=False, kind=None):
    @ophelper
    def operation(self, other):
        if reversed:
            self, other = other, self
        if all_constant([self, other]):
            return ConstantNode(getattr(self.value, "__%s__" % opname)(other.value))
        else:
            return OpNode(opname, (self, other), kind=kind)
    return operation

def func(func, minkind=None):
    @ophelper
    def function(*args):
        if all_constant(args):
            return ConstantNode(func(*[x.value for x in args]))
        kind = common_kind(args)
        if minkind and kind_rank.index(minkind) > kind_rank.index(kind):
            kind = minkind
        return FuncNode(func.__name__, args, kind)
    return function

@ophelper
def where_func(a, b, c):
    if all_constant([a,b,c]):
        return ConstantNode(numpy.where(a, b, c))
    if isinstance(a, ConstantNode):
        raise ValueError("too many dimensions")
    kind = common_kind([b,c])
    return FuncNode('where', [a,b,c])

@ophelper
def div_op(a, b):
    if get_optimization() in ('moderate', 'aggressive'):
        if isinstance(b, ConstantNode) and (a.astKind == b.astKind) and a.astKind in ('float', 'complex'):
            return OpNode('mul', [a, ConstantNode(1./b.value)])
    return OpNode('div', [a,b])


@ophelper
def pow_op(a, b):
    if all_constant([a,b]):
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
                    kind = common_kind([a])
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
            }

class ExpressionNode(object):
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
        return '%s(%s, %s, %s)' % (self.__class__.__name__, self.value, self.astKind,
                               self.children)
    def __repr__(self):
        return self.__str__()

    def __neg__(self):
        return OpNode('neg', (self,))
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

    __gt__ = binop('gt', kind='int')
    __ge__ = binop('ge', kind='int')
    __eq__ = binop('eq', kind='int')
    __ne__ = binop('ne', kind='int')
    __lt__ = binop('gt', reversed=True, kind='int')
    __le__ = binop('ge', reversed=True, kind='int')

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


def normalizeConstant(x):
    for converter in int, float, complex:
        try:
            y = converter(x)
        except StandardError, err:
            continue
        if x == y:
            return y

def getKind(x):
    return {int : 'int',
            float : 'float',
            complex : 'complex'}[type(normalizeConstant(x))]

class ConstantNode(LeafNode):
    astType = 'constant'
    def __init__(self, value=None, children=None):
        kind = getKind(value)
        LeafNode.__init__(self, value=value, kind=kind)
    def __neg__(self):
        return ConstantNode(-self.value)

class OpNode(ExpressionNode):
    astType = 'op'
    def __init__(self, opcode=None, args=None, kind=None):
        if (kind is None) and (args is not None):
            kind = common_kind(args)
        ExpressionNode.__init__(self, value=opcode, kind=kind, children=args)

class FuncNode(OpNode):
    def __init__(self, opcode=None, args=None, kind=None):
        if (kind is None) and (args is not None):
            kind = common_kind(args)
        OpNode.__init__(self, opcode, args, kind)
