import sets, types, re, string, csv, copy
import numpy as N

class Term:

    """
    This class is very simple: it is just a named term in a model formula.

    It is also callable: by default it namespace[self.name], where namespace
    defaults to globals().
    
    """

    def __init__(self, name, _fn=None, termname=None):
        
        self.name = name

        if termname is None:
            self.termname = name
        else:
            self.termname = termname

        if type(self.termname) is not types.StringType:
            raise ValueError, 'expecting a string for termname'
        if _fn:
            self._fn = _fn

    def __str__(self):
        """
        '<term: %s>' % self.termname
        """
        return '<term: %s>' % self.termname

    def __add__(self, other):
        """
        Formula(self) + Formula(other)
        """
        other = Formula(other)
        return other + self

    def __mul__(self, other):
        """
        Formula(self) * Formula(other)
        """

        if other.name is 'intercept':
            return Formula(self)
        elif self.name is 'intercept':
            return Formula(other)
        
        other = Formula(other)
        return other * self

    def names(self):
        """
        Return the names of the columns in design associated to self,
        i.e. len(self.names()) = self().shape[0].
        """
        if type(self.name) is types.StringType:
            return [self.name]
        else:
            return list(self.name)

    def __call__(self, namespace=None, usefn=True, **extra):
        """
        Return the columns associated to self in a design matrix.
        The default behaviour is to return namespace[self.termname]
        where namespace defaults to globals().

        If usefn, and self._fn exists then return
            self._fn(namespace=namspace, **extra).
        """
        
        if namespace is None:
            namespace = globals()
	if not hasattr(self, '_fn') or not usefn:
            val = namespace[self.termname]
            if isinstance(val, Formula):
                val = val(namespace, **extra)
	else:
            val = self._fn(namespace=namespace, **extra)
        val = N.asarray(val)
        return N.squeeze(val)

class Factor(Term):

    """
    A categorical factor.
    """

    def __init__(self, termname, keys, ordinal=False):
        """
        Factor is initialized with keys, representing all valid
        levels of the factor.
        """
        
        self.keys = list(sets.Set(keys))
        self.keys.sort()
        self._name = termname
        self.termname = termname
        self.ordinal = ordinal

        if self.ordinal:
            self._sort = True

            def _fn(namespace=None, key=key):
                v = namespace[self._name]
                col = [float(self.keys.index(v[i])) for i in range(n)]
                return N.array(col)
            Term.__init__(self, self.name, _fn=_fn)

        else:
            def _fn(namespace=None):
                v = namespace[self._name]
                value = []
                for key in self.keys:
                    col = [float((v[i] == key)) for i in range(len(v))]
                    value.append(col)
                return N.array(value)
            Term.__init__(self, ['(%s==%s)' % (self.termname, str(key)) for key in self.keys], _fn=_fn, termname=self.termname)

    def __call__(self, namespace=None, values=False, **extra):
        """
        Return either the columns in the design matrix, or the
        actual values of the factor, if values==True.
        """

        if namespace is None:
            namespace = globals()
        if not values:
            return Term.__call__(self, namespace=namespace, usefn=True, **extra)
        else:
            return Term.__call__(self, namespace=namespace, usefn=False, **extra)

    def verify(self, values):
        """
        Verify that all values correspond to valid keys in self.
        """
        s = sets.Set(values)
        if not s.issubset(self.keys):
            raise ValueError, 'unknown keys in values'

    def __add__(self, other):
        """
        Formula(self) + Formula(other)

        When adding \'intercept\' to a Factor, this just returns self.
        """
        
        if other.name is 'intercept':
            return Formula(self)
        else:
            return Term.__add__(self, other)

    def main_effect(self, reference=None):
        """
        Return the 'main effect' columns of a Factor, choosing
        a reference column number to remove.
        """

        if reference is None:
            reference = 0

        def _fn(namespace=None, reference=reference, names=self.names(), **keywords):
            value = N.asarray(self(namespace=namespace, **keywords))
            rvalue = []
            keep = range(value.shape[0])
            keep.pop(reference)
            for i in range(len(keep)):
                rvalue.append(value[keep[i]] - value[reference])
            return rvalue
        
        keep = range(len(self.names()))
        keep.pop(reference)
        __names = self.names()
        _names = ['%s-%s' % (__names[keep[i]], __names[reference]) for i in range(len(keep))]
        return Term(_names, _fn=_fn, termname='%s:maineffect' % self.termname)

class Quantitative(Term):

    """
    A subclass of Term that presumes namespace[self.termname] is
    an ndarray.

    Basically used for __pow__ method and (looking forward) for splines.

    """

    def __pow__(self, power):
        """
        Raise the quantitative Term's values to an integer power, i.e.
        polynomial.
        """
        if type(power) is not types.IntType:
            raise ValueError, 'expecting an integer'

        name = '%s^%d' % (self.name, power)

        def _fn(namespace=None, power=power):
            x = N.array(namespace[self.name])
            return N.power(x, power)
        return Term(name, _fn=_fn)

class FuncQuant(Quantitative):

    """
    A Term for a quantitative function of a Term.
    """

    counter = 0

    def __init__(self, x, f):
        """
        Return a term whose values are f(x(namespace=namespace)).
        """
        
        self.f = f
        self.x = x
        def _fn(namespace=None, f=self.f):
            x = namespace[x.name]
            return f(x)
        try:
            termname = '%s(%s)' % (f.func_name, quant.name)
        except:
            termname = 'f%d(%s)' % (FuncQuant.counter, quant.name)
            FuncQuant.counter += 1
        Term.__init__(self, termname, _fn=_fn)

class Formula:

    """

    A Formula object for manipulating design matrices in regression models,
    essentially consisting of a list of Term instances.

    The object supports addition and multiplication which correspond
    to concatenation and pairwise multiplication, respectively,
    of the columns of the two Formulas.
    """
    
    def _terms_changed(self):
        self._names = self.names()
        self._termnames = self.termnames()

    def __init__(self, terms):
        """
        Create a Formula from either:

        i) a Formula object
        ii) a sequence of Term instances
        iii) one Term

        """

        if isinstance(terms, Formula):
            self.terms = copy.copy(list(terms.terms))
        elif type(terms) is types.ListType:
            self.terms = terms
        elif isinstance(terms, Term):
            self.terms = [terms]
        else: 
            raise ValueError

        self._terms_changed()

    def __str__(self):
        """
        String representation of list of termnames of a Formula.
        """
        value = []
        for term in self.terms:
            value += [term.termname]
        return '<formula: %s>' % string.join(value, ' + ')  

    def __call__(self, namespace=None, n=-1, **extra):
        """
        Create (transpose) of the design matrix of the Formula within
        namespace. Extra arguments are passed to each Term instance. If
        the Formula just contains an intercept, then the keyword
        argument 'n' indicates the number of rows (observations).
        """
        
        if namespace is None:
            namespace = globals()
        allvals = []
        intercept = False
        for term in self.terms:
            val = term(namespace=namespace, **extra)
            if val.shape == ():
                intercept = True
            elif val.ndim == 1:
                val.shape = (1, val.shape[0])
            allvals.append(val)

        if not intercept:
            allvals = N.concatenate(allvals)
        else:
            if allvals != []:
                n = allvals.shape[1]
                allvals = N.concatenate([N.ones((1,n), N.Float), allvals])
            elif n <= 1:
                raise ValueError, 'with no columns in model, keyword n argument needed for intercept'

        return allvals
    
    def hasterm(self, term):
        """
        Determine whether a given term is in a formula.
        """

        if not isinstance(term, Formula):
            return term.termname in self.termnames()
        elif len(term.terms) == 1:
            term = term.terms[0]
            return term.termname in self.termnames()
        else:
            raise ValueError, 'more than one term passed to hasterm'
        
    def termcolumns(self, term, dict=False):
        """
        Return a list of the indices of all columns associated
        to a given term.
        """

        if self.hasterm(term):
            names = term.names()
            value = {}
            for name in names:
                value[name] = self._names.index(name)
        else:
            raise ValueError, 'term not in formula'
        if dict:
            return value
        else:
            return value.values()

    def names(self):
        """
        Return a list of the names in the Formula. The order of the
        names corresponds to the order of the columns when self
        is evaluated.
        """

        allnames = []
        for term in self.terms:
            allnames += term.names()
        return allnames

    def termnames(self):
        """
        Return a list of the term names in the formula. These
        are the names of each Term instance in self.
        """

        names = []
        for term in self.terms:
            names += [term.termname]
        return names

    def design(self, namespace=None, **keywords):
        """
        transpose(self(namespace=namespace, **keywords))
        """
        return N.transpose(self(namespace=namespace, **keywords))

    def __mul__(self, other, nested=False):
        """
        This returns a Formula whose columns are the pairwise
        product of the columns of self and other.

        TO DO: check for nesting relationship. Should not be too difficult.
        """

        other = Formula(other)

        selftermnames = self.termnames()
        othertermnames = other.termnames()

        I = len(selftermnames)
        J = len(othertermnames)

        terms = []
        termnames = []

        for i in range(I):
            for j in range(J):
                termname = '%s*%s' % (str(selftermnames[i]), str(othertermnames[j]))
                pieces = termname.split('*')
                pieces.sort()
                termname = string.join(pieces, '*')
                termnames.append(termname)

                selfnames = self.terms[i].names()
                othernames = other.terms[j].names()

                if self.terms[i].name is 'intercept':
                    term = other.terms[j]
                elif other.terms[j].name is 'intercept':
                    term = self.terms[i]

                else:
                    names = []
                    for r in range(len(selfnames)):
                        for s in range(len(othernames)):
                            name = '%s*%s' % (str(selfnames[r]), str(othernames[s]))
                            pieces = name.split('*')
                            pieces.sort()
                            name = string.join(pieces, '*')
                            names.append(name)

                    def _fn(namespace=None, selfterm=self.terms[i], otherterm=other.terms[j], **extra):
                        value = []
                        selfval = N.array(selfterm(namespace=namespace, **extra))
                        if len(selfval.shape) == 1:
                            selfval.shape = (1, selfval.shape[0])
                        otherval = N.array(otherterm(namespace=namespace, **extra))
                        if len(otherval.shape) == 1:
                            otherval.shape = (1, otherval.shape[0])

                        for r in range(selfval.shape[0]):
                            for s in range(otherval.shape[0]):
                                value.append(selfval[r] * otherval[s])

                        return N.array(value)
                    term = Term(names, _fn=_fn, termname=termname)
                terms.append(term)

        return Formula(terms)
    
    def __add__(self, other):

        """
        Return a Formula whose columns are the
        concatenation of the columns of self and other.

        Terms in the formula are sorted alphabetically.
        """

        other = Formula(other)
        terms = self.terms + other.terms
        pieces = [(term.name, term) for term in terms]
        pieces.sort()
        terms = [piece[1] for piece in pieces]
        return Formula(terms)

    def __sub__(self, other):

        """
        Return a Formula with all terms in other removed from self.
        If other contains Term instances not in Formula, this
        function does not raise an exception.
        """

        other = Formula(other)
        terms = copy.copy(self.terms)

        for term in other.terms:
            for i in range(len(terms)):
                if terms[i].termname == term.termname:
                    terms.pop(i)
                    break 
        return Formula(terms)

def isnested(A, B, namespace=globals()):
    """
    Is factor B nested within factor A or vice versa: a very crude test
    which depends on the namespace.

    If they are nested, returns (True, F) where F is the finest
    level of the relationship. Otherwise, returns (False, None)

    """

    a = A(namespace, values=True)[0]
    b = B(namespace, values=True)[0]
    
    if len(a) != len(b):
        raise ValueError, 'A() and B() should be sequences of the same length'

    nA = len(sets.Set(a))
    nB = len(sets.Set(b))
    n = max(nA, nB)

    AB = [(a[i],b[i]) for i in range(len(a))]
    nAB = len(sets.Set(AB))

    if nAB == n:
        if nA > nB:
            F = A
        else:
            F = B
        return (True, F)
    else:
        return (False, None)

I = Term('intercept', _fn=lambda x: N.array(1))
