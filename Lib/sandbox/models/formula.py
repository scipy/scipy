import copy
import types
import numpy as N

default_namespace = {}

class term(object):

    """
    This class is very simple: it is just a named term in a model formula.

    It is also callable: by default it namespace[self.name], where namespace
    defaults to formula.default_namespace. 
    When called in an instance of formula, 
    the namespace used is that formula's namespace.
    
    """

    def __pow__(self, power):
        """
        Raise the quantitative term's values to an integer power, i.e.
        polynomial.
        """

        try:
            power = float(power)
        except:
            raise ValueError, 'expecting a float'

        if power == int(power):
            name = '%s^%d' % (self.name, int(power))
        else:
            name = '%s^%0.2f' % (self.name, power)

        value = quantitative(name, func=self, transform=lambda x: N.power(x, power))
        value.power = power
        value.namespace = self.namespace
        return value

    def __init__(self, name, func=None, termname=None):
        
        self.name = name
        self.__namespace = None
        if termname is None:
            self.termname = name
        else:
            self.termname = termname

        if type(self.termname) is not types.StringType:
            raise ValueError, 'expecting a string for termname'
        if func:
            self.func = func

    # Namespace in which self.name will be looked up in, if needed

    def _get_namespace(self): 
	if isinstance(self.__namespace, N.ndarray): 
	    return self.__namespace 
	else: return self.__namespace or default_namespace

    def _set_namespace(self, value):  self.__namespace = value
    def _del_namespace(self): del self.__namespace
    namespace = property(_get_namespace, _set_namespace, _del_namespace)

    def __str__(self):
        """
        '<term: %s>' % self.termname
        """
        return '<term: %s>' % self.termname

    def __add__(self, other):
        """
        formula(self) + formula(other)
        """
        other = formula(other, namespace=self.namespace)
        f = other + self
        f.namespace = self.namespace
        return f

    def __mul__(self, other):
        """
        formula(self) * formula(other)
        """

        if other.name is 'intercept':
            f = formula(self, namespace=self.namespace)
        elif self.name is 'intercept':
            f = formula(other, namespace=other.namespace)
        else:
            other = formula(other, namespace=self.namespace)
            f = other * self
        f.namespace = self.namespace
        return f

    def names(self):
        """
        Return the names of the columns in design associated to the terms,
        i.e. len(self.names()) = self().shape[0].
        """
        if type(self.name) is types.StringType:
            return [self.name]
        else:
            return list(self.name)

    def __call__(self, *args, **kw):
        """
        Return the columns associated to self in a design matrix.
	If the term has no 'func' attribute, it returns
        
	self.namespace[self.termname]

	else, it returns
	
	self.func(*args, **kw)

        """
        
        if not hasattr(self, 'func'):
            val = self.namespace[self.termname]
        else:
            val = self.func
        if callable(val):
            if hasattr(val, "namespace"):
                val.namespace = self.namespace
            val = val(*args, **kw)

        val = N.asarray(val)
        return N.squeeze(val)

class factor(term):

    """
    A categorical factor.
    """



    def __init__(self, termname, keys, ordinal=False):
        """
        factor is initialized with keys, representing all valid
        levels of the factor.
        """
        
        self.keys = list(set(keys))
        self.keys.sort()
        self._name = termname
        self.termname = termname
        self.ordinal = ordinal

        if self.ordinal:
            name = self.name
        else:
            name = ['(%s==%s)' % (self.termname, str(key)) for key in self.keys]

        term.__init__(self, name, termname=self.termname, func=self.get_columns)

    def get_columns(self, *args, **kw):
        """
	Calling function for factor instance.
	"""

        v = self.namespace[self._name]
        while True:
            if callable(v): 
                if hasattr(v, "namespace"):
                    v.namespace = self.namespace
                v = v(*args, **kw)
            else: break 

        if self.ordinal:
            col = [float(self.keys.index(v[i])) for i in range(len(self.keys))]
            return N.array(col)

        else:
            n = len(v)
            value = []
            for key in self.keys:
                col = [float((v[i] == key)) for i in range(n)]
                value.append(col)
            return N.array(value)

    def values(self, *args, **kw):
        """
	Return the keys of the factor, rather than the columns of the design 
	matrix.
	"""

        del(self.func)
        val = self(*args, **kw)
        self.func = self.get_columns
        return val

    def verify(self, values):
        """
        Verify that all values correspond to valid keys in self.
        """
        s = set(values)
        if not s.issubset(self.keys):
            raise ValueError, 'unknown keys in values'

    def __add__(self, other):
        """
        formula(self) + formula(other)

        When adding \'intercept\' to a factor, this just returns 

	formula(self, namespace=self.namespace)

        """
        
        if other.name is 'intercept':
            return formula(self, namespace=self.namespace)
        else:
            return term.__add__(self, other)

    def main_effect(self, reference=None):
        """
        Return the 'main effect' columns of a factor, choosing
        a reference column number to remove.
        """

        if reference is None:
            reference = 0

        names = self.names()

        def maineffect_func(value, reference=reference):
            rvalue = []
            keep = range(value.shape[0])
            keep.pop(reference)
            for i in range(len(keep)):
                rvalue.append(value[keep[i]] - value[reference])
            return N.array(rvalue)
        
        keep = range(len(self.names()))
        keep.pop(reference)
        __names = self.names()
        _names = ['%s-%s' % (__names[keep[i]], __names[reference]) for i in range(len(keep))]
        value = quantitative(_names, func=self, 
		     termname='%s:maineffect' % self.termname,
		     transform=maineffect_func)
        value.namespace = self.namespace
        return value

class quantitative(term):

    """
    A subclass of term that can be used to apply point transformations
    of another term, i.e. to take powers:

    >>> import numpy as N
    >>> from scipy.sandbox.models import formula
    >>> X = N.linspace(0,10,101)
    >>> x = formula.term('X')
    >>> x.namespace={'X':X}
    >>> x2 = x**2
    >>> print N.allclose(x()**2, x2())
    True
    >>> x3 = formula.quantitative('x2', func=x, transform=lambda x: x**2)
    >>> x3.namespace = x.namespace
    >>> print N.allclose(x()**2, x3())
    True

    """

    def __init__(self, name, func=None, termname=None, transform=lambda x: x):
        self.transform = transform
        term.__init__(self, name, func=func, termname=termname)

    def __call__(self, *args, **kw):
        """
	A quantitative is just like term, except there is an additional
	transformation: self.transform.
	"""
        return self.transform(term.__call__(self, *args, **kw))

class formula(object):

    """

    A formula object for manipulating design matrices in regression models,
    essentially consisting of a list of term instances.

    The object supports addition and multiplication which correspond
    to concatenation and pairwise multiplication, respectively,
    of the columns of the two formulas.
    """
    
    def _get_namespace(self): 
	if isinstance(self.__namespace, N.ndarray): 
	    return self.__namespace 
	else: return self.__namespace or default_namespace

    def _set_namespace(self, value):  self.__namespace = value
    def _del_namespace(self): del self.__namespace
    namespace = property(_get_namespace, _set_namespace, _del_namespace)

    def _terms_changed(self):
        self._names = self.names()
        self._termnames = self.termnames()

    def __init__(self, termlist, namespace=default_namespace):
        """
        Create a formula from either:

        i) a formula object
        ii) a sequence of term instances
        iii) one term

        """

        self.__namespace = namespace
        if isinstance(termlist, formula):
            self.terms = copy.copy(list(termlist.terms))
        elif type(termlist) is types.ListType:
            self.terms = termlist
        elif isinstance(termlist, term):
            self.terms = [termlist]
        else: 
            raise ValueError

        self._terms_changed()

    def __str__(self):
        """
        String representation of list of termnames of a formula.
        """
        value = []
        for term in self.terms:
            value += [term.termname]
        return '<formula: %s>' % ' + '.join(value)

    def __call__(self, *args, **kw):

        """
        Create (transpose) of the design matrix of the formula within
        namespace. Extra arguments are passed to each term instance. If
        the formula just contains an intercept, then the keyword
        argument 'n' indicates the number of rows (observations).
        """

        allvals = []
        intercept = False
        iindex = 0
        for t in self.terms:

            t.namespace = self.namespace
            val = t(*args, **kw)

            isintercept = False
            if hasattr(t, "termname"):
                if t.termname == 'intercept': 
                    intercept = True
                    isintercept = True
                    interceptindex = iindex
                    allvals.append(None)

            if val.ndim == 1 and not isintercept:
                val.shape = (1, val.shape[0])
                allvals.append(val)
            elif not isintercept:
                allvals.append(val)
            iindex += 1

        if not intercept:
            try:
                allvals = N.concatenate(allvals)
            except:
                pass
        else:
            if allvals != []:
                if interceptindex > 0:
                    n = allvals[0].shape[1]
                else:
                    n = allvals[1].shape[1]
                allvals[interceptindex] = N.ones((1,n), N.float64) 
                allvals = N.concatenate(allvals)
            elif nrow <= 1: # FIXME: nrow is undefined here
                raise ValueError, 'with only intercept in formula, keyword \'nrow\' argument needed'
            else:
                allvals = I(nrow=nrow) # ... and here
                allvals.shape = (1,) + allvals.shape
        return allvals
    
    def hasterm(self, query_term):
        """
        Determine whether a given term is in a formula.
        """

        if not isinstance(query_term, formula):
	    if type(query_term) == type("name"):
		try: query = self[query_term]
		except: return False
	    elif isinstance(query_term, term):
		return query_term.termname in self.termnames()
        elif len(query_term.terms) == 1:
            query_term = query_term.terms[0]
            return query_term.termname in self.termnames()
        else:
            raise ValueError, 'more than one term passed to hasterm'
        
    def __getitem__(self, name):
        t = self.termnames()
        if name in t:
            return self.terms[t.index(name)]
        else:
            raise KeyError, 'formula has no such term: %s' % repr(name)

    def termcolumns(self, query_term, dict=False):
        """
        Return a list of the indices of all columns associated
        to a given term.
        """

        if self.hasterm(query_term):
            names = query_term.names()
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
        Return a list of the names in the formula. The order of the
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
        are the names of each term instance in self.
        """

        names = []
        for term in self.terms:
            names += [term.termname]
        return names

    def design(self, *args, **kw):
        """
        transpose(self(*args, **kw))
        """
        return self(*args, **kw).T

    def __mul__(self, other, nested=False):
        """
        This returns a formula whose columns are the pairwise
        product of the columns of self and other.

        TO DO: check for nesting relationship. Should not be too difficult.
        """

        other = formula(other, namespace=self.namespace)

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
                termname = '*'.join(pieces)
                termnames.append(termname)

                selfnames = self.terms[i].names()
                othernames = other.terms[j].names()

                if self.terms[i].name is 'intercept':
                    _term = other.terms[j]
                    _term.namespace = other.namespace

                elif other.terms[j].name is 'intercept':
                    _term = self.terms[i]
                    _term.namespace = self.namespace
                else:
                    names = []

                    d1 = len(selfnames) 
                    d2 = len(othernames)

                    for r in range(d1):
                        for s in range(d2):
                            name = '%s*%s' % (str(selfnames[r]), str(othernames[s]))
                            pieces = name.split('*')
                            pieces.sort()
                            name = '*'.join(pieces)
                            names.append(name)

                    def product_func(value, d1=d1, d2=d2):

                        out = []
                        for r in range(d1):
                            for s in range(d2):
                                out.append(value[r] * value[d1+s])
                        return N.array(out)

                    sumterms = self + other
                    sumterms.terms = [self, other] # enforce the order we want
                    sumterms.namespace = self.namespace

                    _term = quantitative(names, func=sumterms, termname=termname,
					 transform=product_func)
                    _term.namespace = self.namespace


                terms.append(_term)

        return formula(terms, namespace=self.namespace)
    
    def __add__(self, other):

        """
        Return a formula whose columns are the
        concatenation of the columns of self and other.

        terms in the formula are sorted alphabetically.
        """

        other = formula(other, namespace=self.namespace)
        terms = self.terms + other.terms
        pieces = [(term.name, term) for term in terms]
        pieces.sort()
        terms = [piece[1] for piece in pieces]
        return formula(terms, namespace=self.namespace)

    def __sub__(self, other):

        """
        Return a formula with all terms in other removed from self.
        If other contains term instances not in formula, this
        function does not raise an exception.
        """

        other = formula(other, namespace=self.namespace)
        terms = copy.copy(self.terms)

        for term in other.terms:
            for i in range(len(terms)):
                if terms[i].termname == term.termname:
                    terms.pop(i)
                    break 
        return formula(terms, namespace=self.namespace)

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

    nA = len(set(a))
    nB = len(set(b))
    n = max(nA, nB)

    AB = [(a[i],b[i]) for i in range(len(a))]
    nAB = len(set(AB))

    if nAB == n:
        if nA > nB:
            F = A
        else:
            F = B
        return (True, F)
    else:
        return (False, None)

def _intercept_fn(nrow=1, **extra):
    return N.ones((1,nrow))
I = term('intercept', func=_intercept_fn)
I.__doc__ = """
Intercept term in a formula. If intercept is the
only term in the formula, then a keywords argument
\'nrow\' is needed.

>>> from formula import *
>>> I()
1
>>> I(nrow=5)
array([1, 1, 1, 1, 1])

>>> f=formula(I)
>>> f(nrow=5)
array([1, 1, 1, 1, 1])

"""
