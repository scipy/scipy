"""
montecarlo.py:

Basic routines for Monte Carlo simulation.  Some of these routines need
to be spun out into separate compiled code in C or Fortran. They might
then belong under stats/ or utils/.  Meanwhile they're useful for small
problems such as the examples for the maximum entropy module.

Author: Ed Schofield
BSD license
"""

__author__ = "Ed Schofield"
__version__ = '2.0-alpha3'

from __future__ import division
import random, math, bisect, string
import scipy


class ordering(list):
    """An ordering or 'hash-indexed list' -- a bijection from a set of
    hashable objects to the positive integers, with constant-time lookups
    in both directions.  This allows efficient sampling from a sequence
    or dictionary containing arbitrary hashable objects.

    >>> a = ordering(['alpha'])
    >>> a.append('blah')
    >>> a
    ['alpha', 'blah']
    >>> a += ['cheese','dodo']
    >>> a.index('cheese')
    2
    >>> a[3]
    'dodo'

    The following should NOT add the word 'cheese' again:
    >>> a
    ['alpha', 'blah', 'cheese', 'dodo']
    >>> a.append('cheese')
    >>> a
    ['alpha', 'blah', 'cheese', 'dodo']
    
    >>> ordering(['blah','blah'])
    ['blah']
    
    >>> a.count('blah')
    1
    >>> a.count('elephant')
    0
    
    Eventually we'd want this:
      del a['blah']
    but I haven't implemented this yet. (Deleting an element might be
    O(n): perhaps all higher elements would need to be reindexed.)

    """
    
    def __add__(self, newlist):
        copy = self.copy()
        copy.extend(newlist)
        return copy
     
    def __iadd__(self, newlist):
        self.extend(newlist)
        return self
    
    def __imul__(self, i):
        raise TypeError, \
                "unsupported operand type(s) for *=: 'ordering' and 'int'"

    def __mul__(self, i):
        raise TypeError, \
                "unsupported operand type(s) for *: 'ordering' and 'int'"
    
    def __init__(self, initList=None):
        list.__init__(self)

        self._index = {}

        if initList is not None:
            self.extend(initList)
    
    def __delitem__(self, index):
        """Deletion of arbitrary elements is not yet supported.  It's not
        apparent to me how to implement this efficiently without
        fragmenting the indices. But perhaps an O(n) operation would be
        fine...
        """

        if index == len(self)-1 or index == -1:
            oldkey = self[index]
            del self._index[oldkey]
            list.__delitem__(self, index)
            
        else:
            print "called with arguments: " + str((self, index))
            raise NotImplementedError, "deletion currently only supported" \
                " for the final element of a doubly-ordered list"
    

    def extend(self, newlist):
        """ This extends an ordering with the values in the argument
        'newlist', which should be unique from each other and from the
        values currently in the ordering.
        
        Examples:
        
        >>> a = ordering()
        >>> a.append('blah')

        This should NOT extend the ordering:
        >>> a
        ['blah']
        >>> a.extend(['blah'])
        >>> a
        ['blah']
        
        #Traceback (most recent call last):
        #  ...
        #ValueError: value blah already in ordering
        
        This should extend it by the one unique element:
        
        >>> a.extend(['cheese','cheese'])
        >>> a
        ['blah', 'cheese']
        >>> a.index('cheese')
        1
        
        #Traceback (most recent call last):
        #  ...
        #ValueError: one or more duplicates in argument
        """

        for x in newlist:
            # We save a little overhead by executing the append() code directly:
            if not self._index.has_key(x):
                self._index[x] = len(self)     # the next available index
                list.append(self, x)
        
   
    def __setitem__(self, index, key):
        try:
            oldkey = self[index]
            del self._index[oldkey]
        except:
            pass
        self._index[key] = index
        list.__setitem__(self, index, key)


    def append(self, key):
        if not self._index.has_key(key):
            self._index[key] = len(self)     # the next available index
            list.append(self, key)
    
   
    def __repr__(self):

        if len(self) == 0:
            return "[]"
        strlist = []
        key = self[0]
        index = 0
        strlist += ['[' + repr(key)]
        for key in self[1:]:
            strlist += [', ' + repr(key)]
        strlist += [']']
        return string.join(strlist,'')
        

    def clear(self):
        self._index = {}
        list.clear(self)

    def key(self, n):
        """ Returns the nth key from the ordered dict. """
        try:
            key = self[n]
        except IndexError:
            raise IndexError('list index out of range')
        return key

    def insert(pos, value):
        raise Exception, 'ordering insertion not implemented'
    
    def count(self,value):
        """ A constant-time implementation of the underlying list.count().
        """
        # There are either 0 or 1 occurrences of 'value' in the ordering
        return 1 * self._index.has_key(value)
        
    def index(self, key):
        return self._index[key]

    def copy(self):
        dup = ordering()
        dup.extend(self)
        return dup

    # Would these functions be useful?
    #def items(self):
    #    return zip(self._index.values(), self._index.keys()).sort()

    #def keys(self):
    #    return list[:]

    def pop(self):
        key = self[-1]
        del self[-1]
        return key

    def reverse(self):
        raise NotImplementedError, 'reversal of an ordering not implemented'

    def sort(self):
        raise NotImplementedError, 'sort of an ordering not implemented'


def tablesampler(table):
    """ A generator that samples items from the distribution whose
    probability mass function is given in the sequence 'table'.  It
    yields the index of each generated element in the table and the log
    of the pmf at the point.

    The function expects the table to be a lookup table stored as a list
    or other sequence like this:
    >>> table = [10,150,20]
    
    The output will be something like:
    >>> sampler = tablesampler(table)
    
    >> [sampler.next()[0] for i in range(20)]
    [1, 1, 1, 0, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1]
    >> sampler.next()
    (1, -0.18232155679395459)
    >> sampler.next()
    (2, -2.1972245773362196)
    >> sampler.next()
    (1, -0.18232155679395459)

    etc., giving a sample from this distribution:
    x       0       1         2
    p(x)   10/180   150/180   20/180

    
    The function tablesampler() uses a binary search and takes time
    proportional to n + klog(n) to sample k points from a distribution on
    n items.

    A test:
    >>> s = [sampler.next()[0] for i in xrange(10000)]
    >>> freq = [s.count(i)*1.0/len(s) for i in range(3)]
    >>> tol = 0.2
    >>> s = sum(table)
    >>> [abs(table[i]*1.0/s - freq[i]) / (table[i]*1.0/s) < tol for i in range(3)]
    [True, True, True]
    """

    # First construct the cumulative distribution F(x) = Prob(X<=x).  This takes
    # linear time.  
    F = []
    sumP = 0.0
    try:
        total = sum(table)
        assert total != 0
    except:
        raise TypeError, "frequency table has invalid format or zero sum"
        
    for x in table:
        sumP += x * 1.0 / total
        F.append(sumP)
    # Now the cumulative dist F should have F[-1] == 1.0000 (perhaps with
    # some roundoff error)

    # Sample the next item from the dist.  This takes O(log(n)) time.
    R = random.Random()
    while True:
        # Use a binary search to look up the appropriate place in F
        w = bisect.bisect(F, R.random())
        yield w, math.log(table[w]*1.0/total)


def dictsampler(table):
    """ A generator that samples objects from the distribution given in
    the table (a dict or sparse vector) and yields the index in the
    dictionary and its log probability under the model.

    The function expects the table to be a lookup table stored as a
    sparse vector or dict like this:
    >>> table = {}
    >>> table['a'] = 10
    >>> table['b'] = 150
    >>> table['c'] = 20
    
    The output will be something like:
    >>> sampler = dictsampler(table)
    
    >> [sampler.next()[0] for i in range(10)]
    ['b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'c', 'b']
    >> sampler.next()
    ('a', -2.890371757896165)
    >> sampler.next()
    ('a', -2.890371757896165)
    >> sampler.next()
    ('b', -0.18232155679395459)

    etc., giving a sample from this distribution:
    x       'a'       'b'       'c'
    p(x)   10/180   150/180   20/180

    For a spmatrix object:
    >>> import spmatrix
    >>> table2 = spmatrix.ll_mat(1,3000)
    >>> table2[0,0] = 10
    >>> table2[0,1000] = 150
    >>> table2[0,2000] = 20
    
    The output will be something like:
    >>> sampler2 = dictsampler(table2)
    
    >> [sampler2.next()[0] for i in range(10)]
    [1000, 1000, 1000, 0, 1000, 1000, 1000, 1000, 1000, 2000]
    
    >> sampler2.next()
    (1000, -0.18232155679395459)
    >> sampler2.next()
    (2000, -2.1972245773362196)
    >> sampler2.next()
    (1000, -0.18232155679395459)

    etc., giving a sample from this distribution:
    x       0       1000      2000
    p(x)   10/180   150/180   20/180

    
    The fn uses a binary search and takes time proportional to 
    n + klog(n) for a table of size n.
    """
    try:
        m,n = table.shape
        # if this works we assume table is a 1-dim spmatrix or scipy array
        # and we need to iterate over it
        assert m == 1
        
        keys = ordering(table.keys())
        newtable = table.values()
        
        sampler =  tablesampler(newtable)

        while True:
            sample, logprob = sampler.next()
            yield keys[sample][1], logprob # the [1] selects the
                                           # column index, since table
                                           # will be a row vector with
                                           # indices (0,j)
    
    except AttributeError:
        # First map the keys to natural numbers
        keys = ordering(table)
        newtable = [table[key] for key in keys]
        sampler = tablesampler(newtable)

        # Sample the next item from the dist.  This takes O(log n) time.
        while True:
            sample, logprob = sampler.next()
            yield keys[sample], logprob


def independencesampler(psampler, log_q_dot):
    """ A generator of samples from q(x) = q_dot(x) / Z, where Z is some
    normalization constant (perhaps unknown), given the generator
    psampler that yields (x, log p(x)) from the distribution p (a
    probability mass or density function), which should be close to q for
    efficiency.  Uses the independent variant of the Metropolis-Hastings
    algorithm.
    
    See, for example, Robert and Casella, "Monte Carlo Statistical
    Methods", Springer 1999.
    """
    
    (x, log_p_x) = psampler.next()
    while True:
        # Candidate point c
        (c, log_p_c) = psampler.next()
        
        # Set acceptance probability
        # aprob=min(1, (q_dot(c)/q_dot(x))      \
        #              / (p_c/p_x))
        # Log version is:
        aprob = min(1.0, math.exp(log_q_dot(c) + log_p_x - log_q_dot(x) \
                - log_p_c))
        
        u=random.random()
        if (u<aprob):
            # Accept the new candidate
            x=c
            log_p_x = log_p_c
        # Otherwise just return the same observation x as before
        yield x


def empsuprejectionsampler(psampler, log_q_dot):
    """ A generator of samples from the density q(x) = q_dot(x) / Z,
    where Z is some normalization constant (perhaps unknown), given the
    generator psampler that yields (x, log p(x)) from the distribution p
    (a probability mass or density function), which should be close to q
    for efficiency.  Uses the Von Neumann rejection method with the
    emprical supremum as an approx to the true supremum.  Brian Caffo et
    al. (2001) show that the resulting sample converges very fast.
    
    See, for example, Robert and Casella, "Monte Carlo Statistical
    Methods", Springer (1999) for the general rejection sampler.
    """

    logT = 0.0

    # Accept X if U <= f(X) / Tg(X)
    # i.e.         if logU <= log q(X) - logT - log p(X)
    # Otherwise loop
    while True:
        # Candidate point x
        (x, log_p_x) = psampler.next()
        
        # Set acceptance probability
        # aprob = q_dot(x) / T p_x
        # Log version is:

        logT = max(logT, log_qdot(x) - log_p_x)    # updating the bound before
                                                   # making the new
                                                   # decision makes no
                                                   # difference; it will
                                                   # be accepted anyway.
        logU = math.log(random.random())
        if (logU < log_qdot(x) - logT - log_p_x):
            # Accept the new candidate
            yield x
        # Otherwise loop
        

def sample_wr(population, k):
    """Chooses k random elements (with replacement) from a population.
    (From the Python Cookbook).

    Example:
    
    >>> a = range(1000)
    
    # Sample 1000 values from 0 to 999 with replacement
    >>> sample = sample_wr(a, 1000)
    
    # The probability that all sample points are identical is negligibly
    # small:
    >>> sample.sort()
    >>> len(sample) > len(set(sample))
    True
    
    """
    n = len(population)
    _random, _int = random.random, int  # speed hack
    return [population[_int(_random() * n)] for i in xrange(k)]


   
def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()

