import numpy as N
from numpy.linalg import pinv
import utils

class ContrastResults:
    """
    Results from looking at a particular contrast of coefficients in
    a parametric model. The class does nothing, it is a container
    for the results from T and F contrasts.
    """

    def __init__(self, t=None, F=None, sd=None, effect=None, df_denom=None,
                 df_num=None):
        if F is not None:
            self.F = F 
            self.df_denom = df_denom
            self.df_num = df_num
        else:
            self.t = t
            self.sd = sd
            self.effect = effect
            self.df_denom = df_denom
            
    def __str__(self):
        if hasattr(self, 'F'):
            return '<F contrast: F=%s, df_denom=%d, df_num=%d>' % \
                   (`self.F`, self.df_denom, self.df_num)
        else:
            return '<T contrast: effect=%s, sd=%s, t=%s, df_denom=%d>' % \
                   (`self.effect`, `self.sd`, `self.t`, self.df_denom)
            

class Contrast:
    """
    This class is used to construct contrast matrices in regression models.
    They are specified by a (term, formula) pair.
    
    The term, T,  is a linear combination of columns of the design
    matrix D=formula(). The getmatrix method constructs
    a contrast matrix C so that

    colspan(dot(D, C)) = colspan(dot(D, dot(pinv(D), T))) 

    where pinv(D) is the generalized inverse of D. Further, the matrix

    Tnew = dot(C, D)
    
    is full rank. The rank attribute is the rank of

    dot(D, dot(pinv(D), T))

    In a regression model, the contrast tests that E(dot(Tnew, Y)) = 0
    for each column of Tnew.

    """

    def __init__(self, term, formula, name=''):
        self.term = term
        self.formula = formula
        if name is '':
            self.name = str(term)
        else:
            self.name = name

    def __str__(self):
        return '<contrast:%s>' % \
               `{'term':str(self.term), 'formula':str(self.formula)}`

    def getmatrix(self, evaldesign=True, **keywords):
        """
        Construct a contrast matrix C so that

        colspan(dot(D, C)) = colspan(dot(D, dot(pinv(D), T))) 

        where pinv(D) is the generalized inverse of D=self.D=self.formula().

        If the design, self.D is already set,
        then evaldesign can be set to False.
        """

        T = N.transpose(N.array(self.term(**keywords)))

        if T.ndim == 1:
            T.shape = (T.shape[0], 1)
        
        T = utils.clean0(T)

        if evaldesign:
            self.D = self.formula.design(**keywords)
            self.pinv = pinv(self.D)

        self.matrix = contrastfromcols(T, self.D)
        try:
            self.rank = self.matrix.shape[1]
        except:
            self.rank = 1

def contrastfromcols(T, D, pseudo=None):
    """
    From an n x p design matrix D and a matrix T, tries
    to determine a p x q contrast matrix C which
    determines a contrast of full rank, i.e. the
    n x q matrix

    dot(transpose(C), pinv(D))

    is full rank.

    T must satisfy either T.shape[0] == n or T.shape[1] == p.

    Note that this always produces a meaningful contrast, not always
    with the intended properties because q is always non-zero unless
    T is identically 0. That is, it produces a contrast that spans
    the column space of T (after projection onto the column space of D).

    """

    n, p = D.shape

    if T.shape[0] != n and T.shape[1] != p:
        raise ValueError, 'shape of T and D mismatched'

    if pseudo is None:
        pseudo = pinv(D)

    if T.shape[0] == n:
        C = N.transpose(N.dot(pseudo, T))
    else:
        C = T

    Tp = N.dot(D, N.transpose(C))

    if utils.rank(Tp) != Tp.shape[1]:
        Tp = utils.fullrank(Tp)
        C = N.transpose(N.dot(pseudo, Tp))

    return N.squeeze(C)

