from cephes import *
import types

class GeneralFunction:
    """
 GeneralFunction(somefunction)  Genearlized Function class.

  Description:
 
    Define a generalized function which takes nested sequence
    objects or Numeric arrays as inputs and returns a
    Numeric array as output, evaluating the function over successive
    tuples of the input arrays like the python map function except it uses
    the broadcasting rules of Numeric Python.

  Input:

    somefunction -- a Python function or method

  Example:

    def myfunc(a,b):
        if a > b:
            return a-b
        else
            return a+b

    gfunc = GeneralFunction(myfunc)

    >>> gfunc([1,2,3,4],2)
    array([3,4,1,2])

    """
    def __init__(self,pyfunc):
        if not callable(pyfunc) or type(pyfunc) is types.ClassType:
            raise TypeError, "Object is not a callable Python object."
        self.thefunc = pyfunc

    def __call__(self,*args):
        return arraymap(self.thefunc,args)





