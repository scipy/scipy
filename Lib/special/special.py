from cephes import *
import types

class GeneralFunction:
    """gfunc = GeneralFunction(somefunc)

    This defines a generalized function gfunc which takes nested sequence
    objects or Numeric arrays as inputs and returns a Numeric array as output,
    evaluating the function over successive tuples of the input arrays like
    the python map function except it uses the broadcasting rules of ufuncs.
    """
    def __init__(self,pyfunc):
        if not callable(pyfunc) or type(pyfunc) is types.ClassType:
            raise TypeError, "Object is not a callable Python object."
        self.thefunc = pyfunc

    def __call__(self,*args):
        return arraymap(self.thefunc,args)





