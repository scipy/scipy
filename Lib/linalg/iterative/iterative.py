
# Iterative methods using reverse-communication raw material
#   These methods solve
#   Ax = b  for x

#   where A must have A.matvec(x,*args) defined
#    or be a numeric array
#    or be a callable object matvec(x,*args)

def type1(self,x):
    return dot(self.obj, x)

def type2(self, x):
    return self.obj(x,*self.args)

class get_matvec:
    def __init__(self, obj, *args):
        self.obj = obj
        self.args = args
        if isinstance(obj, ArrayType):
            self.callfunc = new.instancemethod(type1,self,get_matvec)
            return
        try:
            meth = obj.matvec
        except AttributeError:
            meth = obj
        if not callable(meth):
            raise ValueError, "Object must be an array or a callable object"\
                  "or have a callable matvec attribute."

        self.obj = meth
        self.callfunc = new.instancemethod(type2, self, get_matvec)
        
    def __call__(self, x):
        return self.callfunc(self, x)
        
def bicg(A,b,):
    matvec = get_matvec(A)
    
    


