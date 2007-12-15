import cPickle

# snip on----- DELETE after numpy.deprecate_with_doc is available
import numpy
numpy.deprecate_with_doc = lambda doc: (lambda func: func)
# snip off---- DELETE after numpy.deprecate_with_doc is available

from numpy deprecate_with_doc

@deprecate_with_doc('')
def objsave(file, allglobals, *args):
    """Pickle the part of a dictionary containing the argument list
    into file string.

    Syntax:  objsave(file, globals(), obj1, obj2, ... )
    """
    fid = open(file,'w')
    savedict = {}
    for key in allglobals.keys():
        inarglist = 0
        for obj in args:
            if allglobals[key] is obj:
                inarglist = 1
                break
        if inarglist:
            savedict[key] = obj
    cPickle.dump(savedict,fid,1)
    fid.close()

def objload(file, allglobals):
    """Load a previously pickled dictionary and insert into given dictionary.

    Syntax:  objload(file, globals())
    """
    fid = open(file,'r')
    savedict = cPickle.load(fid)
    allglobals.update(savedict)
    fid.close()
