#At some-point we will have to figure out how to bring in builtin C-functions
#   to this process (big can of worms)

from nevow import rend
from nevow import inevow
from nevow import loaders
from nevow import static
import nevow.tags as T
import re
import types
import pydoc
import inspect

subobj = re.compile("(\s+)(\w+)(\s+)--(\s+)([\w*])")
subobj2 = re.compile("(\s+)(\w+)(\s+)---(\s+)([\w*])")

import scipy
UfuncType = scipy.UfuncType

def makenewparent(name, parent):
    new = name
    if '.' in name:
        return name
    if parent != '':
        new = parent + '.' + new
    return new

# name is child
# mod is module
# parent is old parent
def toPage(obj, name, parent, mod):
    typ = type(obj)
    newparent = makenewparent(mod.__name__, parent)
    #print mod.__name__, parent, newparent, obj
    if typ is types.FunctionType:
        return FuncPage(newparent, obj)
    elif typ is types.BuiltinFunctionType: # maybe different later
        return FuncPage(newparent, obj)
    elif typ is UfuncType:
        return UFuncPage(name, newparent, obj)
    elif typ is types.ClassType:
        return ClassPage(newparent, obj)
    elif typ is types.TypeType:
        return ClassPage(newparent, obj)
    elif typ is types.ModuleType:
        return ModulePage(newparent, obj)    
    elif typ is types.MethodType:
        return FuncPage(newparent, obj)
    elif inspect.ismethoddescriptor(obj):
        return FuncPage(newparent, obj)
    elif typ is types.InstanceType:
        return UFuncPage(name, newparent, obj)
    elif typ is types.DictType:
        return UFuncPage(name, newparent, obj)
    else:
        print "Err 3:", typ
        return None

tohtml = pydoc.HTMLDoc()

class ModulePage(rend.Page):

    addSlash = True
    docFactory = loaders.xmlfile("module.html")
    
    def __init__(self, parent, module, **kw):
        self.parent = parent
        self.mod = module

        # Avoid postponed import modules
        if hasattr(module, '_ppimport_module'):
            self.mod = module._ppimport_module
        
        # If the module.__all__ variable is not defined
        #  then parse the doc-string and pull out all names
        #  where a name occurs on a single line with this pattern:
        #       xxxx --  description
        #  These names are then the children of the module        
        try:
            self.all = module.__all__
        except AttributeError:
            doc = module.__doc__
            if doc is None:
                self.all = []
            else:
                self.all = [x[1] for x in subobj.findall(doc)]
                self.all.extend([x[1] for x in subobj2.findall(doc)])
        if self.all is None:
            self.all = []
        rend.Page.__init__(self, **kw)

    def childFactory(self, context, name):
        if name not in self.all:
            print "Err 1: ", name, self.all
            return None
        child = getattr(self.mod,name,None)

        # special handle postponed import modules
        #  by causing it to import so that a new getattr
        #  returns the actual module
        if hasattr(child, '_ppimport_module'):
            child = getattr(self.mod,name,None)

        if child is None:
            print "Err 2: ", self.mod, name, getattr(self.mod, name, None)
            return None
                    
        return toPage(child, name, self.parent, self.mod)

    def render_title(self, context, data):
        return self.mod.__name__

    def render_docstring(self, context, data):
        doc = inspect.getdoc(self.mod) or "Nothing"
        # put links for names
        doc = subobj.sub("\\1<a href=\\2>\\2</a>\\3--\\4\\5",doc)
        doc = subobj2.sub("\\1<a href=\\2>\\2</a>\\3---\\4\\5",doc)
        return T.xml(doc)

    def render_extra(self, context, data):
        comments = inspect.getcomments(self.mod) or "Nothing"
        return T.xml(comments)

import cStringIO

    
class FuncPage(rend.Page):

    docFactory = loaders.xmlfile("function.html")    

    def __init__(self, parent, func, **kw):
        self.parent = parent
        self.func = func
        rend.Page.__init__(self, **kw)

    def render_title(self, context, data):
        return makenewparent(self.func.__name__, self.parent)

    def render_docstring(self, context, data):
        doc = cStringIO.StringIO()
        scipy.info(self.func, output=doc)
        return T.xml(doc.getvalue())

    def render_extra(self, context, data):
        comments = inspect.getcomments(self.func) or "Nothing"
        return T.xml(comments)

class UFuncPage(FuncPage):

    def __init__(self, name, parent, func, **kw):
        self.name = name
        FuncPage.__init__(self, parent, func, **kw)

    def render_title(self, context, data):
        return makenewparent(self.name, self.parent)


class ClassPage(FuncPage):
    addSlash = True  
    def __init__(self, parent, func, **kw):
        # need to look for methods and add them to all
        FuncPage.__init__(self, parent, func, **kw)
        self.all = []
        for meth in pydoc.allmethods(self.func):
            if meth[0] != '_':
                self.all.append(meth)

    def childFactory(self, context, name):
        if name in self.all:
            child = getattr(self.func,name,None)
            return toPage(child, name, self.parent, self.func)
        print "Err 4: ", name, self.all
        return None
    
    def render_docstring(self, context, data):
        doc = cStringIO.StringIO()
        scipy.info(self.func, output=doc)
        # Now we want to make links in just the Methods: section
        str = doc.getvalue()
        start = str.find("Methods:")
        substr = str[start:]
        substr = subobj.sub("\\1<a href=\\2>\\2</a>\\3--\\4\\5",substr)
        #substr = subobj2.sub("\\1<a href=\\2>\\2</a>\\3---\\4\\5",substr)
        return T.xml(str[:start]+substr)

_packages = ['scipy', 'scipy_base', 'weave', 'scipy_distutils', 'scipy_test']

class LiveDocs(rend.Page):
        
    addSlash = True
    docFactory = loaders.xmlfile("rootpage.html")

    child_styles = static.File('styles')
    
    def childFactory(self, context, name):
        if name not in _packages:
            return None

        # try to import the module.  If impossible return None
        
        try:
            exec("import %s" % name)
        except ImportError:
            return None 

        #  Return the page-renderer for a module.
        return ModulePage('', eval(name))
    
    def render_packageslist(self, context, data):
        packlist = ["<li><a href=%s>%s</a></li>" % (x,x) for x in _packages]
        res = T.xml('\n'.join(packlist))
        return res

