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
import cStringIO

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
# mod is module or Class
# parent is old parent
def toPage(obj, name, parent, mod):
    typ = type(obj)
    newparent = makenewparent(mod.__name__, parent)
    #print mod.__name__, parent, newparent, obj
    if typ is types.FunctionType:
        return BasePage(newparent, obj)
    elif typ is types.BuiltinFunctionType: 
        return BasePage(newparent, obj, canedit=0)
    elif typ is UfuncType:
        return BasePage(newparent, obj, name=name, canedit=0)
    elif typ is types.ClassType:
        return ClassPage(newparent, obj)
    elif typ is types.TypeType:
        return ClassPage(newparent, obj)
    elif typ is types.ModuleType:
        return ModulePage(newparent, obj, canedit=0)    
    elif typ is types.MethodType:
        return BasePage(newparent, obj)
    elif inspect.ismethoddescriptor(obj):
        return BasePage(newparent, obj, canedit=0)
    elif typ is types.InstanceType:
        return BasePage(newparent, obj, name=name, canedit=0)
    elif typ is types.DictType:
        return BasePage(newparent, obj, name=name, canedit=0)
    else:
        print "Err 3:", typ
        return None

class BasePage(rend.Page):
    docFactory = loaders.xmlfile("basepage.html")

    def __init__(self, parent, obj, **kw):
        self.parent = parent
        self.obj = obj
        self.name = kw.pop('name',None)
        if self.name is None:
            self.name = self.obj.__name__
        self.canedit = kw.pop('canedit',1)
        rend.Page.__init__(self, **kw)

    def render_title(self, context, data):
        return makenewparent(self.name, self.parent)

    def dochtml(self, str):
        return str

    def extrahtml(self, str):
        return str

    def getdoc_fromtree(self):
        return None

    def getextra_fromtree(self):
        return None

    def render_docstring(self, context, data):
        # First look for docstring in filetree
        doc = self.getdoc_fromtree()
        if doc is None:
            doc = cStringIO.StringIO()
            scipy.info(self.obj, output=doc)
            doc = doc.getvalue()
        return T.xml(self.dochtml(doc))

    def render_extra(self, context, data):
        # First look for docstring in filetree
        extra = self.getextra_fromtree() or inspect.getcomments(self.obj) or "Nothing"
        return T.xml(self.extrahtml(extra))

    def render_docedit(self, context, data):
        if self.canedit:
            return T.xml("<a href=_doc_edit_%s>edit</a>" % self.name)
        else:
            return ""

    def render_extraedit(self, context, data):
        return T.xml("<a href=_extra_edit_%s>edit</a>" % self.name)

class ModulePage(BasePage):

    addSlash = True
    
    def __init__(self, parent, module, **kw):
        BasePage.__init__(self, parent, module, **kw)

        # Avoid postponed import modules
        if hasattr(module, '_ppimport_module'):
            self.obj = module._ppimport_module
        
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

    def childFactory(self, context, name):
        if name not in self.all:
            print "Err 1: ", name, self.all
            return None
        child = getattr(self.obj,name,None)

        # special handle postponed import modules
        #  by causing it to import so that a new getattr
        #  returns the actual module
        if hasattr(child, '_ppimport_module'):
            child = getattr(self.obj,name,None)

        if child is None:
            print "Err 2: ", self.obj, name, getattr(self.obj, name, None)
            return None
                    
        return toPage(child, name, self.parent, self.obj)

    def render_title(self, context, data):
        return self.obj.__name__

    def dochtml(self, doc):
        # put links for names
        doc = subobj.sub("\\1<a href=\\2>\\2</a>\\3--\\4\\5",doc)
        doc = subobj2.sub("\\1<a href=\\2>\\2</a>\\3---\\4\\5",doc)
        return doc

class ClassPage(BasePage):
    addSlash = True   # I have subchildren
    def __init__(self, parent, func, **kw):
        # need to look for methods and add them to all
        BasePage.__init__(self, parent, func, **kw)
        self.all = []
        for meth in pydoc.allmethods(self.obj):
            if meth[0] != '_':
                self.all.append(meth)

    def childFactory(self, context, name):
        if name in self.all:
            child = getattr(self.obj,name,None)
            return toPage(child, name, self.parent, self.obj)
        print "Err 4: ", name, self.all
        return None

    def dochtml(self, str):
        # Now we want to make links in just the Methods: section
        start = str.find("Methods:")
        substr = str[start:]
        substr = subobj.sub("\\1<a href=\\2>\\2</a>\\3--\\4\\5",substr)
        return str[:start] + substr

_packages = ['scipy', 'scipy_base', 'weave', 'scipy_distutils', 'scipy_test']

class LiveDocs(rend.Page):
        
    addSlash = True
    docFactory = loaders.xmlfile("rootpage.html")

    child_styles = static.File('styles')
    child_images = static.File('images')
    
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

