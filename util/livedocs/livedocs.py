#At some-point we will have to figure out how to bring in builtin C-functions
#   to this process (big can of worms)

from nevow import rend
from nevow import inevow
from nevow import loaders
from nevow import static
import nevow.tags as T
import re
import sys
import types
import pydoc
import inspect
import cStringIO

subobj = re.compile("(\s+)(\w+)(\s+)--(\s+)([\w*])")
subobj2 = re.compile("(\s+)(\w+)(\s+)---(\s+)([\w*])")

import scipy
UfuncType = scipy.UfuncType

try:
    # Disable ppimport hooks.
    # scipy_base has been imported by scipy import
    sys.modules['scipy_base.ppimport'].disable()
except:
    pass

def isfortrantype(obj):
    return repr(type(obj))=="<type 'fortran'>"

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
    modname = getattr(mod,'__name__',name)
    newparent = makenewparent(modname, parent)
    #print mod.__name__, parent, newparent, obj
    if typ is types.FunctionType:
        return BasePage(newparent, obj, havesource=1)
    elif typ is types.BuiltinFunctionType:
        return BasePage(newparent, obj, canedit=0)
    elif typ is UfuncType:
        return BasePage(newparent, obj, name=name, canedit=0)
    elif typ is types.ClassType:
        return ClassPage(newparent, obj, havesource=1)
    elif typ is types.TypeType:
        return ClassPage(newparent, obj, havesource=1)
    elif typ is types.ModuleType:
        if re.match(r'This module \'\w+\' is auto-generated with f2py',
                    getattr(obj,'__doc__','') or ''):
            return F2pyModulePage(newparent, obj, canedit=0)
        return ModulePage(newparent, obj, canedit=0)
    elif typ is types.MethodType:
        return BasePage(newparent, obj, havesource=1)
    elif inspect.ismethoddescriptor(obj):
        return BasePage(newparent, obj, canedit=0)
    elif typ is types.InstanceType:
        return BasePage(newparent, obj, name=name, canedit=0)
    elif typ is types.DictType:
        return BasePage(newparent, obj, name=name, canedit=0)
    elif isfortrantype(obj):
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
        self.havesource = kw.pop('havesource',0)
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
            if type(self.obj) is types.ModuleType:
                reload(self.obj)
            scipy.info(self.obj, output=doc)
            doc = doc.getvalue()
        if self.havesource and hasattr(self.obj,'__name__'):
            if type(self.obj) in [types.TypeType,types.InstanceType,
                                  types.ClassType]:
                doc = re.sub(r'(\b'+self.obj.__name__+r'\b)',
                             r'<a href="__source__">\1</a>',doc,count=1)
            else:
                doc = re.sub(r'(\b'+self.obj.__name__+r'\b)',
                             r'<a href="\1/__source__">\1</a>',doc,count=1)

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

    def childFactory(self, ctx, name):
        if name == '__source__':
            return SourcePage(self.parent, self.obj, name=name)
        return None

class SourcePage(BasePage):
    docFactory = loaders.xmlfile("sourcepage.html")

    def render_docstring(self, context, data):
        from PySourceColor import str2html,lite
        doc = cStringIO.StringIO()
        scipy.source(self.obj, output=doc)            
        doc = doc.getvalue()
        doc = str2html(doc,colors=lite,form='code')
        return T.xml(doc)

class ModulePage(BasePage):

    addSlash = True
    
    def __init__(self, parent, module, **kw):
        BasePage.__init__(self, parent, module, **kw)

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
        if name not in self.all and not hasattr(self.obj,name):
            print "Err 1: ", name, self.all
            return None

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


f2pysubobj = re.compile(r"  ([\w,]+) = (\w+)(.*)")
f2pysubobj2 = re.compile(r"  (\w+)(\(.*)")
class F2pyModulePage(ModulePage):

    def __init__(self, parent, module, **kw):
        BasePage.__init__(self, parent, module, **kw)    
        self.all = [a for a in dir(module) if isfortrantype(getattr(module,a))]

    def dochtml(self, doc):
        doc = f2pysubobj.sub(r"  \1 = <a href=\2>\2</a>\3",doc)
        doc = f2pysubobj2.sub(r"  <a href=\1>\1</a>\2",doc)
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
        r = BasePage.childFactory(self,context, name)
        if r is not None: return r
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

