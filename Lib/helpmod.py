import inspect
import types
import sys

# NOTE:  pydoc defines a help function which works simliarly to this
#  except it uses a pager to take over the screen.

# combine name and arguments and split to multiple lines of
#  width characters.  End lines on a comma and begin argument list
#  indented with the rest of the arguments.
def split_line(name, arguments, width):
    firstwidth = len(name)
    k = firstwidth
    newstr = name
    sepstr = ", "
    arglist = arguments.split(sepstr)
    for argument in arglist:
        if k == firstwidth:
            addstr = ""
        else:
            addstr = sepstr
        k = k + len(argument) + len(addstr)
        if k > width:
            k = firstwidth + 1 + len(argument)
            newstr = newstr + ",\n" + " "*(firstwidth+2) + argument
        else:
            newstr = newstr + addstr + argument
    return newstr

    
def help(object,maxwidth=76,output=sys.stdout):
    """Print help for a Python object.
    """
    if inspect.isfunction(object):
        name = object.func_name
        arguments = apply(inspect.formatargspec, inspect.getargspec(object))

        if len(name+arguments) > maxwidth:
            argstr = split_line(name, arguments, maxwidth)
        else:
            argstr = name + arguments

        print >> output, " " + argstr + "\n"
        print >> output, inspect.getdoc(object)

    elif inspect.isclass(object):  
        name = object.__name__
        if hasattr(object, '__init__'):
            arguments = apply(inspect.formatargspec, inspect.getargspec(object.__init__.im_func))        
            arglist = arguments.split(', ')
            if len(arglist) > 1:
                arglist[1] = "("+arglist[1]
                arguments = ", ".join(arglist[1:])
            else:
                arguments = "()"
        else:
            arguments = "()"

        if len(name+arguments) > maxwidth:
            argstr = split_line(name, arguments, maxwidth)
        else:
            argstr = name + arguments

        print >> output, " " + argstr + "\n"
        doc1 = inspect.getdoc(object) 
        if doc1 is None:
            if hasattr(object,'__init__'):
                print >> output, inspect.getdoc(object.__init__)
        else:
            print >> output, inspect.getdoc(object)

    elif type(object) is types.InstanceType: ## check for __call__ method
        print >> output, "Instance of class: ", object.__class__.__name__
        print >> output
        if hasattr(object, '__call__'):
            arguments = apply(inspect.formatargspec, inspect.getargspec(object.__call__.im_func))
            arglist = arguments.split(', ')
            if len(arglist) > 1:
                arglist[1] = "("+arglist[1]
                arguments = ", ".join(arglist[1:])
            else:
                arguments = "()"

            name = "<name>"
            if len(name+arguments) > maxwidth:
                argstr = split_line(name, arguments, maxwidth)
            else:
                argstr = name + arguments

            print >> output, " " + argstr + "\n"
            doc = inspect.getdoc(object.__call__)
            if doc is not None:
                print >> output, inspect.getdoc(object.__call__)
            print >> output, inspect.getdoc(object)
            
        else:
            print >> output, inspect.getdoc(object)
                                            
    elif inspect.ismethod(object):
        name = object.__name__
        arguments = apply(inspect.formatargspec, inspect.getargspec(object.im_func))
        arglist = arguments.split(', ')
        if len(arglist) > 1:
            arglist[1] = "("+arglist[1]
            arguments = ", ".join(arglist[1:])
        else:
            arguments = "()"

        if len(name+arguments) > maxwidth:
            argstr = split_line(name, arguments, maxwidth)
        else:
            argstr = name + arguments

        print >> output, " " + argstr + "\n"
        print >> output, inspect.getdoc(object)
                
    elif hasattr(object, '__doc__'):
        print >> output, inspect.getdoc(object)
        

def source(object, output=sys.stdout):
    if inspect.isroutine(object):
        print >> output,  inspect.getsource(object)
    else:
        print >> output,  "Not available for this object."
