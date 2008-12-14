import sys

def docformat(docstring, docdict=None):
    ''' Fill a function docstring from variables in dict

    Adapt the indent of the inserted docs
    '''
    if not docstring:
        return docstring
    if docdict is None:
        docdict = {}
    lines = docstring.expandtabs().splitlines()
    # Find the minimum indent of the main docstring, after last line
    indentno = sys.maxint
    for line in lines[1:]:
        stripped = line.lstrip()
        if stripped:
            indentno = min(indentno, len(line) - len(stripped))
    indent = ' ' * indentno
    # Insert this indent to dictionary docstrings
    indented = {}
    for name, dstr in docdict.items():
        lines = dstr.expandtabs().splitlines()
        newlines = [lines[0]]
        for line in lines[1:]:
            newlines.append(indent+line)
        indented[name] = '\n'.join(newlines)
    return docstring % indented


def filldoc(docdict):
    ''' Return docstring decorator using docdict variable dictionary'''
    def decorate(f):
        f.__doc__ = docformat(f.__doc__, docdict)
        return f
    return decorate

