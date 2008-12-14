import sys

def docformat(docstring, docdict=None):
    ''' Fill a function docstring from variables in dictionary

    Adapt the indent of the inserted docs

    Parameters
    ----------
    docstring : string
        docstring from function, possibly with dict formatting strings
    docdict : dict
        dictionary with keys that match the dict formatting strings
        and values that are docstring fragments to be inserted.  The
        indentation of the inserted docstrings is set to match the
        indentation of the ``docstring``.  The string values in the
        docdict are assumed to have no indent in the first line, and
        only indent relative to the first line for following lines.

    Returns
    -------
    outstring : string
        string with any formatted strings inserted
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
    ''' Return docstring decorator using docdict variable dictionary


    '''
    def decorate(f):
        f.__doc__ = docformat(f.__doc__, docdict)
        return f
    return decorate

