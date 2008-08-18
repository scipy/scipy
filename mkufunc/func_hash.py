"""

Author: Ilan Schnell
"""
import md5
import re
from types import CodeType, FunctionType

def md5sum(s):
    return md5.md5(s).hexdigest()


pat_hex = re.compile(r'0x[0-9a-f]{4,}', re.I)


def func_hash(f, salt=None):
    """ Return a MD5 hash for a function or code object as string.
    """
    if type(f) == FunctionType:
        co = f.func_code
    elif type(f) == CodeType:
        co = f
    else:
        raise TypeError("Object %r is not function or code object.")
    
    res = []
    for name in dir(co):
        if not name.startswith('co_'):
            continue
        if name == 'co_consts':
            for c in getattr(co, name):
                if type(c) == CodeType or \
                   type(c) == FunctionType:
                    res.append(func_hash(c))
                else:
                    res.append(repr(c))
        else:
            res.append(repr(getattr(co, name)))
            
    return md5sum(''.join(res) + repr(salt))
