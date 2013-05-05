from __future__ import division, print_function, absolute_import

from shelve import Shelf
try:
    import zlib
except ImportError:
    # Some python installations don't have zlib.
    pass

import pickle


class DbfilenameShelf(Shelf):
    """Shelf implementation using the "anydbm" generic dbm interface.

    This is initialized with the filename for the dbm database.
    See the module's __doc__ string for an overview of the interface.
    """

    def __init__(self, filename, flag='c'):
        from . import _dumbdbm_patched
        Shelf.__init__(self, _dumbdbm_patched.open(filename, flag))

    def __getitem__(self, key):
        compressed = self.dict[key]
        try:
            r = zlib.decompress(compressed)
        except zlib.error:
            r = compressed
        except NameError:
            r = compressed

        return pickle.loads(r)

    def __setitem__(self, key, value):
        s = pickle.dumps(value,1)
        try:
            self.dict[key] = zlib.compress(s)
        except NameError:
            #zlib doesn't exist, leave it uncompressed.
            self.dict[key] = s


def open(filename, flag='c'):
    """Open a persistent dictionary for reading and writing.

    Argument is the filename for the dbm database.
    See the module's __doc__ string for an overview of the interface.
    """

    return DbfilenameShelf(filename, flag)
