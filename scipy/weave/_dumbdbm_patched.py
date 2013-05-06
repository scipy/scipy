"""A dumb and slow but simple dbm clone.

For database spam, spam.dir contains the index (a text file),
spam.bak *may* contain a backup of the index (also a text file),
while spam.dat contains the data (a binary file).

XXX TO DO:

- seems to contain a bug when updating...

- reclaim free space (currently, space once occupied by deleted or expanded
items is never reused)

- support concurrent access (currently, if two processes take turns making
updates, they can mess up the index)

- support efficient access to large databases (currently, the whole index
is read when the database is opened, and some updates rewrite the whole index)

- support opening for read-only (flag = 'm')

"""
from __future__ import division, print_function, absolute_import

_os = __import__('os')

from scipy.lib.six.moves import builtins
from scipy.lib.six import string_types

_open = builtins.open

_BLOCKSIZE = 512

error = IOError             # For anydbm


class _Database(object):

    def __init__(self, file):
        self._dirfile = file + '.dir'
        self._datfile = file + '.dat'
        self._bakfile = file + '.bak'
        # Mod by Jack: create data file if needed
        try:
            f = _open(self._datfile, 'r')
        except IOError:
            f = _open(self._datfile, 'w')
        f.close()
        self._update()

    def _update(self):
        import string
        self._index = {}
        try:
            f = _open(self._dirfile)
        except IOError:
            pass
        else:
            while 1:
                line = string.rstrip(f.readline())
                if not line:
                    break
                key, (pos, siz) = eval(line)
                self._index[key] = (pos, siz)
            f.close()

    def _commit(self):
        try:
            _os.unlink(self._bakfile)
        except _os.error:
            pass
        try:
            _os.rename(self._dirfile, self._bakfile)
        except _os.error:
            pass
        f = _open(self._dirfile, 'w')
        for key, (pos, siz) in self._index.items():
            f.write("%s, (%s, %s)\n" % (repr(key), repr(pos), repr(siz)))
        f.close()

    def __getitem__(self, key):
        pos, siz = self._index[key]  # may raise KeyError
        f = _open(self._datfile, 'rb')
        f.seek(pos)
        dat = f.read(siz)
        f.close()
        return dat

    def __contains__(self, key):
        return key in self._index

    def _addval(self, val):
        f = _open(self._datfile, 'rb+')
        f.seek(0, 2)
        pos = f.tell()
## Does not work under MW compiler
##      pos = ((pos + _BLOCKSIZE - 1) // _BLOCKSIZE) * _BLOCKSIZE
##      f.seek(pos)
        npos = ((pos + _BLOCKSIZE - 1) // _BLOCKSIZE) * _BLOCKSIZE
        f.write('\0'*(npos-pos))
        pos = npos

        f.write(val)
        f.close()
        return (pos, len(val))

    def _setval(self, pos, val):
        f = _open(self._datfile, 'rb+')
        f.seek(pos)
        f.write(val)
        f.close()
        return (pos, len(val))

    def _addkey(self, key, pos_and_siz):
        (pos, siz) = pos_and_siz
        self._index[key] = (pos, siz)
        f = _open(self._dirfile, 'a')
        f.write("%s, (%s, %s)\n" % (repr(key), repr(pos), repr(siz)))
        f.close()

    def __setitem__(self, key, val):
        if not isinstance(key, string_types) or not isinstance(val, string_types):
            raise TypeError("keys and values must be strings")
        if key not in self._index:
            (pos, siz) = self._addval(val)
            self._addkey(key, (pos, siz))
        else:
            pos, siz = self._index[key]
            oldblocks = (siz + _BLOCKSIZE - 1) // _BLOCKSIZE
            newblocks = (len(val) + _BLOCKSIZE - 1) // _BLOCKSIZE
            if newblocks <= oldblocks:
                pos, siz = self._setval(pos, val)
                self._index[key] = pos, siz
            else:
                pos, siz = self._addval(val)
                self._index[key] = pos, siz
            self._addkey(key, (pos, siz))

    def __delitem__(self, key):
        del self._index[key]
        self._commit()

    def keys(self):
        return list(self._index.keys())

    def has_key(self, key):
        return key in self._index

    def __len__(self):
        return len(self._index)

    def close(self):
        self._index = None
        self._datfile = self._dirfile = self._bakfile = None


def open(file, flag=None, mode=None):
    # flag, mode arguments are currently ignored
    return _Database(file)
