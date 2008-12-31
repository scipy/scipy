''' Object for reading from gzipped file-like object

Edited by Matthew Brett, with thanks, from
http://effbot.org/librarybook/zlib-example-4.py

The copyright and license for that code is:

Copyright 1995-2008 by Fredrik Lundh

By obtaining, using, and/or copying this software and/or its
associated documentation, you agree that you have read, understood,
and will comply with the following terms and conditions:

Permission to use, copy, modify, and distribute this software and its
associated documentation for any purpose and without fee is hereby
granted, provided that the above copyright notice appears in all
copies, and that both that copyright notice and this permission notice
appear in supporting documentation, and that the name of Secret Labs
AB or the author not be used in advertising or publicity pertaining to
distribution of the software without specific, written prior
permission.

SECRET LABS AB AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO
THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS. IN NO EVENT SHALL SECRET LABS AB OR THE AUTHOR BE LIABLE FOR
ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT
OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

'''

from zlib import decompressobj


class GzipInputStream(object):
    ''' Fileobject to wrap zlib compressed stream for reading

    >>> from StringIO import StringIO
    >>> from zlib import compress
    >>> S = 'A handy module for reading compressed streams'
    >>> F = StringIO(compress(S))
    >>> ZF = GzipInputStream(F)
    >>> ZF.read()
    'A handy module for reading compressed streams'
    >>> ZF.tell() == len(S)
    True
    >>> F = StringIO(compress(S))
    >>> ZF = GzipInputStream(F)
    >>> ZF.tell()
    0
    >>> ZF.read(6)
    'A hand'
    >>> ZF.tell()
    6
    '''

    blocksize = 16384 # 16K
    def __init__(self, fileobj, zipped_length=None):
        ''' Initialize GzipInputStream

        Parameters
        ----------
        fileobj : file-like object
            Object only need implement ``read`` method
        zipped_length : None or int, optional
            Compressed length of input stream in bytes
        '''
        self.fileobj = fileobj
        self.zipped_length=zipped_length
        self.exhausted = False
        self.unzipped_pos = 0
        self.data = ""
        self._unzipper = decompressobj()
        self._bytes_read = 0

    def __fill(self, bytes):
        ''' Fill self.data with at least *bytes* number of bytes
        If bytes == -1, continue until the end of the stream

        Returns ``None``
        '''
        if self.exhausted:
            return
        # read until we have enough bytes in the buffer
        read_to_end = bytes == -1
        n_to_fetch = self.blocksize
        while read_to_end or len(self.data) < bytes:
            if self.zipped_length: # do not read beyond specified length
                n_to_fetch = min(self.zipped_length-self._bytes_read,
                                 self.blocksize)
                if n_to_fetch == 0:
                    self.exhausted = True
                    break
            data = self.fileobj.read(n_to_fetch)
            self._bytes_read += len(data)
            if data:
                self.data += self._unzipper.decompress(data)
            if len(data) < n_to_fetch: # hit end of file
                self.data += self._unzipper.flush()
                self.exhausted = True
                break

    def seek(self, offset, whence=0):
        ''' Set position in uncompressed stream

        Parameters
        ----------
        offset : int
            byte offset relative to position given by *whence*
            offsets are in terms of uncompressed bytes from stream
        whence : {0,1}
            0 signifies *offset* is relative to beginning of file
            1 means *offset* is relative to current position

        Returns
        -------
        None
        '''
        if whence == 0:
            position = offset
        elif whence == 1:
            position = self.unzipped_pos + offset
        else:
            raise IOError, "Illegal argument"
        if position < self.unzipped_pos:
            raise IOError, "Cannot seek backwards"

        # skip forward, in blocks
        while position > self.unzipped_pos:
            if not self.read(min(position - self.unzipped_pos,
                                 self.blocksize)):
                break

    def tell(self):
        ''' Return current position in terms of uncompressed bytes '''
        return self.unzipped_pos

    def read(self, bytes = -1):
        ''' Read bytes from file

        Parameters
        ----------
        bytes : int, optional
            If *bytes* is a positive integer, read this many bytes
            from file. If *bytes* == -1 (the default), read all bytes
            to the end of file, where the end of the file is detected
            by running out of read data, or by the ``zipped_length``
            attribute.

        Returns
        -------
        data : string
            string containing read data

        '''
        self.__fill(bytes)
        if bytes == -1:
            data = self.data
            self.data = ""
        else:
            data = self.data[:bytes]
            self.data = self.data[bytes:]
        self.unzipped_pos += len(data)
        return data

    def readline(self):
        ''' Read text line from data

        Examples
        --------
        >>> from StringIO import StringIO
        >>> from zlib import compress
        >>> S = 'A handy module\\nfor reading\\ncompressed streams'
        >>> F = StringIO(compress(S))
        >>> ZF = GzipInputStream(F)
        >>> ZF.readline()
        'A handy module\\n'
        >>> ZF.readline()
        'for reading\\n'
        '''
        # make sure we have an entire line
        while not self.exhausted and "\n" not in self.data:
            self.__fill(len(self.data) + 512)
        i = self.data.find("\n") + 1
        if i <= 0:
            return self.read()
        return self.read(i)

    def readlines(self):
        ''' Read all data broken up into list of text lines
        >>> from StringIO import StringIO
        >>> from zlib import compress
        >>> S = 'A handy module\\nfor reading\\ncompressed streams'
        >>> F = StringIO(compress(S))
        >>> ZF = GzipInputStream(F)
        >>> ZF.readlines()
        ['A handy module\\n', 'for reading\\n', 'compressed streams']
        >>>
        '''
        lines = []
        while 1:
            s = self.readline()
            if not s:
                break
            lines.append(s)
        return lines
