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

from StringIO import StringIO
from zlib import decompressobj


class ZlibInputStream(object):
    ''' Fileobject to wrap zlib compressed stream for reading

    >>> from StringIO import StringIO
    >>> from zlib import compress
    >>> s = 'A handy module for reading compressed streams'
    >>> cs = compress(s)
    >>> fobj = StringIO(cs)
    >>> zf = ZlibInputStream(fobj)
    >>> zf.read()
    'A handy module for reading compressed streams'
    >>> zf.tell() == len(s)
    True
    >>> fobj = StringIO(cs)
    >>> zf = ZlibInputStream(fobj)
    >>> zf.tell()
    0
    >>> zf.read(6)
    'A hand'
    >>> zf.tell()
    6

    You can change the blocksize to preserve memory.  Here it is
    ridiculously small for testing.
    
    >>> fobj = StringIO(cs)
    >>> zf = ZlibInputStream(fobj)
    >>> zf.default_blocksize = 3
    >>> zf.read()
    'A handy module for reading compressed streams'

    You can set the known length of the zipped stream.  This is
    normally when the stream is embedded in another stream, so there
    is no end-of-file signal when the zlib stream is finished.

    >>> fobj = StringIO(cs + 'padding')
    >>> zf = ZlibInputStream(fobj, len(cs))
    >>> zf.default_blocksize = 3
    >>> zf.read()
    'A handy module for reading compressed streams'

    >>> fobj = StringIO(cs + 'padding')
    >>> zf = ZlibInputStream(fobj, len(cs))
    >>> zf.default_blocksize = 3
    >>> zf.read(7)
    'A handy'
    >>> zf.tell()
    7
    >>> zf.read(7)
    ' module'
    >>> zf.tell()
    14
    '''

    default_blocksize = 16384 # 16K
    
    def __init__(self, fileobj, zipped_length=None):
        ''' Initialize ZlibInputStream

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
        self.data = StringIO()
        self._unzipper = decompressobj()
        # number of zlib compressed bytes read
        self._z_bytes_read = 0 
        self._blocksize_iterator = self._block_size_generator()
        
    def _block_size_generator(self):
        ''' Generator to give block sizes for reading

        >>> from StringIO import StringIO
        >>> from zlib import compress
        >>> s = 'A handy module\\nfor reading\\ncompressed streams'
        >>> zs = compress(s)
        >>> fobj = StringIO(zs)
        >>> # If length not set, always return default block size
        >>> zis = ZlibInputStream(fobj)
        >>> gen = zis._block_size_generator()
        >>> gen.next() == zis.default_blocksize
        True
        >>> gen.next() == zis.default_blocksize
        True
        >>> # if length is set, return min of remaining length,
        >>> # and default block size
        >>> zis = ZlibInputStream(fobj, len(zs))
        >>> gen = zis._block_size_generator()
        >>> gen.next() == len(zs)
        True
        >>> zis = ZlibInputStream(fobj, len(zs))
        >>> gen = zis._block_size_generator()
        >>> zis.default_blocksize = 5
        >>> gen.next()
        5
        '''
        if self.zipped_length:
            while True:
                yield min(self.zipped_length-self._z_bytes_read,
                          self.default_blocksize)
        else:
            while True:
                yield self.default_blocksize

    def __fill(self, bytes):
        ''' Fill self.data with at least *bytes* number of bytes
        If bytes == -1, continue until the end of the stream

        Parameters
        ----------
        bytes : integer
            Number of bytes to read from zlib stream
            If ``bytes==-1``, read the remaining bytes in stream

        Returns
        -------
        None
        '''
        if self.exhausted:
            return
        # read until we have enough bytes in the buffer
        read_to_end = bytes == -1
        s_data = StringIO(self.data.read())
        s_data.seek(0, 2) # seek to end
        while read_to_end or (bytes - s_data.pos) > 0:
            z_n_to_fetch = self._blocksize_iterator.next()
            if z_n_to_fetch == 0:
                self.exhausted = True
                break
            raw = self.fileobj.read(z_n_to_fetch)
            self._z_bytes_read += len(raw)
            if raw:
                s_data.write(self._unzipper.decompress(raw))
            if len(raw) < z_n_to_fetch: # hit end of file
                s_data.write(self._unzipper.flush())
                self.exhausted = True
                break
        s_data.seek(0)
        self.data = s_data
        
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
                                 self.default_blocksize)):
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
        if (bytes == -1 or
            (self.data.len-self.data.pos) < bytes):
            self.__fill(bytes)
        data = self.data.read(bytes)
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
        >>> zf = ZlibInputStream(F)
        >>> zf.readline()
        'A handy module\\n'
        >>> zf.readline()
        'for reading\\n'

        You can also set the block size
        (here very small for testing)
        
        >>> F = StringIO(compress(S))
        >>> zf = ZlibInputStream(F)
        >>> zf.default_blocksize = 5
        >>> zf.readline()
        'A handy module\\n'
        >>> zf.readline()
        'for reading\\n'
        '''
        # make sure we have an entire line
        data = self.data.read()
        blocks = [data]
        while not self.exhausted and "\n" not in data:
            # fill results in fresh data starting at 0
            data = self.read(512)
            blocks.append(data)
        data = ''.join(blocks)
        i = data.find("\n") + 1
        if i <= 0: # newline at end
            self.unzipped_pos += len(data)
            return data
        # new line not at end
        self.unzipped_pos += i
        self.data = StringIO(data[i:])
        return data[:i]

    def readlines(self):
        ''' Read all data broken up into list of text lines
        >>> from StringIO import StringIO
        >>> from zlib import compress
        >>> S = 'A handy module\\nfor reading\\ncompressed streams'
        >>> F = StringIO(compress(S))
        >>> zf = ZlibInputStream(F)
        >>> zf.readlines()
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


class TwoShotZlibInputStream(ZlibInputStream):
    ''' Class to do one small and a second final read
    of a zlib compressed input stream
    '''
    default_blocksize = 512 # bytes - first read

    def _block_size_generator(self):
        ''' Generator to give block sizes for reading
	'''
	if self.zipped_length:
            # do not read beyond specified length
            yield min(self.zipped_length, self.default_blocksize)
            yield self.zipped_length - self._z_bytes_read
            yield 0
	else:
            while True:
                yield self.default_blocksize


class OneShotZlibInputStream(ZlibInputStream):
    ''' One shot read, for testing '''
    
    def _block_size_generator(self):
        ''' Generator to give block sizes for reading
	'''
	if self.zipped_length:
            yield self.zipped_length
            yield 0
	else:
            while True:
                yield self.default_blocksize
	

class StubbyZlibInputStream(ZlibInputStream):
    ''' One short, then fairly long reads '''

    default_blocksize = 128 * 1024 # 128K
    first_blocksize = 512 # 512 bytes
    
    def _block_size_generator(self):
	if self.zipped_length:
            # do not read beyond specified length
            yield min(self.zipped_length, self.first_blocksize)
            while True:
                yield min(
                    self.zipped_length - self._z_bytes_read,
                    self.default_blocksize)
	else:
            yield self.first_blocksize
            while True:
                yield self.default_blocksize
