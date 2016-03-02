"""
Module to read / write wav files using numpy arrays

Functions
---------
`read`: Return the sample rate (in samples/sec) and data from a WAV file.

`write`: Write a numpy array as a WAV file.

"""
from __future__ import division, print_function, absolute_import

import sys
import numpy
import struct
import warnings


__all__ = [
    'WavFileWarning',
    'read',
    'write'
]


class WavFileWarning(UserWarning):
    pass


WAVE_FORMAT_PCM = 0x0001
WAVE_FORMAT_IEEE_FLOAT = 0x0003
WAVE_FORMAT_EXTENSIBLE = 0xfffe
KNOWN_WAVE_FORMATS = (WAVE_FORMAT_PCM, WAVE_FORMAT_IEEE_FLOAT)

# assumes file pointer is immediately
#  after the 'fmt ' id


def _read_fmt_chunk(fid, is_big_endian):
    if is_big_endian:
        fmt = '>'
    else:
        fmt = '<'

    size = res = struct.unpack(fmt+'I', fid.read(4))[0]
    bytes_read = 0

    if size < 16:
        raise ValueError("Binary structure of wave file is not compliant")

    res = struct.unpack(fmt+'HHIIHH', fid.read(16))
    bytes_read += 16
    comp, noc, rate, sbytes, ba, bits = res

    if comp == WAVE_FORMAT_EXTENSIBLE and size >= (16+2):
        ext_chunk_size = struct.unpack(fmt+'H', fid.read(2))[0]
        bytes_read += 2
        if ext_chunk_size >= 22:
            extensible_chunk_data = fid.read(22)
            bytes_read += 22
            raw_guid = extensible_chunk_data[2+4:2+4+16]
            # GUID template {XXXXXXXX-0000-0010-8000-00AA00389B71} (RFC-2361)
            # MS GUID byte order: first three groups are native byte order, rest is Big Endian
            if is_big_endian:
                tail = b'\x00\x00\x00\x10\x80\x00\x00\xAA\x00\x38\x9B\x71'
            else:
                tail = b'\x00\x00\x10\x00\x80\x00\x00\xAA\x00\x38\x9B\x71'
            if raw_guid.endswith(tail):
                comp = struct.unpack(fmt+'I', raw_guid[:4])[0]
        else:
            raise ValueError("Binary structure of wave file is not compliant")

    if comp not in KNOWN_WAVE_FORMATS:
        raise ValueError("Unknown wave file format")
        
    # move file pointer to next chunk
    if size > (bytes_read):
        fid.read(size - bytes_read)

    return size, comp, noc, rate, sbytes, ba, bits


# assumes file pointer is immediately
#   after the 'data' id
def _read_data_chunk(fid, comp, noc, bits, is_big_endian, mmap=False):
    if is_big_endian:
        fmt = '>I'
    else:
        fmt = '<I'
    size = struct.unpack(fmt, fid.read(4))[0]

    bytes = bits//8
    if bits == 8:
        dtype = 'u1'
    else:
        if is_big_endian:
            dtype = '>'
        else:
            dtype = '<'
        if comp == 1:
            dtype += 'i%d' % bytes
        else:
            dtype += 'f%d' % bytes
    if not mmap:
        data = numpy.fromstring(fid.read(size), dtype=dtype)
    else:
        start = fid.tell()
        data = numpy.memmap(fid, dtype=dtype, mode='c', offset=start,
                            shape=(size//bytes,))
        fid.seek(start + size)

    if noc > 1:
        data = data.reshape(-1, noc)
    return data


def _skip_unknown_chunk(fid, is_big_endian):
    if is_big_endian:
        fmt = '>I'
    else:
        fmt = '<I'

    data = fid.read(4)
    # call unpack() and seek() only if we have really read data from file
    # otherwise empty read at the end of the file would trigger 
    # unnecessary exception at unpack() call
    # in case data equals somehow to 0, there is no need for seek() anyway
    if data:
        size = struct.unpack(fmt, data)[0]
        fid.seek(size, 1)


def _read_riff_chunk(fid):
    str1 = fid.read(4)
    if str1 == b'RIFF':
        is_big_endian = False
    elif str1 == b'RIFX':
        is_big_endian = True
    else:
        raise ValueError("Not a WAV file.")

    if is_big_endian:
        fmt = '>I'
    else:
        fmt = '<I'
    fsize = struct.unpack(fmt, fid.read(4))[0] + 8
    str2 = fid.read(4)
    if (str2 != b'WAVE'):
        raise ValueError("Not a WAV file.")

    return fsize, is_big_endian

# open a wave-file


def read(filename, mmap=False):
    """
    Return the sample rate (in samples/sec) and data from a WAV file

    Parameters
    ----------
    filename : string or open file handle
        Input wav file.
    mmap : bool, optional
        Whether to read data as memory mapped.
        Only to be used on real files (Default: False)

        .. versionadded:: 0.12.0

    Returns
    -------
    rate : int
        Sample rate of wav file
    data : numpy array
        Data read from wav file

    Notes
    -----
    * The file can be an open file or a filename.
    * The returned sample rate is a Python integer.
    * The data is returned as a numpy array with a data-type determined
      from the file.
    * This function cannot read wav files with 24 bit data.

    """
    if hasattr(filename, 'read'):
        fid = filename
        mmap = False
    else:
        fid = open(filename, 'rb')

    try:
        fsize, is_big_endian = _read_riff_chunk(fid)
        fmt_chunk_received = False
        noc = 1
        bits = 8
        comp = WAVE_FORMAT_PCM
        while (fid.tell() < fsize):
            # read the next chunk
            chunk_id = fid.read(4)
            if chunk_id == b'fmt ':
                fmt_chunk_received = True
                size, comp, noc, rate, sbytes, ba, bits = _read_fmt_chunk(fid, is_big_endian=is_big_endian)
                if bits == 24:
                    raise ValueError("Unsupported bit depth: the wav file "
                                     "has 24 bit data.")
            elif chunk_id == b'fact':
                _skip_unknown_chunk(fid, is_big_endian=is_big_endian)
            elif chunk_id == b'data':
                if not fmt_chunk_received:
                    raise ValueError("No fmt chunk before data")
                data = _read_data_chunk(fid, comp, noc, bits, is_big_endian=is_big_endian, mmap=mmap)
            elif chunk_id == b'LIST':
                # Someday this could be handled properly but for now skip it
                _skip_unknown_chunk(fid, is_big_endian=is_big_endian)
            elif chunk_id in (b'JUNK', b'Fake'):
                # Skip alignment chunks without warning
                _skip_unknown_chunk(fid, is_big_endian=is_big_endian)
            else:
                warnings.warn("Chunk (non-data) not understood, skipping it.",
                              WavFileWarning)
                _skip_unknown_chunk(fid, is_big_endian=is_big_endian)
    finally:
        if not hasattr(filename, 'read'):
            fid.close()
        else:
            fid.seek(0)

    return rate, data

# Write a wave-file
# sample rate, data


def write(filename, rate, data):
    """
    Write a numpy array as a WAV file

    Parameters
    ----------
    filename : string or open file handle
        Output wav file
    rate : int
        The sample rate (in samples/sec).
    data : ndarray
        A 1-D or 2-D numpy array of either integer or float data-type.

    Notes
    -----
    * The file can be an open file or a filename.

    * Writes a simple uncompressed WAV file.
    * The bits-per-sample will be determined by the data-type.
    * To write multiple-channels, use a 2-D array of shape
      (Nsamples, Nchannels).

    """
    if hasattr(filename, 'write'):
        fid = filename
    else:
        fid = open(filename, 'wb')

    try:
        dkind = data.dtype.kind
        if not (dkind == 'i' or dkind == 'f' or (dkind == 'u' and
                                                 data.dtype.itemsize == 1)):
            raise ValueError("Unsupported data type '%s'" % data.dtype)
        
        header_data = b''

        header_data += b'RIFF'
        header_data += b'\x00\x00\x00\x00'
        header_data += b'WAVE'

        # fmt chunk
        header_data += b'fmt '
        if dkind == 'f':
            comp = 3
        else:
            comp = 1
        if data.ndim == 1:
            noc = 1
        else:
            noc = data.shape[1]
        bits = data.dtype.itemsize * 8
        sbytes = rate*(bits // 8)*noc
        ba = noc * (bits // 8)
        
        fmt_chunk_data = struct.pack('<HHIIHH', comp, noc, rate, sbytes,
                                     ba, bits)
        if not (dkind == 'i' or dkind == 'u'):
            # add cbSize field for non-PCM files
            fmt_chunk_data += b'\x00\x00'
        
        header_data += struct.pack('<I', len(fmt_chunk_data))
        header_data += fmt_chunk_data

        # fact chunk (non-PCM files)
        if not (dkind == 'i' or dkind == 'u'):
            header_data += b'fact'
            header_data += struct.pack('<II', 4, data.shape[0])

        # check data size (needs to be immediately before the data chunk)
        if ((len(header_data)-4-4) + (4+4+data.nbytes)) > 0xFFFFFFFF:
            raise ValueError("Data exceeds wave file size limit")
        
        fid.write(header_data)
        
        # data chunk
        fid.write(b'data')
        fid.write(struct.pack('<I', data.nbytes))
        if data.dtype.byteorder == '>' or (data.dtype.byteorder == '=' and
                                           sys.byteorder == 'big'):
            data = data.byteswap()
        _array_tofile(fid, data)

        # Determine file size and place it in correct
        #  position at start of the file.
        size = fid.tell()
        fid.seek(4)
        fid.write(struct.pack('<I', size-8))

    finally:
        if not hasattr(filename, 'write'):
            fid.close()
        else:
            fid.seek(0)


if sys.version_info[0] >= 3:
    def _array_tofile(fid, data):
        # ravel gives a c-contiguous buffer
        fid.write(data.ravel().view('b').data)
else:
    def _array_tofile(fid, data):
        fid.write(data.tostring())
