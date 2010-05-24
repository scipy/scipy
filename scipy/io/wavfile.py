"""
Module to read / write wav files using numpy arrays

Functions
---------
read: Return the sample rate (in samples/sec) and data from a WAV file.

write: Write a numpy array as a WAV file.

"""
import numpy
import struct
import warnings

class WavFileWarning(UserWarning):
    pass

_big_endian = False

# assumes file pointer is immediately
#  after the 'fmt ' id
def _read_fmt_chunk(fid):
    if _big_endian:
        fmt = '>'
    else:
        fmt = '<'
    res = struct.unpack(fmt+'ihHIIHH',fid.read(20))
    size, comp, noc, rate, sbytes, ba, bits = res
    if (comp != 1 or size > 16):
        warnings.warn("Unfamiliar format bytes", WavFileWarning)
        if (size>16):
            fid.read(size-16)
    return size, comp, noc, rate, sbytes, ba, bits

# assumes file pointer is immediately
#   after the 'data' id
def _read_data_chunk(fid, noc, bits):
    if _big_endian:
        fmt = '>i'
    else:
        fmt = '<i'
    size = struct.unpack(fmt,fid.read(4))[0]
    if bits == 8:
        data = numpy.fromfile(fid, dtype=numpy.ubyte, count=size)
        if noc > 1:
            data = data.reshape(-1,noc)
    else:
        bytes = bits//8
        if _big_endian:
            dtype = '>i%d' % bytes
        else:
            dtype = '<i%d' % bytes
        data = numpy.fromfile(fid, dtype=dtype, count=size//bytes)
        if noc > 1:
            data = data.reshape(-1,noc)
    return data

def _read_riff_chunk(fid):
    global _big_endian
    str1 = fid.read(4)
    if str1 == 'RIFX':
        _big_endian = True
    elif str1 != 'RIFF':
        raise ValueError("Not a WAV file.")
    if _big_endian:
        fmt = '>I'
    else:
        fmt = '<I'
    fsize = struct.unpack(fmt, fid.read(4))[0] + 8
    str2 = fid.read(4)
    if (str2 != 'WAVE'):
        raise ValueError("Not a WAV file.")
    if str1 == 'RIFX':
        _big_endian = True
    return fsize

# open a wave-file
def read(file):
    """
    Return the sample rate (in samples/sec) and data from a WAV file

    Parameters
    ----------
    file : file
        Input wav file.

    Returns
    -------
    rate : int
        Sample rate of wav file
    data : numpy array
        Data read from wav file

    Notes
    -----

    * The file can be an open file or a filename.

    * The returned sample rate is a Python integer
    * The data is returned as a numpy array with a
      data-type determined from the file.

    """
    if hasattr(file,'read'):
        fid = file
    else:
        fid = open(file, 'rb')

    fsize = _read_riff_chunk(fid)
    noc = 1
    bits = 8
    while (fid.tell() < fsize):
        # read the next chunk
        chunk_id = fid.read(4)
        if chunk_id == 'fmt ':
            size, comp, noc, rate, sbytes, ba, bits = _read_fmt_chunk(fid)
        elif chunk_id == 'data':
            data = _read_data_chunk(fid, noc, bits)
        else:
            warnings.warn("chunk not understood", WavFileWarning)
            size = struct.unpack('I',fid.read(4))[0]
            bytes = fid.read(size)
    fid.close()
    return rate, data

# Write a wave-file
# sample rate, data
def write(filename, rate, data):
    """
    Write a numpy array as a WAV file

    Parameters
    ----------
    filename : file
        The name of the file to write (will be over-written).
    rate : int
        The sample rate (in samples/sec).
    data : ndarray
        A 1-D or 2-D numpy array of integer data-type.

    Notes
    -----
    * Writes a simple uncompressed WAV file.
    * The bits-per-sample will be determined by the data-type.
    * To write multiple-channels, use a 2-D array of shape
      (Nsamples, Nchannels).

    """
    fid = open(filename, 'wb')
    fid.write('RIFF')
    fid.write('\x00\x00\x00\x00')
    fid.write('WAVE')
    # fmt chunk
    fid.write('fmt ')
    if data.ndim == 1:
        noc = 1
    else:
        noc = data.shape[1]
    bits = data.dtype.itemsize * 8
    sbytes = rate*(bits / 8)*noc
    ba = noc * (bits / 8)
    fid.write(struct.pack('<ihHIIHH', 16, 1, noc, rate, sbytes, ba, bits))
    # data chunk
    fid.write('data')
    fid.write(struct.pack('<i', data.nbytes))
    import sys
    if data.dtype.byteorder == '>' or (data.dtype.byteorder == '=' and sys.byteorder == 'big'):
        data = data.byteswap()
    data.tofile(fid)
    # Determine file size and place it in correct
    #  position at start of the file.
    size = fid.tell()
    fid.seek(4)
    fid.write(struct.pack('<i', size-8))
    fid.close()
