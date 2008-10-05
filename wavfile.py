import numpy
import struct

# assumes file pointer is immediately
#  after the 'fmt ' id
def _read_fmt_chunk(fid):
    res = struct.unpack('ihHIIHH',fid.read(20))
    size, comp, noc, rate, sbytes, ba, bits = res
    if (comp != 1 or size > 16):
        print "Warning:  unfamiliar format bytes..."
        if (size>16):
            fid.read(size-16)
    return size, comp, noc, rate, sbytes, ba, bits

# assumes file pointer is immediately
#   after the 'data' id
def _read_data_chunk(fid, noc, bits):
    size = struct.unpack('i',fid.read(4))[0]
    if bits == 8:
        data = numpy.fromfile(fid, dtype=numpy.ubyte, count=size)
        if noc > 1:
            data = data.reshape(-1,noc)
    else:
        bytes = bits//8
        dtype = 'i%d' % bytes
        data = numpy.fromfile(fid, dtype=dtype, count=size//bytes)
        if noc > 1:
            data = data.reshape(-1,noc)
    return data

def _read_riff_chunk(fid):
    str1 = fid.read(4)
    fsize = struct.unpack('I', fid.read(4))[0] + 8
    str2 = fid.read(4)
    if (str1 != 'RIFF' or str2 != 'WAVE'):
        raise ValueError, "Not a WAV file."
    return fsize

# open a wave-file
def read(file):
    """Return the sample rate (in samples/sec) and data from a WAV file

    The file can be an open file or a filename.
    The returned sample rate is a Python integer
    The data is returned as a numpy array with a
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
            print "Reading fmt chunk"
            size, comp, noc, rate, sbytes, ba, bits = _read_fmt_chunk(fid)
        elif chunk_id == 'data':
            print "Reading data chunk"
            data = _read_data_chunk(fid, noc, bits)
        else:
            print "Warning:  %s chunk not understood"
            size = struct.unpack('I',fid.read(4))[0]
            bytes = fid.read(size)
    fid.close()
    return rate, data

# Write a wave-file
# sample rate, data
def write(filename, rate, data):
    """Write a numpy array as a WAV file

    filename -- The name of the file to write (will be over-written)
    rate -- The sample rate (in samples/sec).
    data -- A 1-d or 2-d numpy array of integer data-type.
            The bits-per-sample will be determined by the data-type
            To write multiple-channels, use a 2-d array of shape
            (Nsamples, Nchannels)

    Writes a simple uncompressed WAV file.
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
    fid.write(struct.pack('ihHIIHH', 16, 1, noc, rate, sbytes, ba, bits))
    # data chunk
    fid.write('data')
    fid.write(struct.pack('i', data.nbytes))
    data.tofile(fid)
    # Determine file size and place it in correct
    #  position at start of the file.
    size = fid.tell()
    fid.seek(4)
    fid.write(struct.pack('i', size-8))
    fid.close()
