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
from collections import defaultdict


__all__ = [
    'WavFileWarning',
    'read',
    'write'
]


class WavFileWarning(UserWarning):
    pass


# see the Wave File Format reference for details of chunk layouts:
# https://web.archive.org/web/20141226210234/http://www.sonicspot.com/guide/wavefiles.html

WAVE_FORMAT_PCM = 0x0001
WAVE_FORMAT_IEEE_FLOAT = 0x0003
WAVE_FORMAT_EXTENSIBLE = 0xfffe
KNOWN_WAVE_FORMATS = (WAVE_FORMAT_PCM, WAVE_FORMAT_IEEE_FLOAT)


# assumes file pointer is immediately
#  after the 'fmt ' id
def _read_fmt_chunk(fid, is_big_endian):
    """
    Returns
    -------
    format_tag : int
        PCM, float, or compressed format
    channels : int
        number of channels
    fs : int
        sampling frequency in samples per second
    bytes_per_second : int
        overall byte rate for the file
    block_align : int
        bytes per sample, including all channels
    bit_depth : int
        bits per sample
    """
    if is_big_endian:
        fmt = '>'
    else:
        fmt = '<'

    size = struct.unpack(fmt+'I', fid.read(4))[0]

    if size < 16:
        raise ValueError("Binary structure of wave file is not compliant")

    (format_tag, channels, fs, bytes_per_second, block_align, bit_depth
     ) = struct.unpack(fmt+'HHIIHH', fid.read(16))
    bytes_read = 16

    if format_tag == WAVE_FORMAT_EXTENSIBLE and size >= 18:
        ext_chunk_size = struct.unpack(fmt+'H', fid.read(2))[0]
        bytes_read += 2
        if ext_chunk_size >= 22:
            extensible_chunk_data = fid.read(22)
            bytes_read += 22
            raw_guid = extensible_chunk_data[6:22]
            # GUID template {XXXXXXXX-0000-0010-8000-00AA00389B71} (RFC-2361)
            # MS GUID byte order: first three groups are native byte order,
            # rest is Big Endian
            if is_big_endian:
                tail = b'\x00\x00\x00\x10\x80\x00\x00\xAA\x00\x38\x9B\x71'
            else:
                tail = b'\x00\x00\x10\x00\x80\x00\x00\xAA\x00\x38\x9B\x71'
            if raw_guid.endswith(tail):
                format_tag = struct.unpack(fmt+'I', raw_guid[:4])[0]
        else:
            raise ValueError("Binary structure of wave file is not compliant")

    if format_tag not in KNOWN_WAVE_FORMATS:
        raise ValueError("Unknown wave file format")

    # move file pointer to next chunk
    if size > bytes_read:
        fid.read(size - bytes_read)

    return format_tag, channels, fs, bytes_per_second, block_align, bit_depth


# assumes file pointer is immediately after the 'data' id
def _read_data_chunk(fid, format_tag, channels, bit_depth, is_big_endian,
                     mmap=False):
    if is_big_endian:
        fmt = '>'
    else:
        fmt = '<'

    # Size of the data subchunk in bytes
    size = struct.unpack(fmt+'I', fid.read(4))[0]

    # Number of bytes per sample
    bytes_per_sample = bit_depth//8
    if bit_depth in (8, 24):
        dtype = 'u1'
        bytes_per_sample = 1
    elif format_tag == WAVE_FORMAT_PCM:
        dtype = '%si%d' % (fmt, bytes_per_sample)
    else:
        dtype = '%sf%d' % (fmt, bytes_per_sample)

    if not mmap:
        data = numpy.fromstring(fid.read(size), dtype=dtype)
    else:
        start = fid.tell()
        data = numpy.memmap(fid, dtype=dtype, mode='c', offset=start,
                            shape=(size//bytes_per_sample,))
        fid.seek(start + size)

    # if odd number of bytes, move 1 byte further (data chunk is word-aligned)
    if size % 2 == 1:
        fid.seek(1, 1)

    if bit_depth == 24:
        a = numpy.empty((len(data)/3, 4), dtype='u1')
        a[:, :3] = data.reshape((-1, 3))
        a[:, 3:] = (a[:, 3 - 1:3] >> 7) * 255
        data = a.view('<i4').reshape(a.shape[:-1])

    if channels > 1:
        data = data.reshape(-1, channels)
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
        # if odd number of bytes, move 1 byte further
        # (data chunk is word-aligned)
        if size % 2 == 1:
            size += 1
        fid.seek(size, 1)


def _read_riff_chunk(fid):
    str1 = fid.read(4)  # File signature
    if str1 == b'RIFF':
        is_big_endian = False
        fmt = '<I'
    elif str1 == b'RIFX':
        is_big_endian = True
        fmt = '>I'
    else:
        # There are also .wav files with "FFIR" or "XFIR" signatures?
        raise ValueError("File format %r not understood." % str1)

    # Size of entire file
    file_size = struct.unpack(fmt, fid.read(4))[0] + 8

    str2 = fid.read(4)
    if str2 != b'WAVE':
        raise ValueError("Not a WAV file.")

    return file_size, is_big_endian


def read(filename, mmap=False, return_cues=False, return_pitch=False):
    """
    Open a WAV file

    Return the sample rate (in samples/sec) and data from a WAV file

    Parameters
    ----------
    filename : string or open file handle
        Input wav file.
    mmap : bool, optional
        Whether to read data as memory-mapped.
        Only to be used on real files (Default: False).

        .. versionadded:: 0.12.0

    Returns
    -------
    rate : int
        Sample rate of wav file.
    data : numpy array
        Data read from wav file.  Data-type is determined from the file;
        see Notes.

    Notes
    -----
    The returned sample rate is a Python integer.

    Common data types: [1]_

    =====================  ===========  ===========  =============
         WAV format            Min          Max       NumPy dtype
    =====================  ===========  ===========  =============
    32-bit floating-point  -1.0         +1.0         float32
    32-bit PCM             -2147483648  +2147483647  int32
    16-bit PCM             -32768       +32767       int16
    8-bit PCM              0            255          uint8
    =====================  ===========  ===========  =============

    Note that 8-bit PCM is unsigned.

    References
    ----------
    .. [1] IBM Corporation and Microsoft Corporation, "Multimedia Programming
       Interface and Data Specifications 1.0", section "Data Format of the
       Samples", August 1991
       http://www-mmsp.ece.mcgill.ca/documents/audioformats/wave/Docs/riffmci.pdf
    """
    if hasattr(filename, 'read'):
        fid = filename
        mmap = False
    else:
        fid = open(filename, 'rb')

    try:
        file_size, is_big_endian = _read_riff_chunk(fid)
        fmt_chunk_received = False
        channels = 1
        bit_depth = 8
        format_tag = WAVE_FORMAT_PCM
        cues = defaultdict(dict)
        pitch = 0.0
        while fid.tell() < file_size:
            # read the next chunk
            chunk_id = fid.read(4)

            if not chunk_id:
                raise ValueError("Unexpected end of file.")
            elif len(chunk_id) < 4:
                raise ValueError("Incomplete wav chunk.")

            if chunk_id == b'fmt ':
                fmt_chunk_received = True
                fmt_chunk = _read_fmt_chunk(fid, is_big_endian)
                format_tag, channels, fs, _, __, bit_depth = fmt_chunk
            elif chunk_id == b'data':
                if not fmt_chunk_received:
                    raise ValueError("No fmt chunk before data")
                data = _read_data_chunk(fid, format_tag, channels, bit_depth,
                                        is_big_endian, mmap=mmap)
            elif chunk_id == b'cue ':
                _, num_cues = struct.unpack('<ii', fid.read(8))
                for c in range(num_cues):
                    str1 = fid.read(24)
                    cue_id, position = struct.unpack('<ii', str1)
                    cues[cue_id]['pos'] = position
            elif chunk_id == b'labl':
                str1 = fid.read(8)
                size, cue_id = struct.unpack('<ii', str1)
                size += (size % 2)  # ensure size is even
                cues[cue_id]['label'] = fid.read(size-4).rstrip('\x00')
            elif chunk_id == b'smpl':
                size = struct.unpack('<i', fid.read(4))[0]
                str1 = fid.read(size)
                unity_note, pitch_fraction = struct.unpack('<ii', str1[12:20])
                cents = pitch_fraction / (2**32-1)
                pitch = 440 * 2 ** ((unity_note + cents - 69)/12)
            # see http://www.pjb.com.au/midi/sfspec21.html#i5
            elif chunk_id in [b'ICRD', b'IENG', b'ISFT', b'ISTJ']:
                _skip_unknown_chunk(fid, is_big_endian)
            elif chunk_id in [b'fact', b'LIST', b'JUNK', b'Fake']:
                _skip_unknown_chunk(fid, is_big_endian)
            else:
                warnings.warn("Chunk (non-data) not understood, skipping it.",
                              WavFileWarning)
                _skip_unknown_chunk(fid, is_big_endian)
    finally:
        if not hasattr(filename, 'read'):
            fid.close()
        else:
            fid.seek(0)

    result = [fs, data]
    if return_cues:
        result.append(dict(cues))
    if return_pitch:
        result.append(pitch)
    return tuple(result)


def write(filename, rate, data, cues=None, loops=None, bitrate=None):
    """
    Write a numpy array as a WAV file.

    Parameters
    ----------
    filename : string or open file handle
        Output wav file.
    rate : int
        The sample rate (in samples/sec).
    data : ndarray
        A 1-D or 2-D numpy array of either integer or float data-type.
    cues : sequence of ints, optional
        Play order positions of cues.
    loops : sequence of (int,int) pairs, optional
        Pairs of (unity note, pitch fraction) values.
    bitrate : int, optional
        The number of bits per sample.
        If None, bitrate is determined by the data-type.

    Notes
    -----
    * Writes a simple uncompressed WAV file.
    * To write multiple-channels, use a 2-D array of shape
      (Nsamples, Nchannels).
    * The WAV type (PCM/float) is determined by the input data-type.

    Common data types: [1]_

    =====================  ===========  ===========  =============
         WAV format            Min          Max       NumPy dtype
    =====================  ===========  ===========  =============
    32-bit floating-point  -1.0         +1.0         float32
    32-bit PCM             -2147483648  +2147483647  int32
    16-bit PCM             -32768       +32767       int16
    8-bit PCM              0            255          uint8
    =====================  ===========  ===========  =============

    Note that 8-bit PCM is unsigned.

    References
    ----------
    .. [1] IBM Corporation and Microsoft Corporation, "Multimedia Programming
       Interface and Data Specifications 1.0", section "Data Format of the
       Samples", August 1991
       http://www-mmsp.ece.mcgill.ca/documents/audioformats/wave/Docs/riffmci.pdf

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

        header_data = b'RIFF\x00\x00\x00\x00WAVEfmt '
        # fmt chunk
        if dkind == 'f':
            format_tag = WAVE_FORMAT_IEEE_FLOAT
        else:
            format_tag = WAVE_FORMAT_PCM
        if data.ndim == 1:
            channels = 1
        else:
            channels = data.shape[1]
        if bitrate is None:
            bit_depth = data.dtype.itemsize * 8
        elif bitrate != 24 and bitrate != data.dtype.itemsize * 8:
            raise ValueError("Unsupported bitrate for dtype: %s" % data.dtype)
        else:
            bit_depth = bitrate

        bytes_per_second = rate * (bit_depth // 8) * channels
        block_align = channels * (bit_depth // 8)
        fmt_chunk_data = struct.pack('<HHIIHH', format_tag, channels, rate,
                                     bytes_per_second, block_align, bit_depth)
        if format_tag != WAVE_FORMAT_PCM:
            # add cbSize field for non-PCM files
            fmt_chunk_data += b'\x00\x00'

        header_data += struct.pack('<I', len(fmt_chunk_data))
        header_data += fmt_chunk_data

        # fact chunk (non-PCM files)
        if format_tag != WAVE_FORMAT_PCM:
            header_data += b'fact'
            header_data += struct.pack('<II', 4, data.shape[0])

        # check data size (needs to be immediately before the data chunk)
        if len(header_data) + data.nbytes > 0xFFFFFFFF:
            raise ValueError("Data exceeds wave file size limit")

        fid.write(header_data)

        # data chunk
        if bitrate == 24:
            a32 = numpy.asarray(data, dtype=numpy.int32)
            if a32.ndim == 1:
                # Convert to a 2D array with a single column.
                a32.shape = a32.shape + (1,)
            # By shifting first 0 bits, then 8, then 16,
            # the resulting output is 24 bit little-endian.
            a8 = (a32.reshape(a32.shape + (1,)) >> numpy.array([0,8,16])) & 255
            data = a8.astype(numpy.uint8)

        fid.write(b'data')
        fid.write(struct.pack('<I', data.nbytes))
        if data.dtype.byteorder == '>' or (data.dtype.byteorder == '=' and
                                           sys.byteorder == 'big'):
            data = data.byteswap()
        _array_tofile(fid, data)

        # cue chunk
        if cues:
            fid.write(b'cue ')
            size = 4 + len(cues) * 24
            fid.write(struct.pack('<ii', size, len(cues)))
            for i, c in enumerate(cues):
                # 1635017060 is struct.unpack('<i',b'data')
                fid.write(struct.pack('<iiiiii', i + 1, c, 1635017060, 0, 0, c))

        # smpl chunk
        if loops:
            fid.write(b'smpl')
            size = 36 + len(loops) * 24
            sample_period = int(1000000000 / rate)
            fid.write(struct.pack('<iiiiiiiiii', size, 0, 0, sample_period, 0,
                                  0, 0, 0, len(loops), 0))
            for i, loop in enumerate(loops):
                fid.write(struct.pack('<iiiiii', 0, 0, loop[0], loop[1], 0, 0))

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
