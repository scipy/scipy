''' Classes for read / write of matlab (TM) 5 files

The matfile specification last found here:

http://www.mathworks.com/access/helpdesk/help/pdf_doc/matlab/matfile_format.pdf

(as of December 5 2008)
'''
from __future__ import division, print_function, absolute_import

'''
=================================
 Note on functions and mat files
=================================

The document above does not give any hints as to the storage of matlab
function handles, or anonymous function handles.  I had therefore to
guess the format of matlab arrays of ``mxFUNCTION_CLASS`` and
``mxOPAQUE_CLASS`` by looking at example mat files.

``mxFUNCTION_CLASS`` stores all types of matlab functions.  It seems to
contain a struct matrix with a set pattern of fields.  For anonymous
functions, a sub-fields of one of these fields seems to contain the
well-named ``mxOPAQUE_CLASS``. This seems to cotain:

* array flags as for any matlab matrix
* 3 int8 strings
* a matrix

It seems that, whenever the mat file contains a ``mxOPAQUE_CLASS``
instance, there is also an un-named matrix (name == '') at the end of
the mat file.  I'll call this the ``__function_workspace__`` matrix.

When I saved two anonymous functions in a mat file, or appended another
anonymous function to the mat file, there was still only one
``__function_workspace__`` un-named matrix at the end, but larger than
that for a mat file with a single anonymous function, suggesting that
the workspaces for the two functions had been merged.

The ``__function_workspace__`` matrix appears to be of double class
(``mxCLASS_DOUBLE``), but stored as uint8, the memory for which is in
the format of a mini .mat file, without the first 124 bytes of the file
header (the description and the subsystem_offset), but with the version
U2 bytes, and the S2 endian test bytes.  There follow 4 zero bytes,
presumably for 8 byte padding, and then a series of ``miMATRIX``
entries, as in a standard mat file. The ``miMATRIX`` entries appear to
be series of un-named (name == '') matrices, and may also contain arrays
of this same mini-mat format.

I guess that:

* saving an anonymous function back to a mat file will need the
  associated ``__function_workspace__`` matrix saved as well for the
  anonymous function to work correctly.
* appending to a mat file that has a ``__function_workspace__`` would
  involve first pulling off this workspace, appending, checking whether
  there were any more anonymous functions appended, and then somehow
  merging the relevant workspaces, and saving at the end of the mat
  file.

The mat files I was playing with are in ``tests/data``:

* sqr.mat
* parabola.mat
* some_functions.mat

See ``tests/test_mio.py:test_mio_funcs.py`` for a debugging
script I was working with.

'''

# Small fragments of current code adapted from matfile.py by Heiko
# Henkelmann

from collections import namedtuple
from io import BytesIO
from itertools import islice
import os
import struct
import sys
import time
import warnings
import zlib

import numpy as np
from numpy.compat import asbytes

import scipy.sparse

from scipy.lib.six import string_types, unichr as chr

from . import byteordercodes as boc
from .byteordercodes import native_code, swapped_code

from .miobase import (MatFileReader, docfiller, matdims, arr_to_chars,
                      arr_dtype_number, MatWriteError, MatReadError,
                      MatReadWarning)

# Reader object for matlab 5 format variables
from .mio5_utils import squeeze_element

# Constants and helper objects
from .mio5_params import (
    MatlabObject, MatlabFunction, MatlabOpaque, MDTYPES, NP_TO_MTYPES,
    mclass_dtypes_template, mdtypes_template, NP_TO_MXTYPES, miCOMPRESSED,
    miMATRIX, miINT8, miUTF8, miUINT32, mxCELL_CLASS, mxSTRUCT_CLASS,
    mxOBJECT_CLASS, mxCHAR_CLASS, mxSPARSE_CLASS, mxDOUBLE_CLASS,
    mxFUNCTION_CLASS, mxOPAQUE_CLASS, mclass_info, mat_struct)

from .streams import ZlibInputStream


MatlabArray = namedtuple("MatlabArray", "name data is_global")
MatInfo = namedtuple("MatInfo", "name shape info stream data_position nzmax")


class MatFile5Reader(MatFileReader):
    ''' Reader for Mat 5 mat files
    Adds the following attribute to base class

    uint16_codec - char codec to use for uint16 char arrays
        (defaults to system default codec)

    Uses variable reader that has the following stardard interface (see
    abstract class in ``miobase``::

       __init__(self, file_reader)
       read_header(self)
       array_from_header(self)

    and added interface::
       set_stream(self, stream)
       read_full_tag(self)
    ''' # FIXME

    @docfiller
    def __init__(self,
                 stream,
                 byte_order=None,
                 mat_dtype=False,
                 squeeze_me=False,
                 chars_as_strings=True,
                 matlab_compatible=False,
                 struct_as_record=True,
                 verify_compressed_data_integrity=True,
                 uint16_codec=None): # FIXME
        '''Initializer for matlab 5 file format reader

    %(matstream_arg)s
    %(load_args)s
    %(struct_arg)s
    uint16_codec : {None, string}
        Set codec to use for uint16 char arrays (e.g. 'utf-8').
        Use system default codec if None
        ''' # FIXME

        self._stream = stream
        self._read_header()

        if (byte_order is not None and
            boc.to_numpy_code(byte_order) != self._endian):
            raise ValueError("Incompatible byte order")
        self._mat_dtype = mat_dtype
        self._squeeze_me = squeeze_me
        self._chars_as_strings = chars_as_strings
        if matlab_compatible:
            self.set_matlab_compatible()
        self._struct_as_record = struct_as_record
        self._verify_compressed_data_integrity = (
            verify_compressed_data_integrity)
        self._uint16_codec = uint16_codec or sys.getdefaultencoding()

    def minimat_reader(self, **kwargs):
        self._stream.seek(self._subsys_offset)
        data = next(self._read_iter()).data.tostring()
        if data[4:8] != b"\0" * 4:
            raise ValueError("Invalid padding of function workspace")
        reader = type(self)(BytesIO(b"\0" * 124 + data[:4] + data[8:], **kwargs))
        # The minimat does not always declare sizes properly.
        reader._check_and_pad_stream = (
            lambda stream, _: stream.seek((-stream.tell()) % 8, 1))
        return reader

    def set_matlab_compatible(self):
        ''' Sets options to return arrays as MATLAB loads them '''
        self._mat_dtype = True
        self._squeeze_me = False
        self._chars_as_strings = False

    def _prepare_stream(self):
        self._stream.seek(128)

    def close(self):
        self._stream.close()

    def _read_header(self):
        self._stream.seek(0)
        self._header = self._stream.read(128)
        self._desc = self._header[:116]
        self._endian = {b"IM": "<", b"MI": ">"}[self._header[126:128]]
        self._subsys_offset, ver = self._unpack("QH", self._header[116:126])
        if ver != 0x0100:
            raise ValueError("Unsupported version: {:#04x}".format(ver))
        self._version = "1.0"

    def _unpack(self, fmt, data):
        return struct.unpack(self._endian + fmt, data)

    @staticmethod
    def _as_identifiers(data):
        return [  # Extra call to str to avoid returning unicode on Python 2.
            name for name in str(data.tostring().decode("ascii")).split("\0")
            if name]

    def _read_iter(self, stream=None, info_only=None, load_only=None):
        if stream is None:
            stream = self._stream

        while True:
            entry_start = stream.tell()
            raw_header0 = stream.read(4)
            if not raw_header0:
                return
            header0, = self._unpack("I", raw_header0)
            nbytes, mdtype = divmod(header0, 0x10000)
            if not nbytes:
                mdtype = header0
                nbytes, = self._unpack("I", stream.read(4))
            entry_end = stream.tell() + nbytes

            if mdtype == miCOMPRESSED:
                try:
                    for entry in self._read_iter(
                        ZlibInputStream(stream, nbytes),
                        info_only=info_only, load_only=load_only):
                        yield entry
                except ValueError:
                    if not self._verify_compressed_data_integrity:
                        yield None
                    else:
                        raise ValueError("Invalid compressed data")

            elif mdtype in mdtypes_template:
                dtype = self._endian + mdtypes_template[mdtype]
                data = stream.read(nbytes)
                self._check_and_pad_stream(stream, entry_end)
                yield np.fromstring(data, dtype)

            elif mdtype == miMATRIX:
                reader = self._read_iter(stream)

                flags = next(reader)
                if isinstance(flags, MatlabArray):
                    # This can only occur while reading the function workspace.
                    self._check_and_pad_stream(stream, entry_end)
                    yield flags
                    continue
                else:
                    flags, nzmax = flags
                dims, name = list(islice(reader, 2))
                name, = self._as_identifiers(name) or [""]
                matrix_cls = flags % 0x100
                f_complex = flags & (1 << 11)
                f_global = flags & (1 << 10)
                f_logical = flags & (1 << 9)

                if info_only:
                    if matrix_cls == mxCHAR_CLASS:
                        dims = dims[:-1]
                    class_info = (
                        "logical" if f_logical else mclass_info[matrix_cls])
                    stream.seek(entry_end)
                    self._check_and_pad_stream(stream, entry_end)
                    yield MatInfo(name, dims, class_info, stream,
                                  slice(entry_start, entry_end), nzmax)
                    continue

                if load_only is not None and name not in load_only:
                    stream.seek(entry_end)
                    self._check_and_pad_stream(stream, entry_end)
                    yield None
                    continue

                if matrix_cls == mxCELL_CLASS:
                    dtype = object
                    pr = np.empty(np.product(dims), dtype=dtype)
                    pr[:] = [entry.data for entry in
                             islice(reader, int(np.product(dims)))]
                elif matrix_cls in [mxSTRUCT_CLASS, mxOBJECT_CLASS]:
                    # Struct: field name length, field names
                    # Object: class name, field name length, field names
                    if matrix_cls == mxOBJECT_CLASS:
                        classname, = self._as_identifiers(next(reader))
                    next(reader)  # Drop field name length.
                    fields = self._as_identifiers(next(reader))
                    dtype = ([(field, object) for field in fields]
                             if fields else object)
                    print("DT", dtype)
                    pr = np.empty(np.product(dims), dtype=dtype)
                    for p in pr:
                        for field in fields:
                            p[field] = next(reader).data
                    if matrix_cls == mxOBJECT_CLASS:
                        pr = MatlabObject(pr, classname)
                    elif not self._struct_as_record:
                        # Deprecated
                        pr2 = np.empty_like(pr, dtype=object)
                        dtype = object
                        for i, entry in enumerate(pr):
                            pr2[i] = obj = mat_struct()
                            obj._fieldnames = fields
                            obj.__dict__.update(zip(fields, entry))
                        pr = pr2
                elif matrix_cls == mxSPARSE_CLASS:
                    ir, jc, pr = islice(reader, 3)
                    dtype = np.float64
                elif matrix_cls == mxFUNCTION_CLASS:
                    pr = MatlabFunction(next(reader).data)
                    dtype = object
                elif matrix_cls == mxOPAQUE_CLASS:
                    opaque_components = []
                    while stream.tell() < entry_end:
                        opaque_components.append(next(reader))
                    pr = MatlabOpaque(
                        np.empty(dims, dtype=[
                            ("s{}".format(i), "O")
                            for i in range(len(opaque_components))]))
                    for i, component in enumerate(opaque_components):
                        pr[()]["s{}".format(i)] = component
                    dtype = object
                else:
                    pr = next(reader)
                    dtype = (pr.dtype if not self._mat_dtype else
                             np.bool if f_logical else
                             mclass_dtypes_template[matrix_cls])

                pr = pr.astype(dtype)
                if f_complex:
                    pi = next(reader).astype(dtype)
                    pr = pr + 1j * pi

                if matrix_cls == mxCHAR_CLASS:
                    joiner = "".join if self._chars_as_strings else list
                    # Group according to last dimension, cast to strings,
                    # (join if required), reshape.
                    aux_dims = -1 if all(dims) else 0, dims[-1]
                    final_dims = ((0,) if not all(dims) else
                                  dims[:-1] if self._chars_as_strings else
                                  dims)
                    array = np.array(
                        [joiner(map(chr, line))
                         for line in pr.reshape(aux_dims, order="F").tolist()],
                        dtype="U").reshape(final_dims)
                elif matrix_cls == mxSPARSE_CLASS:
                    array = scipy.sparse.csc_matrix((pr, ir, jc), shape=dims)
                else:
                    array = pr.reshape(dims, order="F")

                self._check_and_pad_stream(stream, entry_end)
                yield MatlabArray(name, array, f_global)

            else:
                raise ValueError("Unsupported mdtype: {}".format(mdtype))

    def _check_and_pad_stream(self, stream, entry_end):
        unread = entry_end - stream.tell()
        if unread > 0:
            raise ValueError("{} bytes not read".format(unread))
        elif unread < 0:
            raise ValueError("Over-read {} bytes".format(-unread))
        stream.seek((-stream.tell()) % 8, 1)  # Padding.

    def list_variables(self):
        '''list variables from stream
        '''
        self._prepare_stream()
        infos = []
        for info in self._read_iter(info_only=True):
            infos.append(info[:3])
        return infos

    def get_variables(self, variable_names=None):
        '''get variables from stream as dictionary

        variable_names   - optional list of variable names to get

        If variable_names is None, then get all variables in file
        '''
        if isinstance(variable_names, string_types):
            variable_names = [variable_names]
        self._prepare_stream()
        variables = {"__header__": self._desc,
                     "__globals__": [],  # FIXME Not covered by tests.
                     "__version__": self._version}
        for entry in self._read_iter(load_only=variable_names):
            if entry is None:
                continue
            if not isinstance(entry, MatlabArray):
                raise ValueError("Expected miMATRIX, got {}".format(entry))
            if entry.is_global:
                variables["__global__"].append(entry.name)
            name = entry.name or "__function_workspace__"
            if name in variables:
                warnings.warn(
                    'Duplicate variable name "{}" in stream - replacing '
                    'previous with new.  Consider mio5.varmats_from_mat to '
                    'split file into single variable files'.format(name),
                    MatReadWarning, stacklevel=2)
            variables[name] = (squeeze_element(entry.data)
                               if self._squeeze_me else entry.data)
        return variables

    def get_varmats(self):
        self._prepare_stream()
        infos = []
        for info in self._read_iter(info_only=True):
            info.stream.seek(info.data_position.start)
            raw = info.stream.read(
                info.data_position.stop - info.data_position.start)
            infos.append((info.name, BytesIO(self._header + raw)))
        return infos


def varmats_from_mat(file_obj):
    """ Pull variables out of mat 5 file as a sequence of mat file objects

    This can be useful with a difficult mat file, containing unreadable
    variables.  This routine pulls the variables out in raw form and puts them,
    unread, back into a file stream for saving or reading.  Another use is the
    pathological case where there is more than one variable of the same name in
    the file; this routine returns the duplicates, whereas the standard reader
    will overwrite duplicates in the returned dictionary.

    The file pointer in `file_obj` will be undefined.  File pointers for the
    returned file-like objects are set at 0.

    Parameters
    ----------
    file_obj : file-like
        file object containing mat file

    Returns
    -------
    named_mats : list
        list contains tuples of (name, BytesIO) where BytesIO is a file-like
        object containing mat file contents as for a single variable.  The
        BytesIO contains a string with the original header and a single var. If
        ``var_file_obj`` is an individual BytesIO instance, then save as a mat
        file with something like ``open('test.mat',
        'wb').write(var_file_obj.read())``

    Examples
    --------
    >>> import scipy.io

    BytesIO is from the ``io`` module in python 3, and is ``cStringIO`` for
    python < 3.

    >>> mat_fileobj = BytesIO()
    >>> scipy.io.savemat(mat_fileobj, {'b': np.arange(10), 'a': 'a string'})
    >>> varmats = varmats_from_mat(mat_fileobj)
    >>> sorted([name for name, str_obj in varmats])
    ['a', 'b']
    """
    return MatFile5Reader(file_obj).get_varmats()


def to_writeable(source):
    ''' Convert input object ``source`` to something we can write

    Parameters
    ----------
    source : object

    Returns
    -------
    arr : ndarray

    Examples
    --------
    >>> to_writeable(np.array([1])) # pass through ndarrays
    array([1])
    >>> expected = np.array([(1, 2)], dtype=[('a', '|O8'), ('b', '|O8')])
    >>> np.all(to_writeable({'a':1,'b':2}) == expected)
    True
    >>> np.all(to_writeable({'a':1,'b':2, '_c':3}) == expected)
    True
    >>> np.all(to_writeable({'a':1,'b':2, 100:3}) == expected)
    True
    >>> np.all(to_writeable({'a':1,'b':2, '99':3}) == expected)
    True
    >>> class klass(object): pass
    >>> c = klass
    >>> c.a = 1
    >>> c.b = 2
    >>> np.all(to_writeable({'a':1,'b':2}) == expected)
    True
    >>> to_writeable([])
    array([], dtype=float64)
    >>> to_writeable(())
    array([], dtype=float64)
    >>> to_writeable(None)

    >>> to_writeable('a string').dtype.type == np.str_
    True
    >>> to_writeable(1)
    array(1)
    >>> to_writeable([1])
    array([1])
    >>> to_writeable([1])
    array([1])
    >>> to_writeable(object()) # not convertable

    dict keys with legal characters are convertible

    >>> to_writeable({'a':1})['a']
    array([1], dtype=object)

    but not with illegal characters

    >>> to_writeable({'1':1}) is None
    True
    >>> to_writeable({'_a':1}) is None
    True
    '''
    if isinstance(source, np.ndarray):
        return source
    if source is None:
        return None
    # Objects that implement mappings
    is_mapping = (hasattr(source, 'keys') and hasattr(source, 'values') and
                  hasattr(source, 'items'))
    # Objects that don't implement mappings, but do have dicts
    if not is_mapping and hasattr(source, '__dict__'):
        source = dict((key, value) for key, value in source.__dict__.items()
                      if not key.startswith('_'))
        is_mapping = True
    if is_mapping:
        dtype = []
        values = []
        for field, value in source.items():
            if (isinstance(field, string_types) and
                    field[0] not in '_0123456789'):
                dtype.append((field,object))
                values.append(value)
        if dtype:
            return np.array([tuple(values)],dtype)
        else:
            return None
    # Next try and convert to an array
    narr = np.asanyarray(source)
    if narr.dtype.type in (np.object, np.object_) and \
       narr.shape == () and narr == source:
        # No interesting conversion possible
        return None
    return narr


# Native byte ordered dtypes for convenience for writers
NDT_FILE_HDR = MDTYPES[native_code]['dtypes']['file_header']
NDT_TAG_FULL = MDTYPES[native_code]['dtypes']['tag_full']
NDT_TAG_SMALL = MDTYPES[native_code]['dtypes']['tag_smalldata']
NDT_ARRAY_FLAGS = MDTYPES[native_code]['dtypes']['array_flags']


class VarWriter5(object):
    ''' Generic matlab matrix writing class '''
    mat_tag = np.zeros((), NDT_TAG_FULL)
    mat_tag['mdtype'] = miMATRIX

    def __init__(self, file_writer):
        self.file_stream = file_writer.file_stream
        self.unicode_strings = file_writer.unicode_strings
        self.long_field_names = file_writer.long_field_names
        self.oned_as = file_writer.oned_as
        # These are used for top level writes, and unset after
        self._var_name = None
        self._var_is_global = False

    def write_bytes(self, arr):
        self.file_stream.write(arr.tostring(order='F'))

    def write_string(self, s):
        self.file_stream.write(s)

    def write_element(self, arr, mdtype=None):
        ''' write tag and data '''
        if mdtype is None:
            mdtype = NP_TO_MTYPES[arr.dtype.str[1:]]
        # Array needs to be in native byte order
        if arr.dtype.byteorder == swapped_code:
            arr = arr.byteswap().newbyteorder()
        byte_count = arr.size*arr.itemsize
        if byte_count <= 4:
            self.write_smalldata_element(arr, mdtype, byte_count)
        else:
            self.write_regular_element(arr, mdtype, byte_count)

    def write_smalldata_element(self, arr, mdtype, byte_count):
        # write tag with embedded data
        tag = np.zeros((), NDT_TAG_SMALL)
        tag['byte_count_mdtype'] = (byte_count << 16) + mdtype
        # if arr.tostring is < 4, the element will be zero-padded as needed.
        tag['data'] = arr.tostring(order='F')
        self.write_bytes(tag)

    def write_regular_element(self, arr, mdtype, byte_count):
        # write tag, data
        tag = np.zeros((), NDT_TAG_FULL)
        tag['mdtype'] = mdtype
        tag['byte_count'] = byte_count
        self.write_bytes(tag)
        self.write_bytes(arr)
        # pad to next 64-bit boundary
        bc_mod_8 = byte_count % 8
        if bc_mod_8:
            self.file_stream.write(b'\x00' * (8-bc_mod_8))

    def write_header(self,
                     shape,
                     mclass,
                     is_complex=False,
                     is_logical=False,
                     nzmax=0):
        ''' Write header for given data options
        shape : sequence
           array shape
        mclass      - mat5 matrix class
        is_complex  - True if matrix is complex
        is_logical  - True if matrix is logical
        nzmax        - max non zero elements for sparse arrays

        We get the name and the global flag from the object, and reset
        them to defaults after we've used them
        '''
        # get name and is_global from one-shot object store
        name = self._var_name
        is_global = self._var_is_global
        # initialize the top-level matrix tag, store position
        self._mat_tag_pos = self.file_stream.tell()
        self.write_bytes(self.mat_tag)
        # write array flags (complex, global, logical, class, nzmax)
        af = np.zeros((), NDT_ARRAY_FLAGS)
        af['data_type'] = miUINT32
        af['byte_count'] = 8
        flags = is_complex << 3 | is_global << 2 | is_logical << 1
        af['flags_class'] = mclass | flags << 8
        af['nzmax'] = nzmax
        self.write_bytes(af)
        # shape
        self.write_element(np.array(shape, dtype='i4'))
        # write name
        name = np.asarray(name)
        if name == '':  # empty string zero-terminated
            self.write_smalldata_element(name, miINT8, 0)
        else:
            self.write_element(name, miINT8)
        # reset the one-shot store to defaults
        self._var_name = ''
        self._var_is_global = False

    def update_matrix_tag(self, start_pos):
        curr_pos = self.file_stream.tell()
        self.file_stream.seek(start_pos)
        byte_count = curr_pos - start_pos - 8
        if byte_count >= 2**32:
            raise MatWriteError("Matrix too large to save with Matlab "
                                "5 format")
        self.mat_tag['byte_count'] = byte_count
        self.write_bytes(self.mat_tag)
        self.file_stream.seek(curr_pos)

    def write_top(self, arr, name, is_global):
        """ Write variable at top level of mat file

        Parameters
        ----------
        arr : array-like
            array-like object to create writer for
        name : str, optional
            name as it will appear in matlab workspace
            default is empty string
        is_global : {False, True}, optional
            whether variable will be global on load into matlab
        """
        # these are set before the top-level header write, and unset at
        # the end of the same write, because they do not apply for lower levels
        self._var_is_global = is_global
        self._var_name = name
        # write the header and data
        self.write(arr)

    def write(self, arr):
        ''' Write `arr` to stream at top and sub levels

        Parameters
        ----------
        arr : array-like
            array-like object to create writer for
        '''
        # store position, so we can update the matrix tag
        mat_tag_pos = self.file_stream.tell()
        # First check if these are sparse
        if scipy.sparse.issparse(arr):
            self.write_sparse(arr)
            self.update_matrix_tag(mat_tag_pos)
            return
        # Try to convert things that aren't arrays
        narr = to_writeable(arr)
        if narr is None:
            raise TypeError('Could not convert %s (type %s) to array'
                            % (arr, type(arr)))
        if isinstance(narr, MatlabObject):
            self.write_object(narr)
        elif isinstance(narr, MatlabFunction):
            raise MatWriteError('Cannot write matlab functions')
        elif narr.dtype.fields:  # struct array
            self.write_struct(narr)
        elif narr.dtype.hasobject:  # cell array
            self.write_cells(narr)
        elif narr.dtype.kind in ('U', 'S'):
            if self.unicode_strings:
                codec = 'UTF8'
            else:
                codec = 'ascii'
            self.write_char(narr, codec)
        else:
            self.write_numeric(narr)
        self.update_matrix_tag(mat_tag_pos)

    def write_numeric(self, arr):
        imagf = arr.dtype.kind == 'c'
        logif = arr.dtype.kind == 'b'
        try:
            mclass = NP_TO_MXTYPES[arr.dtype.str[1:]]
        except KeyError:
            # No matching matlab type, probably complex256 / float128 / float96
            # Cast data to complex128 / float64.
            if imagf:
                arr = arr.astype('c128')
            elif logif:
                arr = arr.astype('i1')  # Should only contain 0/1
            else:
                arr = arr.astype('f8')
            mclass = mxDOUBLE_CLASS
        self.write_header(matdims(arr, self.oned_as),
                          mclass,
                          is_complex=imagf,
                          is_logical=logif)
        if imagf:
            self.write_element(arr.real)
            self.write_element(arr.imag)
        else:
            self.write_element(arr)

    def write_char(self, arr, codec='ascii'):
        ''' Write string array `arr` with given `codec`
        '''
        if arr.size == 0 or np.all(arr == ''):
            # This an empty string array or a string array containing
            # only empty strings.  Matlab cannot distiguish between a
            # string array that is empty, and a string array containing
            # only empty strings, because it stores strings as arrays of
            # char.  There is no way of having an array of char that is
            # not empty, but contains an empty string. We have to
            # special-case the array-with-empty-strings because even
            # empty strings have zero padding, which would otherwise
            # appear in matlab as a string with a space.
            shape = (0,) * np.max([arr.ndim, 2])
            self.write_header(shape, mxCHAR_CLASS)
            self.write_smalldata_element(arr, miUTF8, 0)
            return
        # non-empty string.
        #
        # Convert to char array
        arr = arr_to_chars(arr)
        # We have to write the shape directly, because we are going
        # recode the characters, and the resulting stream of chars
        # may have a different length
        shape = arr.shape
        self.write_header(shape, mxCHAR_CLASS)
        if arr.dtype.kind == 'U' and arr.size:
            # Make one long string from all the characters.  We need to
            # transpose here, because we're flattening the array, before
            # we write the bytes.  The bytes have to be written in
            # Fortran order.
            n_chars = np.product(shape)
            st_arr = np.ndarray(shape=(),
                                dtype=arr_dtype_number(arr, n_chars),
                                buffer=arr.T.copy())  # Fortran order
            # Recode with codec to give byte string
            st = st_arr.item().encode(codec)
            # Reconstruct as one-dimensional byte array
            arr = np.ndarray(shape=(len(st),),
                             dtype='S1',
                             buffer=st)
        self.write_element(arr, mdtype=miUTF8)

    def write_sparse(self, arr):
        ''' Sparse matrices are 2D
        '''
        A = arr.tocsc()  # convert to sparse CSC format
        A.sort_indices()     # MATLAB expects sorted row indices
        is_complex = (A.dtype.kind == 'c')
        is_logical = (A.dtype.kind == 'b')
        nz = A.nnz
        self.write_header(matdims(arr, self.oned_as),
                          mxSPARSE_CLASS,
                          is_complex=is_complex,
                          is_logical=is_logical,
                          # matlab won't load file with 0 nzmax
                          nzmax=1 if nz == 0 else nz)
        self.write_element(A.indices.astype('i4'))
        self.write_element(A.indptr.astype('i4'))
        self.write_element(A.data.real)
        if is_complex:
            self.write_element(A.data.imag)

    def write_cells(self, arr):
        self.write_header(matdims(arr, self.oned_as),
                          mxCELL_CLASS)
        # loop over data, column major
        A = np.atleast_2d(arr).flatten('F')
        for el in A:
            self.write(el)

    def write_struct(self, arr):
        self.write_header(matdims(arr, self.oned_as),
                          mxSTRUCT_CLASS)
        self._write_items(arr)

    def _write_items(self, arr):
        # write fieldnames
        fieldnames = [f[0] for f in arr.dtype.descr]
        length = max([len(fieldname) for fieldname in fieldnames])+1
        max_length = (self.long_field_names and 64) or 32
        if length > max_length:
            raise ValueError("Field names are restricted to %d characters" %
                             (max_length-1))
        self.write_element(np.array([length], dtype='i4'))
        self.write_element(
            np.array(fieldnames, dtype='S%d' % (length)),
            mdtype=miINT8)
        A = np.atleast_2d(arr).flatten('F')
        for el in A:
            for f in fieldnames:
                self.write(el[f])

    def write_object(self, arr):
        '''Same as writing structs, except different mx class, and extra
        classname element after header
        '''
        self.write_header(matdims(arr, self.oned_as),
                          mxOBJECT_CLASS)
        self.write_element(np.array(arr.classname, dtype='S'),
                           mdtype=miINT8)
        self._write_items(arr)


class MatFile5Writer(object):
    ''' Class for writing mat5 files '''

    @docfiller
    def __init__(self, file_stream,
                 do_compression=False,
                 unicode_strings=False,
                 global_vars=None,
                 long_field_names=False,
                 oned_as='row'):
        ''' Initialize writer for matlab 5 format files

        Parameters
        ----------
        %(do_compression)s
        %(unicode_strings)s
        global_vars : None or sequence of strings, optional
            Names of variables to be marked as global for matlab
        %(long_fields)s
        %(oned_as)s
        '''
        self.file_stream = file_stream
        self.do_compression = do_compression
        self.unicode_strings = unicode_strings
        if global_vars:
            self.global_vars = global_vars
        else:
            self.global_vars = []
        self.long_field_names = long_field_names
        self.oned_as = oned_as
        self._matrix_writer = None

    def write_file_header(self):
        # write header
        hdr = np.zeros((), NDT_FILE_HDR)
        hdr['description'] = 'MATLAB 5.0 MAT-file Platform: %s, Created on: %s' \
            % (os.name,time.asctime())
        hdr['version'] = 0x0100
        hdr['endian_test'] = np.ndarray(shape=(),
                                      dtype='S2',
                                      buffer=np.uint16(0x4d49))
        self.file_stream.write(hdr.tostring())

    def put_variables(self, mdict, write_header=None):
        ''' Write variables in `mdict` to stream

        Parameters
        ----------
        mdict : mapping
           mapping with method ``items`` returns name, contents pairs where
           ``name`` which will appear in the matlab workspace in file load, and
           ``contents`` is something writeable to a matlab file, such as a numpy
           array.
        write_header : {None, True, False}
           If True, then write the matlab file header before writing the
           variables.  If None (the default) then write the file header
           if we are at position 0 in the stream.  By setting False
           here, and setting the stream position to the end of the file,
           you can append variables to a matlab file
        '''
        # write header if requested, or None and start of file
        if write_header is None:
            write_header = self.file_stream.tell() == 0
        if write_header:
            self.write_file_header()
        self._matrix_writer = VarWriter5(self)
        for name, var in mdict.items():
            if name[0] == '_':
                continue
            is_global = name in self.global_vars
            if self.do_compression:
                stream = BytesIO()
                self._matrix_writer.file_stream = stream
                self._matrix_writer.write_top(var, asbytes(name), is_global)
                out_str = zlib.compress(stream.getvalue())
                tag = np.empty((), NDT_TAG_FULL)
                tag['mdtype'] = miCOMPRESSED
                tag['byte_count'] = len(out_str)
                self.file_stream.write(tag.tostring())
                self.file_stream.write(out_str)
            else:  # not compressing
                self._matrix_writer.write_top(var, asbytes(name), is_global)
