import numpy as np
import scipy as sp

__all__ = ['save_npz', 'load_npz', 'from_binsparse']


# Make loading safe vs. malicious input
PICKLE_KWARGS = dict(allow_pickle=False)


def save_npz(file, matrix, compressed=True):
    """ Save a sparse matrix or array to a file using ``.npz`` format.

    Parameters
    ----------
    file : str or file-like object
        Either the file name (string) or an open file (file-like object)
        where the data will be saved. If file is a string, the ``.npz``
        extension will be appended to the file name if it is not already
        there.
    matrix: spmatrix or sparray
        The sparse matrix or array to save.
        Supported formats: ``csc``, ``csr``, ``bsr``, ``dia`` or ``coo``.
    compressed : bool, optional
        Allow compressing the file. Default: True

    See Also
    --------
    scipy.sparse.load_npz: Load a sparse matrix from a file using ``.npz`` format.
    numpy.savez: Save several arrays into a ``.npz`` archive.
    numpy.savez_compressed : Save several arrays into a compressed ``.npz`` archive.

    Examples
    --------
    Store sparse matrix to disk, and load it again:

    >>> import numpy as np
    >>> import scipy as sp
    >>> sparse_matrix = sp.sparse.csc_matrix([[0, 0, 3], [4, 0, 0]])
    >>> sparse_matrix
    <Compressed Sparse Column sparse matrix of dtype 'int64'
        with 2 stored elements and shape (2, 3)>
    >>> sparse_matrix.toarray()
    array([[0, 0, 3],
           [4, 0, 0]], dtype=int64)

    >>> sp.sparse.save_npz('/tmp/sparse_matrix.npz', sparse_matrix)
    >>> sparse_matrix = sp.sparse.load_npz('/tmp/sparse_matrix.npz')

    >>> sparse_matrix
    <Compressed Sparse Column sparse matrix of dtype 'int64'
        with 2 stored elements and shape (2, 3)>
    >>> sparse_matrix.toarray()
    array([[0, 0, 3],
           [4, 0, 0]], dtype=int64)
    """
    arrays_dict = {}
    if matrix.format in ('csc', 'csr', 'bsr'):
        arrays_dict.update(indices=matrix.indices, indptr=matrix.indptr)
    elif matrix.format == 'dia':
        arrays_dict.update(offsets=matrix.offsets)
    elif matrix.format == 'coo':
        arrays_dict.update(row=matrix.row, col=matrix.col)
    else:
        msg = f'Save is not implemented for sparse matrix of format {matrix.format}.'
        raise NotImplementedError(msg)
    arrays_dict.update(
        format=matrix.format.encode('ascii'),
        shape=matrix.shape,
        data=matrix.data
    )
    if isinstance(matrix, sp.sparse.sparray):
        arrays_dict.update(_is_array=True)
    if compressed:
        np.savez_compressed(file, **arrays_dict)
    else:
        np.savez(file, **arrays_dict)


def load_npz(file):
    """ Load a sparse array/matrix from a file using ``.npz`` format.

    Parameters
    ----------
    file : str or file-like object
        Either the file name (string) or an open file (file-like object)
        where the data will be loaded.

    Returns
    -------
    result : csc_array, csr_array, bsr_array, dia_array or coo_array
        A sparse array/matrix containing the loaded data.

    Raises
    ------
    OSError
        If the input file does not exist or cannot be read.

    See Also
    --------
    scipy.sparse.save_npz: Save a sparse array/matrix to a file using ``.npz`` format.
    numpy.load: Load several arrays from a ``.npz`` archive.

    Examples
    --------
    Store sparse array/matrix to disk, and load it again:

    >>> import numpy as np
    >>> import scipy as sp
    >>> sparse_array = sp.sparse.csc_array([[0, 0, 3], [4, 0, 0]])
    >>> sparse_array
    <Compressed Sparse Column sparse array of dtype 'int64'
        with 2 stored elements and shape (2, 3)>
    >>> sparse_array.toarray()
    array([[0, 0, 3],
           [4, 0, 0]], dtype=int64)

    >>> sp.sparse.save_npz('/tmp/sparse_array.npz', sparse_array)
    >>> sparse_array = sp.sparse.load_npz('/tmp/sparse_array.npz')

    >>> sparse_array
    <Compressed Sparse Column sparse array of dtype 'int64'
        with 2 stored elements and shape (2, 3)>
    >>> sparse_array.toarray()
    array([[0, 0, 3],
           [4, 0, 0]], dtype=int64)

    In this example we force the result to be csr_array from csr_matrix
    >>> sparse_matrix = sp.sparse.csc_matrix([[0, 0, 3], [4, 0, 0]])
    >>> sp.sparse.save_npz('/tmp/sparse_matrix.npz', sparse_matrix)
    >>> tmp = sp.sparse.load_npz('/tmp/sparse_matrix.npz')
    >>> sparse_array = sp.sparse.csr_array(tmp)
    """
    with np.load(file, **PICKLE_KWARGS) as loaded:
        sparse_format = loaded.get('format')
        if sparse_format is None:
            raise ValueError(f'The file {file} does not contain '
                             f'a sparse array or matrix.')
        sparse_format = sparse_format.item()

        if not isinstance(sparse_format, str):
            # Play safe with Python 2 vs 3 backward compatibility;
            # files saved with SciPy < 1.0.0 may contain unicode or bytes.
            sparse_format = sparse_format.decode('ascii')

        if loaded.get('_is_array'):
            sparse_type = sparse_format + '_array'
        else:
            sparse_type = sparse_format + '_matrix'

        try:
            cls = getattr(sp.sparse, f'{sparse_type}')
        except AttributeError as e:
            raise ValueError(f'Unknown format "{sparse_type}"') from e

        if sparse_format in ('csc', 'csr', 'bsr'):
            return cls((loaded['data'], loaded['indices'], loaded['indptr']),
                       shape=loaded['shape'])
        elif sparse_format == 'dia':
            return cls((loaded['data'], loaded['offsets']),
                       shape=loaded['shape'])
        elif sparse_format == 'coo':
            return cls((loaded['data'], (loaded['row'], loaded['col'])),
                       shape=loaded['shape'])
        else:
            raise NotImplementedError(f'Load is not implemented for '
                                      f'sparse matrix of format {sparse_format}.')

def from_binsparse(arr, /, *, device=None, copy: bool | None = None):
    from scipy.sparse import csr_array, csc_array
    desc = arr.__binsparse_descriptor__()
    arrs = arr.__binsparse__()

    desc = desc["binsparse"]
    version_tuple: tuple[int, ...] = tuple(int(v) for v in desc["version"].split("."))
    if version_tuple != (0, 1):
        raise RuntimeError("Unsupported `__binsparse__` protocol version.")

    format = desc["format"]
    format_err_str = f"Unsupported format: `{format!r}`."

    if isinstance(format, str):
        match format:
            case "CSC" | "CSR":
                desc["format"] = {
                    "custom": {
                        "transpose": [0, 1] if format == "CSR" else [0, 1],
                        "level": {
                            "level_desc": "dense",
                            "level": {
                                "level_desc": "sparse",
                                "level": {
                                    "level_desc": "element",
                                },
                            },
                        },
                    },
                }
            case _:
                raise RuntimeError(format_err_str)

    format = desc["format"]["custom"]
    rank = 0
    level = format
    while "level" in level:
        if "rank" not in level:
            level["rank"] = 1
        rank += level["rank"]
        level = level["level"]
    if "transpose" not in format:
        format["transpose"] = list(range(rank))

    match desc:
        case {
            "format": {
                "custom": {
                    "transpose": transpose,
                    "level": {
                        "level_desc": "dense",
                        "rank": 1,
                        "level": {
                            "level_desc": "sparse",
                            "rank": 1,
                            "level": {
                                "level_desc": "element",
                            },
                        },
                    },
                },
            },
            "shape": shape,
            "data_types": {
                "pointers_to_1": ptr_dtype,
                "indices_1": crd_dtype,
                "values": val_dtype,
            },
            **_kwargs,
        }:
            crd_arr = np.from_dlpack(arrs["pointers_to_1"])
            _check_binsparse_dt(crd_arr, crd_dtype)
            ptr_arr = np.from_dlpack(arrs["indices_1"])
            _check_binsparse_dt(ptr_arr, ptr_dtype)
            val_arr = np.from_dlpack(arrs["values"])
            _check_binsparse_dt(val_arr, val_dtype)

            match transpose:
                case [0, 1]:
                    sparse_type = csr_array
                case [1, 0]:
                    sparse_type = csc_array
                case _:
                    raise RuntimeError(format_err_str)

            return sparse_type((val_arr, ptr_arr, crd_arr), shape=shape)
        case _:
            raise RuntimeError(format_err_str)

def _convert_binsparse_dtype(dt: str) -> np.dtype:
    if dt.startswith("complex[float") and dt.endswith("]"):
        complex_bits = 2 * int(dt[len("complex[float") : -len("]")])
        dt: str = f"complex{complex_bits}"

    return np.dtype(dt)


def _check_binsparse_dt(arr: np.ndarray, dt: str) -> None:
    invalid_dtype_str = "Invalid dtype: `{dtype!s}`, expected `{expected!s}`."
    dt = _convert_binsparse_dtype(dt)
    if dt != arr.dtype:
        raise BufferError(
            invalid_dtype_str.format(
                dtype=arr.dtype,
                expected=dt,
            )
        )
