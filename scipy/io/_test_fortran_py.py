import numpy as np

def read_unformatted_double(m, n, k, filename):
    """
    Read a Fortran-style unformatted binary file written with a single write() call,
    assuming it wraps the data with 4-byte record markers.

    Returns:
        np.ndarray of shape (m, n, k) with dtype float64
    """
    with open(filename.strip(), 'rb') as f:
        f.read(4)  # Skip initial 4-byte record marker
        data = np.fromfile(f, dtype=np.float64, count=m * n * k)
        f.read(4)  # Skip trailing 4-byte record marker

    if data.size != m * n * k:
        raise ValueError(f"Expected {m*n*k} elements, got {data.size}")

    return data.reshape((m, n, k), order='F')  # Fortran column-major order

def read_unformatted_mixed(m, n, k, filename):
    """
    Read a Fortran unformatted binary file that contains a mix of:
    - a double precision array a(m, n)
    - an integer array b(k)

    Assumes a single write(10) a, b was used and file is wrapped
    with Fortran record markers.

    Returns:
        a: np.ndarray of shape (m, n) with dtype float64
        b: np.ndarray of shape (k,) with dtype int32
    """
    with open(filename.strip(), 'rb') as f:
        f.read(4)  # Skip initial 4-byte record marker

        # Read a(m,n): total m*n float64 values
        a_flat = np.fromfile(f, dtype=np.float64, count=m * n)

        # Read b(k): total k int32 values (assuming Fortran default integer*4)
        b = np.fromfile(f, dtype=np.int32, count=k)

        f.read(4)  # Skip trailing 4-byte record marker

    # Reshape a to (m,n) Fortran-style
    a = a_flat.reshape((m, n), order='F')

    return a, b

def read_unformatted_int(m, n, k, filename):
    """
    Read a Fortran unformatted binary file
    containing a 3D integer array (m, n, k).
    Assumes the array is written with a single
    write(10) a and wrapped with record markers.

    Returns:
        np.ndarray: 3D array of shape (m, n, k) with dtype int32
    """
    with open(filename.strip(), 'rb') as f:
        f.read(4)  # Skip Fortran record marker at start

        # Read m*n*k integers (Fortran default = 4 bytes per integer)
        data = np.fromfile(f, dtype=np.int32, count=m * n * k)

        f.read(4)  # Skip Fortran record marker at end

    if data.size != m * n * k:
        raise ValueError(f"Expected {m*n*k} elements, got {data.size}")

    return data.reshape((m, n, k), order='F')  # Fortran-style column-major order
