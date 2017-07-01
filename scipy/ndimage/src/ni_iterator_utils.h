/* Helper functionality used in ni_iterators. */

#ifndef NI_ITERATOR_UTILS_H
#define NI_ITERATOR_UTILS_H


/*
 ***********************************************************************
 ***                         Generic helpers                         ***
 ***********************************************************************
 */

static int
_is_input_array(PyArrayObject *array)
{
    if (array == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Unexpected NULL array.");
        return 0;
    }
    if (!PyArray_ISALIGNED(array)) {
        PyErr_SetString(PyExc_RuntimeError, "Expected an aligned array.");
        return 0;
    }
    if (PyArray_ISBYTESWAPPED(array)) {
        PyErr_SetString(PyExc_RuntimeError, "Unexpected byte swapped array.");
        return 0;
    }
    if (PyArray_SIZE(array) == 0) {
        PyErr_SetString(PyExc_RuntimeError, "Unexpected zero-sized array.");
        return 0;
    }
    return 1;
}


static int
_is_output_array(PyArrayObject *array)
{
    if (!_is_input_array(array)) {
        return 0;
    }
    if (!PyArray_ISWRITEABLE(array)) {
        PyErr_SetString(PyExc_RuntimeError, "Expected a writeable array.");
        return 0;
    }
    return 1;
}


static int
_is_valid_axis(PyArrayObject *array, int *axis)
{
    if (*axis < -PyArray_NDIM(array) || *axis >= PyArray_NDIM(array)) {
        PyErr_Format(PyExc_ValueError,
                     "axis %d out of bounds for %d-dimensional array.",
                     *axis, PyArray_NDIM(array));
        return 0;
    }
    if (*axis < 0) {
        *axis += PyArray_NDIM(array);
    }
    return 1;
}


static int
_intp_arrays_are_equal(npy_intp *one, npy_intp *two, npy_intp length)
{
    while (length--) {
        if (*one++ != *two++) {
            return 0;
        }
    }
    return 1;
}


static int
_arrays_are_compatible(PyArrayObject *one, PyArrayObject *two)
{
    if (PyArray_NDIM(one) != PyArray_NDIM(two)) {
        PyErr_SetString(PyExc_RuntimeError,
                        "Expected arrays with same number of dimensions.");
        return 0;
    }
    if (!_intp_arrays_are_equal(PyArray_DIMS(one), PyArray_DIMS(two),
                                PyArray_NDIM(one))) {
        PyErr_SetString(PyExc_RuntimeError,
                        "Expected arrays with same shape.");
        return 0;
    }
    return 1;
}


static int
_split_filter_size(npy_intp size, npy_intp offset,
                   npy_intp *before, npy_intp *after)
{
    if (size < 1) {
        PyErr_SetString(PyExc_RuntimeError,
                        "Unexpected zero-sized filter.");
        return 0;
    }
    *before = size / 2 + offset;
    *after = size - *before - 1;
    if (*before < 0 || after < 0) {
        PyErr_Format(PyExc_ValueError,
                     "offset %zd out of bounds for filter of size %zd.",
                     offset, size);
        return 0;
    }
    return 1;
}


static int
_is_lbi_supported_type(enum NPY_TYPES numtype)
{
    if (numtype < NPY_BOOL || numtype > NPY_DOUBLE) {
        PyErr_SetString(PyExc_RuntimeError, "Unsupported data type.");
        return 0;
    }
    return 1;
}


/*
 ***********************************************************************
 ***                     Helpers for LBI_ReadLine                    ***
 ***********************************************************************
 */

typedef void (read_line_func)(npy_double*, const char*, npy_intp, npy_intp);

#define READ_LINE_FUNC(_type)                                  \
static void                                                    \
_read_line_##_type(npy_double *buffer, const char *line,       \
                   npy_intp line_length, npy_intp line_stride) \
{                                                              \
    while (line_length--) {                                    \
        *buffer++ = *(_type *)line;                            \
        line += line_stride;                                   \
    }                                                          \
}

READ_LINE_FUNC(npy_bool)
READ_LINE_FUNC(npy_byte)
READ_LINE_FUNC(npy_ubyte)
READ_LINE_FUNC(npy_short)
READ_LINE_FUNC(npy_ushort)
READ_LINE_FUNC(npy_int)
READ_LINE_FUNC(npy_uint)
READ_LINE_FUNC(npy_long)
READ_LINE_FUNC(npy_ulong)
READ_LINE_FUNC(npy_longlong)
READ_LINE_FUNC(npy_ulonglong)
READ_LINE_FUNC(npy_float)
READ_LINE_FUNC(npy_double)

static read_line_func*
_get_read_line_func(enum NPY_TYPES typenum)
{
    /* These functions must be in the same order as enum NPY_TYPES. */
    read_line_func *_lbi_read_funcs[] = {
        &_read_line_npy_bool,
        &_read_line_npy_byte, &_read_line_npy_ubyte,
        &_read_line_npy_short, &_read_line_npy_ushort,
        &_read_line_npy_int, &_read_line_npy_uint,
        &_read_line_npy_long, &_read_line_npy_ulong,
        &_read_line_npy_longlong, &_read_line_npy_ulonglong,
        &_read_line_npy_float, &_read_line_npy_double,
    };
    if (_is_lbi_supported_type(typenum)) {
        return _lbi_read_funcs[typenum];
    }
    return NULL;
}


/*
 ***********************************************************************
 ***                    Helpers for LBI_WriteLine                    ***
 ***********************************************************************
 */


typedef void (write_line_func)(npy_double*, char*, npy_intp, npy_intp);

#define WRITE_LINE_FUNC(_type)                                  \
static void                                                     \
_write_line_##_type(npy_double *buffer, char *line,             \
                    npy_intp line_length, npy_intp line_stride) \
{                                                               \
    while (line_length--) {                                     \
        *(_type *)line = *buffer++;                             \
        line += line_stride;                                    \
    }                                                           \
}

WRITE_LINE_FUNC(npy_bool)
WRITE_LINE_FUNC(npy_ubyte)
WRITE_LINE_FUNC(npy_ushort)
WRITE_LINE_FUNC(npy_uint)
WRITE_LINE_FUNC(npy_ulong)
WRITE_LINE_FUNC(npy_ulonglong)
WRITE_LINE_FUNC(npy_byte)
WRITE_LINE_FUNC(npy_short)
WRITE_LINE_FUNC(npy_int)
WRITE_LINE_FUNC(npy_long)
WRITE_LINE_FUNC(npy_longlong)
WRITE_LINE_FUNC(npy_float)
WRITE_LINE_FUNC(npy_double)


static NPY_INLINE write_line_func*
_get_write_line_func(enum NPY_TYPES typenum)
{
    /* These functions must be in the same order as enum NPY_TYPES. */
    write_line_func *_lbi_write_funcs[] = {
        &_write_line_npy_bool,
        &_write_line_npy_byte, &_write_line_npy_ubyte,
        &_write_line_npy_short, &_write_line_npy_ushort,
        &_write_line_npy_int, &_write_line_npy_uint,
        &_write_line_npy_long, &_write_line_npy_ulong,
        &_write_line_npy_longlong, &_write_line_npy_ulonglong,
        &_write_line_npy_float, &_write_line_npy_double,
    };
    if (_is_lbi_supported_type(typenum)) {
        return _lbi_write_funcs[typenum];
    }
    return NULL;
}


/*
 ***********************************************************************
 ***                    Helpers for LBI_ExtendLine                   ***
 ***********************************************************************
 */

typedef void (extend_line_func)(npy_double*, npy_intp, npy_intp, npy_intp,
                                npy_double);

static void
_extend_line_nearest(npy_double *buffer, npy_intp line_length,
                     npy_intp size_before, npy_intp size_after,
                     npy_double unused_value)
{
    /* aaaaaaaa|abcd|dddddddd */
    npy_double * const first = buffer + size_before;
    npy_double * const last = first + line_length;
    npy_double *dst, val;

    val = *first;
    dst = buffer;
    while (size_before--) {
        *dst++ = val;
    }
    dst = last;
    val = *(last - 1);
    while (size_after--) {
        *dst++ = val;
    }
}

static void
_extend_line_wrap(npy_double *buffer, npy_intp line_length,
                  npy_intp size_before, npy_intp size_after,
                  npy_double unused_value)
{
    /* abcdabcd|abcd|abcdabcd */
    npy_double * const first = buffer + size_before;
    npy_double * const last = first + line_length;
    const npy_double *src;
    npy_double *dst;

    src = last - 1;
    dst = first - 1;
    while (size_before--) {
        *dst-- = *src--;
    }
    src = first;
    dst = last;
    while (size_after--) {
        *dst++ = *src++;
    }
}

static void
_extend_line_reflect(npy_double *buffer, npy_intp line_length,
                     npy_intp size_before, npy_intp size_after,
                     npy_double unused_value)
{
    /* abcddcba|abcd|dcbaabcd */
    npy_double * const first = buffer + size_before;
    npy_double * const last = first + line_length;
    const npy_double *src;
    npy_double *dst;

    src = first;
    dst = first - 1;
    while (src < last && size_before--) {
        *dst-- = *src++;
    }
    src = last - 1;
    while (size_before--) {
        *dst-- = *src--;
    }
    src = last - 1;
    dst = last;
    while (dst >= first && size_after--) {
        *dst++ = *src--;
    }
    src = first;
    while (size_after--) {
        *dst++ = *src++;
    }
}

static void
_extend_line_mirror(npy_double *buffer, npy_intp line_length,
                    npy_intp size_before, npy_intp size_after,
                    npy_double unused_value)
{
    /* abcddcba|abcd|dcbaabcd */
    npy_double * const first = buffer + size_before;
    npy_double * const last = first + line_length;
    const npy_double *src;
    npy_double *dst;

    src = first + 1;
    dst = first - 1;
    while (src < last && size_before--) {
        *dst-- = *src++;
    }
    src = last - 2;
    while (size_before--) {
        *dst-- = *src--;
    }
    src = last - 2;
    dst = last;
    while (dst >= first && size_after--) {
        *dst++ = *src--;
    }
    src = first + 1;
    while (size_after--) {
        *dst++ = *src++;
    }
}

static void
_extend_line_constant(npy_double *buffer, npy_intp line_length,
                      npy_intp size_before, npy_intp size_after,
                      const npy_double value)
{
        /* abcddcba|abcd|dcbaabcd */
    npy_double * const last = buffer + size_before + line_length;
    npy_double *dst;

    dst = buffer;
    while (size_before--) {
        *dst++ = value;
    }
    dst = last;
    while (size_after--) {
        *dst++ = value;
    }
}

static extend_line_func*
_get_extend_line_func(NI_ExtendMode extend_mode)
{
    switch (extend_mode) {
        case NI_EXTEND_NEAREST:
            return &_extend_line_nearest;
        case NI_EXTEND_WRAP:
            return &_extend_line_wrap;
        case NI_EXTEND_REFLECT:
            return &_extend_line_reflect;
        case NI_EXTEND_MIRROR:
            return &_extend_line_mirror;
        case NI_EXTEND_CONSTANT:
            return &_extend_line_constant;
        default:
            PyErr_Format(PyExc_ValueError,
                         "Unexpected extend mode %d.", extend_mode);
            return NULL;
    }
}


#endif /* NI_ITERATOR_UTILS_H */
