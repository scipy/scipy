#ifndef NI_ITERATORS_H
#define NI_ITERATORS_H

#define NO_IMPORT_ARRAY
#include "nd_image.h"
#undef NO_IMPORT_ARRAY


/*
 ***********************************************************************
 ***                           NDIterator                            ***
 ***********************************************************************
 */

/* The struct internals are a private implementation detail. */
typedef struct NDIterator_Internal NDIterator;

/* Creates a new iterator over all axes of the array. */
NDIterator*
NDI_New(PyArrayObject *array);

/* Creates a new iterator over all except one axis of the array. */
NDIterator*
NDI_NewExceptAxis(PyArrayObject *array, int axis);

/* Advances the iterator. Returns 0 if exhausted, else 1. */
int
NDI_Next(NDIterator *nditer);

/* Releases all memory allocated by an iterator */
void
NDI_Delete(NDIterator *nditer);


/*
 ***********************************************************************
 ***                       LineBufferIterator                        ***
 ***********************************************************************
 */

/* The struct internals are a private implementation detail. */
typedef struct LineBufferIterator_Internal LineBufferIterator;

/* Creates a new line buffer iterator and reads the first line in. */
LineBufferIterator*
LBI_New(PyArrayObject *input, int axis, PyArrayObject *output,
        npy_intp filter_size, npy_intp filter_offset,
        NI_ExtendMode extend_mode, npy_double extend_value);

/* Returns a pointer to the iterator's input buffer. */
npy_double*
LBI_GetInputBuffer(LineBufferIterator *lbiter);

/* Returns a pointer to the iterator's output buffer. */
npy_double*
LBI_GetOutputBuffer(LineBufferIterator *lbiter);

/* Returns the length of each of the iterator's lines. */
npy_intp
LBI_GetLineLength(LineBufferIterator *lbiter);

/*
 * Writes the output buffer to the output array, advances the iterator and
 * reads a new line into the input buffer. Returns 0 if exhausted, else 1.
 */
int
LBI_Next(LineBufferIterator *lbiter);

/* Releases all memory allocated by the iterator. Always returns NULL. */
LineBufferIterator*
LBI_Delete(LineBufferIterator *lbiter);


#endif /* NI_ITERATORS_H */
