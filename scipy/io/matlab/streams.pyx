# -*- python -*- or near enough

cimport python_string

cdef extern from "stdlib.h" nogil:
    void *malloc(size_t size)
    void *memcpy(void *str1, void *str2, size_t n)
    void free(void *ptr)

cdef extern from "Python.h":
    void *PyCObject_Import(char *, char *) except NULL
    ctypedef struct PyTypeObject:
        pass
    ctypedef struct PyObject:
        pass
    ctypedef struct FILE
    FILE *PyFile_AsFile(object)
    size_t fread (void *ptr, size_t size, size_t n, FILE* fptr)
    int fseek (FILE * fptr, long int offset, int whence)
    long int ftell (FILE *stream)
    object PyCObject_FromVoidPtr(void *, void (*)(void*))

cdef extern from "fileobject.h":
    ctypedef class __builtin__.file [object PyFileObject]:
        pass

cdef extern from "cStringIO.h":
    # From:
    # http://svn.pyamf.org/pyamf/tags/release-0.4rc2/cpyamf/util.pyx
    # (MIT license) - with thanks
    void PycString_IMPORT()
    object StringIO_NewOutput "PycStringIO->NewOutput" (int)
    object StringIO_NewInput "PycStringIO->NewInput" (object)
    int StringIO_cread "PycStringIO->cread" (object, char **, Py_ssize_t)
    int StringIO_creadline "PycStringIO->creadline" (object, char **)
    int StringIO_cwrite "PycStringIO->cwrite" (object, char *, Py_ssize_t)
    object StringIO_cgetvalue "PycStringIO->cgetvalue" (obj)
    bint PycStringIO_InputCheck(object O)
    bint PycStringIO_OutputCheck(object O)
       
# initialize cStringIO
PycString_IMPORT


cdef class Memholder:
    ''' Object to hold memory pointer, and dealloc on delete '''
    def __dealloc__(self):
       if self.ptr !=NULL:
           free(self.ptr)


cdef class GenericStream:

    def __init__(self, fobj):
        self.fobj = fobj

    cpdef int seek(self, long int offset, int whence=0) except -1:
        self.fobj.seek(offset, whence)
        return 0
        
    cpdef long int tell(self) except -1:
        return self.fobj.tell()

    def read(self, n_bytes):
        return self.fobj.read(n_bytes)

    cdef int read_into(self, void *buf, size_t n) except -1:
        ''' Read n bytes from stream into pre-allocated buffer `buf`
        '''
        cdef char* d_ptr
        data = self.fobj.read(n)
        if python_string.PyString_Size(data) != n:
            raise IOError('could not read bytes')
            return -1
        d_ptr = data
        memcpy(buf, d_ptr, n)
        return 0

    cdef Memholder read_alloc(self, size_t n):
        ''' Make new memory, wrap with object '''
        cdef Memholder mh = Memholder()
        cdef char *d_ptr
        data = self.fobj.read(n)
        if python_string.PyString_Size(data) != n:
            raise IOError('could not read bytes')
        d_ptr = data
        mh.ptr = malloc(n)
        if mh.ptr == NULL:
            raise MemoryError('Could not allocate memory')
        memcpy(mh.ptr, d_ptr, n)
        return mh


cdef class cStringStream(GenericStream):
    
    cpdef int seek(self, long int offset, int whence=0) except -1:
        cdef char *ptr
        if whence == 1 and offset >=0: # forward, from here
            StringIO_cread(self.fobj, &ptr, offset)
            return 0
        else: # use python interface
            return GenericStream.seek(self, offset, whence)

    cdef int read_into(self, void *buf, size_t n) except -1:
        ''' Read n bytes from stream into pre-allocated buffer `buf`
        '''
        cdef:
            size_t n_red
            char* d_ptr
        n_red = StringIO_cread(self.fobj, &d_ptr, n)
        if n_red != n:
            raise IOError('could not read bytes')
        memcpy(buf, <void *>d_ptr, n)
        return 0

    cdef Memholder read_alloc(self, size_t n):
        ''' Make new memory, wrap with object '''
        cdef Memholder mh = Memholder()
        cdef char *d_ptr
        cdef size_t n_red
        n_red = StringIO_cread(self.fobj, &d_ptr, n)
        if n_red != n:
            raise IOError('could not read bytes')
        mh.ptr = malloc(n)
        if mh.ptr == NULL:
            raise MemoryError('Could not allocate memory')
        memcpy(mh.ptr, d_ptr, n)
        return mh
   
cdef class FileStream(GenericStream):
    cdef FILE* file

    def __init__(self, fobj):
        self.fobj = fobj
        self.file = PyFile_AsFile(fobj)
        
    cpdef int seek(self, long int offset, int whence=0) except -1:
        cdef int ret
        ''' move `offset` bytes in stream

        Parameters
        ----------
        offset : long int
           number of bytes to move.  Positive for forward in file,
           negative for backward
        whence : int
           `whence` can be:
           
           * 0 - from beginning of file (`offset` should be >=0)
           * 1 - from current file position
           * 2 - from end of file (`offset nearly always <=0)

        Returns
        -------
        ret : int
        '''
        ret = fseek(self.file, offset, whence)
        if ret:
            raise IOError('Failed seek')
            return -1
        return ret

    cpdef long int tell(self):
        return ftell(self.file)

    cdef int read_into(self, void *buf, size_t n) except -1:
        ''' Read n bytes from stream into pre-allocated buffer `buf`
        '''
        cdef:
            size_t n_red
            char* d_ptr
        n_red = fread(buf, 1, n, self.file)
        if n_red != n:
            raise IOError('Could not read bytes')
            return -1
        return 0

    cdef Memholder read_alloc(self, size_t n):
        ''' Make new memory, wrap with object '''
        cdef Memholder mh = Memholder()
        cdef char *d_ptr
        cdef size_t n_red
        mh.ptr = malloc(n)
        if mh.ptr == NULL:
            raise MemoryError('Could not allocate memory')
        n_red = fread(mh.ptr, 1, n, self.file)
        if n_red != n:
            raise IOError('could not read bytes')
        return mh


def _read_into(GenericStream st, size_t n):
    # for testing only.  Use st.read instead
    cdef char * d_ptr
    my_str = ' ' * n
    d_ptr = my_str
    st.read_into(d_ptr, n)
    return my_str


def _read_alloc(GenericStream st, size_t n):
    # for testing only.  Use st.read instead
    cdef Memholder obj
    cdef char **pp
    cdef char *d_ptr
    obj = st.read_alloc(n)
    my_str = 'A' * n
    d_ptr = my_str
    memcpy(d_ptr, obj.ptr, n)
    return my_str

    
cpdef GenericStream make_stream(object fobj):
    ''' Make stream of correct type for file-like `fobj`
    '''
    if isinstance(fobj, file):
        return FileStream(fobj)
    elif PycStringIO_InputCheck(fobj) or PycStringIO_OutputCheck(fobj):
        return cStringStream(fobj)
    return GenericStream(fobj)


