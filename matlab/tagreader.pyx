# -*- python -*- 
''' Extension to parse matlab 5 tags '''

# Import the pieces of the Python C API we need to use (from c_python.pxd):
cimport c_python as py

def parse(fileobj, int swapf):
    ''' Read in the tag
    The tag can be normal format (mdtype=u4, byte_count=u4)
    or small element format (mdtype=u2, byte_count=u2, data in last 4 bytes)
    Small element format is where mdtype (u4) has non-zero high bytes
    '''
    cdef py.size_t n_out
    cdef char raw_tag[8]
    cdef py.FILE* infile
    infile = py.PyFile_AsFile(fileobj)
    n_out = py.fread(raw_tag, 8, 1, infile)
    # Raise Exception if n_out < 1
    return mdtype, byte_count, buf
