"""
Helper function for checking signatures created by:
 - _cython_blas_signatures.py
 - _cython_lapack_signatures.py

It does so by trying to link the signatures against 
linker options passed.
If the signature succeeds a linking step it is retained, 
otherwise it is removed from the list.
"""
from __future__ import division, print_function, absolute_import

import distutils.ccompiler
import tempfile
import shutil
from os.path import join

__all__ = ['get_linking_signatures']

def get_linking_signatures(signature_file, lib_info):
    """ Reads signature_file and loops on all entries checking for library

    All lines that looks like a function in `signature_file` will be tested
    against the lib_info data to determine whether it exist.
    
    Returns a string containing all the signatures "as-is" which exists.
    """

    # Create compiler
    c = distutils.ccompiler.new_compiler()

    # Create used directories and file names
    tmpdir = tempfile.mkdtemp()
    src = join(tmpdir, 'source.c')
    out = join(tmpdir, 'a.out')
    s = """void {0}();
    int main(int argc, const char *argv[])
    {{
      {0}{1}();
      return 0;
    }}"""
    
    # Add the additional "extra" arguments if possible
    try:
        extra_args = lib_info['extra_link_args']
    except:
        extra_args = []

    # Local function for actual linking
    def check_sig_line(line, info, f_suffix='_'):
        ret = False
        # Retrieve signature routine name
        try:
            sig = line.split('(')[0].split(' ')[-1]
        except:
            return ret
        with open(src, 'wt') as f:
            f.write(s.format(sig, f_suffix))
        obj = c.compile([src], output_dir=tmpdir)
        try:
            c.link_executable(obj, out, libraries=info['libraries'],
                              library_dirs=info['library_dirs'],
                              extra_postargs=extra_args)
            ret = True
        except:
            pass
        # Clean-up the library, this will ensure no overlapping time-stamps
        try:
            os.remove(out)
        except:
            pass
        return ret

    # Returned signature list
    ret_sig = []

    # Loop signatures in file 'signature_file'
    with open(signature_file,'r') as f:
        for line in f:
            if line.startswith('#') or len(line.strip()) == 0:
                # Simply append the comment or empty line
                ret_sig += [line]
            elif check_sig_line(line, lib_info):
                ret_sig += [line]
    
    # Clean-up the temporary folder
    shutil.rmtree(tmpdir)

    # Return list of signatures _as-written_
    return ret_sig
