import re
import sys
import os
import glob


__all__ = ['needs_g77_abi_wrapper', 'split_fortran_files']


def _uses_veclib(info):
    r_accelerate = re.compile("Accelerate|vecLib")

    extra_link_args = info.get('extra_link_args', '')
    for arg in extra_link_args:
        if r_accelerate.search(arg):
            return True

    return False

def uses_mkl(info):
    r_mkl = re.compile("mkl_core")

    libraries = info.get('libraries', '')
    for library in libraries:
        if r_mkl.search(library):
            return True

    return False

def uses_accelerate(info):
    return _uses_veclib(info)

def needs_g77_abi_wrapper(info):
    """Returns true if g77 ABI wrapper must be used."""
    if uses_accelerate(info):
        return True
    # XXX: is this really true only on Mac OS X ?
    elif uses_mkl(info) and sys.platform == "darwin":
        return True
    else:
        return False


def split_fortran_files(source_dir, subroutines=None):
    """Split each file in `source_dir` into separate files per subroutine.

    Parameters
    ----------
    source_dir : str
        Full path to directory in which sources to be split are located.
    subroutines : list of str, optional
        Subroutines to split. (Default: all)

    Returns
    -------
    fnames : list of str
        List of file names (not including any path) that were created
        in `source_dir`.

    Notes
    -----
    This function is useful for code that can't be compiled with g77 because of
    type casting errors which do work with gfortran.

    Created files are named: ``original_name + '_subr_i' + '.f'``, with ``i``
    starting at zero and ending at ``num_subroutines_in_file - 1``.

    """

    if subroutines is not None:
        subroutines = [x.lower() for x in subroutines]

    def split_file(fname):
        with open(fname, 'rb') as f:
            lines = f.readlines()
            subs = []
            need_split_next = True

            # find lines with SUBROUTINE statements
            for ix, line in enumerate(lines):
                m = re.match(b'^\\s+subroutine\\s+([a-z0-9_]+)\s*\\(', line, re.I)
                if m and line[0] not in b'Cc!*':
                    if subroutines is not None:
                        subr_name = m.group(1).decode('ascii').lower()
                        subr_wanted = (subr_name in subroutines)
                    else:
                        subr_wanted = True
                    if subr_wanted or need_split_next:
                        need_split_next = subr_wanted
                        subs.append(ix)

            # check if no split needed
            if len(subs) <= 1:
                return [fname]

            # write out one file per subroutine
            new_fnames = []
            num_files = len(subs)
            for nfile in range(num_files):
                new_fname = fname[:-2] + '_subr_' + str(nfile) + '.f'
                new_fnames.append(new_fname)
                with open(new_fname, 'wb') as fn:
                    if nfile + 1 == num_files:
                        fn.writelines(lines[subs[nfile]:])
                    else:
                        fn.writelines(lines[subs[nfile] : subs[nfile+1]])

        return new_fnames

    exclude_pattern = re.compile('_subr_[0-9]')
    source_fnames = [f for f in glob.glob(os.path.join(source_dir, '*.f'))
                             if not exclude_pattern.search(os.path.basename(f))]
    fnames = []
    for source_fname in source_fnames:
        created_files = split_file(source_fname)
        if created_files is not None:
            for cfile in created_files:
                fnames.append(os.path.basename(cfile))

    return fnames
