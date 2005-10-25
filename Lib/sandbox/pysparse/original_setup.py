#!/usr/bin/env python

from distutils.core import setup, Extension
import glob
import os, socket
import sys

# default settings
library_dirs_list= []
libraries_list = ['lapack', 'blas', 'g2c']
superlu_defs = [('USE_VENDOR_BLAS',1)]
f77_defs = []
linky=[]
compily=[]

umfpack_defs = [('DINT', 1), ('NBLAS', 1)] # most basic configuration, no BLAS
umfpack_libraries = []
umfpack_include_dirs = ['amd', 'umfpack']
umfpack_library_dirs = []

#umfpack_defs = [('DINT', 1), ('CBLAS', 1)] # with atlas c-blas (http://math-atlas.sourceforge.net)
#umfpack_libraries = ['atlas', 'cblas', 'm']
#umfpack_include_dirs = ['amd', 'umfpack'] # you may need to set this to find the atlas
#umfpack_library_dirs = []

hostname = socket.gethostname()
if hostname in ['bree', 'brokoli', 'givens2', 'givens4', 'nedelec', 'gondor']:
    # SuSE Linux 8.x and 9.0
    libraries_list = ['lapack', 'blas', 'g2c']
elif hostname == 'maxwell':
    # AMD Opteron 'x86-64' architecture using ACML
    libraries_list = ['acml', 'g2c']
    library_dirs_list = ['/opt/acml/gnu64/lib']
elif hostname == 'sysiphus':
    # Linux RedHat 7.3 2.4.18-10 i686 with atlas Lapack routines
    library_dirs_list= ['/hg/u/vasseur/Linux/lib/atlas']
    libraries_list = ['lapack', 'f77blas', 'cblas', 'atlas', 'g2c']
elif hostname == 'sophokles':
    #SunOS sophokles 5.8 Generic_108528-15 sun4u sparc SUNW,Sun-Fire-880
    library_dirs_list= ['/hg/s/solaris/8/opt/SUNWspro/lib']
    libraries_list = ['F77', 'sunperf', 'fui', 'fsu', 'sunmath']
elif hostname == 'sim0':
    # Linux sim0 2.2.19-7.0.16enterprise #1 SMP Wed Mar 13 13:23:22 EST 2002 i686 unknown
    library_dirs_list = ['/home/geus/linux/lib']
    libraries_list = ['lapack', 'blas', 'g2c']
elif hostname == 'stardust':
    # HP-UX stardust B.11.11 U 9000/800 3761215035 unlimited-user license
    library_dirs_list = ['/home/infk/geus/lib/pa20_64']
    libraries_list = ['lapack']
elif hostname == 'zuse':
    # SunOS zuse 5.6 Generic_105181-21 sun4u sparc SUNW,Ultra-Enterprise
    library_dirs_list = ['/software/SunOS/5.X/opt/SUNWspro/WS6/lib']
    libraries_list = ['F77', 'sunperf', 'fui', 'fsu', 'sunmath']
elif hostname == 'Rivendell':
    # WinXP, using VC++ and Intel MKL
    #
    # Uses BLAS and LAPACK from the Intel Math Kernel Library 5.2
    # DLL directory must be in the PATH
    library_dirs_list = [r'C:\Program Files\Intel\MKL\ia32\lib']
    libraries_list = ['mkl_c_dll']
    superlu_defs = [("NO_TIMER",1), ("NoChange",1), ('USE_VENDOR_BLAS',1)]
    f77_defs = [("NOF77UNDERSCORE",1)]
elif hostname == 'nedelec-vmware':
    # Win32 using MinGW
    libraries_list = ['lapack', 'blas', 'g2c']
    superlu_defs = [("NO_TIMER",1), ('USE_VENDOR_BLAS',1)]
elif hostname == 'psw283.psi.ch':
    # OSF1 psw283.psi.ch V4.0 1530 alpha
    library_dirs_list = ['/data/geus/lib']
    libraries_list = ['dxml']
elif sys.platform == 'darwin':
    superlu_defs = [('USE_VENDOR_BLAS',1)]
    library_dirs_list = ['/System/Library/Frameworks']
    libraries_list = []
    f77_defs = []

    
    # the following 'linky' arguments must not be concatenated together into a single
    # string, c.f. <http://mail.python.org/pipermail/distutils-sig/2003-December/003532.html>
    
    if sys.exec_prefix == '/sw':
        # fink python
        linky=["-faltivec","-framework","vecLib","-bundle_loader","/sw/bin/python"]
    else:
        # Apple python
        linky=["-faltivec","-framework","vecLib"]
        # The python Framework build is compiled with (and propagates to all library builds) the 
        # '-fno-common' flag. Nobody seems to know why.
        # (c.f. <https://sourceforge.net/tracker/?func=detail&atid=105470&aid=768306&group_id=5470>)
        # This flag wreaks havoc with the nightmarishly circular declarations in the itsolvers module.
        # We reset it by appending this flag:
        compily=["-fcommon"]

elif sys.platform == 'linux2':
    ## fix for Fedora core 4, 'g2c' doesn't exist and isn't required
    if 'redhat-release' in os.listdir('/etc'):
        f = open('/etc/redhat-release', 'r')
        if 'release 4' in f.read():
            libraries_list = ['lapack', 'blas']
        f.close()
        

ext_modules = [Extension('spmatrix', ['Src/spmatrixmodule.c']),
               Extension('itsolvers', ['Src/itsolversmodule.c',
                                       'Src/pcg.c',
                                       'Src/gmres.c',
                                       'Src/minres.c',
                                       'Src/qmrs.c',
                                       'Src/bicgstab.c',
                                       'Src/cgs.c'],
                         library_dirs=library_dirs_list,
                         libraries=libraries_list,
                         define_macros=f77_defs,
                         extra_compile_args=compily,
                         extra_link_args=linky),
               Extension('precon',  [os.path.join('Src', 'preconmodule.c')],
                         library_dirs=library_dirs_list,
                         libraries=libraries_list,
                         define_macros=f77_defs,extra_link_args=linky),
               Extension('superlu', [os.path.join('Src', 'superlumodule.c'),
                                     "superlu/dcolumn_bmod.c",
                                     "superlu/dcolumn_dfs.c",
                                     "superlu/dcomplex.c",
                                     "superlu/scomplex.c",
                                     "superlu/dcopy_to_ucol.c",
                                     "superlu/dgscon.c",
                                     "superlu/dgsequ.c",
                                     "superlu/dgsrfs.c",
                                     "superlu/dgssv.c",
                                     "superlu/dgssvx.c",
                                     "superlu/dgstrf.c",
                                     "superlu/dgstrs.c",
                                     "superlu/dlacon.c",
                                     "superlu/dlamch.c",
                                     "superlu/dlangs.c",
                                     "superlu/dlaqgs.c",
                                     "superlu/dmemory.c",
                                     "superlu/colamd.c",
                                     "superlu/dpanel_bmod.c",
                                     "superlu/dpanel_dfs.c",
                                     "superlu/dpivotL.c",
                                     "superlu/dpivotgrowth.c",
                                     "superlu/dpruneL.c",
                                     "superlu/dreadhb.c",
                                     "superlu/dsnode_bmod.c",
                                     "superlu/dsnode_dfs.c",
                                     "superlu/dsp_blas2.c",
                                     "superlu/dsp_blas3.c",
                                     "superlu/superlu_timer.c",
                                     "superlu/dutil.c",
                                     "superlu/dzsum1.c",
                                     "superlu/get_perm_c.c",
                                     "superlu/icmax1.c",
                                     "superlu/izmax1.c",
                                     "superlu/lsame.c",
                                     "superlu/memory.c",
                                     "superlu/mmd.c",
                                     "superlu/relax_snode.c",
                                     "superlu/sp_coletree.c",
                                     "superlu/sp_ienv.c",
                                     "superlu/sp_preorder.c",
                                     "superlu/util.c",
                                     "superlu/xerbla.c"],
                         define_macros=superlu_defs,
                         include_dirs=["superlu"],
                         library_dirs=library_dirs_list,
                         libraries=libraries_list,extra_link_args=linky),
               Extension('jdsym', [os.path.join('Src', 'jdsymmodule.c')],
                         library_dirs=library_dirs_list,
                         libraries=libraries_list,
                         define_macros=f77_defs,extra_link_args=linky),
               Extension('umfpack', sources=[os.path.join('Src', 'umfpackmodule.c'),
                                     'amd/amd_aat.c',
                                     'amd/amd_1.c',
                                     'amd/amd_2.c',
                                     'amd/amd_dump.c',
                                     'amd/amd_postorder.c',
                                     'amd/amd_post_tree.c',
                                     'amd/amd_defaults.c',
                                     'amd/amd_order.c',
                                     'amd/amd_control.c',
                                     'amd/amd_info.c',
                                     'amd/amd_valid.c',
                                     'umfpack/umf_analyze.c',
                                     'umfpack/umf_apply_order.c',
                                     'umfpack/umf_colamd.c',
                                     'umfpack/umf_free.c',
                                     'umfpack/umf_fsize.c',
                                     'umfpack/umf_is_permutation.c',
                                     'umfpack/umf_malloc.c',
                                     'umfpack/umf_realloc.c',
                                     'umfpack/umf_report_perm.c',
                                     'umfpack/umf_singletons.c',
                                     'umfpack/umfpack_timer.c',
                                     'umfpack/umfpack_tictoc.c',
                                     'umfpack/umf_lhsolve.c',
                                     'umfpack/umf_uhsolve.c',
                                     'umfpack/umf_triplet_map_nox.c',
                                     'umfpack/umf_triplet_nomap_x.c',
                                     'umfpack/umf_triplet_nomap_nox.c',
                                     'umfpack/umf_triplet_map_x.c',
                                     'umfpack/umf_assemble_fixq.c',
                                     'umfpack/umf_assemble.c',
                                     'umfpack/umf_blas3_update.c',
                                     'umfpack/umf_build_tuples.c',
                                     'umfpack/umf_create_element.c',
                                     'umfpack/umf_dump.c',
                                     'umfpack/umf_extend_front.c',
                                     'umfpack/umf_garbage_collection.c',
                                     'umfpack/umf_get_memory.c',
                                     'umfpack/umf_init_front.c',
                                     'umfpack/umf_kernel.c',
                                     'umfpack/umf_kernel_init.c',
                                     'umfpack/umf_kernel_wrapup.c',
                                     'umfpack/umf_local_search.c',
                                     'umfpack/umf_lsolve.c',
                                     'umfpack/umf_ltsolve.c',
                                     'umfpack/umf_mem_alloc_element.c',
                                     'umfpack/umf_mem_alloc_head_block.c',
                                     'umfpack/umf_mem_alloc_tail_block.c',
                                     'umfpack/umf_mem_free_tail_block.c',
                                     'umfpack/umf_mem_init_memoryspace.c',
                                     'umfpack/umf_report_vector.c',
                                     'umfpack/umf_row_search.c',
                                     'umfpack/umf_scale_column.c',
                                     'umfpack/umf_set_stats.c',
                                     'umfpack/umf_solve.c',
                                     'umfpack/umf_symbolic_usage.c',
                                     'umfpack/umf_transpose.c',
                                     'umfpack/umf_tuple_lengths.c',
                                     'umfpack/umf_usolve.c',
                                     'umfpack/umf_utsolve.c',
                                     'umfpack/umf_valid_numeric.c',
                                     'umfpack/umf_valid_symbolic.c',
                                     'umfpack/umf_grow_front.c',
                                     'umfpack/umf_start_front.c',
                                     'umfpack/umf_2by2.c',
                                     'umfpack/umf_store_lu.c',
                                     'umfpack/umf_scale.c',
                                     'umfpack/umfpack_wsolve.c',
                                     'umfpack/umfpack_col_to_triplet.c',
                                     'umfpack/umfpack_defaults.c',
                                     'umfpack/umfpack_free_numeric.c',
                                     'umfpack/umfpack_free_symbolic.c',
                                     'umfpack/umfpack_get_numeric.c',
                                     'umfpack/umfpack_get_lunz.c',
                                     'umfpack/umfpack_get_symbolic.c',
                                     'umfpack/umfpack_numeric.c',
                                     'umfpack/umfpack_qsymbolic.c',
                                     'umfpack/umfpack_report_control.c',
                                     'umfpack/umfpack_report_info.c',
                                     'umfpack/umfpack_report_matrix.c',
                                     'umfpack/umfpack_report_numeric.c',
                                     'umfpack/umfpack_report_perm.c',
                                     'umfpack/umfpack_report_status.c',
                                     'umfpack/umfpack_report_symbolic.c',
                                     'umfpack/umfpack_report_triplet.c',
                                     'umfpack/umfpack_report_vector.c',
                                     'umfpack/umfpack_solve.c',
                                     'umfpack/umfpack_symbolic.c',
                                     'umfpack/umfpack_transpose.c',
                                     'umfpack/umfpack_triplet_to_col.c',
                                     'umfpack/umfpack_scale.c',
                                     'umfpack/umfpack_load_numeric.c',
                                     'umfpack/umfpack_save_numeric.c',
                                     'umfpack/umfpack_load_symbolic.c',
                                     'umfpack/umfpack_save_symbolic.c'],
                         define_macros=umfpack_defs, include_dirs=umfpack_include_dirs,
                         libraries=umfpack_libraries, library_dirs=umfpack_library_dirs)
               
                                     
               ]

execfile(os.path.join('Lib', 'pysparse_version.py'))

setup(name = 'pysparse',
      version = version,
      description = 'Python Sparse Matrix Package',
      author = 'Roman Geus',
      author_email = 'roman@geus.ch',
      license = 'BSD style',
      url = 'http://www.geus.ch',
      packages = [''],
      package_dir = {'': 'Lib'},
      include_dirs = ['Include'],
      headers = glob.glob(os.path.join ("Include","pysparse","*.h")),
      ext_modules = ext_modules
      )
