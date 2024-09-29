fs = import('fs')
if not fs.exists('PROPACK/README')
  error('Missing the `PROPACK` submodule! Run `git submodule update --init` to fix this.')
endif

s_src = [
  'PROPACK/single/printstat.F',
  'PROPACK/single/sblasext.F',
  'PROPACK/single/sbsvd.F',
  'PROPACK/single/sgemm_ovwr.F',
  'PROPACK/single/sgetu0.F',
  'PROPACK/single/slanbpro.F',
  'PROPACK/single/slansvd.F',
  'PROPACK/single/slansvd_irl.F',
  'PROPACK/single/sreorth.F',
  'PROPACK/single/sritzvec.F',
  'PROPACK/single/ssafescal.F'
]

d_src = [
  'PROPACK/double/dblasext.F',
  'PROPACK/double/dbsvd.F',
  'PROPACK/double/dgemm_ovwr.F',
  'PROPACK/double/dgetu0.F',
  'PROPACK/double/dlanbpro.F',
  'PROPACK/double/dlansvd.F',
  'PROPACK/double/dlansvd_irl.F',
  'PROPACK/double/dreorth.F',
  'PROPACK/double/dritzvec.F',
  'PROPACK/double/dsafescal.F',
  'PROPACK/double/printstat.F',
]

c_src = [
  'PROPACK/complex8/ccdotc.F',
  'PROPACK/complex8/cblasext.F',
  'PROPACK/complex8/cgemm_ovwr.F',
  'PROPACK/complex8/cgetu0.F',
  'PROPACK/complex8/clanbpro.F',
  'PROPACK/complex8/clansvd.F',
  'PROPACK/complex8/clansvd_irl.F',
  'PROPACK/complex8/creorth.F',
  'PROPACK/complex8/critzvec.F',
  'PROPACK/complex8/csafescal.F',
  'PROPACK/complex8/printstat.F',
  'PROPACK/complex8/sblasext.F',
  'PROPACK/complex8/sbsvd.F',
  'PROPACK/complex8/sgemm_ovwr.F'
]

z_src = [
  'PROPACK/complex16/dblasext.F',
  'PROPACK/complex16/dbsvd.F',
  'PROPACK/complex16/dgemm_ovwr.F',
  'PROPACK/complex16/printstat.F',
  '../_eigen/arpack/ARPACK/SRC/zzdotc.f',
  'PROPACK/complex16/zblasext.F',
  'PROPACK/complex16/zgemm_ovwr.F',
  'PROPACK/complex16/zgetu0.F',
  'PROPACK/complex16/zlanbpro.F',
  'PROPACK/complex16/zlansvd.F',
  'PROPACK/complex16/zlansvd_irl.F',
  'PROPACK/complex16/zreorth.F',
  'PROPACK/complex16/zritzvec.F',
  'PROPACK/complex16/zsafescal.F'
]

# Use risc msg implementation for 64-bit machines, pentium for 32-bit
if meson.get_compiler('cpp').sizeof('void*') > 4
  s_src += ['PROPACK/single/smgs.risc.F']
  d_src += ['PROPACK/double/dmgs.risc.F']
  c_src += ['PROPACK/complex8/cmgs.risc.F']
  z_src += ['PROPACK/complex16/zmgs.risc.F']
else
  s_src += ['PROPACK/single/smgs.pentium.F']
  d_src += ['PROPACK/double/dmgs.pentium.F']
  c_src += ['PROPACK/complex8/cmgs.pentium.F']
  z_src += ['PROPACK/complex16/zmgs.pentium.F']
endif

elements = [
  ['_spropack', s_src, 'spropack.pyf'],
  ['_dpropack', d_src, 'dpropack.pyf'],
  ['_cpropack', c_src, 'cpropack.pyf'],
  ['_zpropack', z_src, 'zpropack.pyf']
]

fargs_propack = ['-U_OPENMP']
foreach ele: elements
  propack_lib = static_library('lib_' + ele[0],
    ele[1],
    fortran_args: [
      fargs_propack,
      fortran_ignore_warnings,
      _fflag_Wno_intrinsic_shadow,
      _fflag_Wno_uninitialized,
      _fflag_fpp,
    ],
    gnu_symbol_visibility: 'hidden',
  )

  propacklib = py3.extension_module(ele[0],
    f2py_gen.process(ele[2]),
    link_with: propack_lib,
    c_args: ['-U_OPENMP', _cpp_Wno_cpp],
    fortran_args: _fflag_Wno_maybe_uninitialized,
    dependencies: [lapack_dep, blas_dep, fortranobject_dep],
    link_args: version_link_args,
    install: true,
    link_language: 'fortran',
    subdir: 'scipy/sparse/linalg/_propack'
  )
endforeach
