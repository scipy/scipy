'''Generate Cython PYX wrappers for Boost stats distributions.'''

from typing import NamedTuple
from warnings import warn

class _MethodDef(NamedTuple):
    ufunc_name: str
    num_inputs: int
    boost_func_name: str


def _ufunc_gen(scipy_dist: str, types: list, ctor_args: tuple,
              filename: str, boost_dist: str, x_funcs: list, no_x_funcs: list):

    # We need methods defined for each rv_continuous/_discrete internal method:
    #     i.e.: _pdf, _cdf, etc.
    # Some of these methods take constructor arguments and 1 extra argument,
    #     e.g.: _pdf(x, *ctor_args), _ppf(q, *ctor_args)
    # while some of the methods take only constructor arguments:
    #     e.g.: _stats(*ctor_args)
    num_ctor_args = len(ctor_args)
    methods = [_MethodDef(
        ufunc_name=f'_{scipy_dist}_{x_func}',
        num_inputs=num_ctor_args+1,  # +1 for the x argument
        boost_func_name=x_func if boost_dist != 'beta_distribution'
        else 'pdf_beta' if x_func == 'pdf' else x_func,
    ) for x_func in x_funcs]
    methods += [_MethodDef(
        ufunc_name=f'_{scipy_dist}_{func}',
        num_inputs=num_ctor_args,
        boost_func_name=func,
    ) for func in no_x_funcs]

    # Identify potential ufunc issues:
    no_input_methods = [m for m in methods if m.num_inputs == 0]
    if no_input_methods:
        raise ValueError("ufuncs must have >0 arguments! "
                         f"Cannot construct these ufuncs: {no_input_methods}")

    boost_hdr_name = boost_dist.split('_distribution')[0]
    unique_num_inputs = set({m.num_inputs for m in methods})
    has_NPY_FLOAT16 = 'NPY_FLOAT16' in types
    has_NPY_LONGDOUBLE = 'NPY_LONGDOUBLE' in types
    line_joiner = ',\n    '
    Q = '"'
    num_types = len(types)
    loop_fun = 'PyUFunc_T'
    func_defs_cimports = line_joiner.join(
        f"boost_{m.boost_func_name}{num_ctor_args}" for m in methods)
    nontype_params = '\n    '.join(
        f'ctypedef int NINPUTS{n} "{n}"' for n in unique_num_inputs)

    with open(filename, 'w') as fp:
        fp.write('\n'.join([
            f'# distutils: language = c++',
            f'# cython: language_level=3',
            f'',
            f'from numpy cimport (',
            f'    import_array,',
            f'    import_ufunc,',
            f'    PyUFunc_FromFuncAndData,',
            f'    PyUFuncGenericFunction,',
            f'    PyUFunc_None,',
            f'    {line_joiner.join(types)}',
            f')',
            f'from templated_pyufunc cimport PyUFunc_T',
            f'from func_defs cimport (',
            f'    {func_defs_cimports},',
            f')',
            f'cdef extern from "boost/math/distributions/{boost_hdr_name}.hpp"'
            ' namespace "boost::math" nogil:',
            f'    cdef cppclass {boost_dist} nogil:',
            f'        pass',
            f'',
            f'# Workaround for Cython\'s lack of non-type template parameter',
            f'# support',
            f'cdef extern from * nogil:',
            f'    {nontype_params}',
            f'',
            f'_DUMMY = ""',
            f'import_array()',
            f'import_ufunc()',
        ]) + '\n')

        if has_NPY_LONGDOUBLE:
            fp.write('ctypedef long double longdouble\n\n')
        if has_NPY_FLOAT16:
            warn('Boost stats NPY_FLOAT16 ufunc generation not '
                 'currently not supported!')

        # Generate ufuncs for each method
        for ii, m in enumerate(methods):
            fp.write('\n'.join([
                f'cdef PyUFuncGenericFunction loop_func{ii}[{num_types}]',
                f'cdef void* func{ii}[1*{num_types}]',
                #                                  +1 for output arg
                f'cdef char types{ii}[{m.num_inputs+1}*{num_types}]',
            ]) + '\n')

            for jj, T in enumerate(types):
                ctype = {
                    'NPY_LONGDOUBLE': 'longdouble',
                    'NPY_DOUBLE': 'double',
                    'NPY_FLOAT': 'float',
                    'NPY_FLOAT16': 'npy_half',
                }[T]
                boost_fun = f'boost_{m.boost_func_name}{num_ctor_args}'
                boost_tmpl = (f'{boost_dist}, '
                              f'{", ".join([ctype]*(1+num_ctor_args))}')
                fp.write('\n'.join([
                    (f'loop_func{ii}[{jj}] = '
                     f'<PyUFuncGenericFunction>{loop_fun}[{ctype}, '
                     f'NINPUTS{m.num_inputs}]'),
                    f'func{ii}[{jj}] = <void*>{boost_fun}[{boost_tmpl}]',
                ]) + '\n')
                fp.write('\n'.join([
                    f'types{ii}[{tidx}+{jj}*{m.num_inputs+1}] = {T}'
                    for tidx in range(m.num_inputs+1)
                ]) + '\n')
            arg_list_str = ', '.join(ctor_args)
            if m.boost_func_name in x_funcs:
                arg_list_str = 'x, ' + arg_list_str
            fp.write('\n'.join([
                f'{m.ufunc_name} = PyUFunc_FromFuncAndData(',
                f'    loop_func{ii},',
                f'    func{ii},',
                f'    types{ii},',
                f'    {num_types},  # number of supported input types',
                f'    {m.num_inputs},  # number of input args',
                f'    1,  # number of output args',
                f'    PyUFunc_None,  # `identity` element, never mind this',
                f'    "{m.ufunc_name}",  # function name',
                f'    ("{m.ufunc_name}({arg_list_str}) -> computes "',
                f'     "{m.boost_func_name} of {scipy_dist} distribution"),',
                f'    0 # unused',
                f')',
                f'',
            ]) + '\n')
