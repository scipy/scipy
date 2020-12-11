'''First stab at generating PYX files.'''

from typing import NamedTuple

class WrapperDef(NamedTuple):
    ufunc_name: str
    num_inputs: int
    func_name: str

_front_matter = '''# distutils: language = c++
# cython: language_level=3

from numpy cimport (
    import_array,
    import_ufunc,
    PyUFunc_FromFuncAndData,
    PyUFuncGenericFunction,
    PyUFunc_None,
    {NUMPY_CIMPORTS}
)
from templated_pyufunc cimport PyUFunc_T
from func_defs cimport (
    {FUNC_DEFS_CIMPORTS},
)

cdef extern from "boost/math/distributions/{short_boost_name}.hpp" namespace "boost::math" nogil:
    cdef cppclass {boost_name} nogil:
        pass

# Workaround for Cython's lack of non-type template parameter support
cdef extern from * nogil:
    {NINPUT_CTYPEDEFS}

_DUMMY = ""
import_array()
import_ufunc()
'''

_ufunc_front_matter = '''
cdef PyUFuncGenericFunction loop_func{wrapper_idx}[{num_types}]
cdef void* func{wrapper_idx}[1*{num_types}]
cdef char types{wrapper_idx}[4*{num_types}]
'''

_ufunc_type_guts = '''
loop_func{wrapper_idx}[{type_idx}] = <PyUFuncGenericFunction>{loop_fun}[{ctype}, NINPUTS{num_inputs}]
func{wrapper_idx}[{type_idx}] = <void*>boost_{func_name}{num_ctor_args}[{boost_name}, {argtypes}]
types{wrapper_idx}[0+{type_idx}*4] = {NPY_TYPE}
types{wrapper_idx}[1+{type_idx}*4] = {NPY_TYPE}
types{wrapper_idx}[2+{type_idx}*4] = {NPY_TYPE}
types{wrapper_idx}[3+{type_idx}*4] = {NPY_TYPE}
'''

_ufunc_template = '''
{ufunc_name} = PyUFunc_FromFuncAndData(
    loop_func{wrapper_idx},
    func{wrapper_idx},
    types{wrapper_idx},
    {num_types}, # number of supported input types
    {num_inputs}, # number of input args
    1, # number of output args
    PyUFunc_None, # `identity` element, never mind this
    "{ufunc_name}", # function name
    "{ufunc_name}(x) -> computes ... of ... distribution", # docstring
    0 # unused
)
'''

def ufunc_gen(wrapper_prefix: str, types: list, num_ctor_args: int, filename: str, boost_dist: str):

    wrappers = [
        WrapperDef(ufunc_name=f'_{wrapper_prefix}_pdf',
                   num_inputs=1+num_ctor_args,
                   func_name='pdf'),
        WrapperDef(ufunc_name=f'_{wrapper_prefix}_cdf',
                   num_inputs=1+num_ctor_args,
                   func_name='cdf'),
        WrapperDef(ufunc_name=f'_{wrapper_prefix}_icdf',
                   num_inputs=1+num_ctor_args,
                   func_name='icdf'),
        WrapperDef(ufunc_name=f'_{wrapper_prefix}_quantile',
                   num_inputs=1+num_ctor_args,
                   func_name='quantile'),
        WrapperDef(ufunc_name=f'_{wrapper_prefix}_iquantile',
                   num_inputs=1+num_ctor_args,
                   func_name='iquantile'),
        WrapperDef(ufunc_name=f'_{wrapper_prefix}_mean',
                   num_inputs=num_ctor_args,
                   func_name='mean'),
        WrapperDef(ufunc_name=f'_{wrapper_prefix}_variance',
                   num_inputs=num_ctor_args,
                   func_name='variance'),
        WrapperDef(ufunc_name=f'_{wrapper_prefix}_skewness',
                   num_inputs=num_ctor_args,
                   func_name='skewness'),
        WrapperDef(ufunc_name=f'_{wrapper_prefix}_kurtosis_excess',
                   num_inputs=num_ctor_args,
                   func_name='kurtosis_excess'),
    ]

    with open(filename, 'w') as fp:
        short_boost_name = boost_dist.split('_distribution')[0]
        unique_num_inputs = set({w.num_inputs for w in wrappers})
        fp.write(_front_matter.format(
            NUMPY_CIMPORTS=',\n    '.join(types),
            FUNC_DEFS_CIMPORTS=',\n    '.join(f'boost_{w.func_name}{num_ctor_args}' for w in wrappers),
            NINPUT_CTYPEDEFS='\n    '.join(f'ctypedef int NINPUTS{num_inputs} "{num_inputs}"' for num_inputs in unique_num_inputs),
            short_boost_name=short_boost_name,
            boost_name=boost_dist,
        ))
        if 'NPY_LONGDOUBLE' in types:
            fp.write('ctypedef long double longdouble\n')
        for ii, w in enumerate(wrappers):
            if w.num_inputs == 0:
                print(f'skipping {w.ufunc_name} ufunc because it has 0 inputs')
                continue
            fp.write(_ufunc_front_matter.format(
                num_types=len(types),
                wrapper_idx=ii,
                num_inputs=w.num_inputs,
            ))
            for jj, t in enumerate(types):
                ctype = {
                    'NPY_LONGDOUBLE': 'longdouble',
                    'NPY_DOUBLE': 'double',
                    'NPY_FLOAT': 'float',
                    'NPY_FLOAT16': 'float16_t',
                }[t]
                fp.write(_ufunc_type_guts.format(
                    type_idx=jj,
                    NPY_TYPE=t,
                    ctype=ctype,
                    wrapper_idx=ii,
                    loop_fun='PyUFunc_T',
                    num_inputs=w.num_inputs,
                    func_name=w.func_name,
                    argtypes=', '.join([ctype]*(num_ctor_args+1)),
                    boost_name=boost_dist,
                    num_ctor_args=num_ctor_args,
                ))

            fp.write(_ufunc_template.format(
                ufunc_name=w.ufunc_name,
                num_types=len(types),
                wrapper_idx=ii,
                num_inputs=w.num_inputs,
            ))
