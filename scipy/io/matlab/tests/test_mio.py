#!/usr/bin/env python
''' Nose test generators

Need function load / save / roundtrip tests

'''
from os.path import join, dirname
from glob import glob
from StringIO import StringIO
from tempfile import mkdtemp
import warnings
import shutil
import gzip
import copy

from numpy.testing import \
     assert_array_almost_equal, \
     assert_equal, \
     assert_raises, dec

from nose.tools import assert_true

import numpy as np
from numpy import array
import scipy.sparse as SP

from scipy.io.matlab.mio import loadmat, savemat, find_mat_file
from scipy.io.matlab.mio5 import MatlabObject

test_data_path = join(dirname(__file__), 'data')

def mlarr(*args, **kwargs):
    ''' Convenience function to return matlab-compatible 2D array
    Note that matlab writes empty shape as (0,0) - replicated here
    '''
    arr = np.array(*args, **kwargs)
    if arr.size:
        return np.atleast_2d(arr)
    # empty elements return as shape (0,0)
    return arr.reshape((0,0))


# Define cases to test
theta = np.pi/4*np.arange(9,dtype=float).reshape(1,9)
case_table4 = [
    {'name': 'double',
     'expected': {'testdouble': theta}
     }]
case_table4.append(
    {'name': 'string',
     'expected': {'teststring':
                  array([u'"Do nine men interpret?" "Nine men," I nod.'])},
     })
case_table4.append(
    {'name': 'complex',
     'expected': {'testcomplex': np.cos(theta) + 1j*np.sin(theta)}
     })
A = np.zeros((3,5))
A[0] = range(1,6)
A[:,0] = range(1,4)
case_table4.append(
    {'name': 'matrix',
     'expected': {'testmatrix': A},
     })
case_table4.append(
    {'name': 'sparse',
     'expected': {'testsparse': SP.coo_matrix(A)},
     })
B = A.astype(complex)
B[0,0] += 1j
case_table4.append(
    {'name': 'sparsecomplex',
     'expected': {'testsparsecomplex': SP.coo_matrix(B)},
     })
case_table4.append(
    {'name': 'multi',
     'expected': {'theta': theta,
                  'a': A},
     })
case_table4.append(
    {'name': 'minus',
     'expected': {'testminus': mlarr(-1)},
     })
case_table4.append(
    {'name': 'onechar',
     'expected': {'testonechar': array([u'r'])},
     })
# Cell arrays stored as object arrays
CA = mlarr([
    [], # placeholder, object array constructor wierdness otherwise
    mlarr(1),
    mlarr([1,2]),
    mlarr([1,2,3])], dtype=object).reshape(1,-1)
CA[0,0] = array(
    [u'This cell contains this string and 3 arrays of increasing length'])
case_table5 = [
    {'name': 'cell',
     'expected': {'testcell': CA}}]
CAE = mlarr([
    mlarr(1),
    mlarr(2),
    mlarr([]),
    mlarr([]),
    mlarr(3)], dtype=object).reshape(1,-1)
case_table5.append(
    {'name': 'emptycell',
     'expected': {'testemptycell': CAE}})
case_table5.append(
    {'name': 'stringarray',
     'expected': {'teststringarray': array(
    [u'one  ', u'two  ', u'three'])},
     })
case_table5.append(
    {'name': '3dmatrix',
     'expected': {
    'test3dmatrix': np.transpose(np.reshape(range(1,25), (4,3,2)))}
     })
st_sub_arr = array([np.sqrt(2),np.exp(1),np.pi]).reshape(1,3)
dtype = [(n, object) for n in ['stringfield', 'doublefield', 'complexfield']]
st1 = np.zeros((1,1), dtype)
st1['stringfield'][0,0] = array([u'Rats live on no evil star.'])
st1['doublefield'][0,0] = st_sub_arr
st1['complexfield'][0,0] = st_sub_arr * (1 + 1j)
case_table5.append(
    {'name': 'struct',
     'expected': {'teststruct': st1}
     })
CN = np.zeros((1,2), dtype=object)
CN[0,0] = mlarr(1)
CN[0,1] = np.zeros((1,3), dtype=object)
CN[0,1][0,0] = mlarr(2, dtype=np.uint8)
CN[0,1][0,1] = mlarr([[3]], dtype=np.uint8)
CN[0,1][0,2] = np.zeros((1,2), dtype=object)
CN[0,1][0,2][0,0] = mlarr(4, dtype=np.uint8)
CN[0,1][0,2][0,1] = mlarr(5, dtype=np.uint8)
case_table5.append(
    {'name': 'cellnest',
     'expected': {'testcellnest': CN},
     })
st2 = np.empty((1,1), dtype=[(n, object) for n in ['one', 'two']])
st2[0,0]['one'] = mlarr(1)
st2[0,0]['two'] = np.empty((1,1), dtype=[('three', object)])
st2[0,0]['two'][0,0]['three'] = array([u'number 3'])
case_table5.append(
    {'name': 'structnest',
     'expected': {'teststructnest': st2}
     })
a = np.empty((1,2), dtype=[(n, object) for n in ['one', 'two']])
a[0,0]['one'] = mlarr(1)
a[0,0]['two'] = mlarr(2)
a[0,1]['one'] = array([u'number 1'])
a[0,1]['two'] = array([u'number 2'])
case_table5.append(
    {'name': 'structarr',
     'expected': {'teststructarr': a}
     })
ODT = np.dtype([(n, object) for n in
                 ['expr', 'inputExpr', 'args',
                  'isEmpty', 'numArgs', 'version']])
MO = MatlabObject(np.zeros((1,1), dtype=ODT), 'inline')
m0 = MO[0,0]
m0['expr'] = array([u'x'])
m0['inputExpr'] = array([u' x = INLINE_INPUTS_{1};'])
m0['args'] = array([u'x'])
m0['isEmpty'] = mlarr(0)
m0['numArgs'] = mlarr(1)
m0['version'] = mlarr(1)
case_table5.append(
    {'name': 'object',
     'expected': {'testobject': MO}
     })
u_str = file(
    join(test_data_path, 'japanese_utf8.txt'),
    'rb').read().decode('utf-8')
case_table5.append(
    {'name': 'unicode',
    'expected': {'testunicode': array([u_str])}
    })
# These should also have matlab load equivalents,
# but I can't get to matlab at the moment
case_table5_rt = case_table5[:]
case_table5_rt.append(
    {'name': 'sparsefloat',
     'expected': {'testsparsefloat':
                  SP.coo_matrix(array([[1,0,2],[0,-3.5,0]]))},
     })
case_table5_rt.append(
    {'name': 'sparsecomplex',
     'expected': {'testsparsecomplex':
                  SP.coo_matrix(array([[-1+2j,0,2],[0,-3j,0]]))},
     })
case_table5_rt.append(
    {'name': 'objectarray',
     'expected': {'testobjectarray': np.repeat(MO, 2).reshape(1,2)}})

def _check_level(label, expected, actual):
    """ Check one level of a potentially nested array """
    if SP.issparse(expected): # allow different types of sparse matrices
        assert SP.issparse(actual)
        assert_array_almost_equal(actual.todense(),
                                  expected.todense(),
                                  err_msg = label,
                                  decimal = 5)
        return
    # Check types are as expected
    typex = type(expected)
    typac = type(actual)
    assert typex is typac, \
           "Expected type %s, got %s at %s" % (typex, typac, label)
    # A field in a record array may not be an ndarray
    # A scalar from a record array will be type np.void
    if not isinstance(expected,
                      (np.void, np.ndarray, MatlabObject)): 
        assert_equal(expected, actual)
        return
    # This is an ndarray-like thing
    assert_true(expected.shape == actual.shape,
                msg='Expected shape %s, got %s at %s' % (expected.shape,
                                                         actual.shape,
                                                         label)
                )
    ex_dtype = expected.dtype
    if ex_dtype.hasobject: # array of objects
        if isinstance(expected, MatlabObject):
            assert_equal(expected.classname, actual.classname)
        for i, ev in enumerate(expected):
            level_label = "%s, [%d], " % (label, i)
            _check_level(level_label, ev, actual[i])
        return
    if ex_dtype.fields: # probably recarray
        for fn in ex_dtype.fields:
            level_label = "%s, field %s, " % (label, fn)
            _check_level(level_label,
                         expected[fn], actual[fn])
        return
    if ex_dtype.type in (np.unicode, # string
                         np.unicode_):
        assert_equal(actual, expected, err_msg=label)
        return
    # Something numeric
    assert_array_almost_equal(actual, expected, err_msg=label, decimal=5)

def _load_check_case(name, files, case):
    for file_name in files:
        matdict = loadmat(file_name, struct_as_record=True)
        label = "test %s; file %s" % (name, file_name)
        for k, expected in case.items():
            k_label = "%s, variable %s" % (label, k)
            assert k in matdict, "Missing key at %s" % k_label
            _check_level(k_label, expected, matdict[k])

# Round trip tests
def _rt_check_case(name, expected, format):
    mat_stream = StringIO()
    savemat(mat_stream, expected, format=format)
    mat_stream.seek(0)
    _load_check_case(name, [mat_stream], expected)


# generator for load tests
def test_load():
    for case in case_table4 + case_table5:
        name = case['name']
        expected = case['expected']
        filt = join(test_data_path, 'test%s_*.mat' % name)
        files = glob(filt)
        assert files, "No files for test %s using filter %s" % (name, filt)
        yield _load_check_case, name, files, expected

# generator for round trip tests
def test_round_trip():
    for case in case_table4 + case_table5_rt:
        name = case['name'] + '_round_trip'
        expected = case['expected']
        format = case in case_table4 and '4' or '5'
        yield _rt_check_case, name, expected, format

def test_gzip_simple():
    xdense = np.zeros((20,20))
    xdense[2,3]=2.3
    xdense[4,5]=4.5
    x = SP.csc_matrix(xdense)

    name = 'gzip_test'
    expected = {'x':x}
    format='4'

    tmpdir = mkdtemp()
    try:
        fname = join(tmpdir,name)
        mat_stream = gzip.open( fname,mode='wb')
        savemat(mat_stream, expected, format=format)
        mat_stream.close()

        mat_stream = gzip.open( fname,mode='rb')
        actual = loadmat(mat_stream, struct_as_record=True)
        mat_stream.close()
    finally:
        shutil.rmtree(tmpdir)

    assert_array_almost_equal(actual['x'].todense(),
                              expected['x'].todense())


def test_mat73():
    # Check any hdf5 files raise an error
    filenames = glob(
        join(test_data_path, 'testhdf5*.mat'))
    assert len(filenames)
    for filename in filenames:
        assert_raises(NotImplementedError,
                      loadmat,
                      filename,
                      struct_as_record=True)


def test_warnings():
    fname = join(test_data_path, 'testdouble_7.1_GLNX86.mat')
    warnings.simplefilter('error')
    # This should not generate a warning
    mres = loadmat(fname, struct_as_record=True)
    # This neither
    mres = loadmat(fname, struct_as_record=False)
    # This should
    yield assert_raises, FutureWarning, loadmat, fname
    # This too
    yield assert_raises, FutureWarning, find_mat_file, fname
    # we need kwargs for this one
    yield (lambda a, k: assert_raises(*a, **k), 
          (DeprecationWarning, loadmat, fname), 
          {'struct_as_record':True, 'basename':'raw'})
    # Test warning for default format change
    savemat(StringIO(), {}, False, '4')
    savemat(StringIO(), {}, False, '5')
    yield assert_raises, FutureWarning, savemat, StringIO(), {}
    warnings.resetwarnings()


def test_regression_653():
    """Regression test for #653."""
    assert_raises(TypeError, savemat, StringIO(), {'d':{1:2}}, format='5')
