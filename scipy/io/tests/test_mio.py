#!/usr/bin/env python

import os
from glob import glob
from cStringIO import StringIO
from tempfile import mkstemp
from numpy.testing import set_package_path, restore_path, NumpyTestCase, NumpyTest
from numpy.testing import assert_equal, assert_array_almost_equal
from numpy import arange, array, eye, pi, cos, exp, sin, sqrt, ndarray,  \
     zeros, reshape, transpose, empty
import scipy.sparse as SP

set_package_path()
from scipy.io.mio import loadmat, savemat
from scipy.io.mio5 import mat_obj, mat_struct
restore_path()

try:  # Python 2.3 support
    from sets import Set as set
except:
    pass

test_data_path = os.path.join(os.path.dirname(__file__), './data')

class TestMIOArray(NumpyTestCase):
    def __init__(self, *args, **kwargs):
        super(TestMIOArray, self).__init__(*args, **kwargs)

    def _check_level(self, label, expected, actual):
        """ Check one level of a potentially nested object / list """
        # object array is returned from cell array in mat file
        typex = type(expected)
        typac = type(actual)
        if isinstance(expected, ndarray) and expected.dtype.hasobject:
            assert typex is typac, "Different types at %s" % label
            assert len(expected) == len(actual), "Different list lengths at %s" % label
            for i, ev in enumerate(expected):
                level_label = "%s, [%d], " % (label, i)
                self._check_level(level_label, ev, actual[i])
            return
        # object, as container for matlab structs and objects
        elif isinstance(expected, mat_struct) or isinstance(expected, mat_obj):
            assert isinstance(actual, typex), \
                   "Different types %s and %s at %s" % (typex, typac, label)
            ex_fields = dir(expected)
            ac_fields = dir(actual)
            for k in ex_fields:
                if k.startswith('__') and k.endswith('__'):
                    continue
                assert k in ac_fields, "Missing property at %s" % label
                ev = expected.__dict__[k]
                v = actual.__dict__[k]
                level_label = "%s, property %s, " % (label, k)
                self._check_level(level_label, ev, v)
            return
        # hoping this is a single value, which might be an array
        if SP.issparse(expected):
            assert SP.issparse(actual), "Expected sparse at %s" % label
            assert_array_almost_equal(actual.todense(),
                                      expected.todense(),
                                      err_msg = label,
                                      decimal = 5)
        elif isinstance(expected, ndarray):
            if expected.shape: # allow scalar and 0d array comparisons
                assert isinstance(actual, ndarray), "Expected ndarray at %s" % label
            assert_array_almost_equal(actual, expected, err_msg=label, decimal=5)
        else:
            assert isinstance(expected, typac), \
                   "Types %s and %s do not match at %s" % (typex, typac, label)
            assert_equal(actual, expected, err_msg=label)
    
    def _check_case(self, name, files, case):
        for file_name in files:
            matdict = loadmat(file_name)
            label = "test %s; file %s" % (name, file_name)
            for k, expected in case.items():
                k_label = "%s, variable %s" % (label, k)
                assert k in matdict, "Missing key at %s" % k_label
                self._check_level(k_label, expected, matdict[k])

    # Add the load tests dynamically, with given parameters
    def _make_check_case(name, files, expected):
        def cc(self):
            self._check_case(name, files, expected)
        cc.__doc__ = "check loadmat case %s" % name
        return cc

    # Add the round trip tests dynamically, with given parameters
    def _make_rt_check_case(name, expected, format):
        def cc(self):
            mat_stream = StringIO()
            savemat(mat_stream, expected, format=format)
            mat_stream.seek(0)
            self._check_case(name, [mat_stream], expected)
        cc.__doc__ = "check loadmat case %s" % name
        return cc

    # Define cases to test
    theta = pi/4*arange(9,dtype=float)
    case_table4 = [
        {'name': 'double',
         'expected': {'testdouble': theta}
         }]
    case_table4.append(
        {'name': 'string',
         'expected': {'teststring': u'"Do nine men interpret?" "Nine men," I nod.'},
         })
    case_table4.append(
        {'name': 'complex',
         'expected': {'testcomplex': cos(theta) + 1j*sin(theta)}
         })
    A = zeros((3,5))
    A[0] = range(1,6)
    A[:,0] = range(1,4)
    case_table4.append(
        {'name': 'matrix',
         'expected': {'testmatrix': A},
         })
    case_table4.append(
        {'name': 'sparse',
         'expected': {'testsparse': SP.csc_matrix(A)},
         })
    B = A.astype(complex)
    B[0,0] += 1j
    case_table4.append(
        {'name': 'sparsecomplex',
         'expected': {'testsparsecomplex': SP.csc_matrix(B)},
         })
    case_table4.append(
        {'name': 'multi',
         'expected': {'theta': theta,
                      'a': A},
         })
    case_table4.append(
        {'name': 'minus',
         'expected': {'testminus': array(-1)},
         })
    case_table4.append(
        {'name': 'onechar',
         'expected': {'testonechar': u'r'},
         })
    case_table5 = [
        {'name': 'cell',
         'expected': {'testcell':
                      array([u'This cell contains this string and 3 arrays of '+\
                             'increasing length',
                             array(1), array([1,2]), array([1,2,3])], 
                            dtype=object)}
         }]
    case_table5.append(
        {'name': 'emptycell',
         'expected': {'testemptycell':
                      array([array(1), array(2), array([]),
                             array([]), array(3)], dtype=object)}
         })
    case_table5.append(
        {'name': 'stringarray',
         'expected': {'teststringarray': array(
        [u'one  ', u'two  ', u'three'], dtype=object)},
         })
    case_table5.append(
        {'name': '3dmatrix',
         'expected': {'test3dmatrix': transpose(reshape(range(1,25), (4,3,2)))}
         })
    case_table5_rt = [
        {'name': '3dmatrix',
         'expected': {'test3dmatrix': transpose(reshape(range(1,25), (4,3,2)))}
         }]
    st = mat_struct()
    st.stringfield = u'Rats live on no evil star.'
    st.doublefield = array([sqrt(2),exp(1),pi])
    st.complexfield = (1+1j)*array([sqrt(2),exp(1),pi])
    case_table5.append(
        {'name': 'struct', 
         'expected': {'teststruct': st}
         })
    a = array([array(1),
               array([array(2), array(3),
                      array([array(4), array(5)],
                            dtype=object)],
                     dtype=object)],
              dtype=object)
    case_table5.append(
        {'name': 'cellnest',
         'expected': {'testcellnest': a},
         })
    st = mat_struct()
    st.one = array(1)
    st.two = mat_struct()
    st.two.three = u'number 3'
    case_table5.append(
        {'name': 'structnest',
         'expected': {'teststructnest': st}
         })
    a = array([mat_struct(), mat_struct()])
    a[0].one = array(1)
    a[0].two = array(2)
    a[1].one = u'number 1'
    a[1].two = u'number 2'
    case_table5.append(
        {'name': 'structarr',
         'expected': {'teststructarr': a}
         })
    a = mat_obj()
    a._classname = 'inline'
    a.expr = u'x'
    a.inputExpr = u' x = INLINE_INPUTS_{1};'
    a.args = u'x'
    a.isEmpty = array(0)
    a.numArgs = array(1)
    a.version = array(1)
    case_table5.append(
        {'name': 'object',
         'expected': {'testobject': a}
         })
    u_str = file(
        os.path.join(test_data_path, 'japanese_utf8.txt'),
        'rb').read().decode('utf-8')
    case_table5.append(
        {'name': 'unicode',
        'expected': {'testunicode': u_str}
        })
    # add load tests
    for case in case_table4 + case_table5:
        name = case['name']
        expected = case['expected']
        filt = os.path.join(test_data_path, 'test%s_*.mat' % name)
        files = glob(filt)
        assert files, "No files for test %s using filter %s" % (name, filt)
        exec 'check_%s = _make_check_case(name, files, expected)' % name
    # round trip tests
    for case in case_table4 + case_table5_rt:
        name = case['name'] + '_round_trip'
        expected = case['expected']
        format = case in case_table4 and '4' or '5'
        exec 'check_%s = _make_rt_check_case(name, expected, format)' \
             % name


if __name__ == "__main__":
    NumpyTest().run()

