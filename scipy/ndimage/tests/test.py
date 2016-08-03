from __future__ import division, print_function, absolute_import
from scipy import ndimage
import numpy as np
import os.path
import warnings

import numpy as np
from numpy.testing import (assert_, assert_array_almost_equal, assert_equal,
                           assert_almost_equal, assert_array_equal,
                           assert_raises, run_module_suite, TestCase)

import scipy.ndimage as ndimage


def test_find_objects09():
  "test the objects found for a set of labelled images"
  data = np.loadtxt(os.path.join(os.path.dirname(__file__), "data", "label_inputs.txt"))
  strels = np.loadtxt(os.path.join(os.path.dirname(__file__), "data", "label_strels.txt"))
  results = np.loadtxt(os.path.join(os.path.dirname(__file__), "data", "label_results.txt"))
  objects_found = np.loadtxt(os.path.join(os.path.dirname(__file__), "data", "find_objects.txt"))
  data = data.reshape((-1, 7, 7))
  strels = strels.reshape((-1, 3, 3))
  results = results.reshape((-1, 7, 7))
  r = 0
  for i in range(data.shape[0]):
      d = data[i, :, :]
      for j in range(strels.shape[0]):
          s = strels[j, :, :]
          img = ndimage.label(d, s)[0]
          objects = objects.join(str(objects.append(ndimage.find_objects(img))))
          r += 1
  assert_equal(str(objects), str(results))

test_find_objects09()