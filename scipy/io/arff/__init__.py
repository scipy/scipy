"""
Module to read ARFF files, which are the standard data format for WEKA.

ARFF is a text file format which support numerical, string and data values.
The format can also represent missing data and sparse data.

See the `WEKA website
<http://weka.wikispaces.com/ARFF>`_
for more details about arff format and available datasets.

Examples
--------

>>> from scipy.io import arff
>>> from cStringIO import StringIO
>>> content = \"\"\"
... @relation foo
... @attribute width  numeric
... @attribute height numeric
... @attribute color  {red,green,blue,yellow,black}
... @data
... 5.0,3.25,blue
... 4.5,3.75,green
... 3.0,4.00,red
... \"\"\"
>>> f = StringIO(content)
>>> data, meta = arff.loadarff(f)
>>> data
array([(5.0, 3.25, 'blue'), (4.5, 3.75, 'green'), (3.0, 4.0, 'red')],
      dtype=[('width', '<f8'), ('height', '<f8'), ('color', '|S6')])
>>> meta
Dataset: foo
\twidth's type is numeric
\theight's type is numeric
\tcolor's type is nominal, range is ('red', 'green', 'blue', 'yellow', 'black')

"""

from arffread import *
import arffread

__all__ = arffread.__all__

from numpy.testing import Tester
test = Tester().test
