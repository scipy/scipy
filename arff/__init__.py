"""Module to read arff files (weka format).

arff is a simple file format which support numerical, string and data values.
It supports sparse data too.

See http://weka.sourceforge.net/wekadoc/index.php/en:ARFF_(3.4.6) for more
details about arff format and available datasets."""

from arffread import *
import arffread

__all__ = arffread.__all__
