This contains a few hacks to re-use BOOST C++ library
[tests data](https://github.com/boostorg/math/tree/develop/test) for
validating scipy.special functions.

convert.py
----------

This script parses the BOOST data and writes the data into a CVS text file.
For each data file, a subdirectory is created, as one boost data file (.ipp)
can have several test sets.

From the root of this repo
```
# Clone boostorg repo
git clone --depth=1 https://github.com/boostorg/math.git boostmath

# Remove existing data
rm -rf scipy/special/tests/data/boost/*

# Run the coverter script (potentially also update exclude regexes)
python scipy/special/utils/convert.py

# Verify all the new files are used in test_data.py
git diff --stat HEAD | grep "Bin 0"
git diff HEAD -- scipy/special/tests/test_data.py
```

It may be desirable to remove whitespace only changes (see
[#12357](https://github.com/scipy/scipy/pull/12357)) for instructions on that.
