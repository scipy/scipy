This folder contains scripts and utilities to run performance profiles of pyPRIMA Python bindings against the PRIMA Python bindings.

In general we expect the performance profiles to show a straight line for both the output based and history based profiles since the underlying code is the same. However since there are bound to be some differences in how Python and Fortran handle floating point arithmetic we expect to see a few cases where one might perform slightly "better" than the other.

In order to use this you will need both pyPRIMA and the Python bindings available in your environment. See the respective folders for installation instructions. Once they are both available you can simple run `python profiles.py`.
