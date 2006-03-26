#!/bin/bash

# Build and install numpy and scipy
pushd $HOME/Install/numpy
rm -Rf build
python setup.py build 
sudo rm -Rf /usr/lib/python2.4/site-packages/numpy
sudo rm -Rf /usr/lib/python2.4/site-packages/scipy
sudo python setup.py install
popd
rm -Rf build
python setup.py build 
sudo python setup.py install

# Make PySparse symlink
pushd /usr/lib/python2.4/site-packages/scipy/sandbox/pysparse
sudo ln -s pysparse.so spmatrix.so
popd

# Test
pushd $HOME
#source newenv.sh
python -c "import scipy; scipy.test()"

# Build and install matplotlib
pushd $HOME/Install/matplotlibsvn
rm -Rf build
python setup.py build 
sudo rm -Rf /usr/lib/python2.4/site-packages/matplotlib
sudo python setup.py install

