#! /bin/sh
# script to build tarballs, mac os x and windows installers on mac os x

# bootstrap needed to ensure we build the docs from the right scipy version
paver bootstrap
source bootstrap/bin/activate
python setupsconsegg.py install

# we need to copy Sphinx extensions from the numpy source tree
mkdir doc/sphinxext
cp -R ../numpy/doc/sphinxext/ doc/sphinxext/

# build docs
paver pdf

#------------------------------------------------------------------
# Build tarballs, Windows and 64-bit OS X installers (on OS X 10.6)
#------------------------------------------------------------------

paver sdist

export MACOSX_DEPLOYMENT_TARGET=10.6
paver dmg -p 2.7  # 32/64-bit version

paver bdist_superpack -p 3.1
paver bdist_superpack -p 2.7
paver bdist_superpack -p 2.6
paver bdist_superpack -p 2.5


#--------------------------------------------
# Build 32-bit OS X installers (on OS X 10.5)
#--------------------------------------------
#export MACOSX_DEPLOYMENT_TARGET=10.3
#paver dmg -p 2.6
#paver dmg -p 2.7  # 32-bit version
#export CC=/usr/bin/gcc-4.0  # necessary on 10.6, not sure about 10.5
#paver dmg -p 2.5


paver write_release_and_log

