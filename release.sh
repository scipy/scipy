#! /bin/sh
# script to build tarballs, mac os x and windows installers on mac os x

export PYTHONPATH=''

# Check we're using the correct g++/c++ for the 64-bit 2.7 dmg.
# We do this because for Python 2.6 we use a symlink on the PATH to select
# /usr/bin/g++-4.0, while for Python 2.7 we need the default 4.2 version.
export PATH=~/Code/tmp/gpp40temp/:$PATH
gpp="$(g++ --version | grep "4.0")"
if [ -z "$gpp" ]; then
    echo "Wrong g++ version, we need 4.0 to compile scipy with Python 2.6"
    exit 1
fi

# bootstrap needed to ensure we build the docs from the right scipy version
paver bootstrap
source bootstrap/bin/activate
python setupegg.py install

# build docs
paver pdf

#------------------------------------------------------------------
# Build tarballs, Windows and 64-bit OS X installers (on OS X 10.6)
#------------------------------------------------------------------

paver sdist

export MACOSX_DEPLOYMENT_TARGET=10.6
export PATH=~/Code/tmp/gpp42temp/:$PATH
gpp="$(g++ --version | grep "4.2")"
if [ -z "$gpp" ]; then
    echo "Wrong g++ version, we need 4.2 for 64-bit binary for Python 2.7"
    exit 1
fi
paver dmg -p 2.7  # 32/64-bit version
paver dmg -p 3.3  # 32/64-bit version

paver bdist_superpack -p 2.7
paver bdist_superpack -p 2.6
paver bdist_superpack -p 3.2
paver bdist_superpack -p 3.3


#--------------------------------------------
# Build 32-bit OS X installers (on OS X 10.5)
#--------------------------------------------
#export MACOSX_DEPLOYMENT_TARGET=10.3
#paver dmg -p 2.6
#paver dmg -p 2.7  # 32-bit version


paver write_release_and_log
