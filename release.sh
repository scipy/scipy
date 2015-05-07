#! /bin/sh
# script to build tarballs, mac os x and windows installers on mac os x

export PYTHONPATH=''

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
paver dmg -p 3.4  # 32/64-bit version

paver bdist_superpack -p 2.7
paver bdist_superpack -p 3.3
paver bdist_superpack -p 3.4

paver write_release_and_log
