#!/bin/sh
# Make a snapshot of the latest SciPy from the CVS repository.
# Tar file SciPy-$VERSION.tar.gz will be moved into
# the current directory. Usage: just run cvs2tar.sh.
#
# Author: Pearu Peterson, April 2002


PYTHON=python
CVSROOT=:pserver:anonymous@scipy.org:/home/cvsroot
echo
echo "WARNING: this script will retrive ~30MB from scipy.org"
echo "Press ENTER and Ctrl-C to exit the script."
echo
echo "CVS password is blank. Just hit ENTER to carry on:"
cvs -d $CVSROOT login || exit 0
cvs -d $CVSROOT checkout scipy || exit 0
VERSION=`cd scipy && $PYTHON -c 'from scipy_version import *;print scipy_version' && cd -`
TARFILE=SciPy-$VERSION.tar.gz
echo "TARFILE=$TARFILE"

cd scipy && $PYTHON setup.py sdist && cd -

mv -vf scipy/dist/$TARFILE $TARFILE
rm -rf scipy