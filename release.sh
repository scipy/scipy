#! /bin/sh
# script to build tarballs, mac os x and windows installers on mac os x

# copy sphinxext over from numpy, otherwise docs won't build
cp -R ../numpy/doc/sphinxext/ doc/sphinxext/

# launch a virtualenv to install scipy into so we can build docs
paver bootstrap
source bootstrap/bin/activate
python setupsconsegg.py install

paver sdist
paver simple_dmg -p 2.6

# necessary to set CC for Python 2.5 on Snow Leopard
export CC=/usr/bin/gcc-4.0
paver simple_dmg -p 2.5
export CC=''

# XXX: bdist_superpack doesn't take version option.
paver bdist_superpack 
# change PYVER to 2.5, then build another superpack
#paver write_note_changelog
