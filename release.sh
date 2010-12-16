#! /bin/sh
# script to build tarballs, mac os x and windows installers on mac os x

# TODO: fix dmg paver task. Then uncomment the lines below - they are for
# building the docs.
#paver bootstrap
#source bootstrap/bin/activate
#python setupsconsegg.py install
#paver dmg -p 2.6
#paver dmg -p 2.5
#export MACOSX_DEPLOYMENT_TARGET=10.6
#paver dmg -p 2.7


#------------------------------------------------------------------
# Build tarballs, Windows and 64-bit OS X installers (on OS X 10.6)
#------------------------------------------------------------------
paver sdist

paver bdist_superpack -p 3.1
paver bdist_superpack -p 2.7
paver bdist_superpack -p 2.6
paver bdist_superpack -p 2.5

export MACOSX_DEPLOYMENT_TARGET=10.6
paver simple_dmg -p 2.7

#--------------------------------------------
# Build 32-bit OS X installers (on OS X 10.5)
#--------------------------------------------
#export MACOSX_DEPLOYMENT_TARGET=10.3
#paver simple_dmg -p 2.6
#paver simple_dmg -p 2.7
#export CC=/usr/bin/gcc-4.0  # necessary on 10.6, not sure about 10.5
#paver simple_dmg -p 2.5


paver write_release_and_log

