all:

distclean:
	-rm MANIFEST Changelog
	-rm nifti/*.{c,pyc,so} nifti/nifticlib.py
	-rm tests/*.pyc
	-rm -r build
	-rm -r dist

	
orig-src: distclean 
	# clean existing dist dir first to have a single source tarball to process
	-rm -rf dist
	# the debian changelog is also the upstream changelog
	cp debian/changelog Changelog

	# update manpages
	help2man -N -n "compute peristimulus timeseries of fMRI data" \
		bin/pynifti_pst > man/pynifti_pst.1

	if [ ! "$$(dpkg-parsechangelog | egrep ^Version | cut -d ' ' -f 2,2 | cut -d '-' -f 1,1)" == "$$(python setup.py -V)" ]; then \
			printf "WARNING: Changelog version does not match tarball version!\n" ;\
			exit 1; \
	fi
	# let python create the source tarball
	python setup.py sdist --formats=gztar
	# rename to proper Debian orig source tarball and move upwards
	# to keep it out of the Debian diff
	file=$$(ls -1 dist); ver=$${file%*.tar.gz}; ver=$${ver#pynifti-*}; mv dist/$$file ../pynifti_$$ver.orig.tar.gz

bdist_wininst:
	# THIS IS ONLY FOR WINDOWS!
	# Consider this a set of notes on how to build PyNIfTI on win32, rather
	# than an actually working target
	#
	# assumes Dev-Cpp to be installed at C:\Dev-Cpp
	python setup.py build_ext -c mingw32 --swig-opts "-C:\Dev-Cpp\include/nifti -DWIN32" -IC:\Dev-Cpp\include nifti
	
	# for some stupid reason the swig wrapper is in the wrong location
	move /Y nifticlib.py nifti
	
	# now build the installer
	python setup.py bdist_wininst
