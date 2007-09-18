
TAR     ?= /bin/tar
RM      ?= /bin/rm

# get the regression testing data tree
regress_data: nifti_regress_test/nifti_regress_data

# remove any result directories
regress_clean:
	$(RM) -fr nifti_regress_test/results*

# remove any result directories, and remove the data tree
regress_clean_all: regress_clean
	$(RM) -fr nifti_regress_test/nifti_regress_data

# download and unpack the regression data tree
nifti_regress_test/nifti_regress_data:
	wget http://nifti.nimh.nih.gov/pub/dist/data/nifti_regress_data.tgz
	( cd nifti_regress_test; $(TAR) -xzf ../nifti_regress_data.tgz; )
	$(RM) nifti_regress_data.tgz

