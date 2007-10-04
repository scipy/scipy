#!/bin/sh
#
# arg[1] is TESTING_BINARY_DIR
if [ $# -lt 1 ]
then
echo Missing Binary directory name
exit 1
fi


if cd $1
then
echo working in `pwd`
else
echo can\'t cd to $1
exit 1
fi


if wget -nd http://nifti.nimh.nih.gov/pub/dist/data/nifti_regress_data.tgz
then
echo wget succeeded
else
echo wget failed
exit 1
fi

if tar xzvf nifti_regress_data.tgz
then
echo ""
else
echo failed tar xzvf nifti_regress_data.tgz
exit 1
fi

if rm -f nifti_regress_data.tgz
then
echo ""
else
echo can\'t remove ../nifti_regress_data.tgz
exit 1
fi

exit 0

