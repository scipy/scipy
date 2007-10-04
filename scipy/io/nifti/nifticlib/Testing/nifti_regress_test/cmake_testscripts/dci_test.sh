#!/bin/sh

if [ $# -lt 2 ]
then
echo Missing nifti tool and Binary directory name
exit 1
fi

NT=$1
DATA=$2

# get the time series for a slice, get the same sub-bricks,
# compare nifti_images, and display the time series again
#
# this should match the previous

${NT} -keep_hist -cci 3 3 3 -1 -1 -1 -1 -prefix ${DATA}/r.333 -infiles ${DATA}/run.210.nii
if [ $? -ne 0 ] ; then
echo cci failed
fi

${NT} -keep_hist -cci 3 3 3 -1 1 1 1 -prefix ${DATA}/r.333.111 -infiles ${DATA}/run.210.nii
if [ $? -ne 0 ] ; then
echo cci failed
fi

${NT} -disp_ci -1 -1 -1 -1 -1 -1 -1 -infiles ${DATA}/r.333.111.nii -quiet \
        | tee o.10a.dci.1
if [ $? -ne 0 ] ; then
echo disp_ci failed
fi

${NT} -disp_ci  0  0  0 -1  1  1  1 -infiles ${DATA}/r.333.nii -quiet \
        | tee o.10a.dci.2
if [ $? -ne 0 ] ; then
echo disp_ci failed
fi

if diff o.10a.dci.?
then
echo ""
rm -f o.10a.dci* ${DATA}/r.333*
else
echo '** failure, o.10a.dci timeseries files differ'
exit 1
fi
