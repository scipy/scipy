#!/bin/sh

if [ $# -lt 2 ]
then
echo Missing nifti tool and Binary directory name
exit 1
fi

NT=$1
DATA=$2

${NT} -quiet -disp_ts 19 36 11 -infiles ${DATA}/stat0.nii | tee o.08.ts.19.36.11

if [ $? -ne 0 ] ; then
echo disp_ts failed
exit 1
fi
rm -f o.08.ts.19.36.11

${NT} -quiet -disp_ts 19 36 11 -infiles ${DATA}/stat0.nii \
        | awk '{print $(NF-1), $NF, $1, $2}' | tee o.09.ts4.1.awk

if [ $? -ne 0 ] ; then
echo disp_ts failed
exit 1
fi

${NT} -quiet -disp_ts 19 36 11 -infiles ${DATA}/f4.nii | tee o.09.ts4.2.awk

if diff o.09.ts4.1.awk o.09.ts4.2.awk
then
echo ""
else
echo '** failure, ts4 files differ'
exit 1
fi
rm -f 0.09.ts4.1.awk o.09.ts4.2.awk ${DATA}/f4.nii
exit 0

# get the time series for a slice, get the same sub-bricks,
# compare nifti_images, and display the time series again
#
# this should match the previous

if ${NT} -keep_hist -cci 19 36 -1 -1 0 0 0 -prefix ${DATA}f.19.36 \
-infiles ${DATA}/stat0.nii
then
echo ""
else
echo cci failed
exit 1
fi

if ${NT} -keep_hist -cbl -prefix ${DATA}/f.19.36.t4.nii \
           -infiles ${DATA}/f.19.36.nii'[178..$,0,1]'
then
echo ""
else
echo cbi failed
exit 0
fi

${NT} -diff_nim -infiles ${DATA}/f.19.36.nii ${DATA}f.19.36.t4.nii | tee o.10.diff_nim
if [ $? -ne 0 ] ; then
echo f.19.36.nii and f.19.36.t4.nii differ
else
echo f.19.36.nii and f.19.36.t4.nii do not differ
fi

${NT} -quiet -disp_ci 0 0 11 -1 0 0 0 -infiles f.19.36.t4.nii \
          | tee o.10.dci.4
#diff o.09.ts4.1.awk o.10.dci.4
#if( $status ) echo '** failure, o.09 and o.10 ts files differ'

exit 0
