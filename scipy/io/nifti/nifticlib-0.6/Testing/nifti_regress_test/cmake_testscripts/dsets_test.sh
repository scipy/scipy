#!/bin/sh

if [ $# -lt 2 ]
then
echo Missing nifti tool and Binary directory name
exit 1
fi

NT=$1
DATA=$2


# - create datasets out of nothing
#
# - modify some fields and compare against other datasets
rm -f ${DATA}/new_epi.nii
# just test a basic make im, mostly to capture the debug output
${NT} -make_im -debug 3 -new_dim 4 64 64 21 180 0 0 0 -prefix ${DATA}/new_epi.nii
if [ $? -ne 0 ] ; then
echo failed
exit 1
fi

# compare hdr and nim in a fresh image to the existing one
${NT} -diff_hdr -new_dim 4 64 64 21 180 0 0 0      \
           -infiles MAKE_IM ${DATA}/stat0.nii
if [ $? = 0 ] ; then
echo unexpected 0 return code in diff
exit 1
fi

${NT} -diff_nim -new_dim 4 64 64 21 180 0 0 0      \
           -infiles MAKE_IM ${DATA}/stat0.nii
if [ $? = 0 ] ; then
echo unexpected 0 return code in diff
exit 1
fi


rm -f ${DATA}/epi_180_pixdim.nii

# clean up the nim by adjusting pixdim (from empty MAKE_IM)
${NT} -mod_hdr -new_dim 4 64 64 21 180 0 0 0               \
           -mod_field pixdim '0.0 4.0 4.0 6.0 3.0 1.0 1.0 1.0'  \
           -prefix ${DATA}/epi_180_pixdim.nii -infiles MAKE_IM
if [ $? -ne 0 ] ; then
echo mod_hdr failed
exit 1
fi


# and compare again
${NT} -diff_nim -infiles ${DATA}/stat0.nii ${DATA}/epi_180_pixdim.nii
if [ $? -ne 0 ] ; then
echo found changes -- success.
else
echo diff failed to find differences.
exit 1
fi

exit 0
