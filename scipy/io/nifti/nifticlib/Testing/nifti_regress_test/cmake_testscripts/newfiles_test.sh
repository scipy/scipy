#!/bin/sh

if [ $# -lt 2 ]
then
echo Missing nifti tool and Binary directory name
exit 1
fi

NT=$1
DATA=$2

rm -f ${DATA}/new* ${DATA}/ncopy*

# test writing various output file types
${NT} -make_im -prefix ${DATA}/new1.hdr
${NT} -make_im -prefix ${DATA}/new2.hdr.gz
${NT} -make_im -prefix ${DATA}/new3.img.gz
${NT} -make_im -prefix ${DATA}/new4.nii.gz
${NT} -make_im -prefix ${DATA}/new5.nia
${NT} -make_im -prefix ${DATA}/new6.nia.gz

# test reading them
${NT} -copy_im -prefix ${DATA}/ncopy1.nii -infiles ${DATA}/new1.hdr
${NT} -copy_im -prefix ${DATA}/ncopy2.nii -infiles ${DATA}/new2.hdr.gz
${NT} -copy_im -prefix ${DATA}/ncopy3.nii -infiles ${DATA}/new3.img.gz
${NT} -copy_im -prefix ${DATA}/ncopy4.nii -infiles ${DATA}/new4.nii.gz
${NT} -copy_im -prefix ${DATA}/ncopy5.nii -infiles ${DATA}/new5.nia
${NT} -copy_im -prefix ${DATA}/ncopy6.nii -infiles ${DATA}/new6.nia.gz

# verify that they are all the same
set count = 0
for index in 2 3 4 5 6 ; do
    diff ${DATA}/ncopy1.nii ${DATA}/ncopy$index.nii
    if [ $? -ne 0 ] ; then 
        echo "-- failure on test index $index --"
        exit 1
    fi
done
