
# - copy some of the sub-bricks of the default stat file,
#   storing the history (command), and writing to a new 'f4.nii'
# - compare these nifti_image structs

if [ $# -lt 2 ]
then
echo Missing nifti tool and Binary directory name
exit 1
fi

NT=$1
DATA=$2

rm -fr ${DATA}/f4*

if ${NT} -keep_hist -cbl -infiles ${DATA}/stat0.nii'[178..$,0,1]' -prefix ${DATA}/f4
then
echo ""
else
echo "copy sub-bricks failed"
fi

if ${NT} -diff_nim -infiles ${DATA}/stat0.nii ${DATA}/f4.nii
then
echo diff failed '(no diffs found'
exit 1
else
echo diff succeeded
fi

exit 0
