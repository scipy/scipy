#!/bin/sh

if [ $# -lt 2 ]
then
echo Missing nifti tool and Binary directory name
exit 1
fi

NT=$1
DATA=$2


# add some comment and afni extensions, then display them
if \
${NT} -keep_hist -prefix ${DATA}/f4.comment -infiles ${DATA}/anat0.nii   \
           -add_comment '4 slice time series'                      \
           -add_afni_ext 'and an AFNI extension'                   \
           -add_comment 'how about a question AND a comment?'
then
echo ""
else
echo "add comment failed"
exit 1
fi

if ${NT} -disp_ext -infiles ${DATA}/f4.comment.nii
then
echo ""
else
echo "failed"
exit 1
fi

if ${NT} -cbl -infiles ${DATA}/f4.comment.nii -prefix ${DATA}/f4.to.clear.nii
then
echo ""
else
echo "failed"
exit 1
fi

if ${NT} -overwrite -strip -infiles ${DATA}/f4.to.clear.nii
then
echo ""
else
echo "failed"
exit 1
fi

if ${NT} -disp_ext -infiles ${DATA}/f4.to.clear.nii
then
echo ""
else
echo "failed"
exit 1
fi

if ${NT} -diff_nim -infiles ${DATA}/f4.comment.nii ${DATA}/f4.to.clear.nii
then
echo "failed -- no changes found"
exit 1
else
echo "diff succeed -- found changes"
fi

