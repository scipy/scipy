set -xe

WHEEL="$1"
DEST_DIR="$2"
OPENBLAS_DIR=$(python -c"import scipy_openblas32 as sop; print(sop.get_lib_dir())")

delvewheel repair --add-path $OPENBLAS_DIR --no-dll libsf_error_state.dll -w $DEST_DIR $WHEEL
