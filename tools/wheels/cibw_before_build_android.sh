set -xe

PROJECT_DIR="$1"
openblas_dir="$PROJECT_DIR/.openblas"
rm -rf $openblas_dir
mkdir -p $openblas_dir

pip install \
    --target $openblas_dir \
    --extra-index-url https://chaquo.com/pypi-upstream/ \
    chaquopy-openblas==0.3.33
mv $openblas_dir/chaquopy/* $openblas_dir
rmdir $openblas_dir/chaquopy
