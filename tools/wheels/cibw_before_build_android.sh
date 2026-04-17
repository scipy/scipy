set -xe

PROJECT_DIR="$1"
openblas_dir="$PROJECT_DIR/.openblas"
rm -rf $openblas_dir
mkdir -p $openblas_dir

openblas_version=0.2.20
pip install \
    --target $openblas_dir \
    --extra-index-url https://chaquo.com/pypi-13.1/ \
    chaquopy-openblas==$openblas_version
mv $openblas_dir/chaquopy/* $openblas_dir
rmdir $openblas_dir/chaquopy

# There is a pkgconfig file, but it's not usable because it doesn't use the
# `prefix` variable, so overwrite it.
cat <<EOF > $openblas_dir/lib/pkgconfig/openblas.pc
prefix=/usr
Name: openblas
Description: OpenBLAS
Version: ${openblas_version}
Cflags: -I\${prefix}/include
Libs: -L\${prefix}/lib -lopenblas
EOF
