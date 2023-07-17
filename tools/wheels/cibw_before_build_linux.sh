set -xe

NIGHTLY_FLAG=""

usage() {
    echo "Usage: $0 [--nightly] <project_dir>"
    exit 1
}

TEMP=$(getopt -o n -l nightly -- "$@")
if [ $? != 0 ]; then
    usage
fi

eval set -- "$TEMP"

while true ; do
    case "$1" in
        -n|--nightly)
            NIGHTLY_FLAG="--nightly"; shift ;;
        --) shift ; break ;;
    esac
done

if [ -z "$1" ]; then
    usage
else
    PROJECT_DIR="$1"
fi


PLATFORM=$(PYTHONPATH=tools python -c "import openblas_support; print(openblas_support.get_plat())")

printenv
# Update license
cat $PROJECT_DIR/tools/wheels/LICENSE_linux.txt >> $PROJECT_DIR/LICENSE.txt

# Install Openblas
basedir=$(python tools/openblas_support.py $NIGHTLY_FLAG)
cp -r $basedir/lib/* /usr/local/lib
cp $basedir/include/* /usr/local/include
