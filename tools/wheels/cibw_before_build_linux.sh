set -xe


NIGHTLY_FLAG=""

while (( "$#" )); do
  case "$1" in
    --nightly)
      NIGHTLY_FLAG="--nightly"
      shift
      ;;
    --*)
      echo "Error: Unrecognized flag $1" >&2
      exit 1
      ;;
    *)
      break
      ;;
  esac
done

PROJECT_DIR="$1"
PLATFORM=$(PYTHONPATH=tools python -c "import openblas_support; print(openblas_support.get_plat())")

printenv
# Update license
cat $PROJECT_DIR/tools/wheels/LICENSE_linux.txt >> $PROJECT_DIR/LICENSE.txt

# Install Openblas
basedir=$(python tools/openblas_support.py $NIGHTLY_FLAG)
cp -r $basedir/lib/* /usr/local/lib
cp $basedir/include/* /usr/local/include
