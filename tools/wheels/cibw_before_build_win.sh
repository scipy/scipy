set -xe

PROJECT_DIR="$1"
PLATFORM=$(PYTHONPATH=tools python -c "import openblas_support; print(openblas_support.get_plat())")

printenv
# Update license
cat $PROJECT_DIR/tools/wheels/LICENSE_win32.txt >> $PROJECT_DIR/LICENSE.txt

# Install Openblas
PYTHONPATH=tools python -c "import openblas_support; openblas_support.make_init('scipy')"
mkdir -p /c/opt/openblas/if_32/32/lib/pkgconfig
mkdir -p /c/opt/openblas/if_32/64/lib/pkgconfig

# delvewheel is the equivalent of delocate/auditwheel for windows.
python -m pip install delvewheel

# make the DLL available for tools/wheels/repair_windows.sh. If you change
# this location you need to alter that script.
mkdir -p /c/opt/openblas/openblas_dll
which strip

# The 32/64 bit Fortran wheels are currently coming from different locations.
if [[ $PLATFORM == 'win-32' ]]; then
  # 32-bit openBLAS
  # Download 32 bit openBLAS and put it into c/opt/32/lib
  target=$(python -c "import tools.openblas_support as obs; plat=obs.get_plat(); ilp64=obs.get_ilp64(); target=f'openblas_{plat}.zip'; obs.download_openblas(target, plat, ilp64);print(target)")
  unzip $target -d /c/opt/openblas/if_32/
  cp /c/opt/openblas/if_32/32/bin/*.dll /c/opt/openblas/openblas_dll
  # rm /c/opt/openblas/if_32/32/lib/*.dll.a
else
  # 64-bit openBLAS
  curl -L https://github.com/scipy/scipy-ci-artifacts/raw/main/openblas_32_if.zip -o openblas_32_if.zip
  unzip openblas_32_if.zip -d /c
  cp /c/opt/openblas/if_32/64/bin/*.dll /c/opt/openblas/openblas_dll
fi
