lfortran -c scipy/special/specfun/specfun.f -o scipy/special/specfun.o --fixed-form --implicit-typing --implicit-interface --generate-object-code --rtlib
clang -shared -fPIC -o scipy/special/specfun.so scipy/special/specfun.o -L/Users/pgoswami/Desktop/lfortran/src/runtime -llfortran_runtime
cp scipy/special/specfun.so /Users/pgoswami/Desktop/scipy/build-install/lib/python3.10/site-packages/scipy/special/
cp scipy/special/specfun.so $CONDA_PREFIX/lib/scipy/special/
python dev.py test -t scipy.special