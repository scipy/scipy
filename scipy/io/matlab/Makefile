# makefile for cython extensions

clean:
	-rm -rf build *.so

ext: clean
	cython mio5_utils.pyx
	python setup.py build_ext --inplace

