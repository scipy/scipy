# makefile for cython extensions
# Before working, do:
# rm -rf ~/sp/scipy 
# ~/sp links to ~/usr/local ... site-packages
# cd ~/dev_trees/scipy/scipy/io/matlab
# make ext
# cd ~/dev_trees/scipy
# python setup.py install --prefix=~/usr/local
# make development copy be copy on path
# cd ~/sp/scipy/io
# rm -rf matlab
# ln -s ~/dev_trees/scipy/scipy/io/matlab

clean:
	-rm -rf build *.so

ext: clean
	cython mio5_utils.pyx
	cython mio_utils.pyx
	python setup.py build_ext --inplace

