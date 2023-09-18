# Example for use of GNU gettext.
# This file is in the public domain.
#
# Makefile configuration - processed by automake.

# General automake options.
AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4

# The list of subdirectories containing Makefiles.
SUBDIRS = m4 po

# The list of programs that are built.
bin_SCRIPTS = hello

# Making a Python program executable.
# FIXME: How to make a compiled python program (output of automake's
# py-compile command) executable?
hello: hello.py
	(echo '#!@PYTHON@'; cat $<) > $@
	chmod a+x $@
CLEANFILES = hello
#CLEANFILES += *.pyc *.pyo

# Additional files to be distributed.
EXTRA_DIST = autogen.sh autoclean.sh
