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

# Compilation of librep programs.
SUFFIXES = .jl .jlc
.jl.jlc:
	rep --batch --no-rc -l compiler -f compile-batch $<
CLEANFILES = *.jlc

# Making a librep program executable.
hello: hello.jlc
	(echo '#!@REP@'; echo '!#'; cat $<) > $@
	chmod a+x $@
CLEANFILES += hello

# Additional files to be distributed.
EXTRA_DIST = autogen.sh autoclean.sh
