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
all-local: hello
CLEANFILES = hello

# The Smalltalk program is installed as data and invoked through a shell script.
pkgdata_DATA = hello.st
hello:
	{ echo '#!/bin/sh'; \
	  echo "exec '@GST@' -Q '$(pkgdatadir)/hello.st' \"\$$@\""; \
	} > $@

# Additional files to be distributed.
EXTRA_DIST = autogen.sh autoclean.sh
