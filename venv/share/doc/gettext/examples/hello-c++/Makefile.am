# Example for use of GNU gettext.
# This file is in the public domain.
#
# Makefile configuration - processed by automake.

# General automake options.
AUTOMAKE_OPTIONS = foreign no-dependencies
ACLOCAL_AMFLAGS = -I m4

# The list of subdirectories containing Makefiles.
SUBDIRS = m4 po

# The list of programs that are built.
bin_PROGRAMS = hello

# The source files of the 'hello' program.
hello_SOURCES = hello.cc gettext.h

# Define a C macro LOCALEDIR indicating where catalogs will be installed.
DEFS = -DLOCALEDIR=\"$(localedir)\" @DEFS@

# Make sure the gettext.h include file is found.
AM_CPPFLAGS = -I. -I$(srcdir)

# Link time dependencies.
LDADD = @LIBASPRINTF@ @LIBINTL@

# Additional files to be distributed.
EXTRA_DIST = autogen.sh autoclean.sh
