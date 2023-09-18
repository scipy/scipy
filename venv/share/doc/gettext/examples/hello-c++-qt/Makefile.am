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
hello_SOURCES = hello.cc

# Define C macros LOCALEDIR and PKGLOCALEDIR indicating where catalogs will be
# installed.
pkglocaledir = $(pkgdatadir)/locale
DEFS = -DLOCALEDIR=\"$(localedir)\" -DPKGLOCALEDIR=\"$(pkglocaledir)\" @DEFS@

CXXFLAGS += $(QT_CXXFLAGS) $(PTHREAD_CFLAGS)

# Link time dependencies.
LDADD = $(QT_LDFLAGS)
LIBS += $(QT_LIBS)

# Additional files to be distributed.
EXTRA_DIST = autogen.sh autoclean.sh
