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

# Define a C macro LOCALEDIR indicating where catalogs will be installed.
DEFS = -DLOCALEDIR=\"$(localedir)\" @DEFS@

CPPFLAGS += $(WX_CPPFLAGS)
CFLAGS = $(WX_CFLAGS_ONLY)
CXXFLAGS += $(WX_CXXFLAGS_ONLY)

# Link time dependencies.
LDADD = $(WX_LDADD)
LIBS += $(WX_LIBS)

# Additional files to be distributed.
EXTRA_DIST = autogen.sh autoclean.sh
