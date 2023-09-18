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
hello_SOURCES = hello.cc hellowindow.cc

# Use 'moc' for QObject subclasses.
METASOURCES = AUTO

# Define a C macro LOCALEDIR indicating where catalogs will be installed.
DEFS = -DPREFIX=\"$(prefix)\" @DEFS@

# Make sure the KDE, Qt, X include file is found.
AM_CPPFLAGS = $(KDE_INCLUDES) $(QT_INCLUDES) $(X_INCLUDES)

# CXXFLAGS += $(QT_CXXFLAGS) $(PTHREAD_CFLAGS)

# Link time dependencies.
LDADD = $(KDE_LDFLAGS) $(QT_LDFLAGS)
LIBS += $(LIB_KDEUI) $(LIB_KDECORE) $(LIB_QT)

# Generic KDE build mechanism.

$(top_srcdir)/configure.in: configure.in.in $(top_srcdir)/subdirs
	cd $(top_srcdir) && $(MAKE) -f admin/Makefile.common configure.in ;

$(top_srcdir)/subdirs:
	cd $(top_srcdir) && $(MAKE) -f admin/Makefile.common subdirs

$(top_srcdir)/acinclude.m4: $(top_srcdir)/admin/acinclude.m4.in $(top_srcdir)/admin/libtool.m4.in
	@cd $(top_srcdir) && cat admin/acinclude.m4.in admin/libtool.m4.in > acinclude.m4

MAINTAINERCLEANFILES = subdirs configure.in acinclude.m4 configure.files 

package-messages:
	$(MAKE) -f admin/Makefile.common package-messages
	$(MAKE) -C po merge

dist-hook:
	cd $(top_distdir) && perl admin/am_edit -padmin
	cd $(top_distdir) && $(MAKE) -f admin/Makefile.common subdirs

# Additional files to be distributed.
EXTRA_DIST = admin/config.rpath  autogen.sh autoclean.sh
