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
bin_CSHARPPROGRAMS = hello

# The source files of the 'hello' program.
hello_SOURCES = hello.cs

# The link dependencies of the 'hello' program.
hello_CSHARPLIBS = @GNU_GETTEXT_LDADD@ @GNU_GETTEXT_LIBS@ -l System -l System.Drawing -l System.Windows.Forms

# Resources that are generated from PO files.
MAINTAINERCLEANFILES = */*.resources.dll

# Additional files to be distributed.
EXTRA_DIST = autogen.sh autoclean.sh


# ----------------- General rules for compiling C# programs -----------------

pkgdatadir = $(datadir)/$(PACKAGE)
pkglibdir = $(libdir)/$(PACKAGE)

CSHARPCOMP = $(SHELL) csharpcomp.sh
CSHARPCOMPFLAGS = -O -g

EXTRA_DIST += $(hello_SOURCES)
CLEANFILES = hello.net.exe.mdb
DISTCLEANFILES = csharpcomp.sh csharpexec.sh


# Rules for compiling C# programs.

all-local: hello.net.exe hello.sh

hello.net.exe: $(hello_SOURCES)
	$(CSHARPCOMP) $(CSHARPCOMPFLAGS) -o $@ $(hello_CSHARPLIBS) $(srcdir)/hello.cs

hello.sh:
	{ echo '#!/bin/sh'; \
	  echo "exec /bin/sh '$(pkgdatadir)/csharpexec.sh' @GNU_GETTEXT_LDADD@ '$(pkglibdir)/hello.net.exe' \"\$$@\""; \
	} > $@

install-exec-local: all-local
	$(MKDIR_P) $(DESTDIR)$(bindir)
	$(INSTALL_SCRIPT) hello.sh $(DESTDIR)$(bindir)/hello
	$(MKDIR_P) $(DESTDIR)$(pkglibdir)
	$(INSTALL_DATA) hello.net.exe $(DESTDIR)$(pkglibdir)/hello.net.exe

install-data-local: all-local
	$(MKDIR_P) $(DESTDIR)$(pkgdatadir)
	$(INSTALL_DATA) csharpexec.sh $(DESTDIR)$(pkgdatadir)/csharpexec.sh

installdirs-local:
	$(MKDIR_P) $(DESTDIR)$(bindir)
	$(MKDIR_P) $(DESTDIR)$(pkglibdir)
	$(MKDIR_P) $(DESTDIR)$(pkgdatadir)

uninstall-local:
	rm -f $(DESTDIR)$(bindir)/hello
	rm -f $(DESTDIR)$(pkglibdir)/hello.net.exe
	rm -f $(DESTDIR)$(pkgdatadir)/csharpexec.sh

CLEANFILES += hello.net.exe hello.sh
