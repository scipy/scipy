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
bin_JAVAPROGRAMS = hello

# The source files of the 'hello' program.
hello_SOURCES = Hello.java
hello_CLASSES = Hello.class

# The entry point of the 'hello' program.
hello_MAINCLASS = Hello

# The link dependencies of the 'hello' program.
hello_JAVALIBS = @LIBINTL_JAR@

# Additional files to be distributed.
EXTRA_DIST = autogen.sh autoclean.sh


# ----------------- General rules for compiling Java programs -----------------

jardir = $(datadir)/$(PACKAGE)
pkgdatadir = $(datadir)/$(PACKAGE)
pkglibdir = $(libdir)/$(PACKAGE)

JAR = @JAR@
JAVACOMP = $(SHELL) javacomp.sh

EXTRA_DIST += $(hello_SOURCES)
CLEANFILES =
DISTCLEANFILES = javacomp.sh javaexec.sh


# Rules for compiling Java programs as jar libraries.
# This is the preferred mode during development, because you can easily test
# the program without installing it, simply by doing "java -jar hello.jar".

all-local: hello.jar hello.sh

hello.jar: $(hello_CLASSES)
	{ echo "Manifest-Version: 1.0"; echo "Main-Class: $(hello_MAINCLASS)"; echo 'Class-Path: @LIBINTL_JAR@'; } > Manifest.mf
	$(JAR) cfm $@ Manifest.mf Hello*.class
	rm -f Manifest.mf
	abs_jar=`pwd`/$@; (cd po && $(MAKE)) && catalogs=`GNUMAKEFLAGS=--no-print-directory $(MAKE) -s -C po echo-catalogs`; test -n "$$catalogs" && (cd $(srcdir)/po && $(JAR) uf "$$abs_jar" $$catalogs) || { rm -f $@ jartmp*; exit 1; }

Hello.class: $(srcdir)/Hello.java
	CLASSPATH=.@CLASSPATH_SEPARATOR@$(hello_JAVALIBS) $(JAVACOMP) -d . $(srcdir)/Hello.java

hello.sh:
	{ echo '#!/bin/sh'; \
	  echo "CLASSPATH='$(jardir)/hello.jar@CLASSPATH_SEPARATOR@$(hello_JAVALIBS)'\$${CLASSPATH+\"@CLASSPATH_SEPARATOR@\$$CLASSPATH\"}"; \
	  echo "export CLASSPATH"; \
	  echo "exec /bin/sh '$(pkgdatadir)/javaexec.sh' $(hello_MAINCLASS) \"\$$@\""; \
	} > $@

install-exec-local: all-local
	$(MKDIR_P) $(DESTDIR)$(bindir)
	$(INSTALL_SCRIPT) hello.sh $(DESTDIR)$(bindir)/hello

install-data-local: all-local
	$(MKDIR_P) $(DESTDIR)$(jardir)
	$(INSTALL_DATA) hello.jar $(DESTDIR)$(jardir)/hello.jar
	$(MKDIR_P) $(DESTDIR)$(pkgdatadir)
	$(INSTALL_DATA) javaexec.sh $(DESTDIR)$(pkgdatadir)/javaexec.sh

installdirs-local:
	$(MKDIR_P) $(DESTDIR)$(jardir)
	$(MKDIR_P) $(DESTDIR)$(pkgdatadir)

uninstall-local:
	rm -f $(DESTDIR)$(bindir)/hello
	rm -f $(DESTDIR)$(jardir)/hello.jar
	rm -f $(DESTDIR)$(pkgdatadir)/javaexec.sh

CLEANFILES += hello.jar Hello*.class Manifest.mf hello.sh
