# Example for use of GNU gettext.
# This file is in the public domain.
#
# Makefile configuration - processed by automake.

# General automake options.
AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4

# The list of subdirectories containing Makefiles.
SUBDIRS = m4 . po

# The list of programs that are built.
bin_PASCALPROGRAMS = hello

# The source files of the 'hello' program.
hello_SOURCES = $(srcdir)/hello.pas

# Additional files to be distributed.
EXTRA_DIST = autogen.sh autoclean.sh

# ---------------- General rules for compiling Pascal programs ----------------

EXTRA_DIST += $(hello_SOURCES)

# Distribute the RSJ file because it's needed to generate POT files and can
# only be rebuilt on those platforms to which the Pascal compiler is ported.
EXTRA_DIST += hello.rsj
# The GNU Coding Standards say in
# <https://www.gnu.org/prep/standards/html_node/Makefile-Basics.html>:
#   "GNU distributions usually contain some files which are not source files
#    ... . Since these files normally appear in the source directory, they
#    should always appear in the source directory, not in the build directory.
#    So Makefile rules to update them should put the updated files in the
#    source directory."
# Therefore we put this file in the source directory, not the build directory.

# Rules for compiling Pascal programs.

all-local: hello$(EXEEXT)

# How to build the 'hello' program.
hello$(EXEEXT) $(srcdir)/hello.rsj: $(hello_SOURCES)
	LOCALEDIR='@localedir@' $(PPC) -o./hello$(EXEEXT) $(hello_SOURCES)
# Move hello.rsj into $(srcdir). But don't provoke a gratuitous error in a
# VPATH build with read-only $(srcdir).
	if test '$(srcdir)' != .; then \
	  if test -f $(srcdir)/hello.rsj && cmp hello.rsj $(srcdir)/hello.rsj >/dev/null; then \
	    rm -f hello.rsj; \
	  else \
	    mv -f hello.rsj $(srcdir)/hello.rsj; \
	  fi; \
	fi

install-exec-local: all-local
	$(MKDIR_P) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) hello$(EXEEXT) $(DESTDIR)$(bindir)/hello$(EXEEXT)

installdirs-local:
	$(MKDIR_P) $(DESTDIR)$(bindir)

uninstall-local:
	rm -f $(DESTDIR)$(bindir)/hello$(EXEEXT)

# The list of auxiliary files generated during the compilation.
CLEANFILES = hello.o hello$(EXEEXT)

MAINTAINERCLEANFILES = hello.rsj
