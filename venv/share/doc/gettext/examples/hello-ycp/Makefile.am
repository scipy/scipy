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

# The YCP program is installed as data and invoked through a shell script.
pkgdata_DATA = hello.ycp
hello:
	{ echo '#!/bin/sh'; \
	  echo 'if test -n "$$DISPLAY"; then'; \
	  echo "  exec /usr/lib/YaST2/bin/y2base '$(pkgdatadir)/hello.ycp' qt"; \
	  echo 'else'; \
	  echo "  exec /usr/lib/YaST2/bin/y2base '$(pkgdatadir)/hello.ycp' ncurses"; \
	  echo 'fi'; \
	} > $@

# Additional files to be distributed.
EXTRA_DIST = autogen.sh autoclean.sh
