# Example for use of GNU gettext.
# This file is in the public domain.

include $(GNUSTEP_MAKEFILES)/common.make

# Subprojects
SUBPROJECTS = po

# Main application
PACKAGE_NAME = Hello
APP_NAME = Hello
GNUSTEP_INSTALLATION_DIR = $(GNUSTEP_LOCAL_ROOT)
Hello_APPLICATION_ICON =

# Additional libraries
ADDITIONAL_GUI_LIBS +=

# Resource files
Hello_RESOURCE_FILES =
Hello_LOCALIZED_RESOURCE_FILES = Localizable.strings
Hello_LANGUAGES := $(sort English $(shell MAKEFLAGS= $(MAKE) -s -C po echo-languages))

# Header files
Hello_HEADERS = AppController.h Hello.h

# Class files
Hello_OBJC_FILES = main.m AppController.m Hello.m

# C files
Hello_C_FILES =

-include GNUmakefile.preamble
-include GNUmakefile.local
include $(GNUSTEP_MAKEFILES)/aggregate.make
include $(GNUSTEP_MAKEFILES)/application.make
-include GNUmakefile.postamble

internal-distclean::
	rm -rf *.lproj
