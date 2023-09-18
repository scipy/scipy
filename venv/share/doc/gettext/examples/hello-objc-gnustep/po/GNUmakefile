# Example for use of GNU gettext.
# This file is in the public domain.
#
# Makefile - processed by GNU make.

# NLS is enabled by default.
USE_NLS = yes

# List of files which contain translatable strings.
POTFILES = \
  Hello.m

# Usually the message domain is the same as the package name in lower case.
DOMAIN = hello-objc-gnustep

# These options get passed to xgettext.
XGETTEXT_OPTIONS =

# This is the copyright holder that gets inserted into the header of the
# $(DOMAIN).pot file.  Set this to the copyright holder of the surrounding
# package.  (Note that the msgstr strings, extracted from the package's
# sources, belong to the copyright holder of the package.)  Translators are
# expected to transfer the copyright for their translations to this person
# or entity, or to disclaim their copyright.  The empty string stands for
# the public domain; in this case the translators are expected to disclaim
# their copyright.
COPYRIGHT_HOLDER = Yoyodyne, Inc.

# This tells whether or not to prepend "GNU " prefix to the package
# name that gets inserted into the header of the $(DOMAIN).pot file.
# Possible values are "yes", "no", or empty.  If it is empty, try to
# detect it automatically by scanning the files in $(top_srcdir) for
# "GNU packagename" string.
PACKAGE_GNU = no

# This is the email address or URL to which the translators shall report
# bugs in the untranslated strings:
# - Strings which are not entire sentences, see the maintainer guidelines
#   in the GNU gettext documentation, section 'Preparing Strings'.
# - Strings which use unclear terms or require additional context to be
#   understood.
# - Strings which make invalid assumptions about notation of date, time or
#   money.
# - Pluralisation problems.
# - Incorrect English spelling.
# - Incorrect formatting.
# It can be your email address, or a mailing list address where translators
# can write to without being subscribed, or the URL of a web page through
# which the translators can contact you.
MSGID_BUGS_ADDRESS = bug-gettext@gnu.org

# These options get passed to msgmerge.
# Useful options are in particular:
#   --previous            to keep previous msgids of translated messages,
#   --quiet               to reduce the verbosity.
MSGMERGE_OPTIONS = --quiet

# These options get passed to msginit.
# If you want to disable line wrapping when writing PO files, add
# --no-wrap to MSGMERGE_OPTIONS, XGETTEXT_OPTIONS, and
# MSGINIT_OPTIONS.
MSGINIT_OPTIONS =

-include $(GNUSTEP_MAKEFILES)/rules.make

# Main application
PACKAGE_NAME = Hello
PACKAGE_VERSION = 0

XGETTEXT = xgettext
MSGMERGE = msgmerge
MSGMERGE_UPDATE = msgmerge --update
MSGMERGE_FOR_MSGFMT_OPTION = --for-msgfmt
MSGINIT = msginit
MSGCONV = msgconv
MSGFILTER = msgfilter

# Set of available languages.
LINGUAS := $(shell sed -e "/^\#/d" LINGUAS)
# Set of languages which use msgen.
ENLINGUAS := $(filter-out $(LINGUAS), en)

top_srcdir = ..

include LocaleAliases

POTFILES_DEPS = $(foreach file, $(POTFILES), $(top_srcdir)/$(file))
POFILES = $(foreach lang, $(LINGUAS), $(lang).po)
UPDATEPOFILES = $(foreach lang, $(LINGUAS), $(lang).po-update)
DUMMYPOFILES = $(foreach lang, $(LINGUAS), $(lang).nop)
STRINGSFILES = $(foreach lang, $(LINGUAS), $(top_srcdir)/$(NSLOCALE_$(lang)).lproj/Localizable.strings)
ENSTRINGSFILES = $(foreach lang, $(ENLINGUAS), $(top_srcdir)/$(NSLOCALE_$(lang)).lproj/Localizable.strings)
CATALOGS = $(STRINGSFILES) $(ENSTRINGSFILES)

.SUFFIXES: .po .sed .sin .nop .po-create .po-update

.sin.sed:
	sed -e '/^#/d' $< > t-$@
	mv t-$@ $@


internal-all:: all-local-$(USE_NLS)

all-local-yes: $(CATALOGS)
all-local-no:

# Note: Target 'all' must not depend on target '$(DOMAIN).pot-update',
# otherwise packages like GCC can not be built if only parts of the source
# have been downloaded.

# This target rebuilds $(DOMAIN).pot; it is an expensive operation.
# Note that $(DOMAIN).pot is not touched if it doesn't need to be changed.
# The determination of whether the package xyz is a GNU one is based on the
# heuristic whether some file in the top level directory mentions "GNU xyz".
# If GNU 'find' is available, we avoid grepping through monster files.
$(DOMAIN).pot-update: $(POTFILES_DEPS) remove-potcdate.sed
	package_gnu="$(PACKAGE_GNU)"; \
	test -n "$$package_gnu" || { \
	  if { if (LC_ALL=C find --version) 2>/dev/null | grep GNU >/dev/null; then \
	         LC_ALL=C find -L $(top_srcdir) -maxdepth 1 -type f -size -10000000c -exec grep -i 'GNU $(PACKAGE_NAME)' /dev/null '{}' ';' 2>/dev/null; \
	       else \
	         LC_ALL=C grep -i 'GNU $(PACKAGE_NAME)' $(top_srcdir)/* 2>/dev/null; \
	       fi; \
	     } | grep -v 'libtool:' >/dev/null; then \
	     package_gnu=yes; \
	   else \
	     package_gnu=no; \
	   fi; \
	}; \
	if test "$$package_gnu" = "yes"; then \
	  package_prefix='GNU '; \
	else \
	  package_prefix=''; \
	fi; \
	msgid_bugs_address='$(MSGID_BUGS_ADDRESS)'; \
	case `$(XGETTEXT) --version | sed 1q | sed -e 's,^[^0-9]*,,'` in \
	  '' | 0.[0-9] | 0.[0-9].* | 0.1[0-5] | 0.1[0-5].* | 0.16 | 0.16.[0-1]*) \
	    $(XGETTEXT) --default-domain=$(DOMAIN) --directory=$(top_srcdir) \
	      --add-comments=TRANSLATORS: \
	      --copyright-holder='$(COPYRIGHT_HOLDER)' \
	      --msgid-bugs-address="$$msgid_bugs_address" \
	      $(XGETTEXT_OPTIONS) \
	      $(POTFILES) \
	    ;; \
	  *) \
	    $(XGETTEXT) --default-domain=$(DOMAIN) --directory=$(top_srcdir) \
	      --add-comments=TRANSLATORS: \
	      --copyright-holder='$(COPYRIGHT_HOLDER)' \
	      --package-name="$${package_prefix}$(PACKAGE_NAME)" \
	      --package-version='$(VERSION)' \
	      --msgid-bugs-address="$$msgid_bugs_address" \
	      $(XGETTEXT_OPTIONS) \
	      $(POTFILES) \
	    ;; \
	esac
	test ! -f $(DOMAIN).po || { \
	  if test -f $(DOMAIN).pot-header; then \
	    sed -e '1,/^#$$/d' < $(DOMAIN).po > $(DOMAIN).1po && \
	    cat $(DOMAIN).pot-header $(DOMAIN).1po > $(DOMAIN).po && \
	    rm -f $(DOMAIN).1po \
	    || exit 1; \
	  fi; \
	  if test -f $(DOMAIN).pot; then \
	    sed -f remove-potcdate.sed < $(DOMAIN).pot > $(DOMAIN).1po && \
	    sed -f remove-potcdate.sed < $(DOMAIN).po > $(DOMAIN).2po && \
	    if cmp $(DOMAIN).1po $(DOMAIN).2po >/dev/null 2>&1; then \
	      rm -f $(DOMAIN).1po $(DOMAIN).2po $(DOMAIN).po; \
	    else \
	      rm -f $(DOMAIN).1po $(DOMAIN).2po $(DOMAIN).pot && \
	      mv $(DOMAIN).po $(DOMAIN).pot; \
	    fi; \
	  else \
	    mv $(DOMAIN).po $(DOMAIN).pot; \
	  fi; \
	}

# This rule has no dependencies: we don't need to update $(DOMAIN).pot at
# every "make" invocation, only create it when it is missing.
# Only "make $(DOMAIN).pot-update" or "make dist" will force an update.
$(DOMAIN).pot:
	$(MAKE) $(DOMAIN).pot-update

# This target rebuilds a PO file if $(DOMAIN).pot has changed.
# Note that a PO file is not touched if it doesn't need to be changed.
$(POFILES): $(DOMAIN).pot
	@lang=`echo $@ | sed -e 's,.*/,,' -e 's/\.po$$//'`; \
	if test -f $${lang}.po; then \
	  echo "$(MSGMERGE_UPDATE) $(MSGMERGE_OPTIONS) --lang=$${lang} --previous $${lang}.po $(DOMAIN).pot"; \
	  case `$(MSGMERGE_UPDATE) --version | sed 1q | sed -e 's,^[^0-9]*,,'` in \
	    '' | 0.[0-9] | 0.[0-9].* | 0.1[0-5] | 0.1[0-5].*) \
	      $(MSGMERGE_UPDATE) $(MSGMERGE_OPTIONS) $${lang}.po $(DOMAIN).pot;; \
	    0.1[6-7] | 0.1[6-7].*) \
	      $(MSGMERGE_UPDATE) $(MSGMERGE_OPTIONS) --previous $${lang}.po $(DOMAIN).pot;; \
	    *) \
	      $(MSGMERGE_UPDATE) $(MSGMERGE_OPTIONS) --lang=$${lang} --previous $${lang}.po $(DOMAIN).pot;; \
	  esac; \
	else \
	  $(MAKE) $${lang}.po-create; \
	fi

internal-clean::
	rm -f remove-potcdate.sed

update-po:
	$(MAKE) $(DOMAIN).pot-update
	test -z "$(UPDATEPOFILES)" || $(MAKE) $(UPDATEPOFILES)
	$(MAKE) update-strings

# General rule for creating PO files.

.nop.po-create:
	@lang=`echo $@ | sed -e 's/\.po-create$$//'`; \
	echo "File $$lang.po does not exist. If you are a translator, you can create it through 'msginit'." 1>&2; \
	exit 1

# General rule for updating PO files.

.nop.po-update:
	@lang=`echo $@ | sed -e 's/\.po-update$$//'`; \
	tmpdir=`pwd`; \
	echo "$$lang:"; \
	echo "$(MSGMERGE) $(MSGMERGE_OPTIONS) --lang=$$lang $$lang.po $(DOMAIN).pot -o $$lang.new.po"; \
	if { case `$(MSGMERGE) --version | sed 1q | sed -e 's,^[^0-9]*,,'` in \
	       '' | 0.[0-9] | 0.[0-9].* | 0.1[0-7] | 0.1[0-7].*) \
	         $(MSGMERGE) $(MSGMERGE_OPTIONS) -o $$tmpdir/$$lang.new.po $$lang.po $(DOMAIN).pot;; \
	       *) \
	         $(MSGMERGE) $(MSGMERGE_OPTIONS) --lang=$$lang -o $$tmpdir/$$lang.new.po $$lang.po $(DOMAIN).pot;; \
	     esac; \
	   }; then \
	  if cmp $$lang.po $$tmpdir/$$lang.new.po >/dev/null 2>&1; then \
	    rm -f $$tmpdir/$$lang.new.po; \
	  else \
	    if mv -f $$tmpdir/$$lang.new.po $$lang.po; then \
	      :; \
	    else \
	      echo "msgmerge for $$lang.po failed: cannot move $$tmpdir/$$lang.new.po to $$lang.po" 1>&2; \
	      exit 1; \
	    fi; \
	  fi; \
	else \
	  echo "msgmerge for $$lang.po failed!" 1>&2; \
	  rm -f $$tmpdir/$$lang.new.po; \
	fi

$(DUMMYPOFILES):


update-strings: $(CATALOGS)
	@:

# During .po -> .strings conversion, take into account the most recent changes
# to the .pot file. This eliminates the need to update the .po files when the
# .pot file has changed, which would be troublesome if the .po files are put
# under version control.
$(STRINGSFILES): $(DOMAIN).pot $(POFILES)
	mkdir -p $(patsubst %/,%, $(dir $@))
	$(MSGMERGE) $(MSGMERGE_FOR_MSGFMT_OPTION) -o $@.1po $(UNIXLOCALE_$(patsubst %.lproj,%, $(notdir $(patsubst %/,%, $(dir $@))))).po $(DOMAIN).pot
	msgcat --stringtable-output -o $@.tmp $@.1po
	mv $@.tmp $@
	rm -f $@.1po

$(ENSTRINGSFILES): $(DOMAIN).pot
	mkdir -p $(patsubst %/,%, $(dir $@))
	$(MSGINIT) -i $(DOMAIN).pot --no-translator -l $(UNIXLOCALE_$(patsubst %.lproj,%, $(notdir $(patsubst %/,%, $(dir $@))))) --stringtable-output -o $@.tmp
	mv $@.tmp $@

echo-languages:
	@echo $(foreach lang, $(LINGUAS), $(NSLOCALE_$(lang)))

