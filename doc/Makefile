# Makefile for Sphinx documentation

# If `python3` doesn't point to the Python you want to use (with Sphinx installed),
# then you can use, e.g., `make html PYTHON=python3.9` to specify a Python executable
PYTHON = python3

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = LANG=C $(PYTHON) -msphinx

# Environment variables needed for notebooks (see gh-17322):
export SQLALCHEMY_SILENCE_UBER_WARNING?=1
export JUPYTER_PLATFORM_DIRS?=1

FILES=

# Internal variables.
ALLSPHINXOPTS   = -d build/doctrees $(SPHINXOPTS) source

.PHONY: help clean html changes linkcheck dist version-check

#------------------------------------------------------------------------------

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"
	@echo "  pickle    to make pickle files (usable by e.g. sphinx-web)"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  changes   to make an overview over all changed/added/deprecated items"
	@echo "  linkcheck to check all external links for integrity"
	@echo "  clean     to clean all build files/directories"
	@echo "  dist      to make a distribution-ready tree (installing Scipy in venv)"
	@echo "  doc-dist  to make a distribution-ready tree (assuming Scipy is installed)"
	@echo "  upload USERNAME=... RELEASE=... to upload built docs to docs.scipy.org"
	@echo "  show      to show the HTML output"

clean:
	-rm -rf build/* source/reference/generated source/.jupyterlite.doit.db


#------------------------------------------------------------------------------
# Automated generation of all documents
#------------------------------------------------------------------------------

# Build the current scipy version, and extract docs from it.
# We have to be careful of some issues:
#
# - Everything must be done using the same Python version
# - We must use eggs (otherwise they might override PYTHONPATH on import).
# - Different versions of easy_install install to different directories (!)
#

UPLOAD_DIR=/srv/docs_scipy_org/doc/scipy-$(RELEASE)

# `set -o pipefail` is specific to bash
SHELL = /bin/bash
SCIPYVER:=$(shell $(PYTHON) -c "import scipy; print(scipy.version.git_revision[:10])" 2>/dev/null)
GITVER ?= $(shell (cd ..; set -o pipefail && git rev-parse HEAD 2>/dev/null | cut -c1-10) || echo Unknown)

version-check:
ifeq "$(GITVER)" "Unknown"
	# @echo sdist build with unlabeled sources
else ifeq ("", "$(SCIPYVER)")
	@echo scipy not found, cannot build documentation without successful \"import scipy\"
	@exit 1
else ifneq ($(SCIPYVER),$(GITVER))
	@echo installed scipy $(SCIPYVER) != current repo git version \'$(GITVER)\'
	@echo use '"make dist"' or '"GITVER=$(SCIPYVER) make $(MAKECMDGOALS) ..."'
	@exit 1
else
	# for testing
	# @echo installed scipy $(SCIPYVER) matches git version $(GITVER); exit 1
endif

dist:
	install -d build
	rm -rf build/env
	$(PYTHON) -mvenv --system-site-packages build/env
	$(CURDIR)/build/env/bin/python -mpip install -r ../requirements/doc.txt
	$(CURDIR)/build/env/bin/python -mpip install ..
	make PYTHON="$(CURDIR)/build/env/bin/python" doc-dist

doc-dist: VERSIONWARNING=-t versionwarning

doc-dist: html
	-test -d build/htmlhelp || make htmlhelp-build
	-rm -rf build/dist
	mkdir -p build/dist
	cp -r build/html/* build/dist
	(cd build/html && zip -9qr ../dist/scipy-html.zip .)
	chmod ug=rwX,o=rX -R build/dist
	find build/dist -type d -print0 | xargs -0r chmod g+s
	cd build/dist && tar czf ../dist.tar.gz *

upload:
	# SSH must be correctly configured for this to work.
	# Assumes that ``make dist`` was already run
	# Example usage: ``make upload USERNAME=rgommers RELEASE=0.17.0``
	ssh $(USERNAME)@new.scipy.org mkdir -p $(UPLOAD_DIR)
	scp build/dist.tar.gz $(USERNAME)@new.scipy.org:$(UPLOAD_DIR)
	ssh $(USERNAME)@new.scipy.org tar xvC $(UPLOAD_DIR) \
	    -zf $(UPLOAD_DIR)/dist.tar.gz
	ssh $(USERNAME)@new.scipy.org mv $(UPLOAD_DIR)/scipy-html.zip \
	    $(UPLOAD_DIR)/scipy-html-$(RELEASE).zip
	ssh $(USERNAME)@new.scipy.org rm $(UPLOAD_DIR)/dist.tar.gz
	ssh $(USERNAME)@new.scipy.org ln -snf scipy-$(RELEASE) /srv/docs_scipy_org/doc/scipy

#------------------------------------------------------------------------------
# Basic Sphinx generation rules for different formats
#------------------------------------------------------------------------------

html: version-check convert-notebooks html-build
html-build:
	mkdir -p build/html build/doctrees
	$(SPHINXBUILD) -WT --keep-going $(VERSIONWARNING) -b html $(ALLSPHINXOPTS) build/html $(FILES)

convert-notebooks:
	mkdir -p source/_contents
	for file in source/tutorial/stats/*.md; do \
		basename=$$(basename $$file .md); \
		output_name="source/_contents/$${basename}.ipynb"; \
		$(PYTHON) -m jupytext --output "$$output_name"  $$file; \
	done


coverage: build version-check
	mkdir -p build/coverage build/doctrees
	$(SPHINXBUILD) -b coverage $(ALLSPHINXOPTS) build/coverage $(FILES)
	@echo "Coverage finished; see c.txt and python.txt in build/coverage"

changes: version-check
	mkdir -p build/changes build/doctrees
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) build/changes $(FILES)
	@echo
	@echo "The overview file is in build/changes."

linkcheck: version-check
	mkdir -p build/linkcheck build/doctrees
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) build/linkcheck $(FILES)
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in build/linkcheck/output.txt."

show:
	@python -c "import webbrowser; webbrowser.open_new_tab('file://$(PWD)/build/html/index.html')"
