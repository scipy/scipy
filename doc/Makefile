# Makefile for Sphinx documentation

# If `python3` doesn't point to the Python you want to use (with Sphinx installed),
# then you can use, e.g., `make html PYTHON=python3.7` to specify a Python executable
PYTHON = python3

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = LANG=C $(PYTHON) -msphinx
PAPER         =

FILES=

# Internal variables.
PAPEROPT_a4     = -D "latex_elements.papersize=a4paper"
PAPEROPT_letter = -D "latex_elements.papersize=letterpaper"
ALLSPHINXOPTS   = -d build/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) source

.PHONY: help clean html web pickle htmlhelp latex changes linkcheck \
        dist dist-build version-check

#------------------------------------------------------------------------------

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html      to make standalone HTML files"
	@echo "  html-scipyorg  to make standalone HTML files with scipy.org theming"
	@echo "  pickle    to make pickle files (usable by e.g. sphinx-web)"
	@echo "  htmlhelp  to make HTML files and a HTML help project"
	@echo "  latex     to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  changes   to make an overview over all changed/added/deprecated items"
	@echo "  linkcheck to check all external links for integrity"
	@echo "  clean     to clean all build files/directories"
	@echo "  dist      to make a distribution-ready tree (installing Scipy in venv)"
	@echo "  doc-dist  to make a distribution-ready tree (assuming Scipy is installed)"
	@echo "  upload USERNAME=... RELEASE=... to upload built docs to docs.scipy.org"
	@echo "  show      to show the html-scipyorg output"

clean:
	-rm -rf build/* source/generated


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

SCIPYVER:=$(shell $(PYTHON) -c "import scipy; print(scipy.version.git_revision[:10])" 2>/dev/null)
GITVER ?= $(shell cd ..; $(PYTHON) -c "from setup import git_version; \
				print(git_version()[0][:10])")

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
	$(CURDIR)/build/env/bin/python -mpip install Sphinx
	$(CURDIR)/build/env/bin/python -mpip install ..
	make PYTHON="$(CURDIR)/build/env/bin/python" doc-dist

doc-dist: VERSIONWARNING=-t versionwarning

doc-dist: html-scipyorg html
	test -d build/latex || make latex
	make -C build/latex all-pdf LATEXOPTS="-halt-on-error"
	-test -d build/htmlhelp || make htmlhelp-build
	-rm -rf build/dist
	mkdir -p build/dist
	cp -r build/html-scipyorg build/dist
	touch build/dist/index.html
	(cd build/html && zip -9qr ../dist/scipy-html.zip .)
	cp build/latex/scipy-ref.pdf build/dist
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
	ssh $(USERNAME)@new.scipy.org mv $(UPLOAD_DIR)/scipy-ref.pdf \
	    $(UPLOAD_DIR)/scipy-ref-$(RELEASE).pdf
	ssh $(USERNAME)@new.scipy.org mv $(UPLOAD_DIR)/scipy-html.zip \
	    $(UPLOAD_DIR)/scipy-html-$(RELEASE).zip
	ssh $(USERNAME)@new.scipy.org rm $(UPLOAD_DIR)/dist.tar.gz
	ssh $(USERNAME)@new.scipy.org ln -snf scipy-$(RELEASE) /srv/docs_scipy_org/doc/scipy

#------------------------------------------------------------------------------
# Basic Sphinx generation rules for different formats
#------------------------------------------------------------------------------

html: version-check html-build
html-build:
	mkdir -p build/html build/doctrees
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) build/html $(FILES)
	@echo
	@echo "Build finished. The HTML pages are in build/html."

html-scipyorg:
	mkdir -p build/html build/doctrees
	$(SPHINXBUILD) -WT --keep-going -t scipyorg $(VERSIONWARNING) -b html $(ALLSPHINXOPTS) build/html-scipyorg $(FILES)
	@echo
	@echo "Build finished. The HTML pages are in build/html-scipyorg."

latex: version-check latex-build
latex-build:
	mkdir -p build/latex build/doctrees
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) build/latex $(FILES)
	$(PYTHON) postprocess.py tex build/latex/*.tex
	@echo
	@echo "Build finished; the LaTeX files are in build/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

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
	@python -c "import webbrowser; webbrowser.open_new_tab('file://$(PWD)/build/html-scipyorg/index.html')"

