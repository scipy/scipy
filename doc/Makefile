# Makefile for Sphinx documentation
#

# PYVER needs to be major.minor, just "3" doesn't work - it will result in
# issues with the amendments to PYTHONPATH and install paths (see DIST_VARS).

# Use explicit "version_info" indexing since make cannot handle colon characters, and
# evaluate it now to allow easier debugging when printing the varaible

PYVER:=$(shell python3 -c 'from sys import version_info as v; print("{0}.{1}".format(v[0], v[1]))')
PYTHON = python$(PYVER)

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = LANG=C $(PYTHON) $(shell which sphinx-build)
PAPER         =

FILES=

# Internal variables.
PAPEROPT_a4     = -D "latex_elements.papersize=a4paper"
PAPEROPT_letter = -D "latex_elements.papersize=letterpaper"
ALLSPHINXOPTS   = -d build/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) source

.PHONY: help clean html web pickle htmlhelp latex changes linkcheck \
        dist dist-build

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
	@echo "  dist PYVER=... to make a distribution-ready tree"
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

INSTALL_DIR = $(CURDIR)/build/inst-dist/
INSTALL_PPH = $(INSTALL_DIR)/lib/python$(PYVER)/site-packages:$(INSTALL_DIR)/local/lib/python$(PYVER)/site-packages:$(INSTALL_DIR)/lib/python$(PYVER)/dist-packages:$(INSTALL_DIR)/local/lib/python$(PYVER)/dist-packages
UPLOAD_DIR=/srv/docs_scipy_org/doc/scipy-$(RELEASE)

DIST_VARS=PYTHON="PYTHONPATH=$(INSTALL_PPH):$$PYTHONPATH $(PYTHON)" SPHINXBUILD="PYTHONPATH=$(INSTALL_PPH):$$PYTHONPATH $(SPHINXBUILD)"

dist:
	make $(DIST_VARS) real-dist

real-dist: dist-build html html-scipyorg
	test -d build/latex || make latex
	make -C build/latex all-pdf LATEXOPTS="-halt-on-error"
	-test -d build/htmlhelp || make htmlhelp-build
	-rm -rf build/dist
	mkdir -p build/dist
	cp -r build/html-scipyorg build/dist/reference
	touch build/dist/index.html
	(cd build/html && zip -9qr ../dist/scipy-html.zip .)
	cp build/latex/scipy-ref.pdf build/dist
	chmod ug=rwX,o=rX -R build/dist
	find build/dist -type d -print0 | xargs -0r chmod g+s
	cd build/dist && tar czf ../dist.tar.gz *

dist-build:
	rm -f ../dist/*.egg
	cd .. && $(PYTHON) setup.py bdist_egg
	install -d $(subst :, ,$(INSTALL_PPH))
	$(PYTHON) `which easy_install` --prefix=$(INSTALL_DIR) ../dist/*.egg

upload:
	# SSH must be correctly configured for this to work.
	# Assumes that ``make dist`` was already run
	# Example usage: ``make upload USERNAME=rgommers RELEASE=0.17.0``
	ssh $(USERNAME)@new.scipy.org mkdir $(UPLOAD_DIR)
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

html:
	mkdir -p build/html build/doctrees
	$(SPHINXBUILD) -b html $(ALLSPHINXOPTS) build/html $(FILES)
	@echo
	@echo "Build finished. The HTML pages are in build/html."

html-scipyorg:
	mkdir -p build/html build/doctrees
	$(SPHINXBUILD) -WT --keep-going -t scipyorg -b html $(ALLSPHINXOPTS) build/html-scipyorg $(FILES)
	@echo
	@echo "Build finished. The HTML pages are in build/html-scipyorg."

pickle:
	mkdir -p build/pickle build/doctrees
	$(SPHINXBUILD) -b pickle $(ALLSPHINXOPTS) build/pickle $(FILES)
	@echo
	@echo "Build finished; now you can process the pickle files or run"
	@echo "  sphinx-web build/pickle"
	@echo "to start the sphinx-web server."

web: pickle

htmlhelp:
	mkdir -p build/htmlhelp build/doctrees
	$(SPHINXBUILD) -b htmlhelp $(ALLSPHINXOPTS) build/htmlhelp $(FILES)
	@echo
	@echo "Build finished; now you can run HTML Help Workshop with the" \
	      ".hhp project file in build/htmlhelp."

htmlhelp-build: htmlhelp build/htmlhelp/scipy.chm
%.chm: %.hhp
	-hhc.exe $^

latex:
	mkdir -p build/latex build/doctrees
	$(SPHINXBUILD) -b latex $(ALLSPHINXOPTS) build/latex $(FILES)
	$(PYTHON) postprocess.py tex build/latex/*.tex
	@echo
	@echo "Build finished; the LaTeX files are in build/latex."
	@echo "Run \`make all-pdf' or \`make all-ps' in that directory to" \
	      "run these through (pdf)latex."

coverage: build
	mkdir -p build/coverage build/doctrees
	$(SPHINXBUILD) -b coverage $(ALLSPHINXOPTS) build/coverage $(FILES)
	@echo "Coverage finished; see c.txt and python.txt in build/coverage"

changes:
	mkdir -p build/changes build/doctrees
	$(SPHINXBUILD) -b changes $(ALLSPHINXOPTS) build/changes $(FILES)
	@echo
	@echo "The overview file is in build/changes."

linkcheck:
	mkdir -p build/linkcheck build/doctrees
	$(SPHINXBUILD) -b linkcheck $(ALLSPHINXOPTS) build/linkcheck $(FILES)
	@echo
	@echo "Link check complete; look for any errors in the above output " \
	      "or in build/linkcheck/output.txt."

show:
	@python -c "import webbrowser; webbrowser.open_new_tab('file://$(PWD)/build/html-scipyorg/index.html')"

