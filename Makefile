
REVISION="$(shell svnversion ../trunk)"

all: build test

build: build-linux
test: test-linux

test-all: test-linux test-wine
build-all: build-linux build-wine

TEST_STANZA='import sys, os; sys.path.insert(0, os.path.join(os.getcwd(), "site-packages")); import scipy; sys.exit(scipy.test(verbose=2))'

build-linux:
	@echo "version = \"$(REVISION)\"" > scipy/__svn_version__.py
	@echo "--- Building..."
	python2.5 setup.py install --prefix=$$PWD/dist/linux \
		> build.log 2>&1 || { cat build.log; exit 1; }

test-linux:
	@echo "--- Testing in Linux"
	(cd dist/linux/lib/python2.5 && python -c $(TEST_STANZA)) \
		> test.log 2>&1 || { cat test.log; exit 1; }

build-wine:
	@echo "--- Building..."
	wine c:\\Python25\\python.exe setup.py build --compiler=mingw32 install --prefix="dist\\win32" \
		> build.log 2>&1 || { cat build.log; exit 1; }

test-wine:
	@echo "--- Testing in WINE"
	(cd dist/win32/Lib && wine c:\\Python25\\python.exe -c $(TEST_STANZA)) \
		> test.log 2>&1 || { cat test.log; exit 1; }

.PHONY: test build test-linux build-linux test-wine build-wine
