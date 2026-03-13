# Makefile for ARNAUD C library
# Requires BLAS/LAPACK libraries

# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -Wpedantic -O3 -fPIC
CPPFLAGS = -Iinclude -Isrc

# BLAS/LAPACK library configuration
# Override these variables to use different BLAS/LAPACK implementations
# Examples:
#   make BLAS_LIBS="-lopenblas"                    # OpenBLAS
#   make BLAS_LIBS="-lmkl_rt"                      # Intel MKL (single dynamic library)
#   make BLAS_LIBS="-framework Accelerate"         # Apple Accelerate (macOS)
#   make BLAS_LIBS="-lblas -llapack"              # Reference BLAS/LAPACK

# Default: try to detect available BLAS/LAPACK
BLAS_LIBS ?= $(shell pkg-config --libs openblas 2>/dev/null || pkg-config --libs blas lapack 2>/dev/null || echo "-lblas -llapack")

# Additional libraries
LIBS = -lm $(BLAS_LIBS)

# Debug flags (uncomment for debug build)
# CFLAGS += -g -DDEBUG -O0

# Library name (shared library only - static not possible with BLAS/LAPACK dependencies)
LIBNAME = libarnaud
SHARED_LIB = $(LIBNAME).so

# Directories
SRCDIR = src
INCDIR = include
OBJDIR = obj
LIBDIR = lib

# Source files
SOURCES = $(wildcard $(SRCDIR)/*.c)

# Object files
OBJECTS = $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)

# Header dependencies
HEADERS = $(wildcard $(INCDIR)/arnaud/*.h) $(wildcard $(SRCDIR)/*.h)

# Default target - only shared library possible due to BLAS/LAPACK dependencies
all: directories $(SHARED_LIB)

# Create necessary directories
directories:
	@mkdir -p $(OBJDIR) $(LIBDIR)

# Shared library (only viable option with BLAS/LAPACK dependencies)
$(SHARED_LIB): $(OBJECTS)
	@echo "Creating shared library: $(LIBDIR)/$@"
	@echo "Using BLAS/LAPACK libraries: $(BLAS_LIBS)"
	@$(CC) -shared -o $(LIBDIR)/$@ $^ $(LIBS)
	@echo "Shared library created successfully"

# Compile source files to object files
$(OBJDIR)/%.o: $(SRCDIR)/%.c $(HEADERS)
	@echo "Compiling: $<"
	@$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

# Install target (copies library and headers to system directories)
install: all
	@echo "Installing headers to /usr/local/include/arnaud/"
	@sudo mkdir -p /usr/local/include/arnaud
	@sudo cp $(INCDIR)/arnaud/*.h /usr/local/include/arnaud/
	@echo "Installing shared library to /usr/local/lib/"
	@sudo cp $(LIBDIR)/$(SHARED_LIB) /usr/local/lib/
	@sudo ldconfig
	@echo "Installation complete"
	@echo "Note: Ensure BLAS/LAPACK libraries are available in your system"

# Uninstall target
uninstall:
	@echo "Removing installed files..."
	@sudo rm -rf /usr/local/include/arnaud
	@sudo rm -f /usr/local/lib/$(SHARED_LIB)
	@sudo ldconfig
	@echo "Uninstallation complete"

# Clean target
clean:
	@echo "Cleaning build artifacts..."
	@rm -rf $(OBJDIR) $(LIBDIR)
	@echo "Clean complete"

# Clean all (including backup files)
distclean: clean
	@echo "Deep cleaning..."
	@find . -name "*~" -delete
	@find . -name "*.bak" -delete
	@echo "Deep clean complete"

# Debug build
debug: CFLAGS += -g -DDEBUG -O0
debug: CFLAGS := $(filter-out -O3,$(CFLAGS))
debug: all

# Show help
help:
	@echo "Available targets:"
	@echo "  all        - Build shared library (default)"
	@echo "  shared     - Build shared library"
	@echo "  debug      - Build with debug symbols and no optimization"
	@echo "  install    - Install library and headers to system directories"
	@echo "  uninstall  - Remove installed files from system directories"
	@echo "  clean      - Remove build artifacts"
	@echo "  distclean  - Remove all generated files including backups"
	@echo "  help       - Show this help message"
	@echo "  check-blas - Check for available BLAS/LAPACK libraries"
	@echo ""
	@echo "Build configuration:"
	@echo "  Compiler: $(CC)"
	@echo "  CFLAGS: $(CFLAGS)"
	@echo "  CPPFLAGS: $(CPPFLAGS)"
	@echo "  BLAS/LAPACK: $(BLAS_LIBS)"
	@echo ""
	@echo "BLAS/LAPACK Library Examples:"
	@echo "  OpenBLAS:     make BLAS_LIBS='-lopenblas'"
	@echo "  Intel MKL:    make BLAS_LIBS='-lmkl_rt'"
	@echo "  Accelerate:   make BLAS_LIBS='-framework Accelerate'"
	@echo "  Reference:    make BLAS_LIBS='-lblas -llapack'"
	@echo ""
	@echo "Note: Static library compilation not supported due to BLAS/LAPACK dependencies"

# Convenience target (shared library only)
shared: directories $(SHARED_LIB)

# Check for BLAS/LAPACK libraries
check-blas:
	@echo "Checking for BLAS/LAPACK libraries..."
	@echo "Currently configured: $(BLAS_LIBS)"
	@echo ""
	@echo "Testing library availability:"
	@if pkg-config --exists openblas 2>/dev/null; then \
		echo "  OpenBLAS: Available (pkg-config)"; \
		echo "    Use: make BLAS_LIBS=\"$$(pkg-config --libs openblas)\""; \
	else \
		echo "  OpenBLAS: Not found via pkg-config"; \
	fi
	@if pkg-config --exists blas lapack 2>/dev/null; then \
		echo "  Reference BLAS/LAPACK: Available (pkg-config)"; \
		echo "    Use: make BLAS_LIBS=\"$$(pkg-config --libs blas lapack)\""; \
	else \
		echo "  Reference BLAS/LAPACK: Not found via pkg-config"; \
	fi
	@if [ "$$(uname)" = "Darwin" ]; then \
		echo "  Apple Accelerate: Available (macOS detected)"; \
		echo "    Use: make BLAS_LIBS=\"-framework Accelerate\""; \
	fi
	@echo ""
	@echo "Manual library check (common locations):"
	@if ldconfig -p 2>/dev/null | grep -q "libopenblas"; then \
		echo "  libopenblas: Found in system"; \
	else \
		echo "  libopenblas: Not found in system"; \
	fi
	@if ldconfig -p 2>/dev/null | grep -q "libmkl_rt"; then \
		echo "  Intel MKL: Found in system"; \
	else \
		echo "  Intel MKL: Not found in system"; \
	fi

# Check for required dependencies
check-deps:
	@echo "Checking for required dependencies..."
	@which $(CC) > /dev/null || (echo "Error: $(CC) not found" && exit 1)
	@echo "Compiler: $(CC) - OK"
	@echo "All dependencies satisfied"

# Show compiler version and capabilities
compiler-info:
	@echo "Compiler Information:"
	@$(CC) --version
	@echo ""
	@echo "Compiler supports C standards and will use default or system-configured standard"

# Phony targets
.PHONY: all directories install uninstall clean distclean debug help shared check-deps check-blas compiler-info

# Make variables visible
.SUFFIXES:
.SUFFIXES: .c .o

# Pattern rule for object files (explicit for clarity)
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@echo "Compiling: $< -> $@"
	@$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@
