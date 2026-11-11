#ifndef SCIPY_SUPERLU_CONFIG_H
#define SCIPY_SUPERLU_CONFIG_H

/*
 * SciPy's replacement for SuperLU's generated `superlu_config.h`; meson copies
 * this file to `superlu_config.h` in the build directory.
 *
 * Every vendored SuperLU source that needs SciPy's configuration includes
 * "superlu_config.h" before using any of it, either directly or via
 * `slu_[sdcz]defs.h`. Shadowing that header lets us inject our configuration
 * without patching the vendored sources.
 *
 * SciPy's headers must come first: upstream's template only typedefs
 * `slu_blasint` if it is not already defined as a macro.
 */
#include "scipy_slu_blas_config.h"
#include "superlu_config_upstream.h"

#endif /* SCIPY_SUPERLU_CONFIG_H */
