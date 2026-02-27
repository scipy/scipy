#### Reference issue
Closes #22635

#### What does this implement/fix?
This PR adds documentation for default values of optional parameters in `scipy.optimize.root` and its method-specific implementations.

The following changes were made:

**Main `root` function** (`scipy/optimize/_root.py`):
- Documented default values for all optional parameters: `args`, `method`, `jac`, `tol`, `callback`, `options`

**`hybr` method** (`scipy/optimize/_minpack_py.py`):
- Added defaults for: `col_deriv`, `xtol`, `maxfev`, `band`, `eps`, `factor`, `diag`

**`lm` method** (`scipy/optimize/_root.py`):
- Added defaults for: `col_deriv`, `ftol`, `xtol`, `gtol`, `maxiter`, `eps`, `factor`, `diag`

**`df-sane` method** (`scipy/optimize/_spectral.py`):
- Added/standardized defaults for: `ftol`, `fatol`, `fnorm`, `maxfev`, `disp`, `eta_strategy`, `sigma_eps`, `sigma_0`, `M`, `line_search`

#### Additional information
This addresses the issue where users couldn't determine default values like `xtol` and `eps` without looking at the source code. All default values are now explicitly documented in the docstrings following the pattern `parameter : type, optional, default: value`.

This is a documentation-only change with no code modifications.
