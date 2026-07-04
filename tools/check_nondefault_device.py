#!/usr/bin/env python
"""Flag internal array-creation calls that omit an explicit ``device=``.

Under the array API standard, array-*creation* functions (``xp.zeros``,
``xp.full``, ``xp.arange``, ``xp.asarray(<scalar>)``, ...) place the new array on
the library's *default* device when ``device=`` is not passed.  If that array is
later combined with an input array living on a *non-default* device, the
operation raises a device-mismatch error.  The fix is to propagate the input's
device: ``xp.zeros(..., device=xp_device(x))``.

This checker flags a call when ALL of the following hold:

1. It is ``<ns>.<fn>(...)`` where the receiver ``<ns>`` is the array namespace --
   the bare name ``xp`` or an attribute ending in ``xp``/``_xp`` (e.g.
   ``self._xp``, ``self.xp``) -- and ``<fn>`` is a pure-creation function
   (``zeros``, ``ones``, ``empty``, ``full``, ``arange``, ``eye``, ``linspace``),
   OR ``<ns>.fft.rfftfreq`` / ``<ns>.fft.fftfreq``, OR ``<ns>.asarray(<literal>)``
   where the first positional argument is a scalar / list / tuple literal (from
   which no device can be inferred).
2. No ``device=`` keyword is present (and no ``**kwargs``, which might supply
   one -- treated leniently to avoid false positives on wrappers).
3. No line spanned by the call carries the ``# skip device check`` pragma.

Exempt: ``xp.asarray(<array-or-expr>)`` and ``xp.<...>_like(<array>)`` infer the
device from their first argument, so they are fine.

The check is a necessary condition (a ``device=`` is present), not a
sufficient one (that it is actually the right device and propagates end-to-end).
Runtime testing is hence useful too, this is (at the time of writing) still in the
works, with a fixture-based test (see ``devices`` in ``conftest.py``).

The namespace is recognised as the bare name ``xp`` or as an attribute ending in
``xp``/``_xp`` (``self._xp``, ``self.xp``).  Other aliases (e.g. a bare ``mxp``
for a masked-array namespace) are not recognised.

See scipy/scipy#22680 for the canonical issue that motivated the creation of this
linter.

Usage::

    # Integrated as part of the default linting task in spin:
    spin lint

    # Direct invocation of this checker:
    python tools/check_nondefault_device.py            # check the scipy/ tree
    python tools/check_nondefault_device.py scipy/stats/_stats_py.py  # subset

Exit code 1 if any violation is found.  Add ``# skip device check`` to a line to
allow a genuine standalone-scratch array that is never combined with the input.
"""
import ast
import sys
from pathlib import Path

PRAGMA = "# skip device check"

# The array namespace is referenced either as the bare name ``xp`` or as an
# attribute ending in ``xp``/``_xp`` (e.g. ``self._xp``, ``self.xp``).
NS_NAMES = {"xp"}
NS_ATTRS = {"xp", "_xp"}

PURE_CREATION = {
    "zeros", "ones", "empty", "full", "arange", "eye", "linspace",
}
LIKE = {"zeros_like", "ones_like", "empty_like", "full_like"}
FFT_FREQ = {"rfftfreq", "fftfreq"}


def _is_ns(node):
    """True if ``node`` refers to the array namespace (``xp``, ``self._xp``, ...)."""
    if isinstance(node, ast.Name):
        return node.id in NS_NAMES
    if isinstance(node, ast.Attribute):
        return node.attr in NS_ATTRS
    return False


def _is_ns_attr(node, attr_name):
    """True if ``node`` is ``<ns>.<attr_name>`` (``xp.zeros``, ``self._xp.zeros``)."""
    return (
        isinstance(node, ast.Attribute)
        and node.attr == attr_name
        and _is_ns(node.value)
    )


def _is_ns_fft(node):
    """Return the freq fn name if ``node`` is ``<ns>.fft.rfftfreq``/``fftfreq``."""
    if (
        isinstance(node, ast.Attribute)
        and node.attr in FFT_FREQ
        and isinstance(node.value, ast.Attribute)
        and node.value.attr == "fft"
        and _is_ns(node.value.value)
    ):
        return node.attr
    return None


def _is_literal(node):
    """Scalar/list/tuple literal => asarray cannot infer a device."""
    if isinstance(node, ast.Constant):
        return True
    if isinstance(node, (ast.List, ast.Tuple)):
        return True
    if isinstance(node, ast.UnaryOp) and isinstance(node.operand, ast.Constant):
        return True  # e.g. -1
    return False


def _has_device_kwarg(call):
    """True if the call passes ``device=`` or a ``**kwargs`` splat (lenient)."""
    return any(kw.arg == "device" or kw.arg is None for kw in call.keywords)


class DeviceVisitor(ast.NodeVisitor):
    def __init__(self, source):
        self.lines = source.splitlines()
        self.violations = []  # list of (lineno, col, what)

    def _allowed(self, node):
        """True if any line spanned by the call carries the pragma.

        Checking the whole span (not just the first line) keeps the escape
        hatch working when a formatter wraps a call over multiple lines.
        """
        end = getattr(node, "end_lineno", None) or node.lineno
        return any(PRAGMA in line for line in self.lines[node.lineno - 1:end])

    def visit_Call(self, node):
        self.generic_visit(node)
        func = node.func
        name = None
        first = node.args[0] if node.args else None

        # xp.<pure creation>(...)
        for candidate in PURE_CREATION:
            if _is_ns_attr(func, candidate):
                name = candidate
                break

        # xp.<...>_like(<reference array>) infers the device -> OK
        if name is None:
            for candidate in LIKE:
                if _is_ns_attr(func, candidate):
                    if first is not None and not _is_literal(first):
                        return
                    name = candidate
                    break

        # xp.fft.rfftfreq / fftfreq
        if name is None:
            fft = _is_ns_fft(func)
            if fft:
                name = f"fft.{fft}"

        # xp.asarray(<literal>) cannot infer a device;
        # xp.asarray(<array/expr>) can, so it is exempt.
        if name is None and _is_ns_attr(func, "asarray"):
            if first is not None and _is_literal(first):
                name = "asarray"
            else:
                return

        if name is None:
            return
        if _has_device_kwarg(node):
            return
        if self._allowed(node):
            return
        self.violations.append((node.lineno, node.col_offset, ast.unparse(func)))


def check_source(source):
    """Return a list of (lineno, col, what) violations for a source string.

    Exposed for unit testing.
    """
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return []
    visitor = DeviceVisitor(source)
    visitor.visit(tree)
    return visitor.violations


def check_file(path):
    src = path.read_text(encoding="utf-8", errors="replace")
    if "xp." not in src:  # not array-API-aware; skip (matches xp. and self._xp.)
        return []
    try:
        ast.parse(src)
    except SyntaxError as exc:
        # Warn rather than silently skip: usually the interpreter is older than
        # the syntax used in the file (e.g. PEP 695 `type` aliases need 3.12+),
        # which would otherwise hide real violations. Run under the project's
        # minimum supported Python (as CI does).
        print(
            f"{path}:{exc.lineno}: warning: could not parse under Python "
            f"{sys.version_info.major}.{sys.version_info.minor} "
            f"({exc.msg}); file not checked",
            file=sys.stderr,
        )
        return []
    return check_source(src)


def iter_py(root, submodule_paths):
    for path in sorted(root.rglob("*.py")):
        s = str(path.absolute())
        if any(sub in s for sub in submodule_paths):
            continue
        if "/tests/" in path.as_posix() or path.name == "conftest.py":
            continue
        yield path


def main(argv):
    from get_submodule_paths import get_submodule_paths

    submodule_paths = get_submodule_paths()
    # Default to the scipy/ package relative to the repo root, so the checker
    # works regardless of the current working directory (e.g. `cwd = tools`).
    default_target = Path(__file__).resolve().parent.parent / "scipy"
    targets = [Path(a) for a in argv] or [default_target]
    total = 0
    for target in targets:
        files = [target] if target.is_file() else iter_py(target, submodule_paths)
        for path in files:
            for lineno, col, what in check_file(path):
                print(f"{path}:{lineno}:{col}: missing device= in {what}(...)")
                total += 1
    if total:
        print(
            f"\n{total} device-propagation violation(s); add `device=xp_device(x)`"
            f" or a `{PRAGMA}` pragma.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
