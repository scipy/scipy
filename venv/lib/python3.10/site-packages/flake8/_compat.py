"""Expose backports in a single place."""
import sys

if sys.version_info >= (3, 8):  # pragma: no cover (PY38+)
    import importlib.metadata as importlib_metadata
else:  # pragma: no cover (<PY38)
    import importlib_metadata

if sys.version_info[:3] == (3, 8, 0):
    # backported from importlib.metadata in 3.8.1
    importlib_metadata.EntryPoint.__reduce__ = lambda self: (  # type: ignore
        type(self),
        (self.name, self.value, self.group),  # type: ignore
    )

__all__ = ("importlib_metadata",)
