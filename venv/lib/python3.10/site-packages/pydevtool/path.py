
from pathlib import Path

def rglob(paths, pattern='*.py', exclude=()):
    """
    exclude: list of (str, Path) with dir prefix. no GLOB pattern
    """
    for path in paths:
        for match in Path(path).rglob(pattern):
            for exclude_path in exclude:
                if match.is_relative_to(exclude_path):
                    break
            else:
                yield match
