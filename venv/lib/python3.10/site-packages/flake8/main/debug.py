"""Module containing the logic for our debugging logic."""
import platform
from typing import Any
from typing import Dict

from flake8.plugins.finder import Plugins


def information(version: str, plugins: Plugins) -> Dict[str, Any]:
    """Generate the information to be printed for the bug report."""
    versions = sorted(
        {
            (loaded.plugin.package, loaded.plugin.version)
            for loaded in plugins.all_plugins()
            if loaded.plugin.package not in {"flake8", "local"}
        }
    )
    return {
        "version": version,
        "plugins": [
            {"plugin": plugin, "version": version}
            for plugin, version in versions
        ],
        "platform": {
            "python_implementation": platform.python_implementation(),
            "python_version": platform.python_version(),
            "system": platform.system(),
        },
    }
