"""
rich-click is a minimal Python module to combine the efforts of the excellent packages 'rich' and 'click'.

The intention is to provide attractive help output from click, formatted with rich, with minimal
customisation required.
"""

__version__ = "1.6.1"

from typing import TYPE_CHECKING

from click import *  # noqa: F401, F403
from click import command as click_command
from click import group as click_group

from . import rich_click  # noqa: F401

from rich_click.rich_command import RichCommand
from rich_click.rich_group import RichGroup

# MyPy does not like star imports. Therefore when we are type checking, we import each individual module
# from click here. This way MyPy will recognize the import and not throw any errors. Furthermore, because of
# the TYPE_CHECKING check, it does not influence the start routine at all.
if TYPE_CHECKING:
    from click import argument, Choice, option, Path, version_option  # noqa: F401

    __all__ = [
        "argument",
        "Choice",
        "option",
        "Path",
        "version_option",
        "group",
        "command",
    ]


def group(*args, cls=RichGroup, **kwargs):
    """
    Group decorator function.

    Defines the group() function so that it uses the RichGroup class by default.
    """
    return click_group(*args, cls=cls, **kwargs)


def command(name=None, cls=RichCommand, **attrs):
    """
    Command decorator function.

    Defines the command() function so that it uses the RichCommand class by default.
    """
    if callable(name) and cls:
        return click_command(cls=cls, **attrs)(name)

    return click_command(name, cls=cls, **attrs)
