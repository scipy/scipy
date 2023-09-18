"""The command line interface."""

import sys
from importlib import import_module
from textwrap import dedent
from typing import Any, List, Optional

try:
    from importlib.metadata import entry_points
except ImportError:
    # Support Python <3.8
    from importlib_metadata import entry_points

import click
from rich.console import Console
from rich.padding import Padding
from rich.panel import Panel
from rich.text import Text

from rich_click import command as rich_command
from rich_click import group as rich_group
from rich_click import RichCommand, RichGroup
from rich_click.rich_click import (
    ALIGN_ERRORS_PANEL,
    ERRORS_PANEL_TITLE,
    STYLE_ERRORS_PANEL_BORDER,
    STYLE_HELPTEXT,
    STYLE_HELPTEXT_FIRST_LINE,
    STYLE_USAGE,
    STYLE_USAGE_COMMAND,
)

console = Console()


def _print_usage() -> None:
    console.print(
        Padding(
            Text.from_markup(f"[{STYLE_USAGE}]Usage[/]: rich-click [SCRIPT | MODULE:FUNCTION] [-- SCRIPT_ARGS...]"),
            1,
        ),
        style=STYLE_USAGE_COMMAND,
    )


def _print_help() -> None:
    help_paragraphs = dedent(main.__doc__ or "").split("\n\n")
    help_paragraphs = [x.replace("\n", " ").strip() for x in help_paragraphs]
    console.print(
        Padding(
            Text.from_markup(help_paragraphs[0].strip()),
            (0, 1),
        ),
        style=STYLE_HELPTEXT_FIRST_LINE,
    )
    console.print(
        Padding(
            Text.from_markup("\n\n".join(help_paragraphs[1:]).strip()),
            (0, 1),
        ),
        style=STYLE_HELPTEXT,
    )


def patch() -> None:
    """Patch Click internals to use Rich-Click types."""
    click.group = rich_group
    click.command = rich_command
    click.Group = RichGroup
    click.Command = RichCommand


def main(args: Optional[List[str]] = None) -> Any:
    """
    The [link=https://github.com/ewels/rich-click]rich-click[/] CLI provides attractive help output from any
    tool using [link=https://click.palletsprojects.com/]click[/], formatted with
    [link=https://github.com/Textualize/rich]rich[/].

    The rich-click command line tool can be prepended before any Python package
    using native click to provide attractive richified click help output.

    For example, if you have a package called [blue]my_package[/] that uses click,
    you can run:

    [blue]  rich-click my_package --help  [/]

    It only works if the package is using vanilla click without customised [cyan]group()[/]
    or [cyan]command()[/] classes.
    If in doubt, please suggest to the authors that they use rich_click within their
    tool natively - this will always give a better experience.
    """  # noqa: D400, D401
    args = args or sys.argv[1:]
    if not args or args == ["--help"]:
        # Print usage if we got no args, or only --help
        _print_usage()
        _print_help()
        sys.exit(0)
    else:
        script_name = args[0]
    scripts = {script.name: script for script in entry_points().get("console_scripts")}
    if script_name in scripts:
        # a valid script was passed
        script = scripts[script_name]
        module_path, function_name = script.value.split(":", 1)
        prog = script_name
    elif ":" in script_name:
        # the path to a function was passed
        module_path, function_name = args[0].split(":", 1)
        prog = module_path.split(".", 1)[0]
    else:
        _print_usage()
        console.print(
            Panel(
                Text.from_markup(f"No such script: [bold]{script_name}[/]"),
                border_style=STYLE_ERRORS_PANEL_BORDER,
                title=ERRORS_PANEL_TITLE,
                title_align=ALIGN_ERRORS_PANEL,
            )
        )
        console.print(
            Padding(
                "Please run [yellow bold]rich-click --help[/] for usage information.",
                (0, 1),
            ),
            style="dim",
        )
        sys.exit(1)
    if len(args) > 1:
        if args[1] == "--":
            del args[1]
    sys.argv = [prog, *args[1:]]
    # patch click before importing the program function
    patch()
    # import the program function
    module = import_module(module_path)
    function = getattr(module, function_name)
    # simply run it: it should be patched as well
    return function()
