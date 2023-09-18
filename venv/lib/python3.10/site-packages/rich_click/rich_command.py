import sys

import click

from rich_click.rich_click import rich_abort_error, rich_format_error, rich_format_help


class RichCommand(click.Command):
    """Richly formatted click Command.

    Inherits click.Command and overrides help and error methods
    to print richly formatted output.
    """

    standalone_mode = False

    def main(self, *args, standalone_mode: bool = True, **kwargs):
        try:
            rv = super().main(*args, standalone_mode=False, **kwargs)
            if not standalone_mode:
                return rv
        except click.ClickException as e:
            if not standalone_mode:
                raise
            rich_format_error(e)
            sys.exit(e.exit_code)
        except click.exceptions.Abort:
            if not standalone_mode:
                raise
            rich_abort_error()
            sys.exit(1)

    def format_help(self, ctx: click.Context, formatter: click.HelpFormatter):
        rich_format_help(self, ctx, formatter)
