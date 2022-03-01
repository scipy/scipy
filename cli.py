"""
Info: Run tests, builds and other tasks using, typer and system package
-------
    cmd options:
        $ python cli.py --help
        $ python cli.py bench <flag: -t/-s>
        $ python cli.py test <flag: -t/-s> <module>
        $ python cli.py build
"""

from typing import Optional
import typer
import os
app = typer.Typer()


@app.command()
def build():
    """
    Run Scipy build
    """
    os.system("python dev.py --build-only")


@app.command(context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
def bench(ctx: typer.Context):
    """
    Run benchmark options
    """
    for extra_arg in ctx.args:
        os.system(f"python dev.py --bench {extra_arg} integrate.SolveBVP")


@app.command(context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
def test(ctx: typer.Context):
    """
    Run test for a given module
    """
    _arg = ""
    for extra_arg in ctx.args:
        _arg += extra_arg + " "
    os.system(f"python dev.py --no-build {_arg}")


if __name__ == '__main__':
    app()
