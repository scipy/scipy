"""
Info:
-------
    Basic example combining typer and system package
    cmd options:
        $ python cli.py bench -t
        $ python cli.py test -t cluster
        $ python cli.py build
"""

from typing import Optional
import typer
import os
app = typer.Typer()


@app.command()
def build():
    os.system("python dev.py --build-only")


@app.command(context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
def bench(ctx: typer.Context):
    for extra_arg in ctx.args:
        os.system(f"python dev.py --bench {extra_arg} integrate.SolveBVP")


@app.command(context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
def test(ctx: typer.Context):
    _arg = ""
    for extra_arg in ctx.args:
        _arg += extra_arg + " "
    os.system(f"python dev.py --no-build {_arg}")


if __name__ == '__main__':
    app()
