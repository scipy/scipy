"""
Info:
-------
    Basic example combining typer and system package
    cmd to run:  python cli.py build -t
"""

from typing import Optional
import typer
import os
app = typer.Typer()


@app.command(context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
def build(ctx: typer.Context):
    for extra_arg in ctx.args:
        os.system(f'python dev.py --bench {extra_arg} integrate.SolveBVP')

if __name__ == '__main__':
    app()
