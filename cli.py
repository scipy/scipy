"""
Info:
-------
    Basic example combining typer and system package
    cmd to run: typer cli.py run --param -t

Limitations encountered:
------------------------
    - Unable to pass parameters with '-' directly (Error: no such option: -t)
    - Using typer.Options, seems similar to doit params and not solving our use-case
"""

from typing import Optional
import typer
import os
app = typer.Typer()


@app.command()
def build(param: Optional[str] = None):
    os.system(f'python dev.py --bench {param} integrate.SolveBVP')


if __name__ == '__main__':
    typer.run(build)
