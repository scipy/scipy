"""Shared definitions of _UFuncDef to avoid circular dependencies between modules."""

from collections import namedtuple

_UFuncDef = namedtuple("_UFuncDef", ("name", "loop_fun_names", "loop_funs", "types", "num_inputs", "includes"))
