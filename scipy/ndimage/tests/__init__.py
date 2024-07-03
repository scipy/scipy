from __future__ import annotations
import numpy as np

# list of numarray data types
integer_types: list[type] = [
    "int8", "uint8", "int16", "uint16",
    "int32", "uint32", "int64", "uint64"]

float_types: list[type] = ["float32", "float64"]

complex_types: list[type] = ["complex64", "complex128"]

types: list[type] = integer_types + float_types
