#!/usr/bin/env python3
"""Reproducer for issue #23867: special.comb(0, 0, repetition=True) returns 0 instead of 1."""

from scipy.special import comb

print("Testing comb(0, 0) with different parameters:")
print(f"comb(0, 0, exact=True, repetition=False) = {comb(0, 0, exact=True, repetition=False)}")
print(f"comb(0, 0, exact=False, repetition=False) = {comb(0, 0, exact=False, repetition=False)}")
print(f"comb(0, 0, exact=True, repetition=True) = {comb(0, 0, exact=True, repetition=True)}  <- BUG: should be 1")
print(f"comb(0, 0, exact=False, repetition=True) = {comb(0, 0, exact=False, repetition=True)}  <- BUG: should be 1.0")

print("\nAdditional tests - comb(n, 0) should always be 1:")
for n in [0, 1, 2, 5, 10]:
    result_exact = comb(n, 0, exact=True, repetition=True)
    result_float = comb(n, 0, exact=False, repetition=True)
    print(f"comb({n}, 0, repetition=True): exact={result_exact}, float={result_float}")
    if result_exact != 1:
        print(f"  ^^^ BUG: expected 1, got {result_exact}")
