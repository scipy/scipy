'''
This is a module defining exit flags.

Translated from the modern-Fortran reference implementation in PRIMA by Zaikun ZHANG (www.zhangzk.net).

Dedicated to late Professor M. J. D. Powell FRS (1936--2015).

Python implementation by Nickolai Belakovski
'''

INFO_DEFAULT = 0
SMALL_TR_RADIUS = 0
FTARGET_ACHIEVED = 1
TRSUBP_FAILED = 2
MAXFUN_REACHED = 3
MAXTR_REACHED = 20
NAN_INF_X = -1
NAN_INF_F = -2
NAN_INF_MODEL = -3
NO_SPACE_BETWEEN_BOUNDS = 6
DAMAGING_ROUNDING = 7
ZERO_LINEAR_CONSTRAINT = 8
CALLBACK_TERMINATE = 30

# Stop-codes.
# The following codes are used by ERROR STOP as stop-codes, which should be default integers.
INVALID_INPUT = 100
ASSERTION_FAILS = 101
VALIDATION_FAILS = 102
MEMORY_ALLOCATION_FAILS = 103
