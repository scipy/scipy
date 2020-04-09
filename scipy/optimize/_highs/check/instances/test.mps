* Using Xpress-MP extensions
NAME          moselP
ROWS
 N  *OBJ*
 G  r1
 L  r2
 L  r3
 G  r4
 L  r5
 E  r6
 G  r7

COLUMNS
    x(1)      r6        1
    x(1)      r4        1
    x(1)      *OBJ*     1
    x(2)      r7        1
    x(2)      r4        2
    x(2)      *OBJ*     1
    x(3)      r7        -1
    x(3)      r5        -1
    x(3)      r4        -2
    x(3)      r3        2
    x(4)      r7        1
    x(4)      r5        -1
    x(4)      r4        -2
    x(4)      *OBJ*     -1
    x(5)      r2        1
    x(5)      r1        1
    x(5)      *OBJ*     1
    x(6)      r4        1
    x(6)      r6        1
    x(6)      *OBJ*     3
    x(7)      r3        1
    x(7)      *OBJ*     -1

RHS
    *RHS*     r1        1
    *RHS*     r2        3
    *RHS*     r3        10
    *RHS*     r4        1.5
    *RHS*     r5        1
    *RHS*     r6        1
    *RHS*     r7        1

RANGES
    test      r3        6

BOUNDS
 LO BND2      x(3)      1
 UP BND2      x(3)      2

ENDATA
