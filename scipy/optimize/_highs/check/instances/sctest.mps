* Using Xpress-MP extensions
NAME          moselP
ROWS
 N  *OBJ*
 G  r1
 L  r2
 E  r3
 G  r4
 L  r5
 E  r6
 G  r7
 E  r8
 E  r9
 E  r10

COLUMNS
    x(1)      r10       1
    x(1)      r9        1
    x(1)      r8        1
    x(1)      r7        1
    x(1)      r6        1
    x(1)      r4        2
    x(1)      *OBJ*     1
    x(2)      r10       1
    x(2)      r9        1
    x(2)      r7        1
    x(2)      r4        2
    x(2)      *OBJ*     1
    x(3)      r10       -1
    x(3)      r7        -1
    x(3)      r5        -1
    x(3)      r4        -2
    x(3)      r3        2
    x(3)      *OBJ*     1
    x(4)      r10       -1
    x(4)      r8        -1
    x(4)      r7        1
    x(4)      r5        -1
    x(4)      r4        -2
    x(4)      *OBJ*     -1
    x(5)      r2        1
    x(5)      r1        1
    x(5)      *OBJ*     1
    x(6)      r10       1
    x(6)      r9        -1
    x(6)      r8        1
    x(6)      r7        -1
    x(6)      r6        1
    x(6)      r5        1
    x(6)      r4        2

RHS
    *RHS*     r1        3
    *RHS*     r2        3
    *RHS*     r3        1
    *RHS*     r4        1.5
    *RHS*     r5        1
    *RHS*     r6        1
    *RHS*     r7        1
    *RHS*     r8        1
    *RHS*     r9        2
    *RHS*     r10       2

BOUNDS

ENDATA
