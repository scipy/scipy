NAME          AVGAS SIZE: N=8, M=11, NZ=38
ROWS
 G  R1
 G  R2
 G  R3
 G  R4
 G  R5
 G  R6
 G  R7
 G  R8
 G  R9
 G  R10
 N  COST
COLUMNS
    C1        R1           -1.0        R5         -1.0
    C1        R7            2.0        R8          5.0
    C2        R1           -1.0        R6         -1.0
    C2        R9            1.0        R10         1.0
    C2        COST         -2.0
    C3        R2           -1.0        R5         -1.0
    C3        R7            1.0        R8          3.0
    C3        COST         -1.0
    C4        R2           -1.0        R6         -1.0
    C4        R9           -1.0
    C4        COST         -3.0
    C5        R3           -1.0        R5         -1.0
    C5        R8           -3.0
    C5        COST         -2.0
    C6        R3           -1.0        R6         -1.0
    C6        R9           -3.0        R10        -3.0
    C6        COST         -4.0
    C7        R4           -1.0        R5         -1.0
    C7        R7           -1.0        R8         -1.0
    C7        COST         -3.0
    C8        R4           -1.0        R6         -1.0
    C8        R9           -5.0        R10        -2.0
    C8        COST         -5.0
RHS
    DEMANDS   R1           -1.0        R2         -1.0
    DEMANDS   R3           -1.0        R4         -1.0
    DEMANDS   R5           -2.0        R6         -2.0
BOUNDS
 UP SERVINGS  C1            1.0
 UP SERVINGS  C2            1.0
 UP SERVINGS  C3            1.0
 UP SERVINGS  C4            1.0
 UP SERVINGS  C5            1.0
 UP SERVINGS  C6            1.0
 UP SERVINGS  C7            1.0
 UP SERVINGS  C8            1.0
ENDATA
