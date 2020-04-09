************************************************************************
*
*  The data in this file represents the following problem:
*
*  Minimize or maximize Z = x1 + 2x5 - x8
*
*  Subject to:
*
*  2.5 <=   3x1 +  x2         -  2x4 - x5               -    x8
*                 2x2 + 1.1x3                                   <=  2.1
*                          x3              +  x6                 =  4.0
*  1.8 <=                      2.8x4             -1.2x7         <=  5.0
*  3.0 <= 5.6x1                      + x5               + 1.9x8 <= 15.0
*
*  where:
*
*  2.5 <= x1
*    0 <= x2 <= 4.1
*    0 <= x3
*    0 <= x4
*  0.5 <= x5 <= 4.0
*    0 <= x6
*    0 <= x7
*    0 <= x8 <= 4.3
*
*  x3, x4 are 0,1 variables.
*
*************************************************************************
*
*N=8,  M=5,  NZ= 14
*-----cost-----
* 1.0  0.0  0.0  0.0  2.0  0.0  0.0 -1.0 
*------A------
* 3.0  1.0      -2.0 -1.0           -1.0 
*     2.0  1.1                          
*           1.0            1.0           
*                2.8           -1.2      
* 5.6                 1.0            1.9 
*------LB------
* 2.5  -inf 4.0  1.8  3.0 
*------UB------
*  inf 2.1  4.0  5.0 15.0 
*
*
*************************************************************************
NAME          EXAMPLE
ROWS
 N  OBJ
 G  ROW01
 L  ROW02
 E  ROW03
 G  ROW04
 L  ROW05
COLUMNS
    COL01     OBJ                1.0
    COL01     ROW01              3.0   ROW05              5.6
    COL02     ROW01              1.0   ROW02              2.0
*
*  Mark COL03 and COL04 as integer variables.
*
    INT1      'MARKER'                 'INTORG'
    COL03     ROW02              1.1   ROW03              1.0
    COL04     ROW01             -2.0   ROW04              2.8
    INT1END   'MARKER'                 'INTEND'
*
    COL05     OBJ                2.0
    COL05     ROW01             -1.0   ROW05              1.0
    COL06     ROW03              1.0
    COL07     ROW04             -1.2
    COL08     OBJ               -1.0
    COL08     ROW01             -1.0   ROW05              1.9
RHS
    RHS1      ROW01              2.5
    RHS1      ROW02              2.1
    RHS1      ROW03              4.0
    RHS1      ROW04              1.8
    RHS1      ROW05             15.0
RANGES
    RNG1      ROW04              3.2
    RNG1      ROW05             12.0
BOUNDS
 LO BND1      COL01              2.5
 UP BND1      COL02              4.1
 LO BND1      COL05              0.5
 UP BND1      COL05              4.0
 UP BND1      COL08              4.3
ENDATA
