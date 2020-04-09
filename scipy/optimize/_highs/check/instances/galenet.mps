NAME          GALENET
ROWS
 L  S1
 L  S2
 L  S3
 E  NODE4
 E  NODE5
 G  D6
 G  D7
 G  D8
 N  COST
COLUMNS
    T14       S1                  1.   NODE4               1.
    T24       S2                  1.   NODE4               1.
    T25       S2                  1.   NODE5               1.
    T35       S3                  1.   NODE5               1.
    T46       D6                  1.   NODE4              -1.
    T47       D7                  1.   NODE4              -1.
    T57       D7                  1.   NODE5              -1.
    T58       D8                  1.   NODE5              -1.
RHS
    RHS       S1                 20.   S2                 20.
    RHS       S3                 20.   D6                 10.
    RHS       D7                 20.   D8                 30.
BOUNDS
 UP BND       T14                30.
 UP BND       T24                20.
 UP BND       T25                10.
 UP BND       T35                10.
 UP BND       T46                10.
 UP BND       T47                 2.
 UP BND       T57                20.
 UP BND       T58                30.
ENDATA
