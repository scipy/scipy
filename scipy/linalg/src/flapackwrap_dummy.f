      REAL FUNCTION WSLAMCH( CMACH )
      CHARACTER          CMACH
      EXTERNAL           SLAMCH
      REAL               SLAMCH
      WSLAMCH = SLAMCH( CMACH )
      END FUNCTION
