        SUBROUTINE UMD2ER (WHO, ICNTL, INFO, ERROR, S)
        INTEGER WHO, ICNTL (20), INFO (40), ERROR, S
        
C=== UMD2ER ============================================================
C
C  Unsymmetric-pattern MultiFrontal Package (UMFPACK). Version 2.2d
C  Copyright (C) 1997, Timothy A. Davis, University of Florida, USA.
C  ALL RIGHTS RESERVED.
C  Joint work with Iain S. Duff, Rutherford Appleton Laboratory, UK.
C  July 7, 1997. Work supported by the National Science Foundation
C  (DMS-9223088 and DMS-9504974) and the State of Florida; and by CRAY
C  Research Inc. through the allocation of supercomputing resources.

C***********************************************************************
C* NOTICE:  "The UMFPACK Package may be used SOLELY for educational,   *
C* research, and benchmarking purposes by non-profit organizations and *
C* the U.S. government.  Commercial and other organizations may make   *
C* use of UMFPACK SOLELY for benchmarking purposes only.  UMFPACK may  *
C* be modified by or on behalf of the User for such use but at no time *
C* shall UMFPACK or any such modified version of UMFPACK become the    *
C* property of the User.  UMFPACK is provided without warranty of any  *
C* kind, either expressed or implied.  Neither the Authors nor their   *
C* employers shall be liable for any direct or consequential loss or   *
C* damage whatsoever arising out of the use or misuse of UMFPACK by    *
C* the User.  UMFPACK must not be sold.  You may make copies of        *
C* UMFPACK, but this NOTICE and the Copyright notice must appear in    *
C* all copies.  Any other use of UMFPACK requires written permission.  *
C* Your use of UMFPACK is an implicit agreement to these conditions."  *
C*                                                                     *
C* The MA38 Package in Release 12 of the Harwell Subroutine Library    *
C* (HSL) has equivalent functionality (and identical calling interface)*
C* as UMFPACK (the HSL has single and double precision versions only,  *
C* however).  It is available for commercial use.   Technical reports, *
C* information on HSL, and matrices are available via the World Wide   *
C* Web at http://www.cis.rl.ac.uk/struct/ARCD/NUM.html, or by          *
C* anonymous ftp at seamus.cc.rl.ac.uk/pub.  Also contact Dr. Scott    *
C* Roberts, Harwell Subroutine Library, B 552, AEA Technology,         *
C* Harwell, Didcot, Oxon OX11 0RA, England.                            *
C* telephone (44) 1235 434988, fax (44) 1235 434136                    *
C* email Scott.Roberts@aeat.co.uk, who will provide details of price   *
C* and conditions of use.                                              *
C***********************************************************************

C=======================================================================
C  NOT USER-CALLABLE.

C=======================================================================
C  DESCRIPTION:
C=======================================================================
C
C  Print error and warning messages, and set error flags.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       who             which user-callable routine called:
C                       1: UMD2FA, 2: UMD2RF, 3: UMD2SO
C       Icntl (1):      I/O unit for error and warning messages
C       Icntl (3):      printing level
C       Info (1):       the error/warning status
C       error:          the applicable error (<0) or warning (>0).
C                       See UMD2P2 for a description.
C       s:              the relevant offending value

C=======================================================================
C  OUTPUT: 
C=======================================================================
C
C       Info (1):       the error/warning status

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutines:  UMD2CO, UMD2FA, UMD2F0, UMD2F1, UMD2F2,
C                               UMD2RF, UMD2R0, UMD2R2, UMD2SO, UMD2S2
C       subroutines called:     UMD2P2
C       functions called:       MOD
        INTRINSIC MOD

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        LOGICAL BOTH
        DOUBLE PRECISION
     $          IGNORE
        INTEGER IOERR, PRL

C  ioerr:   I/O unit for error and warning messages
C  prl:     printing level
C  both:    if true, then combine errors -3 and -4 into error -5

C=======================================================================
C  EXECUTABLE STATEMENTS:
C=======================================================================

        IGNORE = 0
        IOERR = ICNTL (1)
        PRL = ICNTL (3)
        IF (ERROR .LT. 0) THEN 
C          this is an error message
           BOTH = (INFO (1) .EQ. -3 .AND. ERROR .EQ. -4) .OR.
     $            (INFO (1) .EQ. -4 .AND. ERROR .EQ. -3)
           IF (BOTH) THEN 
C             combine error -3 (out of integer memory) and error -4
C             (out of real memory)
              INFO (1) = -5
           ELSE 
              INFO (1) = ERROR
           ENDIF 
           IF (PRL .GE. 1) THEN 
              CALL UMD2P2 (WHO, ERROR, S, 0, IGNORE, IOERR)
           ENDIF 
        ELSE IF (ERROR .GT. 0) THEN 
C          this is a warning message
           IF (INFO (1) .GE. 0) THEN 
C             do not override a prior error setting, sum up warnings
              IF (MOD (INFO (1) / ERROR, 2) .EQ. 0) THEN 
                 INFO (1) = INFO (1) + ERROR
              ENDIF 
           ENDIF 
           IF (PRL .GE. 2) THEN 
              CALL UMD2P2 (WHO, ERROR, S, 0, IGNORE, IOERR)
           ENDIF 
        ENDIF 
        RETURN
        END 
