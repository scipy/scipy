        SUBROUTINE UMD2RG (XX, XSIZE, XHEAD, XTAIL, XUSE,
     $          LUI, FRDIMC, FRXP, FRNEXT, FRPREV, NLU, LUP,
     $          ICNTL, FFXP, FFSIZE, PFREE, XFREE)
        INTEGER LUI (*), NLU, FRDIMC (NLU+2), FRXP (NLU+2),
     $          FRNEXT (NLU+2), FRPREV (NLU+2), LUP (NLU),
     $          ICNTL (20), XSIZE, XUSE, XHEAD, XTAIL, FFXP, FFSIZE,
     $          PFREE, XFREE
        DOUBLE PRECISION
     $          XX (XSIZE)
        
C=== UMD2RG ============================================================
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
C  Garbage collection for UMD2R2.

C=======================================================================
C  INPUT:
C=======================================================================
C
C       XX:             real workspace, containing matrix being
C                       factorized and partially-computed LU factors
C       xsize:          size of XX
C       xhead:          XX (1..xhead) is in use (matrix, frontal mtc's)
C       xtail:          XX (xtail..xsize) is in use (LU factors)
C       xuse:           memory usage in Value
C       Icntl:          integer control parameters, see UMD21I
C       ffxp:           pointer to current contribution block
C       ffsize:         size of current contribution block
C       nlu:            number of LU arrowheads
C
C       FRdimc (1..nlu+2)       leading dimension of frontal matrices
C       FRxp (1..nlu+2)         pointer to frontal matrices in XX
C       FRnext (1..nlu+2)       pointer to next block in XX
C       FRprev (1..nlu+2)       pointer to previous block in XX
C       LUp (1..nlu)            pointer to LU arrowhead patters in LUi
C       LUi (*)                 pattern of LU factors

C=======================================================================
C  OUTPUT:
C=======================================================================
C
C       XX:             external fragmentation is removed at head 
C       xhead:          XX (1..xhead) is in use, reduced in size
C       xuse:           memory usage in Value, reduced
C       pfree:          pointer to free block in memory list, set to 0
C       xfree:          size of free block in XX, set to -1
C       FRdimc          arrays for frontal matrices are compressed
C       FRxp            frontal matrices have been shifted
C       ffxp            current working array has been shifted

C=======================================================================
C  SUBROUTINES AND FUNCTIONS CALLED / CALLED BY:
C=======================================================================
C
C       called by subroutine:   UMD2R2
C       functions called:       ABS
        INTRINSIC ABS

C=======================================================================
C  LOCAL SCALARS:
C=======================================================================

        INTEGER XDP, I, E, FDIMC, LUDEGR, LUDEGC, J, FLUIP, FXP,
     $          MHEAD, MTAIL

C  xdp:     real destination pointer, current block moved to XX (xdp...)
C  e:       an element
C  fdimc:   column dimension (number of rows) of a frontal matrix
C  ludegr:  row degree (number of columns) of a contribution block
C  ludegc:  column degree (number of rows) of a contribution block
C  fluip:   element is in LUi (fluip...)
C  fxp:     element is in XX (fxp...) prior to compression
C  mhead:   nlu+1, head pointer for contribution block link list
C  mtail:   nlu+2, tail pointer for contribution block link list
C  i:       general loop index
C  j:       general loop index

C=======================================================================
C  EXECUTABLE STATMENTS:
C=======================================================================

C-----------------------------------------------------------------------
C  scan the link list and compress the reals
C-----------------------------------------------------------------------

        MHEAD = NLU+1
        MTAIL = NLU+2
        XDP = FRXP (MHEAD)
        E = FRNEXT (MHEAD)

C       while (e .ne. mtail) do
10      CONTINUE
        IF (E .NE. MTAIL) THEN 

           FDIMC = FRDIMC (E)

C          -------------------------------------------------------------
           IF (FDIMC .EQ. 0) THEN 
C          -------------------------------------------------------------

C             this is a real hole - delete it from the link list

              FRNEXT (FRPREV (E)) = FRNEXT (E)
              FRPREV (FRNEXT (E)) = FRPREV (E)

C          -------------------------------------------------------------
           ELSE 
C          -------------------------------------------------------------

C             this is an unassembled frontal matrix
              FXP = FRXP (E)
              FRXP (E) = XDP
              FLUIP = LUP (E)
              LUDEGR = ABS (LUI (FLUIP+2))
              LUDEGC = ABS (LUI (FLUIP+3))
              IF (FDIMC .EQ. LUDEGC) THEN 
C                contribution block is already compressed
CFPP$ NODEPCHK L
                 DO 20 I = 0, (LUDEGR * LUDEGC) - 1 
                    XX (XDP+I) = XX (FXP+I)
20               CONTINUE 
              ELSE 
C                contribution block is not compressed
C                compress XX (fxp..) to XX (xdp..xdp+(ludegr*ludegc)-1)
                 DO 40 J = 0, LUDEGR - 1 
CFPP$ NODEPCHK L
                    DO 30 I = 0, LUDEGC - 1 
                       XX (XDP + J*LUDEGC + I) = XX (FXP + J*FDIMC + I)
30                  CONTINUE 
40               CONTINUE 
                 FRDIMC (E) = LUDEGC
              ENDIF 
              XDP = XDP + LUDEGR*LUDEGC

           ENDIF 

C          -------------------------------------------------------------
C          get the next item in the link list
C          -------------------------------------------------------------

           E = FRNEXT (E)

C       end while:
        GOTO 10
        ENDIF 

        FRXP (MTAIL) = XDP
        PFREE = 0
        XFREE = -1

C       ----------------------------------------------------------------
C       shift the current working array (if it exists)
C       ----------------------------------------------------------------

        IF (FFXP .NE. 0) THEN 
CFPP$ NODEPCHK L
           DO 50 I = 0, FFSIZE - 1 
              XX (XDP+I) = XX (FFXP+I)
50         CONTINUE 
           FFXP = XDP
           XDP = XDP + FFSIZE
        ENDIF 

C-----------------------------------------------------------------------
C  deallocate the unused space
C-----------------------------------------------------------------------

        XUSE = XUSE - (XHEAD - XDP)
        XHEAD = XDP
        RETURN
        END 
