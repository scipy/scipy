C       COMPUTATION OF SPECIAL FUNCTIONS
C
C          Shanjie Zhang and Jianming Jin
C
C       Copyrighted but permission granted to use code in programs.
C       Buy their book "Computation of Special Functions", 1996, John Wiley & Sons, Inc.
C
C
C      Compiled into a single source file and changed REAL To DBLE throughout.
C
C      Changed according to ERRATA also.
C
C      Changed GAMMA to GAMMA2 and PSI to PSI_SPEC to avoid potential conflicts.
C

        FUNCTION DNAN()
        DOUBLE PRECISION DNAN
        DNAN = 0.0D0
        DNAN = 0.0D0/DNAN
        END

        FUNCTION DINF()
        DOUBLE PRECISION DINF
        DINF = 1.0D300
        DINF = DINF*DINF
        END

        SUBROUTINE CPDSA(N,Z,CDN)
C
C       ===========================================================
C       Purpose: Compute complex parabolic cylinder function Dn(z)
C                for small argument
C       Input:   z   --- complex argument of D(z)
C                n   --- Order of D(z) (n = 0,-1,-2,...)
C       Output:  CDN --- Dn(z)
C       Routine called: GAIH for computing Г(x), x=n/2 (n=1,2,...)
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        EPS=1.0D-15
        PI=3.141592653589793D0
        SQ2=DSQRT(2.0D0)
        CA0=CDEXP(-.25D0*Z*Z)
        VA0=0.5D0*(1.0D0-N)
        IF (N.EQ.0.0) THEN
           CDN=CA0
        ELSE
           IF (CDABS(Z).EQ.0.0) THEN
              IF (VA0.LE.0.0.AND.VA0.EQ.INT(VA0)) THEN
                 CDN=0.0D0
              ELSE
                 CALL GAIH(VA0,GA0)
                 PD=DSQRT(PI)/(2.0D0**(-.5D0*N)*GA0)
                 CDN = DCMPLX(PD, 0.0D0)
              ENDIF
           ELSE
              XN=-N
              CALL GAIH(XN,G1)
              CB0=2.0D0**(-0.5D0*N-1.0D0)*CA0/G1
              VT=-.5D0*N
              CALL GAIH(VT,G0)
              CDN = DCMPLX(G0, 0.0D0)
              CR=(1.0D0,0.0D0)
              DO 10 M=1,250
                 VM=.5D0*(M-N)
                 CALL GAIH(VM,GM)
                 CR=-CR*SQ2*Z/M
                 CDW=GM*CR
                 CDN=CDN+CDW
                 IF (CDABS(CDW).LT.CDABS(CDN)*EPS) GO TO 20
10            CONTINUE
20            CDN=CB0*CDN
           ENDIF
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE CFS(Z,ZF,ZD)
C
C       =========================================================
C       Purpose: Compute complex Fresnel Integral S(z) and S'(z)
C       Input :  z  --- Argument of S(z)
C       Output:  ZF --- S(z)
C                ZD --- S'(z)
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (E,P,W)
        IMPLICIT COMPLEX *16 (C,S,Z)
        EPS=1.0D-14
        PI=3.141592653589793D0
        W0=CDABS(Z)
        ZP=0.5D0*PI*Z*Z
        ZP2=ZP*ZP
        Z0=(0.0D0,0.0D0)
        IF (Z.EQ.Z0) THEN
           S=Z0
        ELSE IF (W0.LE.2.5) THEN
           S=Z*ZP/3.0D0
           CR=S
           WB0=0.0D0
           DO 10 K=1,80
              CR=-.5D0*CR*(4.0D0*K-1.0D0)/K/(2.0D0*K+1.0D0)
     &          /(4.0D0*K+3.0D0)*ZP2
              S=S+CR
              WB=CDABS(S)
              IF (DABS(WB-WB0).LT.EPS.AND.K.GT.10) GO TO 30
10            WB0=WB
        ELSE IF (W0.GT.2.5.AND.W0.LT.4.5) THEN
           M=85
           S=Z0
           CF1=Z0
           CF0=(1.0D-100,0.0D0)
           DO 15 K=M,0,-1
              CF=(2.0D0*K+3.0D0)*CF0/ZP-CF1
              IF (K.NE.INT(K/2)*2) S=S+CF
              CF1=CF0
15            CF0=CF
           S=CDSQRT(2.0D0/(PI*ZP))*CDSIN(ZP)/CF*S
        ELSE
           CR=(1.0D0,0.0D0)
           CF=(1.0D0,0.0D0)
           DO 20 K=1,20
              CR=-.25D0*CR*(4.0D0*K-1.0D0)*(4.0D0*K-3.0D0)/ZP2
20            CF=CF+CR
           CR=1.0D0
           CG=CR
           DO 25 K=1,12
              CR=-.25D0*CR*(4.0D0*K+1.0D0)*(4.0D0*K-1.0D0)/ZP2
25            CG=CG+CR
           CG = CG/(PI*Z*Z)
           S=.5D0-(CF*CDCOS(ZP)+CG*CDSIN(ZP))/(PI*Z)
        ENDIF
30      ZF=S
        ZD=CDSIN(0.5*PI*Z*Z)
        RETURN
        END

C       **********************************

        SUBROUTINE LQMN(MM,M,N,X,QM,QD)
C
C       ==========================================================
C       Purpose: Compute the associated Legendre functions of the
C                second kind, Qmn(x) and Qmn'(x)
C       Input :  x  --- Argument of Qmn(x)
C                m  --- Order of Qmn(x)  ( m = 0,1,2,… )
C                n  --- Degree of Qmn(x) ( n = 0,1,2,… )
C                mm --- Physical dimension of QM and QD
C       Output:  QM(m,n) --- Qmn(x)
C                QD(m,n) --- Qmn'(x)
C       ==========================================================
C
        IMPLICIT DOUBLE PRECISION (Q,X)
        DIMENSION QM(0:MM,0:N),QD(0:MM,0:N)
        IF (DABS(X).EQ.1.0D0) THEN
           DO 10 I=0,M
           DO 10 J=0,N
              QM(I,J)=1.0D+300
              QD(I,J)=1.0D+300
10         CONTINUE
           RETURN
        ENDIF
        LS=1
        IF (DABS(X).GT.1.0D0) LS=-1
        XS=LS*(1.0D0-X*X)
        XQ=DSQRT(XS)
        Q0=0.5D0*DLOG(DABS((X+1.0D0)/(X-1.0D0)))
        IF (DABS(X).LT.1.0001D0) THEN
           QM(0,0)=Q0
           QM(0,1)=X*Q0-1.0D0
           QM(1,0)=-1.0D0/XQ
           QM(1,1)=-LS*XQ*(Q0+X/(1.0D0-X*X))
           DO 15 I=0,1
           DO 15 J=2,N
              QM(I,J)=((2.0D0*J-1.0D0)*X*QM(I,J-1)
     &               -(J+I-1.0D0)*QM(I,J-2))/(J-I)
15         CONTINUE
           DO 20 J=0,N
           DO 20 I=2,M
              QM(I,J)=-2.0D0*(I-1.0D0)*X/XQ*QM(I-1,J)-LS*
     &                (J+I-1.0D0)*(J-I+2.0D0)*QM(I-2,J)
20         CONTINUE
        ELSE
           IF (DABS(X).GT.1.1D0) THEN
              KM=40+M+N
           ELSE
              KM=(40+M+N)*INT(-1.0-1.8*LOG(X-1.0))
           ENDIF
           QF2=0.0D0
           QF1=1.0D0
           QF0=0.0D0
           DO 25 K=KM,0,-1
              QF0=((2*K+3.0D0)*X*QF1-(K+2.0D0)*QF2)/(K+1.0D0)
              IF (K.LE.N) QM(0,K)=QF0
              QF2=QF1
25            QF1=QF0
           DO 30 K=0,N
30            QM(0,K)=Q0*QM(0,K)/QF0
           QF2=0.0D0
           QF1=1.0D0
           DO 35 K=KM,0,-1
              QF0=((2*K+3.0D0)*X*QF1-(K+1.0D0)*QF2)/(K+2.0D0)
              IF (K.LE.N) QM(1,K)=QF0
              QF2=QF1
35            QF1=QF0
           Q10=-1.0D0/XQ
           DO 40 K=0,N
40            QM(1,K)=Q10*QM(1,K)/QF0
           DO 45 J=0,N
              Q0=QM(0,J)
              Q1=QM(1,J)
              DO 45 I=0,M-2
                 QF=-2.0D0*(I+1)*X/XQ*Q1+(J-I)*(J+I+1.0D0)*Q0
                 QM(I+2,J)=QF
                 Q0=Q1
                 Q1=QF
45         CONTINUE
        ENDIF
        QD(0,0)=LS/XS
        DO 50 J=1,N
50         QD(0,J)=LS*J*(QM(0,J-1)-X*QM(0,J))/XS
        DO 55 J=0,N
        DO 55 I=1,M
           QD(I,J)=LS*I*X/XS*QM(I,J)+(I+J)*(J-I+1.0D0)/XQ*QM(I-1,J)
55      CONTINUE
        RETURN
        END

C       **********************************

        SUBROUTINE CLPMN(MM,M,N,X,Y,NTYPE,CPM,CPD)
C
C       =========================================================
C       Purpose: Compute the associated Legendre functions Pmn(z)
C                and their derivatives Pmn'(z) for a complex
C                argument
C       Input :  x     --- Real part of z
C                y     --- Imaginary part of z
C                m     --- Order of Pmn(z),  m = 0,1,2,...,n
C                n     --- Degree of Pmn(z), n = 0,1,2,...,N
C                mm    --- Physical dimension of CPM and CPD
C                ntype --- type of cut, either 2 or 3
C       Output:  CPM(m,n) --- Pmn(z)
C                CPD(m,n) --- Pmn'(z)
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (D,X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CPM(0:MM,0:N),CPD(0:MM,0:N)
        Z = DCMPLX(X, Y)
        DO 10 I=0,N
        DO 10 J=0,M
           CPM(J,I)=(0.0D0,0.0D0)
10         CPD(J,I)=(0.0D0,0.0D0)
        CPM(0,0)=(1.0D0,0.0D0)
        IF (N.EQ.0) RETURN
        IF (DABS(X).EQ.1.0D0.AND.Y.EQ.0.0D0) THEN
           DO 15 I=1,N
              CPM(0,I)=X**I
15            CPD(0,I)=0.5D0*I*(I+1)*X**(I+1)
           DO 20 J=1,N
           DO 20 I=1,M
              IF (I.EQ.1) THEN
                 CPD(I,J)=DINF()
              ELSE IF (I.EQ.2) THEN
                 CPD(I,J)=-0.25D0*(J+2)*(J+1)*J*(J-1)*X**(J+1)
              ENDIF
20         CONTINUE
           RETURN
        ENDIF
        IF (NTYPE.EQ.2) THEN
C       sqrt(1 - z^2) with branch cut on |x|>1
           ZS=(1.0D0-Z*Z)
           ZQ=-CDSQRT(ZS)
           LS=-1
        ELSE
C       sqrt(z^2 - 1) with branch cut between [-1, 1]
           ZS=(Z*Z-1.0D0)
           ZQ=CDSQRT(ZS)
           IF (X.LT.0D0) THEN
              ZQ=-ZQ
           END IF
           LS=1
        END IF
        DO 25 I=1,M
C       DLMF 14.7.15
25         CPM(I,I)=(2.0D0*I-1.0D0)*ZQ*CPM(I-1,I-1)
        DO 30 I=0,MIN(M,N-1)
C       DLMF 14.10.7
30         CPM(I,I+1)=(2.0D0*I+1.0D0)*Z*CPM(I,I)
        DO 35 I=0,M
        DO 35 J=I+2,N
C       DLMF 14.10.3
           CPM(I,J)=((2.0D0*J-1.0D0)*Z*CPM(I,J-1)-(I+J-
     &              1.0D0)*CPM(I,J-2))/(J-I)
35      CONTINUE
        CPD(0,0)=(0.0D0,0.0D0)
        DO 40 J=1,N
C       DLMF 14.10.5
40         CPD(0,J)=LS*J*(Z*CPM(0,J)-CPM(0,J-1))/ZS
        DO 45 I=1,M
        DO 45 J=I,N
C       derivative of DLMF 14.7.11 & DLMF 14.10.6 for type 3
C       derivative of DLMF 14.7.8 & DLMF 14.10.1 for type 2
           CPD(I,J)=LS*(-I*Z*CPM(I,J)/ZS+(J+I)*(J-I+1.0D0)
     &                  /ZQ*CPM(I-1,J))
45      CONTINUE
        RETURN
        END

C       **********************************

        SUBROUTINE VVSA(VA,X,PV)
C
C       ===================================================
C       Purpose: Compute parabolic cylinder function Vv(x)
C                for small argument
C       Input:   x  --- Argument
C                va --- Order
C       Output:  PV --- Vv(x)
C       Routine called : GAMMA2 for computing Г(x)
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EPS=1.0D-15
        PI=3.141592653589793D0
        EP=DEXP(-.25D0*X*X)
        VA0=1.0D0+0.5D0*VA
        IF (X.EQ.0.0) THEN
           IF (VA0.LE.0.0.AND.VA0.EQ.INT(VA0).OR.VA.EQ.0.0) THEN
              PV=0.0D0
           ELSE
              VB0=-0.5D0*VA
              SV0=DSIN(VA0*PI)
              CALL GAMMA2(VA0,GA0)
              PV=2.0D0**VB0*SV0/GA0
           ENDIF
        ELSE
           SQ2=DSQRT(2.0D0)
           A0=2.0D0**(-.5D0*VA)*EP/(2.0D0*PI)
           SV=DSIN(-(VA+.5D0)*PI)
           V1=-.5D0*VA
           CALL GAMMA2(V1,G1)
           PV=(SV+1.0D0)*G1
           R=1.0D0
           FAC=1.0D0
           DO 10 M=1,250
              VM=.5D0*(M-VA)
              CALL GAMMA2(VM,GM)
              R=R*SQ2*X/M
              FAC=-FAC
              GW=FAC*SV+1.0D0
              R1=GW*R*GM
              PV=PV+R1
              IF (DABS(R1/PV).LT.EPS.AND.GW.NE.0.0) GO TO 15
10         CONTINUE
15         PV=A0*PV
        ENDIF
        RETURN
        END



C       **********************************
C       SciPy: Changed P from a character array to an integer array.
        SUBROUTINE JDZO(NT,N,M,P,ZO)
C
C       ===========================================================
C       Purpose: Compute the zeros of Bessel functions Jn(x) and
C                Jn'(x), and arrange them in the order of their
C                magnitudes
C       Input :  NT    --- Number of total zeros ( NT ≤ 1200 )
C       Output:  ZO(L) --- Value of the L-th zero of Jn(x)
C                          and Jn'(x)
C                N(L)  --- n, order of Jn(x) or Jn'(x) associated
C                          with the L-th zero
C                M(L)  --- m, serial number of the zeros of Jn(x)
C                          or Jn'(x) associated with the L-th zero
C                          ( L is the serial number of all the
C                            zeros of Jn(x) and Jn'(x) )
C                P(L)  --- 0 (TM) or 1 (TE), a code for designating the
C                          zeros of Jn(x)  or Jn'(x).
C                          In the waveguide applications, the zeros
C                          of Jn(x) correspond to TM modes and
C                          those of Jn'(x) correspond to TE modes
C       Routine called:    BJNDD for computing Jn(x), Jn'(x) and
C                          Jn''(x)
C       =============================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INTEGER P(1400), P1(70)
        DIMENSION N(1400),M(1400),ZO(0:1400),N1(70),M1(70),
     &            ZOC(0:70),BJ(101),DJ(101),FJ(101)
        X = 0
        ZOC(0) = 0
        IF (NT.LT.600) THEN
           XM=-1.0+2.248485*NT**0.5-.0159382*NT+3.208775E-4
     &        *NT**1.5
           NM=INT(14.5+.05875*NT)
           MM=INT(.02*NT)+6
        ELSE
           XM=5.0+1.445389*NT**.5+.01889876*NT-2.147763E-4
     &        *NT**1.5
           NM=INT(27.8+.0327*NT)
           MM=INT(.01088*NT)+10
        ENDIF
        L0=0
        DO 45 I=1,NM
           X1=.407658+.4795504*(I-1)**.5+.983618*(I-1)
           X2=1.99535+.8333883*(I-1)**.5+.984584*(I-1)
           L1=0
           DO 30 J=1,MM
              IF (I.EQ.1.AND.J.EQ.1) GO TO 15
              X=X1
10            CALL BJNDD(I,X,BJ,DJ,FJ)
              X0=X
              X=X-DJ(I)/FJ(I)
              IF (X1.GT.XM) GO TO 20
              IF (DABS(X-X0).GT.1.0D-10) GO TO 10
15            L1=L1+1
              N1(L1)=I-1
              M1(L1)=J
              IF (I.EQ.1) M1(L1)=J-1
              P1(L1)=1
              ZOC(L1)=X
              IF (I.LE.15) THEN
                 X1=X+3.057+.0122*(I-1)+(1.555+.41575*(I-1))/(J+1)**2
              ELSE
                 X1=X+2.918+.01924*(I-1)+(6.26+.13205*(I-1))/(J+1)**2
              ENDIF
20            X=X2
25            CALL BJNDD(I,X,BJ,DJ,FJ)
              X0=X
              X=X-BJ(I)/DJ(I)
              IF (X.GT.XM) GO TO 30
              IF (DABS(X-X0).GT.1.0D-10) GO TO 25
              L1=L1+1
              N1(L1)=I-1
              M1(L1)=J
              P1(L1)=0
              ZOC(L1)=X
              IF (I.LE.15) THEN
                 X2=X+3.11+.0138*(I-1)+(.04832+.2804*(I-1))/(J+1)**2
              ELSE
                 X2=X+3.001+.0105*(I-1)+(11.52+.48525*(I-1))/(J+3)**2
              ENDIF
30         CONTINUE
           L=L0+L1
           L2=L
35         IF (L0.EQ.0) THEN
              DO 40 K=1,L
                 ZO(K)=ZOC(K)
                 N(K)=N1(K)
                 M(K)=M1(K)
40               P(K)=P1(K)
              L1=0
           ELSE IF (L0.NE.0) THEN
              IF (ZO(L0).GE.ZOC(L1)) THEN
                 ZO(L0+L1)=ZO(L0)
                 N(L0+L1)=N(L0)
                 M(L0+L1)=M(L0)
                 P(L0+L1)=P(L0)
                 L0=L0-1
              ELSE
                 ZO(L0+L1)=ZOC(L1)
                 N(L0+L1)=N1(L1)
                 M(L0+L1)=M1(L1)
                 P(L0+L1)=P1(L1)
                 L1=L1-1
              ENDIF
           ENDIF
           IF (L1.NE.0) GO TO 35
45         L0=L2
        RETURN
        END



C       **********************************

        SUBROUTINE CBK(M,N,C,CV,QT,CK,BK)
C
C       =====================================================
C       Purpose: Compute coefficient Bk's for oblate radial
C                functions with a small argument
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BK(200),CK(200),U(200),V(200),W(200)
        EPS=1.0D-14
        IP=1
        IF (N-M.EQ.2*INT((N-M)/2)) IP=0
        NM=25+INT(0.5*(N-M)+C)
        U(1)=0.0D0
        N2=NM-2
        DO 10 J=2,N2
10         U(J)=C*C
        DO 15 J=1,N2
15         V(J)=(2.0*J-1.0-IP)*(2.0*(J-M)-IP)+M*(M-1.0)-CV
        DO 20 J=1,NM-1
20         W(J)=(2.0*J-IP)*(2.0*J+1.0-IP)
        IF (IP.EQ.0) THEN
           SW=0.0D0
           DO 40 K=0,N2-1
              S1=0.0D0
              I1=K-M+1
              DO 30 I=I1,NM
                 IF (I.LT.0) GO TO 30
                 R1=1.0D0
                 DO 25 J=1,K
25                  R1=R1*(I+M-J)/J
                 S1=S1+CK(I+1)*(2.0*I+M)*R1
                 IF (DABS(S1-SW).LT.DABS(S1)*EPS) GO TO 35
                 SW=S1
30            CONTINUE
35            BK(K+1)=QT*S1
40         CONTINUE
        ELSE IF (IP.EQ.1) THEN
           SW=0.0D0
           DO 60 K=0,N2-1
              S1=0.0D0
              I1=K-M+1
              DO 50 I=I1,NM
                 IF (I.LT.0) GO TO 50
                 R1=1.0D0
                 DO 45 J=1,K
45                  R1=R1*(I+M-J)/J
                 IF (I.GT.0) S1=S1+CK(I)*(2.0*I+M-1)*R1
                 S1=S1-CK(I+1)*(2.0*I+M)*R1
                 IF (DABS(S1-SW).LT.DABS(S1)*EPS) GO TO 55
                 SW=S1
50            CONTINUE
55            BK(K+1)=QT*S1
60         CONTINUE
        ENDIF
        W(1)=W(1)/V(1)
        BK(1)=BK(1)/V(1)
        DO 65 K=2,N2
           T=V(K)-W(K-1)*U(K)
           W(K)=W(K)/T
65         BK(K)=(BK(K)-BK(K-1)*U(K))/T
        DO 70 K=N2-1,1,-1
70         BK(K)=BK(K)-W(K)*BK(K+1)
        RETURN
        END



C       **********************************

        SUBROUTINE CJY01(Z,CBJ0,CDJ0,CBJ1,CDJ1,CBY0,CDY0,CBY1,CDY1)
C
C       =======================================================
C       Purpose: Compute Bessel functions J0(z), J1(z), Y0(z),
C                Y1(z), and their derivatives for a complex
C                argument
C       Input :  z --- Complex argument
C       Output:  CBJ0 --- J0(z)
C                CDJ0 --- J0'(z)
C                CBJ1 --- J1(z)
C                CDJ1 --- J1'(z)
C                CBY0 --- Y0(z)
C                CDY0 --- Y0'(z)
C                CBY1 --- Y1(z)
C                CDY1 --- Y1'(z)
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,E,P,R,W)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION A(12),B(12),A1(12),B1(12)
        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        RP2=2.0D0/PI
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z2=Z*Z
        Z1=Z
        IF (A0.EQ.0.0D0) THEN
           CBJ0=(1.0D0,0.0D0)
           CBJ1=(0.0D0,0.0D0)
           CDJ0=(0.0D0,0.0D0)
           CDJ1=(0.5D0,0.0D0)
           CBY0=-(1.0D300,0.0D0)
           CBY1=-(1.0D300,0.0D0)
           CDY0=(1.0D300,0.0D0)
           CDY1=(1.0D300,0.0D0)
           RETURN
        ENDIF
        IF (DBLE(Z).LT.0.0) Z1=-Z
        IF (A0.LE.12.0) THEN
           CBJ0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 10 K=1,40
              CR=-0.25D0*CR*Z2/(K*K)
              CBJ0=CBJ0+CR
              IF (CDABS(CR).LT.CDABS(CBJ0)*1.0D-15) GO TO 15
10         CONTINUE
15         CBJ1=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 20 K=1,40
              CR=-0.25D0*CR*Z2/(K*(K+1.0D0))
              CBJ1=CBJ1+CR
              IF (CDABS(CR).LT.CDABS(CBJ1)*1.0D-15) GO TO 25
20         CONTINUE
25         CBJ1=0.5D0*Z1*CBJ1
           W0=0.0D0
           CR=(1.0D0,0.0D0)
           CS=(0.0D0,0.0D0)
           DO 30 K=1,40
              W0=W0+1.0D0/K
              CR=-0.25D0*CR/(K*K)*Z2
              CP=CR*W0
              CS=CS+CP
              IF (CDABS(CP).LT.CDABS(CS)*1.0D-15) GO TO 35
30         CONTINUE
35         CBY0=RP2*(CDLOG(Z1/2.0D0)+EL)*CBJ0-RP2*CS
           W1=0.0D0
           CR=(1.0D0,0.0D0)
           CS=(1.0D0,0.0D0)
           DO 40 K=1,40
              W1=W1+1.0D0/K
              CR=-0.25D0*CR/(K*(K+1))*Z2
              CP=CR*(2.0D0*W1+1.0D0/(K+1.0D0))
              CS=CS+CP
              IF (CDABS(CP).LT.CDABS(CS)*1.0D-15) GO TO 45
40         CONTINUE
45         CBY1=RP2*((CDLOG(Z1/2.0D0)+EL)*CBJ1-1.0D0/Z1-.25D0*Z1*CS)
        ELSE
           DATA A/-.703125D-01,.112152099609375D+00,
     &            -.5725014209747314D+00,.6074042001273483D+01,
     &            -.1100171402692467D+03,.3038090510922384D+04,
     &            -.1188384262567832D+06,.6252951493434797D+07,
     &            -.4259392165047669D+09,.3646840080706556D+11,
     &            -.3833534661393944D+13,.4854014686852901D+15/
           DATA B/ .732421875D-01,-.2271080017089844D+00,
     &             .1727727502584457D+01,-.2438052969955606D+02,
     &             .5513358961220206D+03,-.1825775547429318D+05,
     &             .8328593040162893D+06,-.5006958953198893D+08,
     &             .3836255180230433D+10,-.3649010818849833D+12,
     &             .4218971570284096D+14,-.5827244631566907D+16/
           DATA A1/.1171875D+00,-.144195556640625D+00,
     &             .6765925884246826D+00,-.6883914268109947D+01,
     &             .1215978918765359D+03,-.3302272294480852D+04,
     &             .1276412726461746D+06,-.6656367718817688D+07,
     &             .4502786003050393D+09,-.3833857520742790D+11,
     &             .4011838599133198D+13,-.5060568503314727D+15/
           DATA B1/-.1025390625D+00,.2775764465332031D+00,
     &             -.1993531733751297D+01,.2724882731126854D+02,
     &             -.6038440767050702D+03,.1971837591223663D+05,
     &             -.8902978767070678D+06,.5310411010968522D+08,
     &             -.4043620325107754D+10,.3827011346598605D+12,
     &             -.4406481417852278D+14,.6065091351222699D+16/
           K0=12
           IF (A0.GE.35.0) K0=10
           IF (A0.GE.50.0) K0=8
           CT1=Z1-.25D0*PI
           CP0=(1.0D0,0.0D0)
           DO 50 K=1,K0
50            CP0=CP0+A(K)*Z1**(-2*K)
           CQ0=-0.125D0/Z1
           DO 55 K=1,K0
55            CQ0=CQ0+B(K)*Z1**(-2*K-1)
           CU=CDSQRT(RP2/Z1)
           CBJ0=CU*(CP0*CDCOS(CT1)-CQ0*CDSIN(CT1))
           CBY0=CU*(CP0*CDSIN(CT1)+CQ0*CDCOS(CT1))
           CT2=Z1-.75D0*PI
           CP1=(1.0D0,0.0D0)
           DO 60 K=1,K0
60            CP1=CP1+A1(K)*Z1**(-2*K)
           CQ1=0.375D0/Z1
           DO 65 K=1,K0
65            CQ1=CQ1+B1(K)*Z1**(-2*K-1)
           CBJ1=CU*(CP1*CDCOS(CT2)-CQ1*CDSIN(CT2))
           CBY1=CU*(CP1*CDSIN(CT2)+CQ1*CDCOS(CT2))
        ENDIF
        IF (DBLE(Z).LT.0.0) THEN
           IF (DIMAG(Z).LT.0.0) CBY0=CBY0-2.0D0*CI*CBJ0
           IF (DIMAG(Z).GT.0.0) CBY0=CBY0+2.0D0*CI*CBJ0
           IF (DIMAG(Z).LT.0.0) CBY1=-(CBY1-2.0D0*CI*CBJ1)
           IF (DIMAG(Z).GT.0.0) CBY1=-(CBY1+2.0D0*CI*CBJ1)
           CBJ1=-CBJ1
        ENDIF
        CDJ0=-CBJ1
        CDJ1=CBJ0-1.0D0/Z*CBJ1
        CDY0=-CBY1
        CDY1=CBY0-1.0D0/Z*CBY1
        RETURN
        END

C       **********************************

        SUBROUTINE RMN2SP(M,N,C,X,CV,DF,KD,R2F,R2D)
C
C       ======================================================
C       Purpose: Compute prolate spheroidal radial function
C                of the second kind with a small argument
C       Routines called:
C            (1) LPMNS for computing the associated Legendre
C                functions of the first kind
C            (2) LQMNS for computing the associated Legendre
C                functions of the second kind
C            (3) KMN for computing expansion coefficients
C                and joining factors
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION PM(0:251),PD(0:251),QM(0:251),QD(0:251),
     &            DN(200),DF(200)
        IF (DABS(DF(1)).LT.1.0D-280) THEN
           R2F=1.0D+300
           R2D=1.0D+300
           RETURN
        ENDIF
        EPS=1.0D-14
        IP=1
        NM1=INT((N-M)/2)
        IF (N-M.EQ.2*NM1) IP=0
        NM=25+NM1+INT(C)
        NM2=2*NM+M
        CALL KMN(M,N,C,CV,KD,DF,DN,CK1,CK2)
        CALL LPMNS(M,NM2,X,PM,PD)
        CALL LQMNS(M,NM2,X,QM,QD)
        SU0=0.0D0
        SW=0.0D0
        DO 10 K=1,NM
          J=2*K-2+M+IP
          SU0=SU0+DF(K)*QM(J)
          IF (K.GT.NM1.AND.DABS(SU0-SW).LT.DABS(SU0)*EPS) GO TO 15
10        SW=SU0
15      SD0=0.0D0
        DO 20 K=1,NM
          J=2*K-2+M+IP
          SD0=SD0+DF(K)*QD(J)
          IF (K.GT.NM1.AND.DABS(SD0-SW).LT.DABS(SD0)*EPS) GO TO 25
20        SW=SD0
25        SU1=0.0D0
          SD1=0.0D0
          DO 30 K=1,M
             J=M-2*K+IP
             IF (J.LT.0) J=-J-1
             SU1=SU1+DN(K)*QM(J)
30           SD1=SD1+DN(K)*QD(J)
          GA=((X-1.0D0)/(X+1.0D0))**(0.5D0*M)
          DO 55 K=1,M
             J=M-2*K+IP
             IF (J.GE.0) GO TO 55
             IF (J.LT.0) J=-J-1
             R1=1.0D0
             DO 35 J1=1,J
35              R1=(M+J1)*R1
             R2=1.0D0
             DO 40 J2=1,M-J-2
40              R2=J2*R2
             R3=1.0D0
             SF=1.0D0
             DO 45 L1=1,J
                R3=0.5D0*R3*(-J+L1-1.0)*(J+L1)/((M+L1)*L1)*(1.0-X)
45              SF=SF+R3
             IF (M-J.GE.2) GB=(M-J-1.0D0)*R2
             IF (M-J.LE.1) GB=1.0D0
             SPL=R1*GA*GB*SF
             SU1=SU1+(-1)**(J+M)*DN(K)*SPL
             SPD1=M/(X*X-1.0D0)*SPL
             GC=0.5D0*J*(J+1.0)/(M+1.0)
             SD=1.0D0
             R4=1.0D0
             DO 50 L1=1,J-1
                R4=0.5D0*R4*(-J+L1)*(J+L1+1.0)/((M+L1+1.0)*L1)
     &             *(1.0-X)
50              SD=SD+R4
             SPD2=R1*GA*GB*GC*SD
             SD1=SD1+(-1)**(J+M)*DN(K)*(SPD1+SPD2)
55        CONTINUE
          SU2=0.0D0
          KI=(2*M+1+IP)/2
          NM3=NM+KI
          DO 60 K=KI,NM3
             J=2*K-1-M-IP
             SU2=SU2+DN(K)*PM(J)
             IF (J.GT.M.AND.DABS(SU2-SW).LT.DABS(SU2)*EPS) GO TO 65
60           SW=SU2
65        SD2=0.0D0
          DO 70 K=KI,NM3
             J=2*K-1-M-IP
             SD2=SD2+DN(K)*PD(J)
             IF (J.GT.M.AND.DABS(SD2-SW).LT.DABS(SD2)*EPS) GO TO 75
70           SW=SD2
75      SUM=SU0+SU1+SU2
        SDM=SD0+SD1+SD2
        R2F=SUM/CK2
        R2D=SDM/CK2
        RETURN
        END



C       **********************************

        SUBROUTINE BERNOB(N,BN)
C
C       ======================================
C       Purpose: Compute Bernoulli number Bn
C       Input :  n --- Serial number
C       Output:  BN(n) --- Bn
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BN(0:N)
        TPI=6.283185307179586D0
        BN(0)=1.0D0
        BN(1)=-0.5D0
        BN(2)=1.0D0/6.0D0
        R1=(2.0D0/TPI)**2
        DO 20 M=4,N,2
           R1=-R1*(M-1)*M/(TPI*TPI)
           R2=1.0D0
           DO 10 K=2,10000
              S=(1.0D0/K)**M
              R2=R2+S
              IF (S.LT.1.0D-15) GOTO 20
10         CONTINUE
20         BN(M)=R1*R2
        RETURN
        END

C       **********************************

        SUBROUTINE BERNOA(N,BN)
C
C       ======================================
C       Purpose: Compute Bernoulli number Bn
C       Input :  n --- Serial number
C       Output:  BN(n) --- Bn
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BN(0:N)
        BN(0)=1.0D0
        BN(1)=-0.5D0
        DO 30 M=2,N
           S=-(1.0D0/(M+1.0D0)-0.5D0)
           DO 20 K=2,M-1
              R=1.0D0
              DO 10 J=2,K
10               R=R*(J+M-K)/J
20            S=S-R*BN(K)
30         BN(M)=S
        DO 40 M=3,N,2
40         BN(M)=0.0D0
        RETURN
        END

C       **********************************

        SUBROUTINE QSTAR(M,N,C,CK,CK1,QS,QT)
C
C       =========================================================
C       Purpose: Compute Q*mn(-ic) for oblate radial functions
C                with a small argument
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION AP(200),CK(200)
        IP=1
        IF (N-M.EQ.2*INT((N-M)/2)) IP=0
        R=1.0D0/CK(1)**2
        AP(1)=R
        DO 20 I=1,M
           S=0.0D0
           DO 15 L=1,I
              SK=0.0D0
              DO 10 K=0,L
10               SK=SK+CK(K+1)*CK(L-K+1)
15            S=S+SK*AP(I-L+1)
20      AP(I+1)=-R*S
        QS0=AP(M+1)
        DO 30 L=1,M
           R=1.0D0
           DO 25 K=1,L
25            R=R*(2.0D0*K+IP)*(2.0D0*K-1.0D0+IP)/(2.0D0*K)**2
30         QS0=QS0+AP(M-L+1)*R
        QS=(-1)**IP*CK1*(CK1*QS0)/C
        QT=-2.0D0/CK1*QS
        RETURN
        END



C       **********************************

        SUBROUTINE CV0(KD,M,Q,A0)
C
C       =====================================================
C       Purpose: Compute the initial characteristic value of
C                Mathieu functions for m ≤ 12  or q ≤ 300 or
C                q ≥ m*m
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C       Output:  A0 --- Characteristic value
C       Routines called:
C             (1) CVQM for computing initial characteristic
C                 value for q ≤ 3*m
C             (2) CVQL for computing initial characteristic
C                 value for q ≥ m*m
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        Q2=Q*Q
        IF (M.EQ.0) THEN
           IF (Q.LE.1.0) THEN
              A0=(((.0036392*Q2-.0125868)*Q2+.0546875)*Q2-.5)*Q2
           ELSE IF (Q.LE.10.0) THEN
              A0=((3.999267D-3*Q-9.638957D-2)*Q-.88297)*Q
     &           +.5542818
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.1) THEN
           IF (Q.LE.1.0.AND.KD.EQ.2) THEN
              A0=(((-6.51E-4*Q-.015625)*Q-.125)*Q+1.0)*Q+1.0
           ELSE IF (Q.LE.1.0.AND.KD.EQ.3) THEN
              A0=(((-6.51E-4*Q+.015625)*Q-.125)*Q-1.0)*Q+1.0
           ELSE IF (Q.LE.10.0.AND. KD.EQ.2) THEN
              A0=(((-4.94603D-4*Q+1.92917D-2)*Q-.3089229)
     &           *Q+1.33372)*Q+.811752
           ELSE IF (Q.LE.10.0.AND.KD.EQ.3) THEN
              A0=((1.971096D-3*Q-5.482465D-2)*Q-1.152218)
     &           *Q+1.10427
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.2) THEN
           IF (Q.LE.1.0.AND.KD.EQ.1) THEN
              A0=(((-.0036391*Q2+.0125888)*Q2-.0551939)*Q2
     &           +.416667)*Q2+4.0
           ELSE IF (Q.LE.1.0.AND.KD.EQ.4) THEN
              A0=(.0003617*Q2-.0833333)*Q2+4.0
           ELSE IF (Q.LE.15.AND.KD.EQ.1) THEN
              A0=(((3.200972D-4*Q-8.667445D-3)*Q
     &           -1.829032D-4)*Q+.9919999)*Q+3.3290504
           ELSE IF (Q.LE.10.0.AND.KD.EQ.4) THEN
              A0=((2.38446D-3*Q-.08725329)*Q-4.732542D-3)
     &           *Q+4.00909
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.3) THEN
           IF (Q.LE.1.0.AND.KD.EQ.2) THEN
              A0=((6.348E-4*Q+.015625)*Q+.0625)*Q2+9.0
           ELSE IF (Q.LE.1.0.AND.KD.EQ.3) THEN
              A0=((6.348E-4*Q-.015625)*Q+.0625)*Q2+9.0
           ELSE IF (Q.LE.20.0.AND.KD.EQ.2) THEN
              A0=(((3.035731D-4*Q-1.453021D-2)*Q
     &           +.19069602)*Q-.1039356)*Q+8.9449274
           ELSE IF (Q.LE.15.0.AND.KD.EQ.3) THEN
              A0=((9.369364D-5*Q-.03569325)*Q+.2689874)*Q
     &           +8.771735
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.4) THEN
           IF (Q.LE.1.0.AND.KD.EQ.1) THEN
              A0=((-2.1E-6*Q2+5.012E-4)*Q2+.0333333)*Q2+16.0
           ELSE IF (Q.LE.1.0.AND.KD.EQ.4) THEN
              A0=((3.7E-6*Q2-3.669E-4)*Q2+.0333333)*Q2+16.0
           ELSE IF (Q.LE.25.0.AND.KD.EQ.1) THEN
              A0=(((1.076676D-4*Q-7.9684875D-3)*Q
     &           +.17344854)*Q-.5924058)*Q+16.620847
           ELSE IF (Q.LE.20.0.AND.KD.EQ.4) THEN
              A0=((-7.08719D-4*Q+3.8216144D-3)*Q
     &           +.1907493)*Q+15.744
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.5) THEN
           IF (Q.LE.1.0.AND.KD.EQ.2) THEN
              A0=((6.8E-6*Q+1.42E-5)*Q2+.0208333)*Q2+25.0
           ELSE IF (Q.LE.1.0.AND.KD.EQ.3) THEN
              A0=((-6.8E-6*Q+1.42E-5)*Q2+.0208333)*Q2+25.0
           ELSE IF (Q.LE.35.0.AND.KD.EQ.2) THEN
              A0=(((2.238231D-5*Q-2.983416D-3)*Q
     &           +.10706975)*Q-.600205)*Q+25.93515
           ELSE IF (Q.LE.25.0.AND.KD.EQ.3) THEN
              A0=((-7.425364D-4*Q+2.18225D-2)*Q
     &           +4.16399D-2)*Q+24.897
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.6) THEN
           IF (Q.LE.1.0) THEN
              A0=(.4D-6*Q2+.0142857)*Q2+36.0
           ELSE IF (Q.LE.40.0.AND.KD.EQ.1) THEN
              A0=(((-1.66846D-5*Q+4.80263D-4)*Q
     &           +2.53998D-2)*Q-.181233)*Q+36.423
           ELSE IF (Q.LE.35.0.AND.KD.EQ.4) THEN
              A0=((-4.57146D-4*Q+2.16609D-2)*Q-2.349616D-2)*Q
     &           +35.99251
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.EQ.7) THEN
           IF (Q.LE.10.0) THEN
              CALL CVQM(M,Q,A0)
           ELSE IF (Q.LE.50.0.AND.KD.EQ.2) THEN
              A0=(((-1.411114D-5*Q+9.730514D-4)*Q
     &           -3.097887D-3)*Q+3.533597D-2)*Q+49.0547
           ELSE IF (Q.LE.40.0.AND.KD.EQ.3) THEN
              A0=((-3.043872D-4*Q+2.05511D-2)*Q
     &           -9.16292D-2)*Q+49.19035
           ELSE
              CALL CVQL(KD,M,Q,A0)
           ENDIF
        ELSE IF (M.GE.8) THEN
           IF (Q.LE.3.*M) THEN
              CALL CVQM(M,Q,A0)
           ELSE IF (Q.GT.M*M) THEN
              CALL CVQL(KD,M,Q,A0)
           ELSE
              IF (M.EQ.8.AND.KD.EQ.1) THEN
                 A0=(((8.634308D-6*Q-2.100289D-3)*Q+.169072)*Q
     &              -4.64336)*Q+109.4211
              ELSE IF (M.EQ.8.AND.KD.EQ.4) THEN
                 A0=((-6.7842D-5*Q+2.2057D-3)*Q+.48296)*Q+56.59
              ELSE IF (M.EQ.9.AND.KD.EQ.2) THEN
                 A0=(((2.906435D-6*Q-1.019893D-3)*Q+.1101965)*Q
     &              -3.821851)*Q+127.6098
              ELSE IF (M.EQ.9.AND.KD.EQ.3) THEN
                 A0=((-9.577289D-5*Q+.01043839)*Q+.06588934)*Q
     &              +78.0198
              ELSE IF (M.EQ.10.AND.KD.EQ.1) THEN
                 A0=(((5.44927D-7*Q-3.926119D-4)*Q+.0612099)*Q
     &              -2.600805)*Q+138.1923
              ELSE IF (M.EQ.10.AND.KD.EQ.4) THEN
                 A0=((-7.660143D-5*Q+.01132506)*Q-.09746023)*Q
     &              +99.29494
              ELSE IF (M.EQ.11.AND.KD.EQ.2) THEN
                 A0=(((-5.67615D-7*Q+7.152722D-6)*Q+.01920291)*Q
     &              -1.081583)*Q+140.88
              ELSE IF (M.EQ.11.AND.KD.EQ.3) THEN
                 A0=((-6.310551D-5*Q+.0119247)*Q-.2681195)*Q
     &              +123.667
              ELSE IF (M.EQ.12.AND.KD.EQ.1) THEN
                 A0=(((-2.38351D-7*Q-2.90139D-5)*Q+.02023088)*Q
     &              -1.289)*Q+171.2723
              ELSE IF (M.EQ.12.AND.KD.EQ.4) THEN
                 A0=(((3.08902D-7*Q-1.577869D-4)*Q+.0247911)*Q
     &              -1.05454)*Q+161.471
              ENDIF
           ENDIF
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE CVQM(M,Q,A0)
C
C       =====================================================
C       Purpose: Compute the characteristic value of Mathieu
C                functions for q ≤ m*m
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C       Output:  A0 --- Initial characteristic value
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        HM1=.5*Q/(M*M-1.0)
        HM3=.25*HM1**3/(M*M-4.0)
        HM5=HM1*HM3*Q/((M*M-1.0)*(M*M-9.0))
        A0=M*M+Q*(HM1+(5.0*M*M+7.0)*HM3
     &     +(9.0*M**4+58.0*M*M+29.0)*HM5)
        RETURN
        END

C       **********************************

        SUBROUTINE CVQL(KD,M,Q,A0)
C
C       ========================================================
C       Purpose: Compute the characteristic value of Mathieu
C                functions  for q ≥ 3m
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C       Output:  A0 --- Initial characteristic value
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        W=0.0D0
        IF (KD.EQ.1.OR.KD.EQ.2) W=2.0D0*M+1.0D0
        IF (KD.EQ.3.OR.KD.EQ.4) W=2.0D0*M-1.0D0
        W2=W*W
        W3=W*W2
        W4=W2*W2
        W6=W2*W4
        D1=5.0+34.0/W2+9.0/W4
        D2=(33.0+410.0/W2+405.0/W4)/W
        D3=(63.0+1260.0/W2+2943.0/W4+486.0/W6)/W2
        D4=(527.0+15617.0/W2+69001.0/W4+41607.0/W6)/W3
        C1=128.0
        P2=Q/W4
        P1=DSQRT(P2)
        CV1=-2.0*Q+2.0*W*DSQRT(Q)-(W2+1.0)/8.0
        CV2=(W+3.0/W)+D1/(32.0*P1)+D2/(8.0*C1*P2)
        CV2=CV2+D3/(64.0*C1*P1*P2)+D4/(16.0*C1*C1*P2*P2)
        A0=CV1-CV2/(C1*P1)
        RETURN
        END



C       **********************************

        SUBROUTINE CSPHJY(N,Z,NM,CSJ,CDJ,CSY,CDY)
C
C       ==========================================================
C       Purpose: Compute spherical Bessel functions jn(z) & yn(z)
C                and their derivatives for a complex argument
C       Input :  z --- Complex argument
C                n --- Order of jn(z) & yn(z) ( n = 0,1,2,... )
C       Output:  CSJ(n) --- jn(z)
C                CDJ(n) --- jn'(z)
C                CSY(n) --- yn(z)
C                CDY(n) --- yn'(z)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       ==========================================================
C
        IMPLICIT COMPLEX*16 (C,Z)
        DOUBLE PRECISION A0
        DIMENSION CSJ(0:N),CDJ(0:N),CSY(0:N),CDY(0:N)
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-60) THEN
           DO 10 K=0,N
              CSJ(K)=0.0D0
              CDJ(K)=0.0D0
              CSY(K)=-1.0D+300
10            CDY(K)=1.0D+300
           CSJ(0)=(1.0D0,0.0D0)
           IF (N.GT.0) THEN
              CDJ(1)=(.333333333333333D0,0.0D0)
           ENDIF
           RETURN
        ENDIF
        CSJ(0)=CDSIN(Z)/Z
        CDJ(0)=(CDCOS(Z)-CDSIN(Z)/Z)/Z
        CSY(0)=-CDCOS(Z)/Z
        CDY(0)=(CDSIN(Z)+CDCOS(Z)/Z)/Z
        IF (N.LT.1) THEN
           RETURN
        ENDIF
        CSJ(1)=(CSJ(0)-CDCOS(Z))/Z
        IF (N.GE.2) THEN
           CSA=CSJ(0)
           CSB=CSJ(1)
           M=MSTA1(A0,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(A0,N,15)
           ENDIF
           CF0=0.0D0
           CF1=1.0D0-100
           DO 15 K=M,0,-1
              CF=(2.0D0*K+3.0D0)*CF1/Z-CF0
              IF (K.LE.NM) CSJ(K)=CF
              CF0=CF1
15            CF1=CF
           IF (CDABS(CSA).GT.CDABS(CSB)) CS=CSA/CF1
           IF (CDABS(CSA).LE.CDABS(CSB)) CS=CSB/CF0
           DO 20 K=0,NM
20            CSJ(K)=CS*CSJ(K)
        ENDIF
        DO 25 K=1,NM
25         CDJ(K)=CSJ(K-1)-(K+1.0D0)*CSJ(K)/Z
        CSY(1)=(CSY(0)-CDSIN(Z))/Z
        CDY(1)=(2.0D0*CDY(0)-CDCOS(Z))/Z
        DO 30 K=2,NM
           IF (CDABS(CSJ(K-1)).GT.CDABS(CSJ(K-2))) THEN
              CSY(K)=(CSJ(K)*CSY(K-1)-1.0D0/(Z*Z))/CSJ(K-1)
           ELSE
              CSY(K)=(CSJ(K)*CSY(K-2)-(2.0D0*K-1.0D0)/Z**3)/CSJ(K-2)
           ENDIF
30      CONTINUE
        DO 35 K=2,NM
35         CDY(K)=CSY(K-1)-(K+1.0D0)*CSY(K)/Z
        RETURN
        END


        INTEGER FUNCTION MSTA1(X,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that the magnitude of
C                Jn(x) at that point is about 10^(-MP)
C       Input :  x     --- Argument of Jn(x)
C                MP    --- Value of magnitude
C       Output:  MSTA1 --- Starting point
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        N0=INT(1.1D0*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        END


        INTEGER FUNCTION MSTA2(X,N,MP)
C
C       ===================================================
C       Purpose: Determine the starting point for backward
C                recurrence such that all Jn(x) has MP
C                significant digits
C       Input :  x  --- Argument of Jn(x)
C                n  --- Order of Jn(x)
C                MP --- Significant digit
C       Output:  MSTA2 --- Starting point
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)+1
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        END

        REAL*8 FUNCTION ENVJ(N,X)
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END

C       **********************************

        SUBROUTINE ITTJYB(X,TTJ,TTY)
C
C       ==========================================================
C       Purpose: Integrate [1-J0(t)]/t with respect to t from 0
C                to x, and Y0(t)/t with respect to t from x to ∞
C       Input :  x   --- Variable in the limits  ( x ≥ 0 )
C       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
C                TTY --- Integration of Y0(t)/t from x to ∞
C       ==========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        IF (X.EQ.0.0D0) THEN
           TTJ=0.0D0
           TTY=-1.0D+300
        ELSE IF (X.LE.4.0D0) THEN
           X1=X/4.0D0
           T=X1*X1
           TTJ=((((((.35817D-4*T-.639765D-3)*T+.7092535D-2)*T
     &         -.055544803D0)*T+.296292677D0)*T-.999999326D0)
     &         *T+1.999999936D0)*T
           TTY=(((((((-.3546D-5*T+.76217D-4)*T-.1059499D-2)*T
     &         +.010787555D0)*T-.07810271D0)*T+.377255736D0)
     &         *T-1.114084491D0)*T+1.909859297D0)*T
           E0=EL+DLOG(X/2.0D0)
           TTY=PI/6.0D0+E0/PI*(2.0D0*TTJ-E0)-TTY
        ELSE IF (X.LE.8.0D0) THEN
           XT=X+.25D0*PI
           T1=4.0D0/X
           T=T1*T1
           F0=(((((.0145369D0*T-.0666297D0)*T+.1341551D0)*T
     &        -.1647797D0)*T+.1608874D0)*T-.2021547D0)*T
     &        +.7977506D0
           G0=((((((.0160672D0*T-.0759339D0)*T+.1576116D0)*T
     &        -.1960154D0)*T+.1797457D0)*T-.1702778D0)*T
     &        +.3235819D0)*T1
           TTJ=(F0*DCOS(XT)+G0*DSIN(XT))/(DSQRT(X)*X)
           TTJ=TTJ+EL+DLOG(X/2.0D0)
           TTY=(F0*DSIN(XT)-G0*DCOS(XT))/(DSQRT(X)*X)
        ELSE
           T=8.0D0/X
           XT=X+.25D0*PI
           F0=(((((.18118D-2*T-.91909D-2)*T+.017033D0)*T
     &        -.9394D-3)*T-.051445D0)*T-.11D-5)*T+.7978846D0
           G0=(((((-.23731D-2*T+.59842D-2)*T+.24437D-2)*T
     &      -.0233178D0)*T+.595D-4)*T+.1620695D0)*T
           TTJ=(F0*DCOS(XT)+G0*DSIN(XT))/(DSQRT(X)*X)
     &         +EL+DLOG(X/2.0D0)
           TTY=(F0*DSIN(XT)-G0*DCOS(XT))/(DSQRT(X)*X)
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE ITTJYA(X,TTJ,TTY)
C
C       =========================================================
C       Purpose: Integrate [1-J0(t)]/t with respect to t from 0
C                to x, and Y0(t)/t with respect to t from x to ∞
C       Input :  x   --- Variable in the limits  ( x ≥ 0 )
C       Output:  TTJ --- Integration of [1-J0(t)]/t from 0 to x
C                TTY --- Integration of Y0(t)/t from x to ∞
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        IF (X.EQ.0.0D0) THEN
           TTJ=0.0D0
           TTY=-1.0D+300
        ELSE IF (X.LE.20.0D0) THEN
           TTJ=1.0D0
           R=1.0D0
           DO 10 K=2,100
              R=-.25D0*R*(K-1.0D0)/(K*K*K)*X*X
              TTJ=TTJ+R
              IF (DABS(R).LT.DABS(TTJ)*1.0D-12) GO TO 15
10         CONTINUE
15         TTJ=TTJ*.125D0*X*X
           E0=.5D0*(PI*PI/6.0D0-EL*EL)-(.5D0*DLOG(X/2.0D0)+EL)
     &        *DLOG(X/2.0D0)
           B1=EL+DLOG(X/2.0D0)-1.5D0
           RS=1.0D0
           R=-1.0D0
           DO 20 K=2,100
              R=-.25D0*R*(K-1.0D0)/(K*K*K)*X*X
              RS=RS+1.0D0/K
              R2=R*(RS+1.0D0/(2.0D0*K)-(EL+DLOG(X/2.0D0)))
              B1=B1+R2
              IF (DABS(R2).LT.DABS(B1)*1.0D-12) GO TO 25
20         CONTINUE
25         TTY=2.0D0/PI*(E0+.125D0*X*X*B1)
        ELSE
           A0=DSQRT(2.0D0/(PI*X))
           BJ0=0.0D0
           BY0=0.0D0
           BJ1=0.0D0
           DO 50 L=0,1
              VT=4.0D0*L*L
              PX=1.0D0
              R=1.0D0
              DO 30 K=1,14
                 R=-.0078125D0*R*(VT-(4.0D0*K-3.0D0)**2)
     &             /(X*K)*(VT-(4.0D0*K-1.0D0)**2)
     &             /((2.0D0*K-1.0D0)*X)
                 PX=PX+R
                 IF (DABS(R).LT.DABS(PX)*1.0D-12) GO TO 35
30            CONTINUE
35            QX=1.0D0
              R=1.0D0
              DO 40 K=1,14
                 R=-.0078125D0*R*(VT-(4.0D0*K-1.0D0)**2)
     &             /(X*K)*(VT-(4.0D0*K+1.0D0)**2)
     &             /(2.0D0*K+1.0D0)/X
                 QX=QX+R
                 IF (DABS(R).LT.DABS(QX)*1.0D-12) GO TO 45
40            CONTINUE
45            QX=.125D0*(VT-1.0D0)/X*QX
              XK=X-(.25D0+.5D0*L)*PI
              BJ1=A0*(PX*DCOS(XK)-QX*DSIN(XK))
              BY1=A0*(PX*DSIN(XK)+QX*DCOS(XK))
              IF (L.EQ.0) THEN
                 BJ0=BJ1
                 BY0=BY1
              ENDIF
50         CONTINUE
           T=2.0D0/X
           G0=1.0D0
           R0=1.0D0
           DO 55 K=1,10
              R0=-K*K*T*T*R0
55            G0=G0+R0
           G1=1.0D0
           R1=1.0D0
           DO 60 K=1,10
              R1=-K*(K+1.0D0)*T*T*R1
60            G1=G1+R1
           TTJ=2.0D0*G1*BJ0/(X*X)-G0*BJ1/X+EL+DLOG(X/2.0D0)
           TTY=2.0D0*G1*BY0/(X*X)-G0*BY1/X
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE CJYLV(V,Z,CBJV,CDJV,CBYV,CDYV)
C
C       ===================================================
C       Purpose: Compute Bessel functions Jv(z) and Yv(z)
C                and their derivatives with a complex
C                argument and a large order
C       Input:   v --- Order of Jv(z) and Yv(z)
C                z --- Complex argument
C       Output:  CBJV --- Jv(z)
C                CDJV --- Jv'(z)
C                CBYV --- Yv(z)
C                CDYV --- Yv'(z)
C       Routine called:
C                CJK to compute the expansion coefficients
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CF(12),A(91)
        KM=12
        CALL CJK(KM,A)
        PI=3.141592653589793D0
        DO 30 L=1,0,-1
           V0=V-L
           CWS=CDSQRT(1.0D0-(Z/V0)*(Z/V0))
           CETA=CWS+CDLOG(Z/V0/(1.0D0+CWS))
           CT=1.0D0/CWS
           CT2=CT*CT
           DO 15 K=1,KM
              L0=K*(K+1)/2+1
              LF=L0+K
              CF(K)=A(LF)
              DO 10 I=LF-1,L0,-1
10               CF(K)=CF(K)*CT2+A(I)
15            CF(K)=CF(K)*CT**K
           VR=1.0D0/V0
           CSJ=(1.0D0,0.0D0)
           DO 20 K=1,KM
20            CSJ=CSJ+CF(K)*VR**K
           CBJV=CDSQRT(CT/(2.0D0*PI*V0))*CDEXP(V0*CETA)*CSJ
           IF (L.EQ.1) CFJ=CBJV
           CSY=(1.0D0,0.0D0)
           DO 25 K=1,KM
25            CSY=CSY+(-1)**K*CF(K)*VR**K
           CBYV=-CDSQRT(2.0D0*CT/(PI*V0))*CDEXP(-V0*CETA)*CSY
           IF (L.EQ.1) CFY=CBYV
30      CONTINUE
        CDJV=-V/Z*CBJV+CFJ
        CDYV=-V/Z*CBYV+CFY
        RETURN
        END



C       **********************************

        SUBROUTINE RMN2L(M,N,C,X,DF,KD,R2F,R2D,ID)
C
C       ========================================================
C       Purpose: Compute prolate and oblate spheroidal radial
C                functions of the second kind for given m, n,
C                c and a large cx
C       Routine called:
C                SPHY for computing the spherical Bessel
C                functions of the second kind
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION DF(200),SY(0:251),DY(0:251)
        EPS=1.0D-14
        IP=1
        NM1=INT((N-M)/2)
        IF (N-M.EQ.2*NM1) IP=0
        NM=25+NM1+INT(C)
        REG=1.0D0
        IF (M+NM.GT.80) REG=1.0D-200
        NM2=2*NM+M
        CX=C*X
        CALL SPHY(NM2,CX,NM2,SY,DY)
        R0=REG
        DO 10 J=1,2*M+IP
10         R0=R0*J
        R=R0
        SUC=R*DF(1)
        SW=0.0D0
        DO 15 K=2,NM
           R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
           SUC=SUC+R*DF(K)
           IF (K.GT.NM1.AND.DABS(SUC-SW).LT.DABS(SUC)*EPS) GO TO 20
15         SW=SUC
20      A0=(1.0D0-KD/(X*X))**(0.5D0*M)/SUC
        R2F=0.0D0
        EPS1=0.0D0
        NP=0
        DO 50 K=1,NM
           L=2*K+M-N-2+IP
           LG=1
           IF (L.NE.4*INT(L/4)) LG=-1
           IF (K.EQ.1) THEN
              R=R0
           ELSE
              R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
           ENDIF
           NP=M+2*K-2+IP
           R2F=R2F+LG*R*(DF(K)*SY(NP))
           EPS1=DABS(R2F-SW)
           IF (K.GT.NM1.AND.EPS1.LT.DABS(R2F)*EPS) GO TO 55
50         SW=R2F
55      ID1=INT(LOG10(EPS1/DABS(R2F)+EPS))
        R2F=R2F*A0
        IF (NP.GE.NM2) THEN
           ID=10
           RETURN
        ENDIF
        B0=KD*M/X**3.0D0/(1.0-KD/(X*X))*R2F
        SUD=0.0D0
        EPS2=0.0D0
        DO 60 K=1,NM
           L=2*K+M-N-2+IP
           LG=1
           IF (L.NE.4*INT(L/4)) LG=-1
           IF (K.EQ.1) THEN
              R=R0
           ELSE
              R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
           ENDIF
           NP=M+2*K-2+IP
           SUD=SUD+LG*R*(DF(K)*DY(NP))
           EPS2=DABS(SUD-SW)
           IF (K.GT.NM1.AND.EPS2.LT.DABS(SUD)*EPS) GO TO 65
60         SW=SUD
65      R2D=B0+A0*C*SUD
        ID2=INT(LOG10(EPS2/DABS(SUD)+EPS))
        ID=MAX(ID1,ID2)
        RETURN
        END



C       **********************************

        SUBROUTINE PSI_SPEC(X,PS)
C
C       ======================================
C       Purpose: Compute Psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=XA
           DO 10 K=1 ,N-1
10            S=S+1.0D0/K
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5)) THEN
           N=XA-.5
           DO 20 K=1,N
20            S=S+1.0/(2.0D0*K-1.0D0)
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
30               S=S+1.0D0/(XA+K)
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END

C       **********************************

        SUBROUTINE CVA2(KD,M,Q,A)
C
C       ======================================================
C       Purpose: Calculate a specific characteristic value of
C                Mathieu functions
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C                KD --- Case code
C                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
C                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
C                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
C                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
C       Output:  A  --- Characteristic value
C       Routines called:
C             (1) REFINE for finding accurate characteristic
C                 value using an iteration method
C             (2) CV0 for finding initial characteristic
C                 values using polynomial approximation
C             (3) CVQM for computing initial characteristic
C                 values for q ≤ 3*m
C             (3) CVQL for computing initial characteristic
C                 values for q ≥ m*m
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (M.LE.12.OR.Q.LE.3.0*M.OR.Q.GT.M*M) THEN
            CALL CV0(KD,M,Q,A)
            IF (Q.NE.0.0D0.AND.M.NE.2) CALL REFINE(KD,M,Q,A)
            IF (Q.GT.2.0D-3.AND.M.EQ.2) CALL REFINE(KD,M,Q,A)
        ELSE
           NDIV=10
           DELTA=(M-3.0)*M/NDIV
           IF ((Q-3.0*M).LE.(M*M-Q)) THEN
5             NN=INT((Q-3.0*M)/DELTA)+1
              DELTA=(Q-3.0*M)/NN
              Q1=2.0*M
              CALL CVQM(M,Q1,A1)
              Q2=3.0*M
              CALL CVQM(M,Q2,A2)
              QQ=3.0*M
              DO 10 I=1,NN
                 QQ=QQ+DELTA
                 A=(A1*Q2-A2*Q1+(A2-A1)*QQ)/(Q2-Q1)
                 IFLAG=1
                 IF (I.EQ.NN) IFLAG=-1
                 CALL REFINE(KD,M,QQ,A)
                 Q1=Q2
                 Q2=QQ
                 A1=A2
                 A2=A
10            CONTINUE
              IF (IFLAG.EQ.-10) THEN
                 NDIV=NDIV*2
                 DELTA=(M-3.0)*M/NDIV
                 GO TO 5
              ENDIF
           ELSE
15            NN=INT((M*M-Q)/DELTA)+1
              DELTA=(M*M-Q)/NN
              Q1=M*(M-1.0)
              CALL CVQL(KD,M,Q1,A1)
              Q2=M*M
              CALL CVQL(KD,M,Q2,A2)
              QQ=M*M
              DO 20 I=1,NN
                 QQ=QQ-DELTA
                 A=(A1*Q2-A2*Q1+(A2-A1)*QQ)/(Q2-Q1)
                 IFLAG=1
                 IF (I.EQ.NN) IFLAG=-1
                 CALL REFINE(KD,M,QQ,A)
                 Q1=Q2
                 Q2=QQ
                 A1=A2
                 A2=A
20            CONTINUE
              IF (IFLAG.EQ.-10) THEN
                 NDIV=NDIV*2
                 DELTA=(M-3.0)*M/NDIV
                 GO TO 15
              ENDIF
           ENDIF
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE LPMNS(M,N,X,PM,PD)
C
C       ========================================================
C       Purpose: Compute associated Legendre functions Pmn(x)
C                and Pmn'(x) for a given order
C       Input :  x --- Argument of Pmn(x)
C                m --- Order of Pmn(x),  m = 0,1,2,...,n
C                n --- Degree of Pmn(x), n = 0,1,2,...,N
C       Output:  PM(n) --- Pmn(x)
C                PD(n) --- Pmn'(x)
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION PM(0:N),PD(0:N)
        DO 10 K=0,N
           PM(K)=0.0D0
10         PD(K)=0.0D0
        IF (DABS(X).EQ.1.0D0) THEN
           DO 15 K=0,N
              IF (M.EQ.0) THEN
                 PM(K)=1.0D0
                 PD(K)=0.5D0*K*(K+1.0)
                 IF (X.LT.0.0) THEN
                    PM(K)=(-1)**K*PM(K)
                    PD(K)=(-1)**(K+1)*PD(K)
                 ENDIF
              ELSE IF (M.EQ.1) THEN
                 PD(K)=1.0D+300
              ELSE IF (M.EQ.2) THEN
                 PD(K)=-0.25D0*(K+2.0)*(K+1.0)*K*(K-1.0)
                 IF (X.LT.0.0) PD(K)=(-1)**(K+1)*PD(K)
              ENDIF
15         CONTINUE
           RETURN
        ENDIF
        X0=DABS(1.0D0-X*X)
        PM0=1.0D0
        PMK=PM0
        DO 20 K=1,M
           PMK=(2.0D0*K-1.0D0)*DSQRT(X0)*PM0
20         PM0=PMK
        PM1=(2.0D0*M+1.0D0)*X*PM0
        PM(M)=PMK
        PM(M+1)=PM1
        DO 25 K=M+2,N
           PM2=((2.0D0*K-1.0D0)*X*PM1-(K+M-1.0D0)*PMK)/(K-M)
           PM(K)=PM2
           PMK=PM1
25         PM1=PM2
        PD(0)=((1.0D0-M)*PM(1)-X*PM(0))/(X*X-1.0)
        DO 30 K=1,N
30         PD(K)=(K*X*PM(K)-(K+M)*PM(K-1))/(X*X-1.0D0)
        DO 35 K=1,N
           PM(K)=(-1)**M*PM(K)
35         PD(K)=(-1)**M*PD(K)
        RETURN
        END

C       **********************************

        SUBROUTINE CERF(Z,CER,CDER)
C
C       ==========================================================
C       Purpose: Compute complex Error function erf(z) & erf'(z)
C       Input:   z   --- Complex argument of erf(z)
C                x   --- Real part of z
C                y   --- Imaginary part of z
C       Output:  CER --- erf(z)
C                CDER --- erf'(z)
C       ==========================================================
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        COMPLEX *16 Z,CER,CDER
        EPS=1.0D-12
        PI=3.141592653589793D0
        X=DBLE(Z)
        Y=DIMAG(Z)
        X2=X*X
        IF (X.LE.3.5D0) THEN
           ER=1.0D0
           R=1.0D0
           W=0.0D0
           DO 10 K=1,100
              R=R*X2/(K+0.5D0)
              ER=ER+R
              IF (DABS(ER-W).LE.EPS*DABS(ER)) GO TO 15
10            W=ER
15         C0=2.0D0/DSQRT(PI)*X*DEXP(-X2)
           ER0=C0*ER
        ELSE
           ER=1.0D0
           R=1.0D0
           DO 20 K=1,12
              R=-R*(K-0.5D0)/X2
20            ER=ER+R
           C0=DEXP(-X2)/(X*DSQRT(PI))
           ER0=1.0D0-C0*ER
        ENDIF
        IF (Y.EQ.0.0D0) THEN
           ERR=ER0
           ERI=0.0D0
        ELSE
           CS=DCOS(2.0D0*X*Y)
           SS=DSIN(2.0D0*X*Y)
           ER1=DEXP(-X2)*(1.0D0-CS)/(2.0D0*PI*X)
           EI1=DEXP(-X2)*SS/(2.0D0*PI*X)
           ER2=0.0D0
           W1=0.0D0
           DO 25 N=1,100
              ER2=ER2+DEXP(-.25D0*N*N)/(N*N+4.0D0*X2)*(2.0D0*X
     &            -2.0D0*X*DCOSH(N*Y)*CS+N*DSINH(N*Y)*SS)
              IF (DABS((ER2-W1)/ER2).LT.EPS) GO TO 30
25            W1=ER2
30         C0=2.0D0*DEXP(-X2)/PI
           ERR=ER0+ER1+C0*ER2
           EI2=0.0D0
           W2=0.0D0
           DO 35 N=1,100
              EI2=EI2+DEXP(-.25D0*N*N)/(N*N+4.0D0*X2)*(2.0D0*X
     &            *DCOSH(N*Y)*SS+N*DSINH(N*Y)*CS)
              IF (DABS((EI2-W2)/EI2).LT.EPS) GO TO 40
35            W2=EI2
40         ERI=EI1+C0*EI2
        ENDIF
        CER = DCMPLX(ERR, ERI)
        CDER=2.0D0/DSQRT(PI)*CDEXP(-Z*Z)
        RETURN
        END

C       **********************************

        SUBROUTINE RSWFP(M,N,C,X,CV,KF,R1F,R1D,R2F,R2D)
C
C       ==============================================================
C       Purpose: Compute prolate spheriodal radial functions of the
C                first and second kinds, and their derivatives
C       Input :  m  --- Mode parameter, m = 0,1,2,...
C                n  --- Mode parameter, n = m,m+1,m+2,...
C                c  --- Spheroidal parameter
C                x  --- Argument of radial function ( x > 1.0 )
C                cv --- Characteristic value
C                KF --- Function code
C                       KF=1 for the first kind
C                       KF=2 for the second kind
C                       KF=3 for both the first and second kinds
C       Output:  R1F --- Radial function of the first kind
C                R1D --- Derivative of the radial function of
C                        the first kind
C                R2F --- Radial function of the second kind
C                R2D --- Derivative of the radial function of
C                        the second kind
C       Routines called:
C            (1) SDMN for computing expansion coefficients dk
C            (2) RMN1 for computing prolate and oblate radial
C                functions of the first kind
C            (3) RMN2L for computing prolate and oblate radial
C                functions of the second kind for a large argument
C            (4) RMN2SP for computing the prolate radial function
C                of the second kind for a small argument
C       ==============================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION DF(200)
        KD=1
        CALL SDMN(M,N,C,CV,KD,DF)
        IF (KF.NE.2) THEN
           CALL RMN1(M,N,C,X,DF,KD,R1F,R1D)
        ENDIF
        IF (KF.GT.1) THEN
           CALL RMN2L(M,N,C,X,DF,KD,R2F,R2D,ID)
           IF (ID.GT.-8) THEN
              CALL RMN2SP(M,N,C,X,CV,DF,KD,R2F,R2D)
           ENDIF
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
C
C       ===========================================================
C       Purpose: Compute Bessel functions Jn(x) and Yn(x), and
C                their first and second derivatives
C       Input:   x   ---  Argument of Jn(x) and Yn(x) ( x > 0 )
C                n   ---  Order of Jn(x) and Yn(x)
C       Output:  BJN ---  Jn(x)
C                DJN ---  Jn'(x)
C                FJN ---  Jn"(x)
C                BYN ---  Yn(x)
C                DYN ---  Yn'(x)
C                FYN ---  Yn"(x)
C       Routines called:
C                JYNBH to compute Jn and Yn
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BJ(2),BY(2)
        CALL JYNBH(N+1,N,X,NM,BJ,BY)
C       Compute derivatives by differentiation formulas
        BJN=BJ(1)
        BYN=BY(1)
        DJN=-BJ(2)+N*BJ(1)/X
        DYN=-BY(2)+N*BY(1)/X
        FJN=(N*N/(X*X)-1.0D0)*BJN-DJN/X
        FYN=(N*N/(X*X)-1.0D0)*BYN-DYN/X
        RETURN
        END


C       **********************************

        SUBROUTINE GAM0 (X,GA)
C
C       ================================================
C       Purpose: Compute gamma function Г(x)
C       Input :  x  --- Argument of Г(x)  ( |x| ≤ 1 )
C       Output:  GA --- Г(x)
C       ================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(25)
        DATA G/1.0D0,0.5772156649015329D0,
     &       -0.6558780715202538D0, -0.420026350340952D-1,
     &        0.1665386113822915D0, -.421977345555443D-1,
     &        -.96219715278770D-2, .72189432466630D-2,
     &        -.11651675918591D-2, -.2152416741149D-3,
     &         .1280502823882D-3, -.201348547807D-4,
     &        -.12504934821D-5, .11330272320D-5,
     &        -.2056338417D-6, .61160950D-8,
     &         .50020075D-8, -.11812746D-8,
     &         .1043427D-9, .77823D-11,
     &        -.36968D-11, .51D-12,
     &        -.206D-13, -.54D-14, .14D-14/
        GR=(25)
        DO 20 K=24,1,-1
20         GR=GR*X+G(K)
        GA=1.0D0/(GR*X)
        RETURN
        END


C       **********************************

        SUBROUTINE CISIB(X,CI,SI)
C
C       =============================================
C       Purpose: Compute cosine and sine integrals
C                Si(x) and Ci(x) ( x ≥ 0 )
C       Input :  x  --- Argument of Ci(x) and Si(x)
C       Output:  CI --- Ci(x)
C                SI --- Si(x)
C       =============================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        X2=X*X
        IF (X.EQ.0.0) THEN
           CI=-1.0D+300
           SI=0.0D0
        ELSE IF (X.LE.1.0D0) THEN
           CI=((((-3.0D-8*X2+3.10D-6)*X2-2.3148D-4)
     &        *X2+1.041667D-2)*X2-0.25)*X2+0.577215665D0+LOG(X)
           SI=((((3.1D-7*X2-2.834D-5)*X2+1.66667D-003)
     &        *X2-5.555556D-002)*X2+1.0)*X
        ELSE
           FX=((((X2+38.027264D0)*X2+265.187033D0)*X2
     &        +335.67732D0)*X2+38.102495D0)/((((X2
     &        +40.021433D0)*X2+322.624911D0)*X2
     &        +570.23628D0)*X2+157.105423D0)
           GX=((((X2+42.242855D0)*X2+302.757865D0)*X2
     &        +352.018498D0)*X2+21.821899D0)/((((X2
     &        +48.196927D0)*X2+482.485984D0)*X2
     &        +1114.978885D0)*X2+449.690326D0)/X
           CI=FX*SIN(X)/X-GX*COS(X)/X
           SI=1.570796327D0-FX*COS(X)/X-GX*SIN(X)/X
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE EULERA(N,EN)
C
C       ======================================
C       Purpose: Compute Euler number En
C       Input :  n --- Serial number
C       Output:  EN(n) --- En
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION EN(0:N)
        EN(0)=1.0D0
        DO 30 M=1,N/2
           S=1.0D0
           DO 20 K=1,M-1
              R=1.0D0
              DO 10 J=1,2*K
10               R=R*(2.0D0*M-2.0D0*K+J)/J
20            S=S+R*EN(2*K)
30         EN(2*M)=-S
        RETURN
        END

C       **********************************

        SUBROUTINE REFINE(KD,M,Q,A)
C
C       =====================================================
C       Purpose: calculate the accurate characteristic value
C                by the secant method
C       Input :  m --- Order of Mathieu functions
C                q --- Parameter of Mathieu functions
C                A --- Initial characteristic value
C       Output:  A --- Refineed characteristic value
C       Routine called:  CVF for computing the value of F for
C                        characteristic equation
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EPS=1.0D-14
        MJ=10+M
        CA=A
        DELTA=0.0D0
        X0=A
        CALL CVF(KD,M,Q,X0,MJ,F0)
        X1=1.002*A
        CALL CVF(KD,M,Q,X1,MJ,F1)
        DO 10 IT=1,100
           MJ=MJ+1
           X=X1-(X1-X0)/(1.0D0-F0/F1)
           CALL CVF(KD,M,Q,X,MJ,F)
           IF (ABS(1.0-X1/X).LT.EPS.OR.F.EQ.0.0) GO TO 15
           X0=X1
           F0=F1
           X1=X
10         F1=F
15      A=X
        RETURN
        END



C       **********************************

        SUBROUTINE CISIA(X,CI,SI)
C
C       =============================================
C       Purpose: Compute cosine and sine integrals
C                Si(x) and Ci(x)  ( x ≥ 0 )
C       Input :  x  --- Argument of Ci(x) and Si(x)
C       Output:  CI --- Ci(x)
C                SI --- Si(x)
C       =============================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BJ(101)
        P2=1.570796326794897D0
        EL=.5772156649015329D0
        EPS=1.0D-15
        X2=X*X
        IF (X.EQ.0.0D0) THEN
           CI=-1.0D+300
           SI=0.0D0
        ELSE IF (X.LE.16.0D0) THEN
           XR=-.25D0*X2
           CI=EL+DLOG(X)+XR
           DO 10 K=2,40
              XR=-.5D0*XR*(K-1)/(K*K*(2*K-1))*X2
              CI=CI+XR
              IF (DABS(XR).LT.DABS(CI)*EPS) GO TO 15
10         CONTINUE
15         XR=X
           SI=X
           DO 20 K=1,40
              XR=-.5D0*XR*(2*K-1)/K/(4*K*K+4*K+1)*X2
              SI=SI+XR
              IF (DABS(XR).LT.DABS(SI)*EPS) RETURN
20         CONTINUE
        ELSE IF (X.LE.32.0D0) THEN
           M=INT(47.2+.82*X)
           XA1=0.0D0
           XA0=1.0D-100
           DO 25 K=M,1,-1
              XA=4.0D0*K*XA0/X-XA1
              BJ(K)=XA
              XA1=XA0
25            XA0=XA
           XS=BJ(1)
           DO 30 K=3,M,2
30            XS=XS+2.0D0*BJ(K)
           BJ(1)=BJ(1)/XS
           DO 35 K=2,M
35            BJ(K)=BJ(K)/XS
           XR=1.0D0
           XG1=BJ(1)
           DO 40 K=2,M
              XR=.25D0*XR*(2.0*K-3.0)**2/((K-1.0)*(2.0*K-1.0)**2)*X
40            XG1=XG1+BJ(K)*XR
           XR=1.0D0
           XG2=BJ(1)
           DO 45 K=2,M
              XR=.25D0*XR*(2.0*K-5.0)**2/((K-1.0)*(2.0*K-3.0)**2)*X
45            XG2=XG2+BJ(K)*XR
           XCS=DCOS(X/2.0D0)
           XSS=DSIN(X/2.0D0)
           CI=EL+DLOG(X)-X*XSS*XG1+2*XCS*XG2-2*XCS*XCS
           SI=X*XCS*XG1+2*XSS*XG2-DSIN(X)
        ELSE
           XR=1.0D0
           XF=1.0D0
           DO 50 K=1,9
              XR=-2.0D0*XR*K*(2*K-1)/X2
50            XF=XF+XR
           XR=1.0D0/X
           XG=XR
           DO 55 K=1,8
              XR=-2.0D0*XR*(2*K+1)*K/X2
55            XG=XG+XR
           CI=XF*DSIN(X)/X-XG*DCOS(X)/X
           SI=P2-XF*DCOS(X)/X-XG*DSIN(X)/X
        ENDIF
        RETURN
        END


C       **********************************

        SUBROUTINE ITSL0(X,TL0)
C
C       ===========================================================
C       Purpose: Evaluate the integral of modified Struve function
C                L0(t) with respect to t from 0 to x
C       Input :  x   --- Upper limit  ( x ≥ 0 )
C       Output:  TL0 --- Integration of L0(t) from 0 to x
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(18)
        PI=3.141592653589793D0
        R=1.0D0
        IF (X.LE.20.0) THEN
           S=0.5D0
           DO 10 K=1,100
              RD=1.0D0
              IF (K.EQ.1) RD=0.5D0
              R=R*RD*K/(K+1.0D0)*(X/(2.0D0*K+1.0D0))**2
              S=S+R
              IF (DABS(R/S).LT.1.0D-12) GO TO 15
10         CONTINUE
15         TL0=2.0D0/PI*X*X*S
        ELSE
           S=1.0D0
           DO 20 K=1,10
              R=R*K/(K+1.0D0)*((2.0D0*K+1.0D0)/X)**2
              S=S+R
              IF (DABS(R/S).LT.1.0D-12) GO TO 25
20            CONTINUE
25         EL=.57721566490153D0
           S0=-S/(PI*X*X)+2.0D0/PI*(DLOG(2.0D0*X)+EL)
           A0=1.0D0
           A1=5.0D0/8.0D0
           A(1)=A1
           DO 30 K=1,10
              AF=((1.5D0*(K+.50D0)*(K+5.0D0/6.0D0)*A1-.5D0*
     &            (K+.5D0)**2*(K-.5D0)*A0))/(K+1.0D0)
              A(K+1)=AF
              A0=A1
30            A1=AF
           TI=1.0D0
           R=1.0D0
           DO 35 K=1,11
              R=R/X
35            TI=TI+A(K)*R
           TL0=TI/DSQRT(2*PI*X)*DEXP(X)+S0
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE STVL1(X,SL1)
C
C       ================================================
C       Purpose: Compute modified Struve function L1(x)
C       Input :  x   --- Argument of L1(x) ( x ≥ 0 )
C       Output:  SL1 --- L1(x)
C       ================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        R=1.0D0
        IF (X.LE.20.0D0) THEN
           S=0.0D0
           DO 10 K=1,60
              R=R*X*X/(4.0D0*K*K-1.0D0)
              S=S+R
              IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 15
10         CONTINUE
15         SL1=2.0D0/PI*S
        ELSE
           S=1.0D0
           KM=INT(.50*X)
           IF (X.GT.50) KM=25
           DO 20 K=1,KM
              R=R*(2.0D0*K+3.0D0)*(2.0D0*K+1.0D0)/(X*X)
              S=S+R
              IF (DABS(R/S).LT.1.0D-12) GO TO 25
20            CONTINUE
25         SL1=2.0D0/PI*(-1.0D0+1.0D0/(X*X)+3.0D0*S/X**4)
           A1=DEXP(X)/DSQRT(2.0D0*PI*X)
           R=1.0D0
           BI1=1.0D0
           DO 30 K=1,16
              R=-0.125D0*R*(4.0D0-(2.0D0*K-1.0D0)**2)/(K*X)
              BI1=BI1+R
              IF (DABS(R/BI1).LT.1.0D-12) GO TO 35
30         CONTINUE
35         SL1=SL1+A1*BI1
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE CLQN(N,X,Y,CQN,CQD)
C
C       ==================================================
C       Purpose: Compute the Legendre functions Qn(z) and
C                their derivatives Qn'(z) for a complex
C                argument
C       Input :  x --- Real part of z
C                y --- Imaginary part of z
C                n --- Degree of Qn(z), n = 0,1,2,...
C       Output:  CQN(n) --- Qn(z)
C                CQD(n) --- Qn'(z)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CQN(0:N),CQD(0:N)
        Z = DCMPLX(X, Y)
        IF (Z.EQ.1.0D0) THEN
           DO 10 K=0,N
              CQN(K)=(1.0D+300,0.0D0)
10            CQD(K)=(1.0D+300,0.0D0)
           RETURN
        ENDIF
        LS=1
        IF (CDABS(Z).GT.1.0D0) LS=-1
        CQ0=0.5D0*CDLOG(LS*(1.0D0+Z)/(1.0D0-Z))
        CQ1=Z*CQ0-1.0D0
        CQN(0)=CQ0
        CQN(1)=CQ1
        IF (CDABS(Z).LT.1.0001D0) THEN
           CQF0=CQ0
           CQF1=CQ1
           DO 15 K=2,N
              CQF2=((2.0D0*K-1.0D0)*Z*CQF1-(K-1.0D0)*CQF0)/K
              CQN(K)=CQF2
              CQF0=CQF1
15            CQF1=CQF2
        ELSE
           IF (CDABS(Z).GT.1.1D0) THEN
              KM=40+N
           ELSE
              KM=(40+N)*INT(-1.0-1.8*LOG(CDABS(Z-1.0)))
           ENDIF
           CQF2=0.0D0
           CQF1=1.0D0
           DO 20 K=KM,0,-1
              CQF0=((2*K+3.0D0)*Z*CQF1-(K+2.0D0)*CQF2)/(K+1.0D0)
              IF (K.LE.N) CQN(K)=CQF0
              CQF2=CQF1
20            CQF1=CQF0
           DO 25 K=0,N
25            CQN(K)=CQN(K)*CQ0/CQF0
        ENDIF
        CQD(0)=(CQN(1)-Z*CQN(0))/(Z*Z-1.0D0)
        DO 30 K=1,N
30         CQD(K)=(K*Z*CQN(K)-K*CQN(K-1))/(Z*Z-1.0D0)
        RETURN
        END

C       **********************************

        SUBROUTINE STVL0(X,SL0)
C
C       ================================================
C       Purpose: Compute modified Struve function L0(x)
C       Input :  x   --- Argument of L0(x) ( x ≥ 0 )
C       Output:  SL0 --- L0(x)
C       ================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        S=1.0D0
        R=1.0D0
        IF (X.LE.20.0D0) THEN
           A0=2.0D0*X/PI
           DO 10 K=1,60
              R=R*(X/(2.0D0*K+1.0D0))**2
              S=S+R
              IF (DABS(R/S).LT.1.0D-12) GO TO 15
10         CONTINUE
15         SL0=A0*S
        ELSE
           KM=INT(.5*(X+1.0))
           IF (X.GE.50.0) KM=25
           DO 20 K=1,KM
              R=R*((2.0D0*K-1.0D0)/X)**2
              S=S+R
              IF (DABS(R/S).LT.1.0D-12) GO TO 25
20         CONTINUE
25         A1=DEXP(X)/DSQRT(2.0D0*PI*X)
           R=1.0D0
           BI0=1.0D0
           DO 30 K=1,16
              R=0.125D0*R*(2.0D0*K-1.0D0)**2/(K*X)
              BI0=BI0+R
              IF (DABS(R/BI0).LT.1.0D-12) GO TO 35
30         CONTINUE
35         BI0=A1*BI0
           SL0=-2.0D0/(PI*X)*S+BI0
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE AIRYZO(NT,KF,XA,XB,XC,XD)
C
C       ========================================================
C       Purpose: Compute the first NT zeros of Airy functions
C                Ai(x) and Ai'(x), a and a', and the associated
C                values of Ai(a') and Ai'(a); and the first NT
C                zeros of Airy functions Bi(x) and Bi'(x), b and
C                b', and the associated values of Bi(b') and
C                Bi'(b)
C       Input :  NT    --- Total number of zeros
C                KF    --- Function code
C                          KF=1 for Ai(x) and Ai'(x)
C                          KF=2 for Bi(x) and Bi'(x)
C       Output:  XA(m) --- a, the m-th zero of Ai(x) or
C                          b, the m-th zero of Bi(x)
C                XB(m) --- a', the m-th zero of Ai'(x) or
C                          b', the m-th zero of Bi'(x)
C                XC(m) --- Ai(a') or Bi(b')
C                XD(m) --- Ai'(a) or Bi'(b)
C                          ( m --- Serial number of zeros )
C       Routine called: AIRYB for computing Airy functions and
C                       their derivatives
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION XA(NT),XB(NT),XC(NT),XD(NT)
        PI=3.141592653589793D0
        RT0=0.0D0
        RT=0.0D0
        DO 15 I=1,NT
           IF (KF.EQ.1) THEN
              U=3.0*PI*(4.0*I-1)/8.0D0
              U1=1/(U*U)
              RT0=-(U*U)**(1.0/3.0)*((((-15.5902*U1+.929844)*U1
     &            -.138889)*U1+.10416667D0)*U1+1.0D0)
           ELSE IF (KF.EQ.2) THEN
              IF (I.EQ.1) THEN
                 RT0=-1.17371
              ELSE
                 U=3.0*PI*(4.0*I-3.0)/8.0
                 U1=1.0D0/(U*U)
                 RT0=-(U*U)**(1.0/3.0)*((((-15.5902*U1+.929844)*U1
     &               -.138889)*U1+.10416667)*U1+1.0)
              ENDIF
           ENDIF
10         X=RT0
           CALL AIRYB(X,AI,BI,AD,BD)
           IF (KF.EQ.1) RT=RT0-AI/AD
           IF (KF.EQ.2) RT=RT0-BI/BD
           IF (DABS((RT-RT0)/RT).GT.1.D-9) THEN
              RT0=RT
              GOTO 10
           ELSE
              XA(I)=RT
              IF (KF.EQ.1) XD(I)=AD
              IF (KF.EQ.2) XD(I)=BD
           ENDIF
15      CONTINUE
        DO 25 I=1,NT
           IF (KF.EQ.1) THEN
              IF (I.EQ.1) THEN
                 RT0=-1.01879
              ELSE
                 U=3.0*PI*(4.0*I-3.0)/8.0
                 U1=1/(U*U)
                 RT0=-(U*U)**(1.0/3.0)*((((15.0168*U1-.873954)
     &            *U1+.121528)*U1-.145833D0)*U1+1.0D0)
              ENDIF
           ELSE IF (KF.EQ.2) THEN
              IF (I.EQ.1) THEN
                 RT0=-2.29444
              ELSE
                 U=3.0*PI*(4.0*I-1.0)/8.0
                 U1=1.0/(U*U)
                 RT0=-(U*U)**(1.0/3.0)*((((15.0168*U1-.873954)
     &               *U1+.121528)*U1-.145833)*U1+1.0)
              ENDIF
           ENDIF
20         X=RT0
           CALL AIRYB(X,AI,BI,AD,BD)
           IF (KF.EQ.1) RT=RT0-AD/(AI*X)
           IF (KF.EQ.2) RT=RT0-BD/(BI*X)
           IF (DABS((RT-RT0)/RT).GT.1.0D-9) THEN
              RT0=RT
              GOTO 20
           ELSE
              XB(I)=RT
              IF (KF.EQ.1) XC(I)=AI
              IF (KF.EQ.2) XC(I)=BI
           ENDIF
25      CONTINUE
        RETURN
        END



C       **********************************

        SUBROUTINE ERROR(X,ERR)
C
C       =========================================
C       Purpose: Compute error function erf(x)
C       Input:   x   --- Argument of erf(x)
C       Output:  ERR --- erf(x)
C       =========================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EPS=1.0D-15
        PI=3.141592653589793D0
        X2=X*X
        IF (DABS(X).LT.3.5D0) THEN
           ER=1.0D0
           R=1.0D0
           DO 10 K=1,50
              R=R*X2/(K+0.5D0)
              ER=ER+R
              IF (DABS(R).LE.DABS(ER)*EPS) GO TO 15
10         CONTINUE
15         C0=2.0D0/DSQRT(PI)*X*DEXP(-X2)
           ERR=C0*ER
        ELSE
           ER=1.0D0
           R=1.0D0
           DO 20 K=1,12
              R=-R*(K-0.5D0)/X2
20            ER=ER+R
           C0=DEXP(-X2)/(DABS(X)*DSQRT(PI))
           ERR=1.0D0-C0*ER
           IF (X.LT.0.0) ERR=-ERR
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE CERROR(Z,CER)
C
C       ====================================================
C       Purpose: Compute error function erf(z) for a complex
C                argument (z=x+iy)
C       Input :  z   --- Complex argument
C       Output:  CER --- erf(z)
C       ====================================================
C
        IMPLICIT COMPLEX *16 (C,Z)
        DOUBLE PRECISION A0,PI
        A0=CDABS(Z)
        C0=CDEXP(-Z*Z)
        PI=3.141592653589793D0
        Z1=Z
        IF (DBLE(Z).LT.0.0) THEN
           Z1=-Z
        ENDIF
C
C       Cutoff radius R = 4.36; determined by balancing rounding error
C       and asymptotic expansion error, see below.
C
C       The resulting maximum global accuracy expected is around 1e-8
C
        IF (A0.LE.4.36D0) THEN
C
C          Rounding error in the Taylor expansion is roughly
C
C          ~ R*R * EPSILON * R**(2 R**2) / (2 R**2 Gamma(R**2 + 1/2))
C
           CS=Z1
           CR=Z1
           DO 10 K=1,120
              CR=CR*Z1*Z1/(K+0.5D0)
              CS=CS+CR
              IF (CDABS(CR/CS).LT.1.0D-15) GO TO 15
10         CONTINUE
15         CER=2.0D0*C0*CS/DSQRT(PI)
        ELSE
           CL=1.0D0/Z1
           CR=CL
C
C          Asymptotic series; maximum K must be at most ~ R^2.
C
C          The maximum accuracy obtainable from this expansion is roughly
C
C          ~ Gamma(2R**2 + 2) / (
C                   (2 R**2)**(R**2 + 1/2) Gamma(R**2 + 3/2) 2**(R**2 + 1/2))
C
           DO 20 K=1,20
              CR=-CR*(K-0.5D0)/(Z1*Z1)
              CL=CL+CR
              IF (CDABS(CR/CL).LT.1.0D-15) GO TO 25
20         CONTINUE
25         CER=1.0D0-C0*CL/DSQRT(PI)
        ENDIF
        IF (DBLE(Z).LT.0.0) THEN
           CER=-CER
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE EULERB(N,EN)
C
C       ======================================
C       Purpose: Compute Euler number En
C       Input :  n --- Serial number
C       Output:  EN(n) --- En
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION EN(0:N)
        HPI=2.0D0/3.141592653589793D0
        EN(0)=1.0D0
        EN(2)=-1.0D0
        R1=-4.0D0*HPI**3
        DO 20 M=4,N,2
           R1=-R1*(M-1)*M*HPI*HPI
           R2=1.0D0
           ISGN=1.0D0
           DO 10 K=3,1000,2
              ISGN=-ISGN
              S=(1.0D0/K)**(M+1)
              R2=R2+ISGN*S
              IF (S.LT.1.0D-15) GOTO 20
10         CONTINUE
20         EN(M)=R1*R2
        RETURN
        END

C       **********************************

        SUBROUTINE CVA1(KD,M,Q,CV)
C
C       ============================================================
C       Purpose: Compute a sequence of characteristic values of
C                Mathieu functions
C       Input :  M  --- Maximum order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C                KD --- Case code
C                       KD=1 for cem(x,q)  ( m = 0,2,4,… )
C                       KD=2 for cem(x,q)  ( m = 1,3,5,… )
C                       KD=3 for sem(x,q)  ( m = 1,3,5,… )
C                       KD=4 for sem(x,q)  ( m = 2,4,6,… )
C       Output:  CV(I) --- Characteristic values; I = 1,2,3,...
C                For KD=1, CV(1), CV(2), CV(3),..., correspond to
C                the characteristic values of cem for m = 0,2,4,...
C                For KD=2, CV(1), CV(2), CV(3),..., correspond to
C                the characteristic values of cem for m = 1,3,5,...
C                For KD=3, CV(1), CV(2), CV(3),..., correspond to
C                the characteristic values of sem for m = 1,3,5,...
C                For KD=4, CV(1), CV(2), CV(3),..., correspond to
C                the characteristic values of sem for m = 0,2,4,...
C       ============================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(200),H(200),D(500),E(500),F(500),CV(200)
        EPS=1.0D-14
        ICM=INT(M/2)+1
        IF (KD.EQ.4) ICM=M/2
        IF (Q.EQ.0.0D0) THEN
           IF (KD.EQ.1) THEN
              DO 10 IC=1,ICM
10               CV(IC)=4.0D0*(IC-1.0D0)**2
           ELSE IF (KD.NE.4) THEN
              DO 15 IC=1,ICM
15               CV(IC)=(2.0D0*IC-1.0D0)**2
           ELSE
              DO 20 IC=1,ICM
20               CV(IC)=4.0D0*IC*IC
           ENDIF
        ELSE
           NM=INT(10+1.5*M+0.5*Q)
           E(1)=0.0D0
           F(1)=0.0D0
           IF (KD.EQ.1) THEN
              D(1)=0.0D0
              DO 25 I=2,NM
                 D(I)=4.0D0*(I-1.0D0)**2
                 E(I)=Q
25               F(I)=Q*Q
              E(2)=DSQRT(2.0D0)*Q
              F(2)=2.0D0*Q*Q
           ELSE IF (KD.NE.4) THEN
              D(1)=1.0D0+(-1)**KD*Q
              DO 30 I=2,NM
                 D(I)=(2.0D0*I-1.0D0)**2
                 E(I)=Q
30               F(I)=Q*Q
           ELSE
              D(1)=4.0D0
              DO 35 I=2,NM
                 D(I)=4.0D0*I*I
                 E(I)=Q
35               F(I)=Q*Q
           ENDIF
           XA=D(NM)+DABS(E(NM))
           XB=D(NM)-DABS(E(NM))
           NM1=NM-1
           DO 40 I=1,NM1
              T=DABS(E(I))+DABS(E(I+1))
              T1=D(I)+T
              IF (XA.LT.T1) XA=T1
              T1=D(I)-T
              IF (T1.LT.XB) XB=T1
40         CONTINUE
           DO 45 I=1,ICM
              G(I)=XA
45            H(I)=XB
           DO 75 K=1,ICM
              DO 50 K1=K,ICM
                 IF (G(K1).LT.G(K)) THEN
                    G(K)=G(K1)
                    GO TO 55
                 ENDIF
50            CONTINUE
55            IF (K.NE.1.AND.H(K).LT.H(K-1)) H(K)=H(K-1)
60            X1=(G(K)+H(K))/2.0D0
              CV(K)=X1
              IF (DABS((G(K)-H(K))/X1).LT.EPS) GO TO 70
              J=0
              S=1.0D0
              DO 65 I=1,NM
                 IF (S.EQ.0.0D0) S=S+1.0D-30
                 T=F(I)/S
                 S=D(I)-T-X1
                 IF (S.LT.0.0) J=J+1
65            CONTINUE
              IF (J.LT.K) THEN
                 H(K)=X1
              ELSE
                 G(K)=X1
                 IF (J.GE.ICM) THEN
                    G(ICM)=X1
                 ELSE
                    IF (H(J+1).LT.X1) H(J+1)=X1
                    IF (X1.LT.G(J)) G(J)=X1
                 ENDIF
              ENDIF
              GO TO 60
70            CV(K)=X1
75         CONTINUE
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE ITTIKB(X,TTI,TTK)
C
C       =========================================================
C       Purpose: Integrate [I0(t)-1]/t with respect to t from 0
C                to x, and K0(t)/t with respect to t from x to ∞
C       Input :  x   --- Variable in the limits  ( x ≥ 0 )
C       Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
C                TTK --- Integration of K0(t)/t from x to ∞
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        IF (X.EQ.0.0D0) THEN
           TTI=0.0D0
        ELSE IF (X.LE.5.0D0) THEN
           X1=X/5.0D0
           T=X1*X1
           TTI=(((((((.1263D-3*T+.96442D-3)*T+.968217D-2)*T
     &         +.06615507D0)*T+.33116853D0)*T+1.13027241D0)
     &         *T+2.44140746D0)*T+3.12499991D0)*T
        ELSE
           T=5.0D0/X
           TTI=(((((((((2.1945464D0*T-3.5195009D0)*T
     &         -11.9094395D0)*T+40.394734D0)*T-48.0524115D0)
     &         *T+28.1221478D0)*T-8.6556013D0)*T+1.4780044D0)
     &         *T-.0493843D0)*T+.1332055D0)*T+.3989314D0
           TTI=TTI*DEXP(X)/(DSQRT(X)*X)
        ENDIF
        IF (X.EQ.0.0D0) THEN
           TTK=1.0D+300
        ELSE IF (X.LE.2.0D0) THEN
           T1=X/2.0D0
           T=T1*T1
           TTK=(((((.77D-6*T+.1544D-4)*T+.48077D-3)*T
     &         +.925821D-2)*T+.10937537D0)*T+.74999993D0)*T
           E0=EL+DLOG(X/2.0D0)
           TTK=PI*PI/24.0D0+E0*(.5D0*E0+TTI)-TTK
        ELSE IF (X.LE.4.0D0) THEN
           T=2.0D0/X
           TTK=(((.06084D0*T-.280367D0)*T+.590944D0)*T
     &         -.850013D0)*T+1.234684D0
           TTK=TTK*DEXP(-X)/(DSQRT(X)*X)
        ELSE
           T=4.0D0/X
           TTK=(((((.02724D0*T-.1110396D0)*T+.2060126D0)*T
     &         -.2621446D0)*T+.3219184D0)*T-.5091339D0)*T
     &         +1.2533141D0
           TTK=TTK*DEXP(-X)/(DSQRT(X)*X)
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE LQNB(N,X,QN,QD)
C
C       ====================================================
C       Purpose: Compute Legendre functions Qn(x) & Qn'(x)
C       Input :  x  --- Argument of Qn(x)
C                n  --- Degree of Qn(x)  ( n = 0,1,2,…)
C       Output:  QN(n) --- Qn(x)
C                QD(n) --- Qn'(x)
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION QN(0:N),QD(0:N)
        EPS=1.0D-14
        IF (DABS(X).EQ.1.0D0) THEN
           DO 10 K=0,N
              QN(K)=1.0D+300
10            QD(K)=1.0D+300
           RETURN
        ENDIF
        IF (X.LE.1.021D0) THEN
           X2=DABS((1.0D0+X)/(1.0D0-X))
           Q0=0.5D0*DLOG(X2)
           Q1=X*Q0-1.0D0
           QN(0)=Q0
           QN(1)=Q1
           QD(0)=1.0D0/(1.0D0-X*X)
           QD(1)=QN(0)+X*QD(0)
           DO 15 K=2,N
              QF=((2.0D0*K-1.0D0)*X*Q1-(K-1.0D0)*Q0)/K
              QN(K)=QF
              QD(K)=(QN(K-1)-X*QF)*K/(1.0D0-X*X)
              Q0=Q1
15            Q1=QF
        ELSE
           QC1=0.0D0
           QC2=1.0D0/X
           DO 20 J=1,N
              QC2=QC2*J/((2.0*J+1.0D0)*X)
              IF (J.EQ.N-1) QC1=QC2
20         CONTINUE
           DO 35 L=0,1
              NL=N+L
              QF=1.0D0
              QR=1.0D0
              DO 25 K=1,500
                 QR=QR*(0.5D0*NL+K-1.0D0)*(0.5D0*(NL-1)+K)
     &              /((NL+K-0.5D0)*K*X*X)
                 QF=QF+QR
                 IF (DABS(QR/QF).LT.EPS) GO TO 30
25            CONTINUE
30            IF (L.EQ.0) THEN
                 QN(N-1)=QF*QC1
              ELSE
                 QN(N)=QF*QC2
              ENDIF
35         CONTINUE
           QF2=QN(N)
           QF1=QN(N-1)
           DO 40 K=N,2,-1
              QF0=((2*K-1.0D0)*X*QF1-K*QF2)/(K-1.0D0)
              QN(K-2)=QF0
              QF2=QF1
40            QF1=QF0
           QD(0)=1.0D0/(1.0D0-X*X)
           DO 45 K=1,N
45            QD(K)=K*(QN(K-1)-X*QN(K))/(1.0D0-X*X)
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE CJK(KM,A)
C
C       ========================================================
C       Purpose: Compute the expansion coefficients for the
C                asymptotic expansion of Bessel functions
C                with large orders
C       Input :  Km   --- Maximum k
C       Output:  A(L) --- Cj(k) where j and k are related to L
C                         by L=j+1+[k*(k+1)]/2; j,k=0,1,...,Km
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(*)
        A(1)=1.0D0
        F0=1.0D0
        G0=1.0D0
        DO 10 K=0,KM-1
           L1=(K+1)*(K+2)/2+1
           L2=(K+1)*(K+2)/2+K+2
           F=(0.5D0*K+0.125D0/(K+1))*F0
           G=-(1.5D0*K+0.625D0/(3.0*(K+1.0D0)))*G0
           A(L1)=F
           A(L2)=G
           F0=F
10         G0=G
        DO 15 K=1,KM-1
           DO 15 J=1,K
              L3=K*(K+1)/2+J+1
              L4=(K+1)*(K+2)/2+J+1
              A(L4)=(J+0.5D0*K+0.125D0/(2.0*J+K+1.0))*A(L3)
     &             -(J+0.5D0*K-1.0+0.625D0/(2.0*J+K+1.0))*A(L3-1)
15         CONTINUE
        RETURN
        END


C       **********************************

        SUBROUTINE ITTIKA(X,TTI,TTK)
C
C       =========================================================
C       Purpose: Integrate [I0(t)-1]/t with respect to t from 0
C                to x, and K0(t)/t with respect to t from x to ∞
C       Input :  x   --- Variable in the limits  ( x ≥ 0 )
C       Output:  TTI --- Integration of [I0(t)-1]/t from 0 to x
C                TTK --- Integration of K0(t)/t from x to ∞
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION C(8)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        DATA C/1.625D0,4.1328125D0,
     &       1.45380859375D+1,6.553353881835D+1,
     &       3.6066157150269D+2,2.3448727161884D+3,
     &       1.7588273098916D+4,1.4950639538279D+5/
        IF (X.EQ.0.0D0) THEN
           TTI=0.0D0
           TTK=1.0D+300
           RETURN
        ENDIF
        IF (X.LT.40.0D0) THEN
           TTI=1.0D0
           R=1.0D0
           DO 10 K=2,50
              R=.25D0*R*(K-1.0D0)/(K*K*K)*X*X
              TTI=TTI+R
              IF (DABS(R/TTI).LT.1.0D-12) GO TO 15
10         CONTINUE
15         TTI=TTI*.125D0*X*X
        ELSE
           TTI=1.0D0
           R=1.0D0
           DO 20 K=1,8
              R=R/X
20            TTI=TTI+C(K)*R
           RC=X*DSQRT(2.0D0*PI*X)
           TTI=TTI*DEXP(X)/RC
        ENDIF
        IF (X.LE.12.0D0) THEN
           E0=(.5D0*DLOG(X/2.0D0)+EL)*DLOG(X/2.0D0)
     &        +PI*PI/24.0D0+.5D0*EL*EL
           B1=1.5D0-(EL+DLOG(X/2.0D0))
           RS=1.0D0
           R=1.0D0
           DO 25 K=2,50
              R=.25D0*R*(K-1.0D0)/(K*K*K)*X*X
              RS=RS+1.0D0/K
              R2=R*(RS+1.0D0/(2.0D0*K)-(EL+DLOG(X/2.0D0)))
              B1=B1+R2
              IF (DABS(R2/B1).LT.1.0D-12) GO TO 30
25         CONTINUE
30         TTK=E0-.125D0*X*X*B1
        ELSE
           TTK=1.0D0
           R=1.0D0
           DO 35 K=1,8
              R=-R/X
35            TTK=TTK+C(K)*R
           RC=X*DSQRT(2.0D0/PI*X)
           TTK=TTK*DEXP(-X)/RC
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE LAMV(V,X,VM,VL,DL)
C
C       =========================================================
C       Purpose: Compute lambda function with arbitrary order v,
C                and their derivative
C       Input :  x --- Argument of lambda function
C                v --- Order of lambda function
C       Output:  VL(n) --- Lambda function of order n+v0
C                DL(n) --- Derivative of lambda function
C                VM --- Highest order computed
C       Routines called:
C            (1) MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C            (2) GAM0 for computing gamma function (|x| ≤ 1)
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION VL(0:*),DL(0:*)
        PI=3.141592653589793D0
        RP2=0.63661977236758D0
        X=DABS(X)
        X2=X*X
        N=INT(V)
        V0=V-N
        VM=V
        IF (X.LE.12.0D0) THEN
           DO 25 K=0,N
              VK=V0+K
              BK=1.0D0
              R=1.0D0
              DO 10 I=1,50
                 R=-0.25D0*R*X2/(I*(I+VK))
                 BK=BK+R
                 IF (DABS(R).LT.DABS(BK)*1.0D-15) GO TO 15
10            CONTINUE
15            VL(K)=BK
              UK=1.0D0
              R=1.0D0
              DO 20 I=1,50
                 R=-0.25D0*R*X2/(I*(I+VK+1.0D0))
                 UK=UK+R
                 IF (DABS(R).LT.DABS(UK)*1.0D-15) GO TO 25
20            CONTINUE
25            DL(K)=-0.5D0*X/(VK+1.0D0)*UK
           RETURN
        ENDIF
        K0=11
        IF (X.GE.35.0D0) K0=10
        IF (X.GE.50.0D0) K0=8
        BJV0=0.0D0
        BJV1=0.0D0
        DO 40 J=0,1
           VV=4.0D0*(J+V0)*(J+V0)
           PX=1.0D0
           RP=1.0D0
           DO 30 K=1,K0
              RP=-0.78125D-2*RP*(VV-(4.0*K-3.0)**2.0)*(VV-
     &            (4.0*K-1.0)**2.0)/(K*(2.0*K-1.0)*X2)
30            PX=PX+RP
           QX=1.0D0
           RQ=1.0D0
           DO 35 K=1,K0
              RQ=-0.78125D-2*RQ*(VV-(4.0*K-1.0)**2.0)*(VV-
     &            (4.0*K+1.0)**2.0)/(K*(2.0*K+1.0)*X2)
35            QX=QX+RQ
           QX=0.125D0*(VV-1.0D0)*QX/X
           XK=X-(0.5D0*(J+V0)+0.25D0)*PI
           A0=DSQRT(RP2/X)
           CK=DCOS(XK)
           SK=DSIN(XK)
           IF (J.EQ.0) BJV0=A0*(PX*CK-QX*SK)
           IF (J.EQ.1) BJV1=A0*(PX*CK-QX*SK)
40      CONTINUE
        IF (V0.EQ.0.0D0) THEN
           GA=1.0D0
        ELSE
           CALL GAM0(V0,GA)
           GA=V0*GA
        ENDIF
        FAC=(2.0D0/X)**V0*GA
        VL(0)=BJV0
        DL(0)=-BJV1+V0/X*BJV0
        VL(1)=BJV1
        DL(1)=BJV0-(1.0D0+V0)/X*BJV1
        R0=2.0D0*(1.0D0+V0)/X
        IF (N.LE.1) THEN
           VL(0)=FAC*VL(0)
           DL(0)=FAC*DL(0)-V0/X*VL(0)
           VL(1)=FAC*R0*VL(1)
           DL(1)=FAC*R0*DL(1)-(1.0D0+V0)/X*VL(1)
           RETURN
        ENDIF
        IF (N.GE.2.AND.N.LE.INT(0.9*X)) THEN
           F0=BJV0
           F1=BJV1
           DO 45 K=2,N
              F=2.0D0*(K+V0-1.0D0)/X*F1-F0
              F0=F1
              F1=F
45            VL(K)=F
        ELSE IF (N.GE.2) THEN
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              N=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F=0.0D0
           F2=0.0D0
           F1=1.0D-100
           DO 50 K=M,0,-1
              F=2.0D0*(V0+K+1.0D0)/X*F1-F2
              IF (K.LE.N) VL(K)=F
              F2=F1
50            F1=F
           CS=0.0D0
           IF (DABS(BJV0).GT.DABS(BJV1)) CS=BJV0/F
           ELSE CS=BJV1/F2
           DO 55 K=0,N
55            VL(K)=CS*VL(K)
        ENDIF
        VL(0)=FAC*VL(0)
        DO 65 J=1,N
           RC=FAC*R0
           VL(J)=RC*VL(J)
           DL(J-1)=-0.5D0*X/(J+V0)*VL(J)
65         R0=2.0D0*(J+V0+1)/X*R0
        DL(N)=2.0D0*(V0+N)*(VL(N-1)-VL(N))/X
        VM=N+V0
        RETURN
        END



C       **********************************

        SUBROUTINE CHGUIT(A,B,X,HU,ID)
C
C       ======================================================
C       Purpose: Compute hypergeometric function U(a,b,x) by
C                using Gaussian-Legendre integration (n=60)
C       Input  : a  --- Parameter ( a > 0 )
C                b  --- Parameter
C                x  --- Argument ( x > 0 )
C       Output:  HU --- U(a,b,z)
C                ID --- Estimated number of significant digits
C       Routine called: GAMMA2 for computing Г(x)
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION T(30),W(30)
        DATA T/ .259597723012478D-01, .778093339495366D-01,
     &          .129449135396945D+00, .180739964873425D+00,
     &          .231543551376029D+00, .281722937423262D+00,
     &          .331142848268448D+00, .379670056576798D+00,
     &          .427173741583078D+00, .473525841761707D+00,
     &          .518601400058570D+00, .562278900753945D+00,
     &          .604440597048510D+00, .644972828489477D+00,
     &          .683766327381356D+00, .720716513355730D+00,
     &          .755723775306586D+00, .788693739932264D+00,
     &          .819537526162146D+00, .848171984785930D+00,
     &          .874519922646898D+00, .898510310810046D+00,
     &          .920078476177628D+00, .939166276116423D+00,
     &          .955722255839996D+00, .969701788765053D+00,
     &          .981067201752598D+00, .989787895222222D+00,
     &          .995840525118838D+00, .999210123227436D+00/
        DATA W/ .519078776312206D-01, .517679431749102D-01,
     &          .514884515009810D-01, .510701560698557D-01,
     &          .505141845325094D-01, .498220356905502D-01,
     &          .489955754557568D-01, .480370318199712D-01,
     &          .469489888489122D-01, .457343797161145D-01,
     &          .443964787957872D-01, .429388928359356D-01,
     &          .413655512355848D-01, .396806954523808D-01,
     &          .378888675692434D-01, .359948980510845D-01,
     &          .340038927249464D-01, .319212190192963D-01,
     &          .297524915007890D-01, .275035567499248D-01,
     &          .251804776215213D-01, .227895169439978D-01,
     &          .203371207294572D-01, .178299010142074D-01,
     &          .152746185967848D-01, .126781664768159D-01,
     &          .100475571822880D-01, .738993116334531D-02,
     &          .471272992695363D-02, .202681196887362D-02/
        ID=9
C       DLMF 13.4.4, integration up to C=12/X
        A1=A-1.0D0
        B1=B-A-1.0D0
        C=12.0D0/X
        HU0=0.0D0
        DO 20 M=10,100,5
           HU1=0.0D0
           G=0.5D0*C/M
           D=G
           DO 15 J=1,M
              S=0.0D0
              DO 10 K=1,30
                 T1=D+G*T(K)
                 T2=D-G*T(K)
                 F1=DEXP(-X*T1)*T1**A1*(1.0D0+T1)**B1
                 F2=DEXP(-X*T2)*T2**A1*(1.0D0+T2)**B1
                 S=S+W(K)*(F1+F2)
10            CONTINUE
              HU1=HU1+S*G
              D=D+2.0D0*G
15         CONTINUE
           IF (DABS(1.0D0-HU0/HU1).LT.1.0D-9) GO TO 25
           HU0=HU1
20      CONTINUE
25      CALL GAMMA2(A,GA)
        HU1=HU1/GA
C       DLMF 13.4.4 with substitution t=C/(1-u)
C       integration u from 0 to 1, i.e. t from C=12/X to infinity
        DO 40 M=2,10,2
           HU2=0.0D0
           G=0.5D0/M
           D=G
           DO 35 J=1,M
              S=0.0D0
              DO 30 K=1,30
                 T1=D+G*T(K)
                 T2=D-G*T(K)
                 T3=C/(1.0D0-T1)
                 T4=C/(1.0D0-T2)
                 F1=T3*T3/C*DEXP(-X*T3)*T3**A1*(1.0D0+T3)**B1
                 F2=T4*T4/C*DEXP(-X*T4)*T4**A1*(1.0D0+T4)**B1
                 S=S+W(K)*(F1+F2)
30            CONTINUE
              HU2=HU2+S*G
              D=D+2.0D0*G
35         CONTINUE
           IF (DABS(1.0D0-HU0/HU2).LT.1.0D-9) GO TO 45
           HU0=HU2
40      CONTINUE
45      CALL GAMMA2(A,GA)
        HU2=HU2/GA
        HU=HU1+HU2
        RETURN
        END



C       **********************************

        SUBROUTINE KMN(M,N,C,CV,KD,DF,DN,CK1,CK2)
C
C       ===================================================
C       Purpose: Compute the expansion coefficients of the
C                prolate and oblate spheroidal functions
C                and joining factors
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION U(200),V(200),W(200),DF(200),DN(200),
     &            TP(200),RK(200)
        NM=25+INT(0.5*(N-M)+C)
        NN=NM+M
        CS=C*C*KD
        IP=1
        IF (N-M.EQ.2*INT((N-M)/2)) IP=0
        K=0
        DO 10 I=1,NN+3
           IF (IP.EQ.0) K=-2*(I-1)
           IF (IP.EQ.1) K=-(2*I-3)
           GK0=2.0D0*M+K
           GK1=(M+K)*(M+K+1.0D0)
           GK2=2.0D0*(M+K)-1.0D0
           GK3=2.0D0*(M+K)+3.0D0
           U(I)=GK0*(GK0-1.0D0)*CS/(GK2*(GK2+2.0D0))
           V(I)=GK1-CV+(2.0D0*(GK1-M*M)-1.0D0)*CS/(GK2*GK3)
10         W(I)=(K+1.0D0)*(K+2.0D0)*CS/((GK2+2.0D0)*GK3)
        DO 20 K=1,M
           T=V(M+1)
           DO 15 L=0,M-K-1
15            T=V(M-L)-W(M-L+1)*U(M-L)/T
20         RK(K)=-U(K)/T
        R=1.0D0
        DO 25 K=1,M
           R=R*RK(K)
25         DN(K)=DF(1)*R
        TP(NN)=V(NN+1)
        DO 30 K=NN-1,M+1,-1
           TP(K)=V(K+1)-W(K+2)*U(K+1)/TP(K+1)
           IF (K.GT.M+1) RK(K)=-U(K)/TP(K)
30      CONTINUE
        IF (M.EQ.0) DNP=DF(1)
        IF (M.NE.0) DNP=DN(M)
        DN(M+1)=(-1)**IP*DNP*CS/((2.0*M-1.0)*(2.0*M+1.0-4.0*IP)
     &          *TP(M+1))
        DO 35 K=M+2,NN
35         DN(K)=RK(K)*DN(K-1)
        R1=1.0D0
        DO 40 J=1,(N+M+IP)/2
40         R1=R1*(J+0.5D0*(N+M+IP))
        NM1=(N-M)/2
        R=1.0D0
        DO 45 J=1,2*M+IP
45         R=R*J
        SU0=R*DF(1)
        SW=0.0D0
        DO 50 K=2,NM
           R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
           SU0=SU0+R*DF(K)
           IF (K.GT.NM1.AND.DABS((SU0-SW)/SU0).LT.1.0D-14) GO TO 55
50         SW=SU0
55      IF (KD.EQ.1) GOTO 70
        R2=1.0D0
        DO 60 J=1,M
60         R2=2.0D0*C*R2*J
        R3=1.0D0
        DO 65 J=1,(N-M-IP)/2
65         R3=R3*J
        SA0=(2.0*(M+IP)+1.0)*R1/(2.0**N*C**IP*R2*R3*DF(1))
        CK1=SA0*SU0
        IF (KD.EQ.-1) RETURN
70      R4=1.0D0
        DO 75 J=1,(N-M-IP)/2
75         R4=4.0D0*R4*J
        R5=1.0D0
        DO 80 J=1,M
80         R5=R5*(J+M)/C
        G0=DN(M)
        IF (M.EQ.0) G0=DF(1)
        SB0=(IP+1.0)*C**(IP+1)/(2.0*IP*(M-2.0)+1.0)/(2.0*M-1.0)
        CK2=(-1)**IP*SB0*R4*R5*G0/R1*SU0
        RETURN
        END



C       **********************************

        SUBROUTINE LAGZO(N,X,W)
C
C       =========================================================
C       Purpose : Compute the zeros of Laguerre polynomial Ln(x)
C                 in the interval [0,∞], and the corresponding
C                 weighting coefficients for Gauss-Laguerre
C                 integration
C       Input :   n    --- Order of the Laguerre polynomial
C                 X(n) --- Zeros of the Laguerre polynomial
C                 W(n) --- Corresponding weighting coefficients
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(N),W(N)
        HN=1.0D0/N
        PF=0.0D0
        PD=0.0D0
        DO 35 NR=1,N
           Z=HN
           IF (NR.GT.1) Z=X(NR-1)+HN*NR**1.27
           IT=0
10         IT=IT+1
           Z0=Z
           P=1.0D0
           DO 15 I=1,NR-1
15            P=P*(Z-X(I))
           F0=1.0D0
           F1=1.0D0-Z
           DO 20 K=2,N
              PF=((2.0D0*K-1.0D0-Z)*F1-(K-1.0D0)*F0)/K
              PD=K/Z*(PF-F1)
              F0=F1
20            F1=PF
           FD=PF/P
           Q=0.0D0
           DO 30 I=1,NR-1
              WP=1.0D0
              DO 25 J=1,NR-1
                 IF (J.EQ.I) GO TO 25
                 WP=WP*(Z-X(J))
25            CONTINUE
              Q=Q+WP
30         CONTINUE
           GD=(PD-Q*FD)/P
           Z=Z-FD/GD
           IF (IT.LE.40.AND.DABS((Z-Z0)/Z).GT.1.0D-15) GO TO 10
           X(NR)=Z
           W(NR)=1.0D0/(Z*PD*PD)
35      CONTINUE
        RETURN
        END

C       **********************************

        SUBROUTINE VVLA(VA,X,PV)
C
C       ===================================================
C       Purpose: Compute parabolic cylinder function Vv(x)
C                for large argument
C       Input:   x  --- Argument
C                va --- Order
C       Output:  PV --- Vv(x)
C       Routines called:
C             (1) DVLA for computing Dv(x) for large |x|
C             (2) GAMMA2 for computing Г(x)
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EPS=1.0D-12
        QE=DEXP(0.25*X*X)
        A0=DABS(X)**(-VA-1.0D0)*DSQRT(2.0D0/PI)*QE
        R=1.0D0
        PV=1.0D0
        DO 10 K=1,18
           R=0.5D0*R*(2.0*K+VA-1.0)*(2.0*K+VA)/(K*X*X)
           PV=PV+R
           IF (DABS(R/PV).LT.EPS) GO TO 15
10      CONTINUE
15      PV=A0*PV
        IF (X.LT.0.0D0) THEN
           X1=-X
           CALL DVLA(VA,X1,PDL)
           CALL GAMMA2(-VA,GL)
           DSL=DSIN(PI*VA)*DSIN(PI*VA)
           PV=DSL*GL/PI*PDL-DCOS(PI*VA)*PV
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE CJYVA(V,Z,VM,CBJ,CDJ,CBY,CDY)
C
C       ===========================================================
C       Purpose: Compute Bessel functions Jv(z), Yv(z) and their
C                derivatives for a complex argument
C       Input :  z --- Complex argument
C                v --- Order of Jv(z) and Yv(z)
C                      ( v = n+v0, n = 0,1,2,..., 0 ≤ v0 < 1 )
C       Output:  CBJ(n) --- Jn+v0(z)
C                CDJ(n) --- Jn+v0'(z)
C                CBY(n) --- Yn+v0(z)
C                CDY(n) --- Yn+v0'(z)
C                VM --- Highest order computed
C       Routines called:
C            (1) GAMMA2 for computing the gamma function
C            (2) MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,G,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:*),CDJ(0:*),CBY(0:*),CDY(0:*)
        PI=3.141592653589793D0
        RP2=.63661977236758D0
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z1=Z
        Z2=Z*Z
        N=INT(V)
        V0=V-N
        PV0=PI*V0
        PV1=PI*(1.0D0+V0)
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBJ(K)=(0.0D0,0.0D0)
              CDJ(K)=(0.0D0,0.0D0)
              CBY(K)=-(1.0D+300,0.0D0)
10            CDY(K)=(1.0D+300,0.0D0)
           IF (V0.EQ.0.0) THEN
              CBJ(0)=(1.0D0,0.0D0)
              CDJ(1)=(0.5D0,0.0D0)
           ELSE
              CDJ(0)=(1.0D+300,0.0D0)
           ENDIF
           VM=V
           RETURN
        ENDIF
        LB0=0.0D0
        IF (DBLE(Z).LT.0.0) Z1=-Z
        IF (A0.LE.12.0) THEN
           DO 25 L=0,1
              VL=V0+L
              CJVL=(1.0D0,0.0D0)
              CR=(1.0D0,0.0D0)
              DO 15 K=1,40
                 CR=-0.25D0*CR*Z2/(K*(K+VL))
                 CJVL=CJVL+CR
                 IF (CDABS(CR).LT.CDABS(CJVL)*1.0D-15) GO TO 20
15            CONTINUE
20            VG=1.0D0+VL
              CALL GAMMA2(VG,GA)
              CA=(0.5D0*Z1)**VL/GA
              IF (L.EQ.0) CJV0=CJVL*CA
              IF (L.EQ.1) CJV1=CJVL*CA
25         CONTINUE
        ELSE
           K0=11
           IF (A0.GE.35.0) K0=10
           IF (A0.GE.50.0) K0=8
           DO 40 J=0,1
              VV=4.0D0*(J+V0)*(J+V0)
              CPZ=(1.0D0,0.0D0)
              CRP=(1.0D0,0.0D0)
              DO 30 K=1,K0
                 CRP=-0.78125D-2*CRP*(VV-(4.0*K-3.0)**2.0)*(VV-
     &               (4.0*K-1.0)**2.0)/(K*(2.0*K-1.0)*Z2)
30               CPZ=CPZ+CRP
              CQZ=(1.0D0,0.0D0)
              CRQ=(1.0D0,0.0D0)
              DO 35 K=1,K0
                 CRQ=-0.78125D-2*CRQ*(VV-(4.0*K-1.0)**2.0)*(VV-
     &               (4.0*K+1.0)**2.0)/(K*(2.0*K+1.0)*Z2)
35               CQZ=CQZ+CRQ
              CQZ=0.125D0*(VV-1.0)*CQZ/Z1
              ZK=Z1-(0.5D0*(J+V0)+0.25D0)*PI
              CA0=CDSQRT(RP2/Z1)
              CCK=CDCOS(ZK)
              CSK=CDSIN(ZK)
              IF (J.EQ.0) THEN
                 CJV0=CA0*(CPZ*CCK-CQZ*CSK)
                 CYV0=CA0*(CPZ*CSK+CQZ*CCK)
              ELSE IF (J.EQ.1) THEN
                 CJV1=CA0*(CPZ*CCK-CQZ*CSK)
                 CYV1=CA0*(CPZ*CSK+CQZ*CCK)
              ENDIF
40         CONTINUE
        ENDIF
        IF (A0.LE.12.0) THEN
           IF (V0.NE.0.0) THEN
              DO 55 L=0,1
                 VL=V0+L
                 CJVL=(1.0D0,0.0D0)
                 CR=(1.0D0,0.0D0)
                 DO 45 K=1,40
                    CR=-0.25D0*CR*Z2/(K*(K-VL))
                    CJVL=CJVL+CR
                    IF (CDABS(CR).LT.CDABS(CJVL)*1.0D-15) GO TO 50
45               CONTINUE
50               VG=1.0D0-VL
                 CALL GAMMA2(VG,GB)
                 CB=(2.0D0/Z1)**VL/GB
                 IF (L.EQ.0) CJU0=CJVL*CB
                 IF (L.EQ.1) CJU1=CJVL*CB
55            CONTINUE
              CYV0=(CJV0*DCOS(PV0)-CJU0)/DSIN(PV0)
              CYV1=(CJV1*DCOS(PV1)-CJU1)/DSIN(PV1)
           ELSE
              CEC=CDLOG(Z1/2.0D0)+.5772156649015329D0
              CS0=(0.0D0,0.0D0)
              W0=0.0D0
              CR0=(1.0D0,0.0D0)
              DO 60 K=1,30
                 W0=W0+1.0D0/K
                 CR0=-0.25D0*CR0/(K*K)*Z2
60               CS0=CS0+CR0*W0
              CYV0=RP2*(CEC*CJV0-CS0)
              CS1=(1.0D0,0.0D0)
              W1=0.0D0
              CR1=(1.0D0,0.0D0)
              DO 65 K=1,30
                 W1=W1+1.0D0/K
                 CR1=-0.25D0*CR1/(K*(K+1))*Z2
65               CS1=CS1+CR1*(2.0D0*W1+1.0D0/(K+1.0D0))
              CYV1=RP2*(CEC*CJV1-1.0D0/Z1-0.25D0*Z1*CS1)
           ENDIF
        ENDIF
        IF (DBLE(Z).LT.0.0D0) THEN
           CFAC0=CDEXP(PV0*CI)
           CFAC1=CDEXP(PV1*CI)
           IF (DIMAG(Z).LT.0.0D0) THEN
              CYV0=CFAC0*CYV0-2.0D0*CI*DCOS(PV0)*CJV0
              CYV1=CFAC1*CYV1-2.0D0*CI*DCOS(PV1)*CJV1
              CJV0=CJV0/CFAC0
              CJV1=CJV1/CFAC1
           ELSE IF (DIMAG(Z).GT.0.0D0) THEN
              CYV0=CYV0/CFAC0+2.0D0*CI*DCOS(PV0)*CJV0
              CYV1=CYV1/CFAC1+2.0D0*CI*DCOS(PV1)*CJV1
              CJV0=CFAC0*CJV0
              CJV1=CFAC1*CJV1
           ENDIF
        ENDIF
        CBJ(0)=CJV0
        CBJ(1)=CJV1
        IF (N.GE.2.AND.N.LE.INT(0.25*A0)) THEN
           CF0=CJV0
           CF1=CJV1
           DO 70 K=2,N
              CF=2.0D0*(K+V0-1.0D0)/Z*CF1-CF0
              CBJ(K)=CF
              CF0=CF1
70            CF1=CF
        ELSE IF (N.GE.2) THEN
           M=MSTA1(A0,200)
           IF (M.LT.N) THEN
              N=M
           ELSE
              M=MSTA2(A0,N,15)
           ENDIF
           CF2=(0.0D0,0.0D0)
           CF1=(1.0D-100,0.0D0)
           DO 75 K=M,0,-1
              CF=2.0D0*(V0+K+1.0D0)/Z*CF1-CF2
              IF (K.LE.N) CBJ(K)=CF
              CF2=CF1
75            CF1=CF
           IF (CDABS(CJV0).GT.CDABS(CJV1)) CS=CJV0/CF
           IF (CDABS(CJV0).LE.CDABS(CJV1)) CS=CJV1/CF2
           DO 80 K=0,N
80            CBJ(K)=CS*CBJ(K)
        ENDIF
        CDJ(0)=V0/Z*CBJ(0)-CBJ(1)
        DO 85 K=1,N
85         CDJ(K)=-(K+V0)/Z*CBJ(K)+CBJ(K-1)
        CBY(0)=CYV0
        CBY(1)=CYV1
        YA0=CDABS(CYV0)
        LB=0
        CG0=CYV0
        CG1=CYV1
        DO 90 K=2,N
           CYK=2.0D0*(V0+K-1.0D0)/Z*CG1-CG0
           IF (CDABS(CYK).GT.1.0D+290) GO TO 90
           YAK=CDABS(CYK)
           YA1=CDABS(CG0)
           IF (YAK.LT.YA0.AND.YAK.LT.YA1) LB=K
           CBY(K)=CYK
           CG0=CG1
           CG1=CYK
90      CONTINUE
        IF (LB.LE.4.OR.DIMAG(Z).EQ.0.0D0) GO TO 125
95      IF (LB.EQ.LB0) GO TO 125
        CH2=(1.0D0,0.0D0)
        CH1=(0.0D0,0.0D0)
        LB0=LB
        DO 100 K=LB,1,-1
           CH0=2.0D0*(K+V0)/Z*CH1-CH2
           CH2=CH1
100        CH1=CH0
        CP12=CH0
        CP22=CH2
        CH2=(0.0D0,0.0D0)
        CH1=(1.0D0,0.0D0)
        DO 105 K=LB,1,-1
           CH0=2.0D0*(K+V0)/Z*CH1-CH2
           CH2=CH1
105        CH1=CH0
        CP11=CH0
        CP21=CH2
        IF (LB.EQ.N) CBJ(LB+1)=2.0D0*(LB+V0)/Z*CBJ(LB)-CBJ(LB-1)
        IF (CDABS(CBJ(0)).GT.CDABS(CBJ(1))) THEN
           CBY(LB+1)=(CBJ(LB+1)*CYV0-2.0D0*CP11/(PI*Z))/CBJ(0)
           CBY(LB)=(CBJ(LB)*CYV0+2.0D0*CP12/(PI*Z))/CBJ(0)
        ELSE
           CBY(LB+1)=(CBJ(LB+1)*CYV1-2.0D0*CP21/(PI*Z))/CBJ(1)
           CBY(LB)=(CBJ(LB)*CYV1+2.0D0*CP22/(PI*Z))/CBJ(1)
        ENDIF
        CYL2=CBY(LB+1)
        CYL1=CBY(LB)
        DO 110 K=LB-1,0,-1
           CYLK=2.0D0*(K+V0+1.0D0)/Z*CYL1-CYL2
           CBY(K)=CYLK
           CYL2=CYL1
110        CYL1=CYLK
        CYL1=CBY(LB)
        CYL2=CBY(LB+1)
        DO 115 K=LB+1,N-1
           CYLK=2.0D0*(K+V0)/Z*CYL2-CYL1
           CBY(K+1)=CYLK
           CYL1=CYL2
115        CYL2=CYLK
        DO 120 K=2,N
           WA=CDABS(CBY(K))
           IF (WA.LT.CDABS(CBY(K-1))) LB=K
120     CONTINUE
        GO TO 95
125     CDY(0)=V0/Z*CBY(0)-CBY(1)
        DO 130 K=1,N
130        CDY(K)=CBY(K-1)-(K+V0)/Z*CBY(K)
        VM=N+V0
        RETURN
        END



C       **********************************

        SUBROUTINE CJYVB(V,Z,VM,CBJ,CDJ,CBY,CDY)
C
C       ===========================================================
C       Purpose: Compute Bessel functions Jv(z), Yv(z) and their
C                derivatives for a complex argument
C       Input :  z --- Complex argument
C                v --- Order of Jv(z) and Yv(z)
C                      ( v = n+v0, n = 0,1,2,..., 0 ≤ v0 < 1 )
C       Output:  CBJ(n) --- Jn+v0(z)
C                CDJ(n) --- Jn+v0'(z)
C                CBY(n) --- Yn+v0(z)
C                CDY(n) --- Yn+v0'(z)
C                VM --- Highest order computed
C       Routines called:
C            (1) GAMMA2 for computing the gamma function
C            (2) MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,G,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:*),CDJ(0:*),CBY(0:*),CDY(0:*)
        PI=3.141592653589793D0
        RP2=.63661977236758D0
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z1=Z
        Z2=Z*Z
        N=INT(V)
        V0=V-N
        PV0=PI*V0
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBJ(K)=(0.0D0,0.0D0)
              CDJ(K)=(0.0D0,0.0D0)
              CBY(K)=-(1.0D+300,0.0D0)
10            CDY(K)=(1.0D+300,0.0D0)
           IF (V0.EQ.0.0) THEN
              CBJ(0)=(1.0D0,0.0D0)
              CDJ(1)=(0.5D0,0.0D0)
           ELSE
              CDJ(0)=(1.0D+300,0.0D0)
           ENDIF
           VM=V
           RETURN
        ENDIF
        IF (DBLE(Z).LT.0.0D0) Z1=-Z
        IF (A0.LE.12.0) THEN
           CJV0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 15 K=1,40
              CR=-0.25D0*CR*Z2/(K*(K+V0))
              CJV0=CJV0+CR
              IF (CDABS(CR).LT.CDABS(CJV0)*1.0D-15) GO TO 20
15         CONTINUE
20         VG=1.0D0+V0
           CALL GAMMA2(VG,GA)
           CA=(0.5D0*Z1)**V0/GA
           CJV0=CJV0*CA
        ELSE
           K0=11
           IF (A0.GE.35.0) K0=10
           IF (A0.GE.50.0) K0=8
           VV=4.0D0*V0*V0
           CPZ=(1.0D0,0.0D0)
           CRP=(1.0D0,0.0D0)
           DO 25 K=1,K0
              CRP=-0.78125D-2*CRP*(VV-(4.0*K-3.0)**2.0)*(VV-
     &            (4.0*K-1.0)**2.0)/(K*(2.0*K-1.0)*Z2)
25            CPZ=CPZ+CRP
           CQZ=(1.0D0,0.0D0)
           CRQ=(1.0D0,0.0D0)
           DO 30 K=1,K0
              CRQ=-0.78125D-2*CRQ*(VV-(4.0*K-1.0)**2.0)*(VV-
     &            (4.0*K+1.0)**2.0)/(K*(2.0*K+1.0)*Z2)
30            CQZ=CQZ+CRQ
           CQZ=0.125D0*(VV-1.0)*CQZ/Z1
           ZK=Z1-(0.5D0*V0+0.25D0)*PI
           CA0=CDSQRT(RP2/Z1)
           CCK=CDCOS(ZK)
           CSK=CDSIN(ZK)
           CJV0=CA0*(CPZ*CCK-CQZ*CSK)
           CYV0=CA0*(CPZ*CSK+CQZ*CCK)
        ENDIF
        IF (A0.LE.12.0) THEN
           IF (V0.NE.0.0) THEN
              CJVN=(1.0D0,0.0D0)
              CR=(1.0D0,0.0D0)
              DO 35 K=1,40
                 CR=-0.25D0*CR*Z2/(K*(K-V0))
                 CJVN=CJVN+CR
                 IF (CDABS(CR).LT.CDABS(CJVN)*1.0D-15) GO TO 40
35            CONTINUE
40            VG=1.0D0-V0
              CALL GAMMA2(VG,GB)
              CB=(2.0D0/Z1)**V0/GB
              CJU0=CJVN*CB
              CYV0=(CJV0*DCOS(PV0)-CJU0)/DSIN(PV0)
           ELSE
              CEC=CDLOG(Z1/2.0D0)+.5772156649015329D0
              CS0=(0.0D0,0.0D0)
              W0=0.0D0
              CR0=(1.0D0,0.0D0)
              DO 45 K=1,30
                 W0=W0+1.0D0/K
                 CR0=-0.25D0*CR0/(K*K)*Z2
45               CS0=CS0+CR0*W0
              CYV0=RP2*(CEC*CJV0-CS0)
           ENDIF
        ENDIF
        IF (N.EQ.0) N=1
        M=MSTA1(A0,200)
        IF (M.LT.N) THEN
           N=M
        ELSE
           M=MSTA2(A0,N,15)
        ENDIF
        CF2=(0.0D0,0.0D0)
        CF1=(1.0D-100,0.0D0)
        DO 50 K=M,0,-1
           CF=2.0D0*(V0+K+1.0D0)/Z1*CF1-CF2
           IF (K.LE.N) CBJ(K)=CF
           CF2=CF1
50         CF1=CF
        CS=CJV0/CF
        DO 55 K=0,N
55         CBJ(K)=CS*CBJ(K)
        IF (DBLE(Z).LT.0.0D0) THEN
           CFAC0=CDEXP(PV0*CI)
           IF (DIMAG(Z).LT.0.0D0) THEN
              CYV0=CFAC0*CYV0-2.0D0*CI*DCOS(PV0)*CJV0
           ELSE IF (DIMAG(Z).GT.0.0D0) THEN
              CYV0=CYV0/CFAC0+2.0D0*CI*DCOS(PV0)*CJV0
           ENDIF
           DO 60 K=0,N
              IF (DIMAG(Z).LT.0.0D0) THEN
                 CBJ(K)=CDEXP(-PI*(K+V0)*CI)*CBJ(K)
              ELSE IF (DIMAG(Z).GT.0.0D0) THEN
                 CBJ(K)=CDEXP(PI*(K+V0)*CI)*CBJ(K)
              ENDIF
60         CONTINUE
           Z1=Z1
        ENDIF
        CBY(0)=CYV0
        DO 65 K=1,N
           CYY=(CBJ(K)*CBY(K-1)-2.0D0/(PI*Z))/CBJ(K-1)
           CBY(K)=CYY
65      CONTINUE
        CDJ(0)=V0/Z*CBJ(0)-CBJ(1)
        DO 70 K=1,N
70         CDJ(K)=-(K+V0)/Z*CBJ(K)+CBJ(K-1)
        CDY(0)=V0/Z*CBY(0)-CBY(1)
        DO 75 K=1,N
75         CDY(K)=CBY(K-1)-(K+V0)/Z*CBY(K)
        VM=N+V0
        RETURN
        END



C       **********************************

        SUBROUTINE JY01A(X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
C
C       =======================================================
C       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
C                Y1(x), and their derivatives
C       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ≥ 0 )
C       Output:  BJ0 --- J0(x)
C                DJ0 --- J0'(x)
C                BJ1 --- J1(x)
C                DJ1 --- J1'(x)
C                BY0 --- Y0(x)
C                DY0 --- Y0'(x)
C                BY1 --- Y1(x)
C                DY1 --- Y1'(x)
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(12),B(12),A1(12),B1(12)
        PI=3.141592653589793D0
        RP2=0.63661977236758D0
        X2=X*X
        IF (X.EQ.0.0D0) THEN
           BJ0=1.0D0
           BJ1=0.0D0
           DJ0=0.0D0
           DJ1=0.5D0
           BY0=-1.0D+300
           BY1=-1.0D+300
           DY0=1.0D+300
           DY1=1.0D+300
           RETURN
        ENDIF
        IF (X.LE.12.0D0) THEN
           BJ0=1.0D0
           R=1.0D0
           DO 5 K=1,30
              R=-0.25D0*R*X2/(K*K)
              BJ0=BJ0+R
              IF (DABS(R).LT.DABS(BJ0)*1.0D-15) GO TO 10
5          CONTINUE
10         BJ1=1.0D0
           R=1.0D0
           DO 15 K=1,30
              R=-0.25D0*R*X2/(K*(K+1.0D0))
              BJ1=BJ1+R
              IF (DABS(R).LT.DABS(BJ1)*1.0D-15) GO TO 20
15         CONTINUE
20         BJ1=0.5D0*X*BJ1
           EC=DLOG(X/2.0D0)+0.5772156649015329D0
           CS0=0.0D0
           W0=0.0D0
           R0=1.0D0
           DO 25 K=1,30
              W0=W0+1.0D0/K
              R0=-0.25D0*R0/(K*K)*X2
              R=R0*W0
              CS0=CS0+R
              IF (DABS(R).LT.DABS(CS0)*1.0D-15) GO TO 30
25         CONTINUE
30         BY0=RP2*(EC*BJ0-CS0)
           CS1=1.0D0
           W1=0.0D0
           R1=1.0D0
           DO 35 K=1,30
              W1=W1+1.0D0/K
              R1=-0.25D0*R1/(K*(K+1))*X2
              R=R1*(2.0D0*W1+1.0D0/(K+1.0D0))
              CS1=CS1+R
              IF (DABS(R).LT.DABS(CS1)*1.0D-15) GO TO 40
35         CONTINUE
40         BY1=RP2*(EC*BJ1-1.0D0/X-0.25D0*X*CS1)
        ELSE
           DATA A/-.7031250000000000D-01,.1121520996093750D+00,
     &            -.5725014209747314D+00,.6074042001273483D+01,
     &            -.1100171402692467D+03,.3038090510922384D+04,
     &            -.1188384262567832D+06,.6252951493434797D+07,
     &            -.4259392165047669D+09,.3646840080706556D+11,
     &            -.3833534661393944D+13,.4854014686852901D+15/
           DATA B/ .7324218750000000D-01,-.2271080017089844D+00,
     &             .1727727502584457D+01,-.2438052969955606D+02,
     &             .5513358961220206D+03,-.1825775547429318D+05,
     &             .8328593040162893D+06,-.5006958953198893D+08,
     &             .3836255180230433D+10,-.3649010818849833D+12,
     &             .4218971570284096D+14,-.5827244631566907D+16/
           DATA A1/.1171875000000000D+00,-.1441955566406250D+00,
     &             .6765925884246826D+00,-.6883914268109947D+01,
     &             .1215978918765359D+03,-.3302272294480852D+04,
     &             .1276412726461746D+06,-.6656367718817688D+07,
     &             .4502786003050393D+09,-.3833857520742790D+11,
     &             .4011838599133198D+13,-.5060568503314727D+15/
           DATA B1/-.1025390625000000D+00,.2775764465332031D+00,
     &             -.1993531733751297D+01,.2724882731126854D+02,
     &             -.6038440767050702D+03,.1971837591223663D+05,
     &             -.8902978767070678D+06,.5310411010968522D+08,
     &             -.4043620325107754D+10,.3827011346598605D+12,
     &             -.4406481417852278D+14,.6065091351222699D+16/
           K0=12
           IF (X.GE.35.0) K0=10
           IF (X.GE.50.0) K0=8
           T1=X-0.25D0*PI
           P0=1.0D0
           Q0=-0.125D0/X
           DO 45 K=1,K0
              P0=P0+A(K)*X**(-2*K)
45            Q0=Q0+B(K)*X**(-2*K-1)
           CU=DSQRT(RP2/X)
           BJ0=CU*(P0*DCOS(T1)-Q0*DSIN(T1))
           BY0=CU*(P0*DSIN(T1)+Q0*DCOS(T1))
           T2=X-0.75D0*PI
           P1=1.0D0
           Q1=0.375D0/X
           DO 50 K=1,K0
              P1=P1+A1(K)*X**(-2*K)
50            Q1=Q1+B1(K)*X**(-2*K-1)
           CU=DSQRT(RP2/X)
           BJ1=CU*(P1*DCOS(T2)-Q1*DSIN(T2))
           BY1=CU*(P1*DSIN(T2)+Q1*DCOS(T2))
        ENDIF
        DJ0=-BJ1
        DJ1=BJ0-BJ1/X
        DY0=-BY1
        DY1=BY0-BY1/X
        RETURN
        END

C       **********************************

        SUBROUTINE INCOG(A,X,GIN,GIM,GIP)
C
C       ===================================================
C       Purpose: Compute the incomplete gamma function
C                r(a,x), Г(a,x) and P(a,x)
C       Input :  a   --- Parameter ( a ≤ 170 )
C                x   --- Argument
C       Output:  GIN --- r(a,x)
C                GIM --- Г(a,x)
C                GIP --- P(a,x)
C       Routine called: GAMMA2 for computing Г(x)
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XAM=-X+A*DLOG(X)
        IF (XAM.GT.700.0.OR.A.GT.170.0) THEN
           WRITE(*,*)'a and/or x too large'
           STOP
        ENDIF
        IF (X.EQ.0.0) THEN
           GIN=0.0
           CALL GAMMA2(A,GA)
           GIM=GA
           GIP=0.0
        ELSE IF (X.LE.1.0+A) THEN
           S=1.0D0/A
           R=S
           DO 10 K=1,60
              R=R*X/(A+K)
              S=S+R
              IF (DABS(R/S).LT.1.0D-15) GO TO 15
10         CONTINUE
15         GIN=DEXP(XAM)*S
           CALL GAMMA2(A,GA)
           GIP=GIN/GA
           GIM=GA-GIN
        ELSE IF (X.GT.1.0+A) THEN
           T0=0.0D0
           DO 20 K=60,1,-1
              T0=(K-A)/(1.0D0+K/(X+T0))
20         CONTINUE
           GIM=DEXP(XAM)/(X+T0)
           CALL GAMMA2(A,GA)
           GIN=GA-GIM
           GIP=1.0D0-GIM/GA
        ENDIF
        END



C       **********************************

        SUBROUTINE ITIKB(X,TI,TK)
C
C       =======================================================
C       Purpose: Integrate Bessel functions I0(t) and K0(t)
C                with respect to t from 0 to x
C       Input :  x  --- Upper limit of the integral ( x ≥ 0 )
C       Output:  TI --- Integration of I0(t) from 0 to x
C                TK --- Integration of K0(t) from 0 to x
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           TI=0.0D0
        ELSE IF (X.LT.5.0D0) THEN
           T1=X/5.0D0
           T=T1*T1
           TI=((((((((.59434D-3*T+.4500642D-2)*T
     &        +.044686921D0)*T+.300704878D0)*T+1.471860153D0)
     &        *T+4.844024624D0)*T+9.765629849D0)*T
     &        +10.416666367D0)*T+5.0D0)*T1
        ELSE IF (X.GE.5.0.AND.X.LE.8.0D0) THEN
           T=5.0D0/X
           TI=(((-.015166D0*T-.0202292D0)*T+.1294122D0)*T
     &        -.0302912D0)*T+.4161224D0
           TI=TI*DEXP(X)/DSQRT(X)
        ELSE
           T=8.0D0/X
           TI=(((((-.0073995D0*T+.017744D0)*T-.0114858D0)*T
     &        +.55956D-2)*T+.59191D-2)*T+.0311734D0)*T
     &        +.3989423D0
           TI=TI*DEXP(X)/DSQRT(X)
        ENDIF
        IF (X.EQ.0.0D0) THEN
           TK=0.0D0
        ELSE IF (X.LE.2.0D0) THEN
           T1=X/2.0D0
           T=T1*T1
           TK=((((((.116D-5*T+.2069D-4)*T+.62664D-3)*T
     &        +.01110118D0)*T+.11227902D0)*T+.50407836D0)*T
     &        +.84556868D0)*T1
              TK=TK-DLOG(X/2.0D0)*TI
        ELSE IF (X.GT.2.0.AND.X.LE.4.0D0) THEN
           T=2.0D0/X
           TK=(((.0160395D0*T-.0781715D0)*T+.185984D0)*T
     &        -.3584641D0)*T+1.2494934D0
           TK=PI/2.0D0-TK*DEXP(-X)/DSQRT(X)
        ELSE IF (X.GT.4.0.AND.X.LE.7.0D0) THEN
           T=4.0D0/X
           TK=(((((.37128D-2*T-.0158449D0)*T+.0320504D0)*T
     &        -.0481455D0)*T+.0787284D0)*T-.1958273D0)*T
     &        +1.2533141D0
           TK=PI/2.0D0-TK*DEXP(-X)/DSQRT(X)
        ELSE
           T=7.0D0/X
           TK=(((((.33934D-3*T-.163271D-2)*T+.417454D-2)*T
     &        -.933944D-2)*T+.02576646D0)*T-.11190289D0)*T
     &        +1.25331414D0
           TK=PI/2.0D0-TK*DEXP(-X)/DSQRT(X)
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE ITIKA(X,TI,TK)
C
C       =======================================================
C       Purpose: Integrate modified Bessel functions I0(t) and
C                K0(t) with respect to t from 0 to x
C       Input :  x  --- Upper limit of the integral  ( x ≥ 0 )
C       Output:  TI --- Integration of I0(t) from 0 to x
C                TK --- Integration of K0(t) from 0 to x
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(10)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        DATA A/.625D0,1.0078125D0,
     &       2.5927734375D0,9.1868591308594D0,
     &       4.1567974090576D+1,2.2919635891914D+2,
     &       1.491504060477D+3,1.1192354495579D+4,
     &       9.515939374212D+4,9.0412425769041D+5/
        IF (X.EQ.0.0D0) THEN
           TI=0.0D0
           TK=0.0D0
           RETURN
        ELSE IF (X.LT.20.0D0) THEN
           X2=X*X
           TI=1.0D0
           R=1.0D0
           DO 10 K=1,50
              R=.25D0*R*(2*K-1.0D0)/(2*K+1.0D0)/(K*K)*X2
              TI=TI+R
              IF (DABS(R/TI).LT.1.0D-12) GO TO 15
10         CONTINUE
15         TI=TI*X
        ELSE
           X2=0.0D0
           TI=1.0D0
           R=1.0D0
           DO 20 K=1,10
              R=R/X
20            TI=TI+A(K)*R
           RC1=1.0D0/DSQRT(2.0D0*PI*X)
           TI=RC1*DEXP(X)*TI
        ENDIF
        IF (X.LT.12.0D0) THEN
           E0=EL+DLOG(X/2.0D0)
           B1=1.0D0-E0
           B2=0.0D0
           RS=0.0D0
           R=1.0D0
           TW=0.0D0
           DO 25 K=1,50
              R=.25D0*R*(2*K-1.0D0)/(2*K+1.0D0)/(K*K)*X2
              B1=B1+R*(1.0D0/(2*K+1)-E0)
              RS=RS+1.0D0/K
              B2=B2+R*RS
              TK=B1+B2
              IF (DABS((TK-TW)/TK).LT.1.0D-12) GO TO 30
25            TW=TK
30         TK=TK*X
        ELSE
           TK=1.0D0
           R=1.0D0
           DO 35 K=1,10
              R=-R/X
35            TK=TK+A(K)*R
           RC2=DSQRT(PI/(2.0D0*X))
           TK=PI/2.0D0-RC2*TK*DEXP(-X)
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE JYV(V,X,VM,BJ,DJ,BY,DY)
C
C       =======================================================
C       Purpose: Compute Bessel functions Jv(x) and Yv(x)
C                and their derivatives
C       Input :  x --- Argument of Jv(x) and Yv(x)
C                v --- Order of Jv(x) and Yv(x)
C                      ( v = n+v0, 0 ≤ v0 < 1, n = 0,1,2,... )
C       Output:  BJ(n) --- Jn+v0(x)
C                DJ(n) --- Jn+v0'(x)
C                BY(n) --- Yn+v0(x)
C                DY(n) --- Yn+v0'(x)
C                VM --- Highest order computed
C       Routines called:
C            (1) GAMMA2 for computing gamma function
C            (2) MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BJ(0:*),DJ(0:*),BY(0:*),DY(0:*)
        EL=.5772156649015329D0
        PI=3.141592653589793D0
        RP2=.63661977236758D0
        X2=X*X
        N=INT(V)
        V0=V-N
        IF (X.LT.1.0D-100) THEN
           DO 10 K=0,N
              BJ(K)=0.0D0
              DJ(K)=0.0D0
              BY(K)=-1.0D+300
10            DY(K)=1.0D+300
           IF (V0.EQ.0.0) THEN
              BJ(0)=1.0D0
              DJ(1)=0.5D0
           ELSE
              DJ(0)=1.0D+300
           ENDIF
           VM=V
           RETURN
        ENDIF
        BJV0=0.0D0
        BJV1=0.0D0
        BYV0=0.0D0
        BYV1=0.0D0
        IF (X.LE.12.0) THEN
           DO 25 L=0,1
              VL=V0+L
              BJVL=1.0D0
              R=1.0D0
              DO 15 K=1,40
                 R=-0.25D0*R*X2/(K*(K+VL))
                 BJVL=BJVL+R
                 IF (DABS(R).LT.DABS(BJVL)*1.0D-15) GO TO 20
15            CONTINUE
20            VG=1.0D0+VL
              CALL GAMMA2(VG,GA)
              A=(0.5D0*X)**VL/GA
              IF (L.EQ.0) BJV0=BJVL*A
              IF (L.EQ.1) BJV1=BJVL*A
25         CONTINUE
        ELSE
           K0=11
           IF (X.GE.35.0) K0=10
           IF (X.GE.50.0) K0=8
           DO 40 J=0,1
              VV=4.0D0*(J+V0)*(J+V0)
              PX=1.0D0
              RP=1.0D0
              DO 30 K=1,K0
                 RP=-0.78125D-2*RP*(VV-(4.0*K-3.0)**2.0)*(VV-
     &              (4.0*K-1.0)**2.0)/(K*(2.0*K-1.0)*X2)
30               PX=PX+RP
              QX=1.0D0
              RQ=1.0D0
              DO 35 K=1,K0
                 RQ=-0.78125D-2*RQ*(VV-(4.0*K-1.0)**2.0)*(VV-
     &              (4.0*K+1.0)**2.0)/(K*(2.0*K+1.0)*X2)
35               QX=QX+RQ
              QX=0.125D0*(VV-1.0)*QX/X
              XK=X-(0.5D0*(J+V0)+0.25D0)*PI
              A0=DSQRT(RP2/X)
              CK=DCOS(XK)
              SK=DSIN(XK)
              IF (J.EQ.0) THEN
                 BJV0=A0*(PX*CK-QX*SK)
                 BYV0=A0*(PX*SK+QX*CK)
              ELSE IF (J.EQ.1) THEN
                 BJV1=A0*(PX*CK-QX*SK)
                 BYV1=A0*(PX*SK+QX*CK)
              ENDIF
40         CONTINUE
        ENDIF
        BJ(0)=BJV0
        BJ(1)=BJV1
        DJ(0)=V0/X*BJ(0)-BJ(1)
        DJ(1)=-(1.0D0+V0)/X*BJ(1)+BJ(0)
        IF (N.GE.2.AND.N.LE.INT(0.9*X)) THEN
           F0=BJV0
           F1=BJV1
           DO 45 K=2,N
              F=2.0D0*(K+V0-1.0D0)/X*F1-F0
              BJ(K)=F
              F0=F1
45            F1=F
        ELSE IF (N.GE.2) THEN
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              N=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F=0.0D0
           F2=0.0D0
           F1=1.0D-100
           DO 50 K=M,0,-1
              F=2.0D0*(V0+K+1.0D0)/X*F1-F2
              IF (K.LE.N) BJ(K)=F
              F2=F1
50            F1=F
           IF (DABS(BJV0).GT.DABS(BJV1)) THEN
               CS=BJV0/F
           ELSE
               CS=BJV1/F2
           ENDIF
           DO 55 K=0,N
55            BJ(K)=CS*BJ(K)
        ENDIF
        DO 60 K=2,N
60         DJ(K)=-(K+V0)/X*BJ(K)+BJ(K-1)
        IF (X.LE.12.0D0) THEN
           IF (V0.NE.0.0) THEN
              BJU0=0.0D0
              BJU1=0.0D0
              DO 75 L=0,1
                 VL=V0+L
                 BJVL=1.0D0
                 R=1.0D0
                 DO 65 K=1,40
                    R=-0.25D0*R*X2/(K*(K-VL))
                    BJVL=BJVL+R
                    IF (DABS(R).LT.DABS(BJVL)*1.0D-15) GO TO 70
65               CONTINUE
70               VG=1.0D0-VL
                 CALL GAMMA2(VG,GB)
                 B=(2.0D0/X)**VL/GB
                 IF (L.EQ.0) BJU0=BJVL*B
                 IF (L.EQ.1) BJU1=BJVL*B
75            CONTINUE
              PV0=PI*V0
              PV1=PI*(1.0D0+V0)
              BYV0=(BJV0*DCOS(PV0)-BJU0)/DSIN(PV0)
              BYV1=(BJV1*DCOS(PV1)-BJU1)/DSIN(PV1)
           ELSE
              EC=DLOG(X/2.0D0)+EL
              CS0=0.0D0
              W0=0.0D0
              R0=1.0D0
              DO 80 K=1,30
                 W0=W0+1.0D0/K
                 R0=-0.25D0*R0/(K*K)*X2
80               CS0=CS0+R0*W0
              BYV0=RP2*(EC*BJV0-CS0)
              CS1=1.0D0
              W1=0.0D0
              R1=1.0D0
              DO 85 K=1,30
                 W1=W1+1.0D0/K
                 R1=-0.25D0*R1/(K*(K+1))*X2
85               CS1=CS1+R1*(2.0D0*W1+1.0D0/(K+1.0D0))
              BYV1=RP2*(EC*BJV1-1.0D0/X-0.25D0*X*CS1)
           ENDIF
        ENDIF
        BY(0)=BYV0
        BY(1)=BYV1
        DO 90 K=2,N
           BYVK=2.0D0*(V0+K-1.0D0)/X*BYV1-BYV0
           BY(K)=BYVK
           BYV0=BYV1
90         BYV1=BYVK
        DY(0)=V0/X*BY(0)-BY(1)
        DO 95 K=1,N
95         DY(K)=-(K+V0)/X*BY(K)+BY(K-1)
        VM=N+V0
        RETURN
        END



C       **********************************

        SUBROUTINE JYNB(N,X,NM,BJ,DJ,BY,DY)
C
C       =====================================================
C       Purpose: Compute Bessel functions Jn(x), Yn(x) and
C                their derivatives
C       Input :  x --- Argument of Jn(x) and Yn(x) ( x ≥ 0 )
C                n --- Order of Jn(x) and Yn(x)
C       Output:  BJ(n) --- Jn(x)
C                DJ(n) --- Jn'(x)
C                BY(n) --- Yn(x)
C                DY(n) --- Yn'(x)
C                NM --- Highest order computed
C       Routines called:
C                JYNBH to calculate the Jn and Yn
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BJ(0:N),DJ(0:N),BY(0:N),DY(0:N)
        CALL JYNBH(N,0,X,NM,BJ,BY)
C       Compute derivatives by differentiation formulas
        IF (X.LT.1.0D-100) THEN
           DO 10 K=0,N
              DJ(K) = 0.0D0
 10           DY(K) = 1.0D+300
           DJ(1)=0.5D0
        ELSE
           DJ(0)=-BJ(1)
           DO 40 K=1,NM
 40           DJ(K)=BJ(K-1)-K/X*BJ(K)
           DY(0)=-BY(1)
           DO 50 K=1,NM
 50           DY(K)=BY(K-1)-K*BY(K)/X
        END IF
        RETURN
        END


C       **********************************

        SUBROUTINE JYNBH(N,NMIN,X,NM,BJ,BY)
C
C       =====================================================
C       Purpose: Compute Bessel functions Jn(x), Yn(x)
C       Input :  x --- Argument of Jn(x) and Yn(x) ( x ≥ 0 )
C                n --- Highest order of Jn(x) and Yn(x) computed  ( n ≥ 0 )
C                nmin -- Lowest order computed  ( nmin ≥ 0 )
C       Output:  BJ(n-NMIN) --- Jn(x)   ; if indexing starts at 0
C                BY(n-NMIN) --- Yn(x)   ; if indexing starts at 0
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 to calculate the starting
C                point for backward recurrence
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BJ(0:N-NMIN),BY(0:N-NMIN),A(4),B(4),A1(4),B1(4)
        PI=3.141592653589793D0
        R2P=.63661977236758D0
        NM=N
        IF (X.LT.1.0D-100) THEN
           DO 10 K=NMIN,N
              BJ(K-NMIN)=0.0D0
10            BY(K-NMIN)=-1.0D+300
           IF (NMIN.EQ.0) BJ(0)=1.0D0
           RETURN
        ENDIF
        IF (X.LE.300.0.OR.N.GT.INT(0.9*X)) THEN
C          Backward recurrence for Jn
           IF (N.EQ.0) NM=1
           M=MSTA1(X,200)
           IF (M.LT.NM) THEN
              NM=M
           ELSE
              M=MSTA2(X,NM,15)
           ENDIF
           BS=0.0D0
           SU=0.0D0
           SV=0.0D0
           F2=0.0D0
           F1=1.0D-100
           F=0.0D0
           DO 15 K=M,0,-1
              F=2.0D0*(K+1.0D0)/X*F1-F2
              IF (K.LE.NM .AND. K.GE.NMIN) BJ(K-NMIN)=F
              IF (K.EQ.2*INT(K/2).AND.K.NE.0) THEN
                 BS=BS+2.0D0*F
                 SU=SU+(-1)**(K/2)*F/K
              ELSE IF (K.GT.1) THEN
                 SV=SV+(-1)**(K/2)*K/(K*K-1.0D0)*F
              ENDIF
              F2=F1
15            F1=F
           S0=BS+F
           DO 20 K=NMIN,NM
20            BJ(K-NMIN)=BJ(K-NMIN)/S0
C          Estimates for Yn at start of recurrence
           BJ0 = F1 / S0
           BJ1 = F2 / S0
           EC=DLOG(X/2.0D0)+0.5772156649015329D0
           BY0=R2P*(EC*BJ0-4.0D0*SU/S0)
           BY1=R2P*((EC-1.0D0)*BJ1-BJ0/X-4.0D0*SV/S0)
           IF (0.GE.NMIN) BY(0-NMIN)=BY0
           IF (1.GE.NMIN) BY(1-NMIN)=BY1
           KY=2
        ELSE
C          Hankel expansion
           DATA A/-.7031250000000000D-01,.1121520996093750D+00,
     &            -.5725014209747314D+00,.6074042001273483D+01/
           DATA B/ .7324218750000000D-01,-.2271080017089844D+00,
     &             .1727727502584457D+01,-.2438052969955606D+02/
           DATA A1/.1171875000000000D+00,-.1441955566406250D+00,
     &             .6765925884246826D+00,-.6883914268109947D+01/
           DATA B1/-.1025390625000000D+00,.2775764465332031D+00,
     &             -.1993531733751297D+01,.2724882731126854D+02/
           T1=X-0.25D0*PI
           P0=1.0D0
           Q0=-0.125D0/X
           DO 25 K=1,4
              P0=P0+A(K)*X**(-2*K)
25            Q0=Q0+B(K)*X**(-2*K-1)
           CU=DSQRT(R2P/X)
           BJ0=CU*(P0*DCOS(T1)-Q0*DSIN(T1))
           BY0=CU*(P0*DSIN(T1)+Q0*DCOS(T1))
           IF (0.GE.NMIN) BJ(0-NMIN)=BJ0
           IF (0.GE.NMIN) BY(0-NMIN)=BY0
           T2=X-0.75D0*PI
           P1=1.0D0
           Q1=0.375D0/X
           DO 30 K=1,4
              P1=P1+A1(K)*X**(-2*K)
30            Q1=Q1+B1(K)*X**(-2*K-1)
           BJ1=CU*(P1*DCOS(T2)-Q1*DSIN(T2))
           BY1=CU*(P1*DSIN(T2)+Q1*DCOS(T2))
           IF (1.GE.NMIN) BJ(1-NMIN)=BJ1
           IF (1.GE.NMIN) BY(1-NMIN)=BY1
           DO 35 K=2,NM
              BJK=2.0D0*(K-1.0D0)/X*BJ1-BJ0
              IF (K.GE.NMIN) BJ(K-NMIN)=BJK
              BJ0=BJ1
35            BJ1=BJK
           KY=2
        ENDIF
C       Forward recurrence for Yn
        DO 45 K=KY,NM
           BYK=2.0D0*(K-1.0D0)*BY1/X-BY0
           IF (K.GE.NMIN) BY(K-NMIN)=BYK
           BY0=BY1
45         BY1=BYK
        RETURN
        END


C       **********************************

        SUBROUTINE STVH1(X,SH1)
C
C       =============================================
C       Purpose: Compute Struve function H1(x)
C       Input :  x   --- Argument of H1(x) ( x ≥ 0 )
C       Output:  SH1 --- H1(x)
C       =============================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        R=1.0D0
        IF (X.LE.20.0D0) THEN
           S=0.0D0
           A0=-2.0D0/PI
           DO 10 K=1,60
              R=-R*X*X/(4.0D0*K*K-1.0D0)
              S=S+R
              IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 15
10         CONTINUE
15         SH1=A0*S
        ELSE
           S=1.0D0
           KM=INT(.5*X)
           IF (X.GT.50.D0) KM=25
           DO 20 K=1,KM
              R=-R*(4.0D0*K*K-1.0D0)/(X*X)
              S=S+R
              IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 25
20         CONTINUE
25         T=4.0D0/X
           T2=T*T
           P1=((((.42414D-5*T2-.20092D-4)*T2+.580759D-4)*T2
     &        -.223203D-3)*T2+.29218256D-2)*T2+.3989422819D0
           Q1=T*(((((-.36594D-5*T2+.1622D-4)*T2-.398708D-4)*
     &        T2+.1064741D-3)*T2-.63904D-3)*T2+.0374008364D0)
           TA1=X-.75D0*PI
           BY1=2.0D0/DSQRT(X)*(P1*DSIN(TA1)+Q1*DCOS(TA1))
           SH1=2.0/PI*(1.0D0+S/(X*X))+BY1
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE LEGZO(N,X,W)
C
C       =========================================================
C       Purpose : Compute the zeros of Legendre polynomial Pn(x)
C                 in the interval [-1,1], and the corresponding
C                 weighting coefficients for Gauss-Legendre
C                 integration
C       Input :   n    --- Order of the Legendre polynomial
C       Output:   X(n) --- Zeros of the Legendre polynomial
C                 W(n) --- Corresponding weighting coefficients
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(N),W(N)
        N0=(N+1)/2
        PF=0.0D0
        PD=0.0D0
        DO 45 NR=1,N0
           Z=DCOS(3.1415926D0*(NR-0.25D0)/N)
10         Z0=Z
           P=1.0D0
           DO 15 I=1,NR-1
15            P=P*(Z-X(I))
           F0=1.0D0
           IF (NR.EQ.N0.AND.N.NE.2*INT(N/2)) Z=0.0D0
           F1=Z
           DO 20 K=2,N
              PF=(2.0D0-1.0D0/K)*Z*F1-(1.0D0-1.0D0/K)*F0
              PD=K*(F1-Z*PF)/(1.0D0-Z*Z)
              F0=F1
20            F1=PF
           IF (Z.EQ.0.0) GO TO 40
           FD=PF/P
           Q=0.0D0
           DO 35 I=1,NR
              WP=1.0D0
              DO 30 J=1,NR
                 IF (J.NE.I) WP=WP*(Z-X(J))
30            CONTINUE
35            Q=Q+WP
           GD=(PD-Q*FD)/P
           Z=Z-FD/GD
           IF (DABS(Z-Z0).GT.DABS(Z)*1.0D-15) GO TO 10
40         X(NR)=Z
           X(N+1-NR)=-Z
           W(NR)=2.0D0/((1.0D0-Z*Z)*PD*PD)
45         W(N+1-NR)=W(NR)
        RETURN
        END

C       **********************************

        SUBROUTINE ASWFA(M,N,C,X,KD,CV,S1F,S1D)
C
C       ===========================================================
C       Purpose: Compute the prolate and oblate spheroidal angular
C                functions of the first kind and their derivatives
C       Input :  m  --- Mode parameter,  m = 0,1,2,...
C                n  --- Mode parameter,  n = m,m+1,...
C                c  --- Spheroidal parameter
C                x  --- Argument of angular function, |x| < 1.0
C                KD --- Function code
C                       KD=1 for prolate;  KD=-1 for oblate
C                cv --- Characteristic value
C       Output:  S1F --- Angular function of the first kind
C                S1D --- Derivative of the angular function of
C                        the first kind
C       Routine called:
C                SCKB for computing expansion coefficients ck
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION CK(200),DF(200)
        EPS=1.0D-14
        X0=X
        X=DABS(X)
        IP=1
        IF (N-M.EQ.2*INT((N-M)/2)) IP=0
        NM=40+INT((N-M)/2+C)
        NM2=NM/2-2
        CALL SDMN(M,N,C,CV,KD,DF)
        CALL SCKB(M,N,C,DF,CK)
        X1=1.0D0-X*X
        IF (M.EQ.0.AND.X1.EQ.0.0D0) THEN
           A0=1.0D0
        ELSE
           A0=X1**(0.5D0*M)
        ENDIF
        SU1=CK(1)
        DO 10 K=1,NM2
           R=CK(K+1)*X1**K
           SU1=SU1+R
           IF (K.GE.10.AND.DABS(R/SU1).LT.EPS) GO TO 15
10         CONTINUE
15      S1F=A0*X**IP*SU1
        IF (X.EQ.1.0D0) THEN
           IF (M.EQ.0) S1D=IP*CK(1)-2.0D0*CK(2)
           IF (M.EQ.1) S1D=-1.0D+100
           IF (M.EQ.2) S1D=-2.0D0*CK(1)
           IF (M.GE.3) S1D=0.0D0
        ELSE
           D0=IP-M/X1*X**(IP+1.0D0)
           D1=-2.0D0*A0*X**(IP+1.0D0)
           SU2=CK(2)
           DO 20 K=2,NM2
              R=K*CK(K+1)*X1**(K-1.0D0)
              SU2=SU2+R
              IF (K.GE.10.AND.DABS(R/SU2).LT.EPS) GO TO 25
20            CONTINUE
25         S1D=D0*A0*SU1+D1*SU2
        ENDIF
        IF (X0.LT.0.0D0.AND.IP.EQ.0) S1D=-S1D
        IF (X0.LT.0.0D0.AND.IP.EQ.1) S1F=-S1F
        X=X0
        RETURN
        END



C       **********************************

        SUBROUTINE JYNA(N,X,NM,BJ,DJ,BY,DY)
C
C       ==========================================================
C       Purpose: Compute Bessel functions Jn(x) & Yn(x) and
C                their derivatives
C       Input :  x --- Argument of Jn(x) & Yn(x)  ( x ≥ 0 )
C                n --- Order of Jn(x) & Yn(x)
C       Output:  BJ(n) --- Jn(x)
C                DJ(n) --- Jn'(x)
C                BY(n) --- Yn(x)
C                DY(n) --- Yn'(x)
C                NM --- Highest order computed
C       Routines called:
C            (1) JY01B to calculate J0(x), J1(x), Y0(x) & Y1(x)
C            (2) MSTA1 and MSTA2 to calculate the starting
C                point for backward recurrence
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BJ(0:N),BY(0:N),DJ(0:N),DY(0:N)
        NM=N
        IF (X.LT.1.0D-100) THEN
           DO 10 K=0,N
              BJ(K)=0.0D0
              DJ(K)=0.0D0
              BY(K)=-1.0D+300
10            DY(K)=1.0D+300
           BJ(0)=1.0D0
           DJ(1)=0.5D0
           RETURN
        ENDIF
        CALL JY01B(X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
        BJ(0)=BJ0
        BJ(1)=BJ1
        BY(0)=BY0
        BY(1)=BY1
        DJ(0)=DJ0
        DJ(1)=DJ1
        DY(0)=DY0
        DY(1)=DY1
        IF (N.LE.1) RETURN
        IF (N.LT.INT(0.9*X)) THEN
           DO 20 K=2,N
              BJK=2.0D0*(K-1.0D0)/X*BJ1-BJ0
              BJ(K)=BJK
              BJ0=BJ1
20            BJ1=BJK
        ELSE
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F2=0.0D0
           F1=1.0D-100
           F=0.0D0
           DO 30 K=M,0,-1
              F=2.0D0*(K+1.0D0)/X*F1-F2
              IF (K.LE.NM) BJ(K)=F
              F2=F1
30            F1=F
           IF (DABS(BJ0).GT.DABS(BJ1)) THEN
              CS=BJ0/F
           ELSE
              CS=BJ1/F2
           ENDIF
           DO 40 K=0,NM
40            BJ(K)=CS*BJ(K)
        ENDIF
        DO 50 K=2,NM
50         DJ(K)=BJ(K-1)-K/X*BJ(K)
        F0=BY(0)
        F1=BY(1)
        DO 60 K=2,NM
           F=2.0D0*(K-1.0D0)/X*F1-F0
           BY(K)=F
           F0=F1
60         F1=F
        DO 70 K=2,NM
70         DY(K)=BY(K-1)-K*BY(K)/X
        RETURN
        END



C       **********************************

        SUBROUTINE PBDV(V,X,DV,DP,PDF,PDD)
C
C       ====================================================
C       Purpose: Compute parabolic cylinder functions Dv(x)
C                and their derivatives
C       Input:   x --- Argument of Dv(x)
C                v --- Order of Dv(x)
C       Output:  DV(na) --- Dn+v0(x)
C                DP(na) --- Dn+v0'(x)
C                ( na = |n|, v0 = v-n, |v0| < 1,
C                  n = 0,±1,±2,… )
C                PDF --- Dv(x)
C                PDD --- Dv'(x)
C       Routines called:
C             (1) DVSA for computing Dv(x) for small |x|
C             (2) DVLA for computing Dv(x) for large |x|
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION DV(0:*),DP(0:*)
        XA=DABS(X)
        VH=V
        V=V+DSIGN(1.0D0,V)
        NV=INT(V)
        V0=V-NV
        NA=ABS(NV)
        EP=DEXP(-.25D0*X*X)
        JA=0
        IF (NA.GE.1) JA=1
        IF (V.GE.0.0) THEN
           IF (V0.EQ.0.0) THEN
              PD0=EP
              PD1=X*EP
           ELSE
              DO 10 L=0,JA
                 V1=V0+L
                 IF (XA.LE.5.8) CALL DVSA(V1,X,PD1)
                 IF (XA.GT.5.8) CALL DVLA(V1,X,PD1)
                 IF (L.EQ.0) PD0=PD1
10            CONTINUE
           ENDIF
           DV(0)=PD0
           DV(1)=PD1
           DO 15 K=2,NA
              PDF=X*PD1-(K+V0-1.0D0)*PD0
              DV(K)=PDF
              PD0=PD1
15            PD1=PDF
        ELSE
           IF (X.LE.0.0) THEN
              IF (XA.LE.5.8D0)  THEN
                 CALL DVSA(V0,X,PD0)
                 V1=V0-1.0D0
                 CALL DVSA(V1,X,PD1)
              ELSE
                 CALL DVLA(V0,X,PD0)
                 V1=V0-1.0D0
                 CALL DVLA(V1,X,PD1)
              ENDIF
              DV(0)=PD0
              DV(1)=PD1
              DO 20 K=2,NA
                 PD=(-X*PD1+PD0)/(K-1.0D0-V0)
                 DV(K)=PD
                 PD0=PD1
20               PD1=PD
           ELSE IF (X.LE.2.0) THEN
              V2=NV+V0
              IF (NV.EQ.0) V2=V2-1.0D0
              NK=INT(-V2)
              CALL DVSA(V2,X,F1)
              V1=V2+1.0D0
              CALL DVSA(V1,X,F0)
              DV(NK)=F1
              DV(NK-1)=F0
              DO 25 K=NK-2,0,-1
                 F=X*F0+(K-V0+1.0D0)*F1
                 DV(K)=F
                 F1=F0
25               F0=F
           ELSE
              IF (XA.LE.5.8) CALL DVSA(V0,X,PD0)
              IF (XA.GT.5.8) CALL DVLA(V0,X,PD0)
              DV(0)=PD0
              M=100+NA
              F1=0.0D0
              F0=1.0D-30
              F=0.0D0
              DO 30 K=M,0,-1
                 F=X*F0+(K-V0+1.0D0)*F1
                 IF (K.LE.NA) DV(K)=F
                 F1=F0
30               F0=F
              S0=PD0/F
              DO 35 K=0,NA
35               DV(K)=S0*DV(K)
           ENDIF
        ENDIF
        DO 40 K=0,NA-1
           V1=ABS(V0)+K
           IF (V.GE.0.0D0) THEN
              DP(K)=0.5D0*X*DV(K)-DV(K+1)
           ELSE
              DP(K)=-0.5D0*X*DV(K)-V1*DV(K+1)
           ENDIF
40      CONTINUE
        PDF=DV(NA-1)
        PDD=DP(NA-1)
        V=VH
        RETURN
        END



C       **********************************

        SUBROUTINE ITSH0(X,TH0)
C
C       ===================================================
C       Purpose: Evaluate the integral of Struve function
C                H0(t) with respect to t from 0 and x
C       Input :  x   --- Upper limit  ( x ≥ 0 )
C       Output:  TH0 --- Integration of H0(t) from 0 and x
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(25)
        PI=3.141592653589793D0
        R=1.0D0
        IF (X.LE.30.0) THEN
           S=0.5D0
           DO 10 K=1,100
              RD=1.0D0
              IF (K.EQ.1) RD=0.5D0
              R=-R*RD*K/(K+1.0D0)*(X/(2.0D0*K+1.0D0))**2
              S=S+R
              IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 15
10         CONTINUE
15         TH0=2.0D0/PI*X*X*S
        ELSE
           S=1.0D0
           DO 20 K=1,12
              R=-R*K/(K+1.0D0)*((2.0D0*K+1.0D0)/X)**2
              S=S+R
              IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 25
20         CONTINUE
25         EL=.57721566490153D0
           S0=S/(PI*X*X)+2.0D0/PI*(DLOG(2.0D0*X)+EL)
           A0=1.0D0
           A1=5.0D0/8.0D0
           A(1)=A1
           DO 30 K=1,20
              AF=((1.5D0*(K+.5D0)*(K+5.0D0/6.0D0)*A1-.5D0
     &           *(K+.5D0)*(K+.5D0)*(K-.5D0)*A0))/(K+1.0D0)
              A(K+1)=AF
              A0=A1
30            A1=AF
           BF=1.0D0
           R=1.0D0
           DO 35 K=1,10
              R=-R/(X*X)
35            BF=BF+A(2*K)*R
           BG=A(1)/X
           R=1.0D0/X
           DO 40 K=1,10
              R=-R/(X*X)
40            BG=BG+A(2*K+1)*R
           XP=X+.25D0*PI
           TY=DSQRT(2.0D0/(PI*X))*(BG*DCOS(XP)-BF*DSIN(XP))
           TH0=TY+S0
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE CERZO(NT,ZO)
C
C       ===============================================================
C       Purpose : Evaluate the complex zeros of error function erf(z)
C                 using the modified Newton's iteration method
C       Input :   NT --- Total number of zeros
C       Output:   ZO(L) --- L-th zero of erf(z), L=1,2,...,NT
C       Routine called: CERF for computing erf(z) and erf'(z)
C       ===============================================================
C
        IMPLICIT DOUBLE PRECISION (E,P,W)
        IMPLICIT COMPLEX *16 (C,Z)
        DIMENSION ZO(NT)
        PI=3.141592653589793D0
        W=0.0D0
        DO 35 NR=1,NT
           PU=DSQRT(PI*(4.0D0*NR-0.5D0))
           PV=PI*DSQRT(2.0D0*NR-0.25D0)
           PX=0.5*PU-0.5*DLOG(PV)/PU
           PY=0.5*PU+0.5*DLOG(PV)/PU
           Z = DCMPLX(PX, PY)
           IT=0
15         IT=IT+1
           CALL CERF(Z,ZF,ZD)
           ZP=(1.0D0,0.0D0)
           DO 20 I=1,NR-1
20            ZP=ZP*(Z-ZO(I))
           ZFD=ZF/ZP
           ZQ=(0.0D0,0.0D0)
           DO 30 I=1,NR-1
              ZW=(1.0D0,0.0D0)
              DO 25 J=1,NR-1
                 IF (J.EQ.I) GO TO 25
                 ZW=ZW*(Z-ZO(J))
25            CONTINUE
30            ZQ=ZQ+ZW
           ZGD=(ZD-ZQ*ZFD)/ZP
           Z=Z-ZFD/ZGD
           W0=W
           W=CDABS(Z)
           IF (IT.LE.50.AND.DABS((W-W0)/W).GT.1.0D-11) GO TO 15
35         ZO(NR)=Z
        RETURN
        END



C       **********************************

        SUBROUTINE GAMMA2(X,GA)
C
C       ==================================================
C       Purpose: Compute gamma function Г(x)
C       Input :  x  --- Argument of Г(x)
C                       ( x is not equal to 0,-1,-2,…)
C       Output:  GA --- Г(x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           R=1.0D0
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE CHGU(A,B,X,HU,MD)
C
C       =======================================================
C       Purpose: Compute the confluent hypergeometric function
C                U(a,b,x)
C       Input  : a  --- Parameter
C                b  --- Parameter
C                x  --- Argument  ( x > 0 )
C       Output:  HU --- U(a,b,x)
C                MD --- Method code
C       Routines called:
C            (1) CHGUS for small x ( MD=1 )
C            (2) CHGUL for large x ( MD=2 )
C            (3) CHGUBI for integer b ( MD=3 )
C            (4) CHGUIT for numerical integration ( MD=4 )
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL IL1,IL2,IL3,BL1,BL2,BL3,BN
        AA=A-B+1.0D0
        IL1=A.EQ.INT(A).AND.A.LE.0.0
        IL2=AA.EQ.INT(AA).AND.AA.LE.0.0
        IL3=ABS(A*(A-B+1.0))/X.LE.2.0
        BL1=X.LE.5.0.OR.(X.LE.10.0.AND.A.LE.2.0)
        BL2=(X.GT.5.0.AND.X.LE.12.5).AND.(A.GE.1.0.AND.B.GE.A+4.0)
        BL3=X.GT.12.5.AND.A.GE.5.0.AND.B.GE.A+5.0
        BN=B.EQ.INT(B).AND.B.NE.0.0
        ID1=-100
        HU1=0.0D0
        IF (B.NE.INT(B)) THEN
           CALL CHGUS(A,B,X,HU,ID1)
           MD=1
           IF (ID1.GE.9) RETURN
           HU1=HU
        ENDIF
        IF (IL1.OR.IL2.OR.IL3) THEN
           CALL CHGUL(A,B,X,HU,ID)
           MD=2
           IF (ID.GE.9) RETURN
           IF (ID1.GT.ID) THEN
              MD=1
              ID=ID1
              HU=HU1
           ENDIF
        ENDIF
        IF (A.GE.1.0) THEN
           IF (BN.AND.(BL1.OR.BL2.OR.BL3)) THEN
              CALL CHGUBI(A,B,X,HU,ID)
              MD=3
           ELSE
              CALL CHGUIT(A,B,X,HU,ID)
              MD=4
           ENDIF
        ELSE
           IF (B.LE.A) THEN
              A00=A
              B00=B
              A=A-B+1.0D0
              B=2.0D0-B
              CALL CHGUIT(A,B,X,HU,ID)
              HU=X**(1.0D0-B00)*HU
              A=A00
              B=B00
              MD=4
           ELSE IF (BN.AND.(.NOT.IL1)) THEN
              CALL CHGUBI(A,B,X,HU,ID)
              MD=3
           ENDIF
        ENDIF
        IF (ID.LT.6) WRITE(*,*)'No accurate result obtained'
        RETURN
        END



C       **********************************

        SUBROUTINE LAMN(N,X,NM,BL,DL)
C
C       =========================================================
C       Purpose: Compute lambda functions and their derivatives
C       Input:   x --- Argument of lambda function
C                n --- Order of lambda function
C       Output:  BL(n) --- Lambda function of order n
C                DL(n) --- Derivative of lambda function
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the start
C                point for backward recurrence
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BL(0:N),DL(0:N)
        NM=N
        IF (DABS(X).LT.1.0D-100) THEN
           DO 10 K=0,N
              BL(K)=0.0D0
10            DL(K)=0.0D0
           BL(0)=1.0D0
           DL(1)=0.5D0
           RETURN
        ENDIF
        IF (X.LE.12.0D0) THEN
           X2=X*X
           DO 25 K=0,N
              BK=1.0D0
              R=1.0D0
              DO 15 I=1,50
                 R=-0.25D0*R*X2/(I*(I+K))
                 BK=BK+R
                 IF (DABS(R).LT.DABS(BK)*1.0D-15) GO TO 20
15            CONTINUE
20            BL(K)=BK
25            IF (K.GE.1) DL(K-1)=-0.5D0*X/K*BK
           UK=1.0D0
           R=1.0D0
           DO 30 I=1,50
              R=-0.25D0*R*X2/(I*(I+N+1.0D0))
              UK=UK+R
              IF (DABS(R).LT.DABS(UK)*1.0D-15) GO TO 35
30            CONTINUE
35         DL(N)=-0.5D0*X/(N+1.0D0)*UK
           RETURN
        ENDIF
        IF (N.EQ.0) NM=1
        M=MSTA1(X,200)
        IF (M.LT.NM) THEN
           NM=M
        ELSE
           M=MSTA2(X,NM,15)
        ENDIF
        BS=0.0D0
        F=0.0D0
        F0=0.0D0
        F1=1.0D-100
        DO 40 K=M,0,-1
           F=2.0D0*(K+1.0D0)*F1/X-F0
           IF (K.LE.NM) BL(K)=F
           IF (K.EQ.2*INT(K/2)) BS=BS+2.0D0*F
           F0=F1
40         F1=F
        BG=BS-F
        DO 45 K=0,NM
45         BL(K)=BL(K)/BG
        R0=1.0D0
        DO 50 K=1,NM
           R0=2.0D0*R0*K/X
50         BL(K)=R0*BL(K)
        DL(0)=-0.5D0*X*BL(1)
        DO 55 K=1,NM
55         DL(K)=2.0D0*K/X*(BL(K-1)-BL(K))
        RETURN
        END



C       **********************************

        SUBROUTINE COMELP(HK,CK,CE)
C
C       ==================================================
C       Purpose: Compute complete elliptic integrals K(k)
C                and E(k)
C       Input  : K  --- Modulus k ( 0 ≤ k ≤ 1 )
C       Output : CK --- K(k)
C                CE --- E(k)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PK=1.0D0-HK*HK
        IF (HK.EQ.1.0) THEN
           CK=1.0D+300
           CE=1.0D0
        ELSE
           AK=(((.01451196212D0*PK+.03742563713D0)*PK
     &        +.03590092383D0)*PK+.09666344259D0)*PK+
     &        1.38629436112D0
           BK=(((.00441787012D0*PK+.03328355346D0)*PK+
     &        .06880248576D0)*PK+.12498593597D0)*PK+.5D0
           CK=AK-BK*DLOG(PK)
           AE=(((.01736506451D0*PK+.04757383546D0)*PK+
     &        .0626060122D0)*PK+.44325141463D0)*PK+1.0D0
           BE=(((.00526449639D0*PK+.04069697526D0)*PK+
     &        .09200180037D0)*PK+.2499836831D0)*PK
           CE=AE-BE*DLOG(PK)
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE INCOB(A,B,X,BIX)
C
C       ========================================================
C       Purpose: Compute the incomplete beta function Ix(a,b)
C       Input :  a --- Parameter
C                b --- Parameter
C                x --- Argument ( 0 ≤ x ≤ 1 )
C       Output:  BIX --- Ix(a,b)
C       Routine called: BETA for computing beta function B(p,q)
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION DK(51),FK(51)
        S0=(A+1.0D0)/(A+B+2.0D0)
        CALL BETA(A,B,BT)
        IF (X.LE.S0) THEN
           DO 10 K=1,20
10            DK(2*K)=K*(B-K)*X/(A+2.0D0*K-1.0D0)/(A+2.0D0*K)
           DO 15 K=0,20
15            DK(2*K+1)=-(A+K)*(A+B+K)*X/(A+2.D0*K)/(A+2.0*K+1.0)
           T1=0.0D0
           DO 20 K=20,1,-1
20            T1=DK(K)/(1.0D0+T1)
           TA=1.0D0/(1.0D0+T1)
           BIX=X**A*(1.0D0-X)**B/(A*BT)*TA
        ELSE
           DO 25 K=1,20
25            FK(2*K)=K*(A-K)*(1.0D0-X)/(B+2.*K-1.0)/(B+2.0*K)
           DO 30 K=0,20
30            FK(2*K+1)=-(B+K)*(A+B+K)*(1.D0-X)/
     &                   (B+2.D0*K)/(B+2.D0*K+1.D0)
           T2=0.0D0
           DO 35 K=20,1,-1
35            T2=FK(K)/(1.0D0+T2)
           TB=1.0D0/(1.0D0+T2)
           BIX=1.0D0-X**A*(1.0D0-X)**B/(B*BT)*TB
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE CVF(KD,M,Q,A,MJ,F)
C
C       ======================================================
C       Purpose: Compute the value of F for characteristic
C                equation of Mathieu functions
C       Input :  m --- Order of Mathieu functions
C                q --- Parameter of Mathieu functions
C                A --- Characteristic value
C       Output:  F --- Value of F for characteristic equation
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        B=A
        IC=INT(M/2)
        L=0
        L0=0
        J0=2
        JF=IC
        IF (KD.EQ.1) L0=2
        IF (KD.EQ.1) J0=3
        IF (KD.EQ.2.OR.KD.EQ.3) L=1
        IF (KD.EQ.4) JF=IC-1
        T1=0.0D0
        DO 10 J=MJ,IC+1,-1
10         T1=-Q*Q/((2.0D0*J+L)**2-B+T1)
        IF (M.LE.2) THEN
           T2=0.0D0
           IF (KD.EQ.1.AND.M.EQ.0) T1=T1+T1
           IF (KD.EQ.1.AND.M.EQ.2) T1=-2.0D0*Q*Q/(4.0D0-B+T1)-4.0D0
           IF (KD.EQ.2.AND.M.EQ.1) T1=T1+Q
           IF (KD.EQ.3.AND.M.EQ.1) T1=T1-Q
        ELSE
           T0=0.0D0
           IF (KD.EQ.1) T0=4.0D0-B+2.0D0*Q*Q/B
           IF (KD.EQ.2) T0=1.0D0-B+Q
           IF (KD.EQ.3) T0=1.0D0-B-Q
           IF (KD.EQ.4) T0=4.0D0-B
           T2=-Q*Q/T0
           DO 15 J=J0,JF
15            T2=-Q*Q/((2.0D0*J-L-L0)**2-B+T2)
        ENDIF
        F=(2.0D0*IC+L)**2+T1+T2-B
        RETURN
        END



C       **********************************

        SUBROUTINE CLPN(N,X,Y,CPN,CPD)
C
C       ==================================================
C       Purpose: Compute Legendre polynomials Pn(z) and
C                their derivatives Pn'(z) for a complex
C                argument
C       Input :  x --- Real part of z
C                y --- Imaginary part of z
C                n --- Degree of Pn(z), n = 0,1,2,...
C       Output:  CPN(n) --- Pn(z)
C                CPD(n) --- Pn'(z)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX *16 (C,Z)
        DIMENSION CPN(0:N),CPD(0:N)
        Z = DCMPLX(X, Y)
        CPN(0)=(1.0D0,0.0D0)
        CPN(1)=Z
        CPD(0)=(0.0D0,0.0D0)
        CPD(1)=(1.0D0,0.0D0)
        CP0=(1.0D0,0.0D0)
        CP1=Z
        DO 10 K=2,N
           CPF=(2.0D0*K-1.0D0)/K*Z*CP1-(K-1.0D0)/K*CP0
           CPN(K)=CPF
           IF (DABS(X).EQ.1.0D0.AND.Y.EQ.0.0D0) THEN
              CPD(K)=0.5D0*X**(K+1)*K*(K+1.0D0)
           ELSE
              CPD(K)=K*(CP1-Z*CPF)/(1.0D0-Z*Z)
           ENDIF
           CP0=CP1
10         CP1=CPF
        RETURN
        END

C       **********************************

        SUBROUTINE LQMNS(M,N,X,QM,QD)
C
C       ========================================================
C       Purpose: Compute associated Legendre functions Qmn(x)
C                and Qmn'(x) for a given order
C       Input :  x --- Argument of Qmn(x)
C                m --- Order of Qmn(x),  m = 0,1,2,...
C                n --- Degree of Qmn(x), n = 0,1,2,...
C       Output:  QM(n) --- Qmn(x)
C                QD(n) --- Qmn'(x)
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION QM(0:N),QD(0:N)
        DO 10 K=0,N
           QM(K)=0.0D0
10         QD(K)=0.0D0
        IF (DABS(X).EQ.1.0D0) THEN
           DO 15 K=0,N
              QM(K)=1.0D+300
15            QD(K)=1.0D+300
           RETURN
        ENDIF
        LS=1
        IF (DABS(X).GT.1.0D0) LS=-1
        XQ=DSQRT(LS*(1.0D0-X*X))
        Q0=0.5D0*DLOG(DABS((X+1.0)/(X-1.0)))
        Q00=Q0
        Q10=-1.0D0/XQ
        Q01=X*Q0-1.0D0
        Q11=-LS*XQ*(Q0+X/(1.0D0-X*X))
        QF0=Q00
        QF1=Q10
        QM0=0.0D0
        QM1=0.0D0
        DO 20 K=2,M
           QM0=-2.0D0*(K-1.0)/XQ*X*QF1-LS*(K-1.0)*(2.0-K)*QF0
           QF0=QF1
20         QF1=QM0
        IF (M.EQ.0) QM0=Q00
        IF (M.EQ.1) QM0=Q10
        QM(0)=QM0
        IF (DABS(X).LT.1.0001D0) THEN
           IF (M.EQ.0.AND.N.GT.0) THEN
              QF0=Q00
              QF1=Q01
              DO 25 K=2,N
                 QF2=((2.0*K-1.0D0)*X*QF1-(K-1.0)*QF0)/K
                 QM(K)=QF2
                 QF0=QF1
25               QF1=QF2
           ENDIF
           QG0=Q01
           QG1=Q11
           DO 30 K=2,M
              QM1=-2.0D0*(K-1.0)/XQ*X*QG1-LS*K*(3.0-K)*QG0
              QG0=QG1
30            QG1=QM1
           IF (M.EQ.0) QM1=Q01
           IF (M.EQ.1) QM1=Q11
           QM(1)=QM1
           IF (M.EQ.1.AND.N.GT.1) THEN
              QH0=Q10
              QH1=Q11
              DO 35 K=2,N
                 QH2=((2.0*K-1.0D0)*X*QH1-K*QH0)/(K-1.0)
                 QM(K)=QH2
                 QH0=QH1
35               QH1=QH2
           ELSE IF (M.GE.2) THEN
              QG0=Q00
              QG1=Q01
              QH0=Q10
              QH1=Q11
              QMK=0.0D0
              DO 45 L=2,N
                 Q0L=((2.0D0*L-1.0D0)*X*QG1-(L-1.0D0)*QG0)/L
                 Q1L=((2.0*L-1.0D0)*X*QH1-L*QH0)/(L-1.0D0)
                 QF0=Q0L
                 QF1=Q1L
                 DO 40 K=2,M
                    QMK=-2.0D0*(K-1.0)/XQ*X*QF1-LS*(K+L-1.0)*
     &                  (L+2.0-K)*QF0
                    QF0=QF1
40                  QF1=QMK
                 QM(L)=QMK
                 QG0=QG1
                 QG1=Q0L
                 QH0=QH1
45               QH1=Q1L
           ENDIF
        ELSE
           IF (DABS(X).GT.1.1) THEN
              KM=40+M+N
           ELSE
              KM=(40+M+N)*INT(-1.0-1.8*LOG(X-1.0))
           ENDIF
           QF2=0.0D0
           QF1=1.0D0
           DO 50 K=KM,0,-1
              QF0=((2.0*K+3.0D0)*X*QF1-(K+2.0-M)*QF2)/(K+M+1.0)
              IF (K.LE.N) QM(K)=QF0
              QF2=QF1
50            QF1=QF0
           DO 55 K=0,N
55            QM(K)=QM(K)*QM0/QF0
        ENDIF
        IF (DABS(X).LT.1.0D0) THEN
           DO 60 K=0,N
60            QM(K)=(-1)**M*QM(K)
        ENDIF
        QD(0)=((1.0D0-M)*QM(1)-X*QM(0))/(X*X-1.0)
        DO 65 K=1,N
65         QD(K)=(K*X*QM(K)-(K+M)*QM(K-1))/(X*X-1.0)
        RETURN
        END

C       **********************************

        SUBROUTINE CIKLV(V,Z,CBIV,CDIV,CBKV,CDKV)
C
C       =====================================================
C       Purpose: Compute modified Bessel functions Iv(z) and
C                Kv(z) and their derivatives with a complex
C                argument and a large order
C       Input:   v --- Order of Iv(z) and Kv(z)
C                z --- Complex argument
C       Output:  CBIV --- Iv(z)
C                CDIV --- Iv'(z)
C                CBKV --- Kv(z)
C                CDKV --- Kv'(z)
C       Routine called:
C                CJK to compute the expansion coefficients
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CF(12),A(91)
        PI=3.141592653589793D0
        KM=12
        CALL CJK(KM,A)
        DO 30 L=1,0,-1
           V0=V-L
           CWS=CDSQRT(1.0D0+(Z/V0)*(Z/V0))
           CETA=CWS+CDLOG(Z/V0/(1.0D0+CWS))
           CT=1.0D0/CWS
           CT2=CT*CT
           DO 15 K=1,KM
              L0=K*(K+1)/2+1
              LF=L0+K
              CF(K)=A(LF)
              DO 10 I=LF-1,L0,-1
10               CF(K)=CF(K)*CT2+A(I)
15            CF(K)=CF(K)*CT**K
           VR=1.0D0/V0
           CSI=(1.0D0,0.0D0)
           DO 20 K=1,KM
20            CSI=CSI+CF(K)*VR**K
           CBIV=CDSQRT(CT/(2.0D0*PI*V0))*CDEXP(V0*CETA)*CSI
           IF (L.EQ.1) CFI=CBIV
           CSK=(1.0D0,0.0D0)
           DO 25 K=1,KM
25            CSK=CSK+(-1)**K*CF(K)*VR**K
           CBKV=CDSQRT(PI*CT/(2.0D0*V0))*CDEXP(-V0*CETA)*CSK
           IF (L.EQ.1) CFK=CBKV
30      CONTINUE
        CDIV=CFI-V/Z*CBIV
        CDKV=-CFK-V/Z*CBKV
        RETURN
        END



C       **********************************

        SUBROUTINE ELIT(HK,PHI,FE,EE)
C
C       ==================================================
C       Purpose: Compute complete and incomplete elliptic
C                integrals F(k,phi) and E(k,phi)
C       Input  : HK  --- Modulus k ( 0 ≤ k ≤ 1 )
C                Phi --- Argument ( in degrees )
C       Output : FE  --- F(k,phi)
C                EE  --- E(k,phi)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        G=0.0D0
        PI=3.14159265358979D0
        A0=1.0D0
        B0=DSQRT(1.0D0-HK*HK)
        D0=(PI/180.0D0)*PHI
        R=HK*HK
        IF (HK.EQ.1.0D0.AND.PHI.EQ.90.0D0) THEN
           FE=1.0D+300
           EE=1.0D0
        ELSE IF (HK.EQ.1.0D0) THEN
           FE=DLOG((1.0D0+DSIN(D0))/DCOS(D0))
           EE=DSIN(D0)
        ELSE
           FAC=1.0D0
           D=0.0D0
           DO 10 N=1,40
              A=(A0+B0)/2.0D0
              B=DSQRT(A0*B0)
              C=(A0-B0)/2.0D0
              FAC=2.0D0*FAC
              R=R+FAC*C*C
              IF (PHI.NE.90.0D0) THEN
                 D=D0+DATAN((B0/A0)*DTAN(D0))
                 G=G+C*DSIN(D)
                 D0=D+PI*INT(D/PI+.5D0)
              ENDIF
              A0=A
              B0=B
              IF (C.LT.1.0D-7) GO TO 15
10         CONTINUE
15         CK=PI/(2.0D0*A)
           CE=PI*(2.0D0-R)/(4.0D0*A)
           IF (PHI.EQ.90.0D0) THEN
              FE=CK
              EE=CE
           ELSE
              FE=D/(FAC*A)
              EE=FE*CE/CK+G
           ENDIF
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE ELIT3(PHI,HK,C,EL3)
C
C       =========================================================
C       Purpose: Compute the elliptic integral of the third kind
C                using Gauss-Legendre quadrature
C       Input :  Phi --- Argument ( in degrees )
C                 k  --- Modulus   ( 0 ≤ k ≤ 1.0 )
C                 c  --- Parameter ( 0 ≤ c ≤ 1.0 )
C       Output:  EL3 --- Value of the elliptic integral of the
C                        third kind
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION T(10),W(10)
        LOGICAL LB1,LB2
        DATA T/.9931285991850949D0,.9639719272779138D0,
     &         .9122344282513259D0,.8391169718222188D0,
     &         .7463319064601508D0,.6360536807265150D0,
     &         .5108670019508271D0,.3737060887154195D0,
     &         .2277858511416451D0,.7652652113349734D-1/
        DATA W/.1761400713915212D-1,.4060142980038694D-1,
     &         .6267204833410907D-1,.8327674157670475D-1,
     &         .1019301198172404D0,.1181945319615184D0,
     &         .1316886384491766D0,.1420961093183820D0,
     &         .1491729864726037D0,.1527533871307258D0/
        LB1=HK.EQ.1.0D0.AND.DABS(PHI-90.0).LE.1.0D-8
        LB2=C.EQ.1.0D0.AND.DABS(PHI-90.0).LE.1.0D-8
        IF (LB1.OR.LB2) THEN
            EL3=1.0D+300
            RETURN
        ENDIF
        C1=0.87266462599716D-2*PHI
        C2=C1
        EL3=0.0D0
        DO 10 I=1,10
           C0=C2*T(I)
           T1=C1+C0
           T2=C1-C0
           F1=1.0D0/((1.0D0-C*DSIN(T1)*DSIN(T1))*
     &              DSQRT(1.0D0-HK*HK*DSIN(T1)*DSIN(T1)))
           F2=1.0D0/((1.0D0-C*DSIN(T2)*DSIN(T2))*
     &              DSQRT(1.0D0-HK*HK*DSIN(T2)*DSIN(T2)))
10         EL3=EL3+W(I)*(F1+F2)
        EL3=C1*EL3
        RETURN
        END

C       **********************************

        SUBROUTINE EIX(X,EI)
C
C       ============================================
C       Purpose: Compute exponential integral Ei(x)
C       Input :  x  --- Argument of Ei(x)
C       Output:  EI --- Ei(x)
C       ============================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (X.EQ.0.0) THEN
           EI=-1.0D+300
        ELSE IF (X .LT. 0) THEN
           CALL E1XB(-X, EI)
           EI = -EI
        ELSE IF (DABS(X).LE.40.0) THEN
C          Power series around x=0
           EI=1.0D0
           R=1.0D0
           DO 15 K=1,100
              R=R*K*X/(K+1.0D0)**2
              EI=EI+R
              IF (DABS(R/EI).LE.1.0D-15) GO TO 20
15         CONTINUE
20         GA=0.5772156649015328D0
           EI=GA+DLOG(X)+X*EI
        ELSE
C          Asymptotic expansion (the series is not convergent)
           EI=1.0D0
           R=1.0D0
           DO 25 K=1,20
              R=R*K/X
25            EI=EI+R
           EI=DEXP(X)/X*EI
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE EIXZ(Z,CEI)
C
C       ============================================
C       Purpose: Compute exponential integral Ei(x)
C       Input :  x  --- Complex argument of Ei(x)
C       Output:  EI --- Ei(x)
C       ============================================
C
        IMPLICIT NONE
        DOUBLE COMPLEX Z, CEI, IMF
        DOUBLE PRECISION PI
        PI=3.141592653589793D0
        CALL E1Z(-Z, CEI)
        CEI = -CEI
        IF (DIMAG(Z).GT.0) THEN
           CEI = CEI + (0d0,1d0)*PI
        ELSE IF (DIMAG(Z).LT.0) THEN
           CEI = CEI - (0d0,1d0)*PI
        ELSE IF (DIMAG(Z).EQ.0) THEN
           IF (DBLE(Z).GT.0) THEN
              CEI = CEI - (0d0,1d0)*PI
           ENDIF
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE E1XB(X,E1)
C
C       ============================================
C       Purpose: Compute exponential integral E1(x)
C       Input :  x  --- Argument of E1(x)
C       Output:  E1 --- E1(x)  ( x > 0 )
C       ============================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (X.EQ.0.0) THEN
           E1=1.0D+300
        ELSE IF (X.LE.1.0) THEN
           E1=1.0D0
           R=1.0D0
           DO 10 K=1,25
              R=-R*K*X/(K+1.0D0)**2
              E1=E1+R
              IF (DABS(R).LE.DABS(E1)*1.0D-15) GO TO 15
10         CONTINUE
15         GA=0.5772156649015328D0
           E1=-GA-DLOG(X)+X*E1
        ELSE
           M=20+INT(80.0/X)
           T0=0.0D0
           DO 20 K=M,1,-1
              T0=K/(1.0D0+K/(X+T0))
20         CONTINUE
           T=1.0D0/(X+T0)
           E1=DEXP(-X)*T
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE CHGM(A,B,X,HG)
C
C       ===================================================
C       Purpose: Compute confluent hypergeometric function
C                M(a,b,x)
C       Input  : a  --- Parameter
C                b  --- Parameter ( b <> 0,-1,-2,... )
C                x  --- Argument
C       Output:  HG --- M(a,b,x)
C       Routine called: CGAMA for computing complex ln[Г(x)]
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z)
        IMPLICIT COMPLEX*16 (C)
        PI=3.141592653589793D0
        A0=A
        A1=A
        X0=X
        HG=0.0D0
        IF (B.EQ.0.0D0.OR.B.EQ.-ABS(INT(B))) THEN
           HG=1.0D+300
        ELSE IF (A.EQ.0.0D0.OR.X.EQ.0.0D0) THEN
           HG=1.0D0
        ELSE IF (A.EQ.-1.0D0) THEN
           HG=1.0D0-X/B
        ELSE IF (A.EQ.B) THEN
           HG=DEXP(X)
        ELSE IF (A-B.EQ.1.0D0) THEN
           HG=(1.0D0+X/B)*DEXP(X)
        ELSE IF (A.EQ.1.0D0.AND.B.EQ.2.0D0) THEN
           HG=(DEXP(X)-1.0D0)/X
        ELSE IF (A.EQ.INT(A).AND.A.LT.0.0D0) THEN
           M=INT(-A)
           R=1.0D0
           HG=1.0D0
           DO 10 K=1,M
              R=R*(A+K-1.0D0)/K/(B+K-1.0D0)*X
10            HG=HG+R
        ENDIF
        IF (HG.NE.0.0D0) RETURN
        IF (X.LT.0.0D0) THEN
           A=B-A
           A0=A
           X=DABS(X)
        ENDIF
        NL=0
        LA=0
        IF (A.GE.2.0D0) THEN
           NL=1
           LA=INT(A)
           A=A-LA-1.0D0
        ENDIF
        Y0=0.0D0
        Y1=0.0D0
        DO 30 N=0,NL
           IF (A0.GE.2.0D0) A=A+1.0D0
           IF (X.LE.30.0D0+DABS(B).OR.A.LT.0.0D0) THEN
              HG=1.0D0
              RG=1.0D0
              DO 15 J=1,500
                 RG=RG*(A+J-1.0D0)/(J*(B+J-1.0D0))*X
                 HG=HG+RG
                 IF (HG.NE.0D0.AND.DABS(RG/HG).LT.1.0D-15) GO TO 25
15            CONTINUE
           ELSE
              Y=0.0D0
              CALL CGAMA(A,Y,0,TAR,TAI)
              CTA = DCMPLX(TAR, TAI)
              Y=0.0D0
              CALL CGAMA(B,Y,0,TBR,TBI)
              CTB = DCMPLX(TBR, TBI)
              XG=B-A
              Y=0.0D0
              CALL CGAMA(XG,Y,0,TBAR,TBAI)
              CTBA = DCMPLX(TBAR, TBAI)
              SUM1=1.0D0
              SUM2=1.0D0
              R1=1.0D0
              R2=1.0D0
              DO 20 I=1,8
                 R1=-R1*(A+I-1.0D0)*(A-B+I)/(X*I)
                 R2=-R2*(B-A+I-1.0D0)*(A-I)/(X*I)
                 SUM1=SUM1+R1
20               SUM2=SUM2+R2
              HG1=DBLE(CDEXP(CTB-CTBA))*X**(-A)*DCOS(PI*A)*SUM1
              HG2=DBLE(CDEXP(CTB-CTA+X))*X**(A-B)*SUM2
              HG=HG1+HG2
           ENDIF
25         IF (N.EQ.0) Y0=HG
           IF (N.EQ.1) Y1=HG
30      CONTINUE
        IF (A0.GE.2.0D0) THEN
           DO 35 I=1,LA-1
              HG=((2.0D0*A-B+X)*Y1+(B-A)*Y0)/A
              Y0=Y1
              Y1=HG
35            A=A+1.0D0
        ENDIF
        IF (X0.LT.0.0D0) HG=HG*DEXP(X0)
        A=A1
        X=X0
        RETURN
        END



C       **********************************

        SUBROUTINE STVH0(X,SH0)
C
C       =============================================
C       Purpose: Compute Struve function H0(x)
C       Input :  x   --- Argument of H0(x) ( x ≥ 0 )
C       Output:  SH0 --- H0(x)
C       =============================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        S=1.0D0
        R=1.0D0
        IF (X.LE.20.0D0) THEN
           A0=2.0*X/PI
           DO 10 K=1,60
              R=-R*X/(2.0D0*K+1.0D0)*X/(2.0D0*K+1.0D0)
              S=S+R
              IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 15
10         CONTINUE
15         SH0=A0*S
        ELSE
           KM=INT(.5*(X+1.0))
           IF (X.GE.50.0) KM=25
           DO 20 K=1,KM
              R=-R*((2.0D0*K-1.0D0)/X)**2
              S=S+R
              IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 25
20         CONTINUE
25         T=4.0D0/X
           T2=T*T
           P0=((((-.37043D-5*T2+.173565D-4)*T2-.487613D-4)
     &        *T2+.17343D-3)*T2-.1753062D-2)*T2+.3989422793D0
           Q0=T*(((((.32312D-5*T2-.142078D-4)*T2+.342468D-4)*
     &        T2-.869791D-4)*T2+.4564324D-3)*T2-.0124669441D0)
           TA0=X-.25D0*PI
           BY0=2.0D0/DSQRT(X)*(P0*DSIN(TA0)+Q0*DCOS(TA0))
           SH0=2.0D0/(PI*X)*S+BY0
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE HYGFX(A,B,C,X,HF)
C
C       ====================================================
C       Purpose: Compute hypergeometric function F(a,b,c,x)
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter, c <> 0,-1,-2,...
C                x --- Argument   ( x < 1 )
C       Output:  HF --- F(a,b,c,x)
C       Routines called:
C            (1) GAMMA2 for computing gamma function
C            (2) PSI_SPEC for computing psi function
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL L0,L1,L2,L3,L4,L5
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        L0=C.EQ.INT(C).AND.C.LT.0.0
        L1=1.0D0-X.LT.1.0D-15.AND.C-A-B.LE.0.0
        L2=A.EQ.INT(A).AND.A.LT.0.0
        L3=B.EQ.INT(B).AND.B.LT.0.0
        L4=C-A.EQ.INT(C-A).AND.C-A.LE.0.0
        L5=C-B.EQ.INT(C-B).AND.C-B.LE.0.0
        IF (L0.OR.L1) THEN
           WRITE(*,*)'The hypergeometric series is divergent'
           RETURN
        ENDIF
        EPS=1.0D-15
        IF (X.GT.0.95) EPS=1.0D-8
        IF (X.EQ.0.0.OR.A.EQ.0.0.OR.B.EQ.0.0) THEN
           HF=1.0D0
           RETURN
        ELSE IF (1.0D0-X.EQ.EPS.AND.C-A-B.GT.0.0) THEN
           CALL GAMMA2(C,GC)
           CALL GAMMA2(C-A-B,GCAB)
           CALL GAMMA2(C-A,GCA)
           CALL GAMMA2(C-B,GCB)
           HF=GC*GCAB/(GCA*GCB)
           RETURN
        ELSE IF (1.0D0+X.LE.EPS.AND.DABS(C-A+B-1.0).LE.EPS) THEN
           G0=DSQRT(PI)*2.0D0**(-A)
           CALL GAMMA2(C,G1)
           CALL GAMMA2(1.0D0+A/2.0-B,G2)
           CALL GAMMA2(0.5D0+0.5*A,G3)
           HF=G0*G1/(G2*G3)
           RETURN
        ELSE IF (L2.OR.L3) THEN
           IF (L2) NM=INT(ABS(A))
           IF (L3) NM=INT(ABS(B))
           HF=1.0D0
           R=1.0D0
           DO 10 K=1,NM
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
10            HF=HF+R
           RETURN
        ELSE IF (L4.OR.L5) THEN
           IF (L4) NM=INT(ABS(C-A))
           IF (L5) NM=INT(ABS(C-B))
           HF=1.0D0
           R=1.0D0
           DO 15 K=1,NM
              R=R*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*X
15            HF=HF+R
           HF=(1.0D0-X)**(C-A-B)*HF
           RETURN
        ENDIF
        AA=A
        BB=B
        X1=X
        IF (X.LT.0.0D0) THEN
           X=X/(X-1.0D0)
           IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0) THEN
              A=BB
              B=AA
           ENDIF
           B=C-B
        ENDIF
        HW=0.0D0
        IF (X.GE.0.75D0) THEN
           GM=0.0D0
           IF (DABS(C-A-B-INT(C-A-B)).LT.1.0D-15) THEN
              M=INT(C-A-B)
              CALL GAMMA2(A,GA)
              CALL GAMMA2(B,GB)
              CALL GAMMA2(C,GC)
              CALL GAMMA2(A+M,GAM)
              CALL GAMMA2(B+M,GBM)
              CALL PSI_SPEC(A,PA)
              CALL PSI_SPEC(B,PB)
              IF (M.NE.0) GM=1.0D0
              DO 30 J=1,ABS(M)-1
30               GM=GM*J
              RM=1.0D0
              DO 35 J=1,ABS(M)
35               RM=RM*J
              F0=1.0D0
              R0=1.0D0
              R1=1.0D0
              SP0=0.D0
              SP=0.0D0
              IF (M.GE.0) THEN
                 C0=GM*GC/(GAM*GBM)
                 C1=-GC*(X-1.0D0)**M/(GA*GB*RM)
                 DO 40 K=1,M-1
                    R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(K-M))*(1.0-X)
40                  F0=F0+R0
                 DO 45 K=1,M
45                  SP0=SP0+1.0D0/(A+K-1.0)+1.0/(B+K-1.0)-1.0/K
                 F1=PA+PB+SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 55 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 50 J=1,M
50                     SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0))+1.0/
     &                    (B+J+K-1.0)
                    RP=PA+PB+2.0D0*EL+SP+SM+DLOG(1.0D0-X)
                    R1=R1*(A+M+K-1.0D0)*(B+M+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 60
55                  HW=F1
60               HF=F0*C0+F1*C1
              ELSE IF (M.LT.0) THEN
                 M=-M
                 C0=GM*GC/(GA*GB*(1.0D0-X)**M)
                 C1=-(-1)**M*GC/(GAM*GBM*RM)
                 DO 65 K=1,M-1
                    R0=R0*(A-M+K-1.0D0)*(B-M+K-1.0)/(K*(K-M))*(1.0-X)
65                  F0=F0+R0
                 DO 70 K=1,M
70                  SP0=SP0+1.0D0/K
                 F1=PA+PB-SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 80 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 75 J=1,M
75                     SM=SM+1.0D0/(J+K)
                    RP=PA+PB+2.0D0*EL+SP-SM+DLOG(1.0D0-X)
                    R1=R1*(A+K-1.0D0)*(B+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 85
80                  HW=F1
85               HF=F0*C0+F1*C1
              ENDIF
           ELSE
              CALL GAMMA2(A,GA)
              CALL GAMMA2(B,GB)
              CALL GAMMA2(C,GC)
              CALL GAMMA2(C-A,GCA)
              CALL GAMMA2(C-B,GCB)
              CALL GAMMA2(C-A-B,GCAB)
              CALL GAMMA2(A+B-C,GABC)
              C0=GC*GCAB/(GCA*GCB)
              C1=GC*GABC/(GA*GB)*(1.0D0-X)**(C-A-B)
              HF=0.0D0
              R0=C0
              R1=C1
              DO 90 K=1,250
                 R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(A+B-C+K))*(1.0-X)
                 R1=R1*(C-A+K-1.0D0)*(C-B+K-1.0)/(K*(C-A-B+K))
     &              *(1.0-X)
                 HF=HF+R0+R1
                 IF (DABS(HF-HW).LT.DABS(HF)*EPS) GO TO 95
90               HW=HF
95            HF=HF+C0+C1
           ENDIF
        ELSE
           A0=1.0D0
           IF (C.GT.A.AND.C.LT.2.0D0*A.AND.
     &         C.GT.B.AND.C.LT.2.0D0*B) THEN
              A0=(1.0D0-X)**(C-A-B)
              A=C-A
              B=C-B
           ENDIF
           HF=1.0D0
           R=1.0D0
           DO 100 K=1,250
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
              HF=HF+R
              IF (DABS(HF-HW).LE.DABS(HF)*EPS) GO TO 105
100           HW=HF
105        HF=A0*HF
        ENDIF
        IF (X1.LT.0.0D0) THEN
           X=X1
           C0=1.0D0/(1.0D0-X)**AA
           HF=C0*HF
        ENDIF
        A=AA
        B=BB
        IF (K.GT.120) WRITE(*,115)
115     FORMAT(1X,'Warning! You should check the accuracy')
        RETURN
        END



C       **********************************

        SUBROUTINE CCHG(A,B,Z,CHG)
C
C       ===================================================
C       Purpose: Compute confluent hypergeometric function
C                M(a,b,z) with real parameters a, b and a
C                complex argument z
C       Input :  a --- Parameter
C                b --- Parameter
C                z --- Complex argument
C       Output:  CHG --- M(a,b,z)
C       Routine called: CGAMA for computing complex ln[Г(x)]
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX *16 (C,Z)
        PI=3.141592653589793D0
        CI=(0.0D0,1.0D0)
        A0=A
        A1=A
        Z0=Z
        IF (B.EQ.0.0.OR.B.EQ.-INT(ABS(B))) THEN
           CHG=(1.0D+300,0.0D0)
        ELSE IF (A.EQ.0.0D0.OR.Z.EQ.0.0D0) THEN
           CHG=(1.0D0,0.0D0)
        ELSE IF (A.EQ.-1.0D0) THEN
           CHG=1.0D0-Z/B
        ELSE IF (A.EQ.B) THEN
           CHG=CDEXP(Z)
        ELSE IF (A-B.EQ.1.0D0) THEN
           CHG=(1.0D0+Z/B)*CDEXP(Z)
        ELSE IF (A.EQ.1.0D0.AND.B.EQ.2.0D0) THEN
           CHG=(CDEXP(Z)-1.0D0)/Z
        ELSE IF (A.EQ.INT(A).AND.A.LT.0.0D0) THEN
           M=INT(-A)
           CR=(1.0D0,0.0D0)
           CHG=(1.0D0,0.0D0)
           DO 10 K=1,M
              CR=CR*(A+K-1.0D0)/K/(B+K-1.0D0)*Z
10            CHG=CHG+CR
        ELSE
           X0=DBLE(Z)
           IF (X0.LT.0.0D0) THEN
              A=B-A
              A0=A
              Z=-Z
           ENDIF
           NL=0
           LA=0
           IF (A.GE.2.0D0) THEN
              NL=1
              LA=INT(A)
              A=A-LA-1.0D0
           ENDIF
           NS=0
           DO 30 N=0,NL
              IF (A0.GE.2.0D0) A=A+1.0D0
              IF (CDABS(Z).LT.20.0D0+ABS(B).OR.A.LT.0.0D0) THEN
                 CHG=(1.0D0,0.0D0)
                 CRG=(1.0D0,0.0D0)
                 DO 15 J=1,500
                    CRG=CRG*(A+J-1.0D0)/(J*(B+J-1.0D0))*Z
                    CHG=CHG+CRG
                    IF (CDABS((CHG-CHW)/CHG).LT.1.D-15) GO TO 25
                    CHW=CHG
15               CONTINUE
              ELSE
                 Y=0.0D0
                 CALL CGAMA(A,Y,0,G1R,G1I)
                 CG1 = DCMPLX(G1R, G1I)
                 Y=0.0D0
                 CALL CGAMA(B,Y,0,G2R,G2I)
                 CG2 = DCMPLX(G2R,G2I)
                 BA=B-A
                 Y=0.0D0
                 CALL CGAMA(BA,Y,0,G3R,G3I)
                 CG3 = DCMPLX(G3R, G3I)
                 CS1=(1.0D0,0.0D0)
                 CS2=(1.0D0,0.0D0)
                 CR1=(1.0D0,0.0D0)
                 CR2=(1.0D0,0.0D0)
                 DO 20 I=1,8
                    CR1=-CR1*(A+I-1.0D0)*(A-B+I)/(Z*I)
                    CR2=CR2*(B-A+I-1.0D0)*(I-A)/(Z*I)
                    CS1=CS1+CR1
20                  CS2=CS2+CR2
                 X=DBLE(Z)
                 Y=DIMAG(Z)
                 IF (X.EQ.0.0.AND.Y.GE.0.0) THEN
                    PHI=0.5D0*PI
                 ELSE IF (X.EQ.0.0.AND.Y.LE.0.0) THEN
                    PHI=-0.5D0*PI
                 ELSE
                    PHI=DATAN(Y/X)
                 ENDIF
                 IF (PHI.GT.-0.5*PI.AND.PHI.LT.1.5*PI) NS=1
                 IF (PHI.GT.-1.5*PI.AND.PHI.LE.-0.5*PI) NS=-1
                 CFAC=CDEXP(NS*CI*PI*A)
                 IF (Y.EQ.0.0D0) CFAC=DCOS(PI*A)
                 CHG1=CDEXP(CG2-CG3)*Z**(-A)*CFAC*CS1
                 CHG2=CDEXP(CG2-CG1+Z)*Z**(A-B)*CS2
                 CHG=CHG1+CHG2
              ENDIF
25            IF (N.EQ.0) CY0=CHG
              IF (N.EQ.1) CY1=CHG
30         CONTINUE
           IF (A0.GE.2.0D0) THEN
              DO 35 I=1,LA-1
                 CHG=((2.0D0*A-B+Z)*CY1+(B-A)*CY0)/A
                 CY0=CY1
                 CY1=CHG
35               A=A+1.0D0
           ENDIF
           IF (X0.LT.0.0D0) CHG=CHG*CDEXP(-Z)
        ENDIF
        A=A1
        Z=Z0
        RETURN
        END



C       **********************************

        SUBROUTINE HYGFZ(A,B,C,Z,ZHF)
C
C       ======================================================
C       Purpose: Compute the hypergeometric function for a
C                complex argument, F(a,b,c,z)
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter,  c <> 0,-1,-2,...
C                z --- Complex argument
C       Output:  ZHF --- F(a,b,c,z)
C       Routines called:
C            (1) GAMMA2 for computing gamma function
C            (2) PSI_SPEC for computing psi function
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Y)
        IMPLICIT COMPLEX *16 (Z)
        LOGICAL L0,L1,L2,L3,L4,L5,L6
        X=DBLE(Z)
        Y=DIMAG(Z)
        EPS=1.0D-15
        L0=C.EQ.INT(C).AND.C.LT.0.0D0
        L1=DABS(1.0D0-X).LT.EPS.AND.Y.EQ.0.0D0.AND.C-A-B.LE.0.0D0
        L2=CDABS(Z+1.0D0).LT.EPS.AND.DABS(C-A+B-1.0D0).LT.EPS
        L3=A.EQ.INT(A).AND.A.LT.0.0D0
        L4=B.EQ.INT(B).AND.B.LT.0.0D0
        L5=C-A.EQ.INT(C-A).AND.C-A.LE.0.0D0
        L6=C-B.EQ.INT(C-B).AND.C-B.LE.0.0D0
        AA=A
        BB=B
        A0=CDABS(Z)
        IF (A0.GT.0.95D0) EPS=1.0D-8
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        IF (L0.OR.L1) THEN
C           WRITE(*,*)'The hypergeometric series is divergent'
           ZHF = 1.0D300
           RETURN
        ENDIF
        NM=0
        IF (A0.EQ.0.0D0.OR.A.EQ.0.0D0.OR.B.EQ.0.0D0) THEN
           ZHF=(1.0D0,0.0D0)
        ELSE IF (Z.EQ.1.0D0.AND.C-A-B.GT.0.0D0) THEN
           CALL GAMMA2(C,GC)
           CALL GAMMA2(C-A-B,GCAB)
           CALL GAMMA2(C-A,GCA)
           CALL GAMMA2(C-B,GCB)
           ZHF=GC*GCAB/(GCA*GCB)
        ELSE IF (L2) THEN
           G0=DSQRT(PI)*2.0D0**(-A)
           CALL GAMMA2(C,G1)
           CALL GAMMA2(1.0D0+A/2.0D0-B,G2)
           CALL GAMMA2(0.5D0+0.5D0*A,G3)
           ZHF=G0*G1/(G2*G3)
        ELSE IF (L3.OR.L4) THEN
           IF (L3) NM=INT(ABS(A))
           IF (L4) NM=INT(ABS(B))
           ZHF=(1.0D0,0.0D0)
           ZR=(1.0D0,0.0D0)
           DO 10 K=1,NM
              ZR=ZR*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*Z
10            ZHF=ZHF+ZR
        ELSE IF (L5.OR.L6) THEN
           IF (L5) NM=INT(ABS(C-A))
           IF (L6) NM=INT(ABS(C-B))
           ZHF=(1.0D0,0.0D0)
           ZR=(1.0D0,0.0D0)
           DO 15 K=1,NM
              ZR=ZR*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*Z
15            ZHF=ZHF+ZR
           ZHF=(1.0D0-Z)**(C-A-B)*ZHF
        ELSE IF (A0.LE.1.0D0) THEN
           IF (X.LT.0.0D0) THEN
              Z1=Z/(Z-1.0D0)
              IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0) THEN
                 A=BB
                 B=AA
              ENDIF
              ZC0=1.0D0/((1.0D0-Z)**A)
              ZHF=(1.0D0,0.0D0)
              ZR0=(1.0D0,0.0D0)
              DO 20 K=1,500
                 ZR0=ZR0*(A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*Z1
                 ZHF=ZHF+ZR0
                 IF (CDABS(ZHF-ZW).LT.CDABS(ZHF)*EPS) GO TO 25
20               ZW=ZHF
25            ZHF=ZC0*ZHF
           ELSE IF (A0.GE.0.90D0) THEN
              GM=0.0D0
              MCAB=INT(C-A-B+EPS*DSIGN(1.0D0,C-A-B))
              IF (DABS(C-A-B-MCAB).LT.EPS) THEN
                 M=INT(C-A-B)
                 CALL GAMMA2(A,GA)
                 CALL GAMMA2(B,GB)
                 CALL GAMMA2(C,GC)
                 CALL GAMMA2(A+M,GAM)
                 CALL GAMMA2(B+M,GBM)
                 CALL PSI_SPEC(A,PA)
                 CALL PSI_SPEC(B,PB)
                 IF (M.NE.0) GM=1.0D0
                 DO 30 J=1,ABS(M)-1
30                  GM=GM*J
                 RM=1.0D0
                 DO 35 J=1,ABS(M)
35                  RM=RM*J
                 ZF0=(1.0D0,0.0D0)
                 ZR0=(1.0D0,0.0D0)
                 ZR1=(1.0D0,0.0D0)
                 SP0=0.D0
                 SP=0.0D0
                 IF (M.GE.0) THEN
                    ZC0=GM*GC/(GAM*GBM)
                    ZC1=-GC*(Z-1.0D0)**M/(GA*GB*RM)
                    DO 40 K=1,M-1
                       ZR0=ZR0*(A+K-1.D0)*(B+K-1.D0)/(K*(K-M))*(1.D0-Z)
40                     ZF0=ZF0+ZR0
                    DO 45 K=1,M
45                     SP0=SP0+1.0D0/(A+K-1.0D0)+1.0/(B+K-1.0D0)-1.D0/K
                    ZF1=PA+PB+SP0+2.0D0*EL+CDLOG(1.0D0-Z)
                    DO 55 K=1,500
                       SP=SP+(1.0D0-A)/(K*(A+K-1.0D0))+(1.0D0-B)/
     &                    (K*(B+K-1.0D0))
                       SM=0.0D0
                       DO 50 J=1,M
                          SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0D0))
     &                       +1.0D0/(B+J+K-1.0D0)
50                     CONTINUE
                       ZP=PA+PB+2.0D0*EL+SP+SM+CDLOG(1.0D0-Z)
                       ZR1=ZR1*(A+M+K-1.0D0)*(B+M+K-1.0D0)/(K*(M+K))
     &                     *(1.0D0-Z)
                       ZF1=ZF1+ZR1*ZP
                       IF (CDABS(ZF1-ZW).LT.CDABS(ZF1)*EPS) GO TO 60
55                     ZW=ZF1
60                  ZHF=ZF0*ZC0+ZF1*ZC1
                 ELSE IF (M.LT.0) THEN
                    M=-M
                    ZC0=GM*GC/(GA*GB*(1.0D0-Z)**M)
                    ZC1=-(-1)**M*GC/(GAM*GBM*RM)
                    DO 65 K=1,M-1
                       ZR0=ZR0*(A-M+K-1.0D0)*(B-M+K-1.0D0)/(K*(K-M))
     &                     *(1.0D0-Z)
65                     ZF0=ZF0+ZR0
                    DO 70 K=1,M
70                     SP0=SP0+1.0D0/K
                    ZF1=PA+PB-SP0+2.0D0*EL+CDLOG(1.0D0-Z)
                    DO 80 K=1,500
                       SP=SP+(1.0D0-A)/(K*(A+K-1.0D0))+(1.0D0-B)/(K*
     &                    (B+K-1.0D0))
                       SM=0.0D0
                       DO 75 J=1,M
75                        SM=SM+1.0D0/(J+K)
                       ZP=PA+PB+2.0D0*EL+SP-SM+CDLOG(1.0D0-Z)
                       ZR1=ZR1*(A+K-1.D0)*(B+K-1.D0)/(K*(M+K))*(1.D0-Z)
                       ZF1=ZF1+ZR1*ZP
                       IF (CDABS(ZF1-ZW).LT.CDABS(ZF1)*EPS) GO TO 85
80                     ZW=ZF1
85                  ZHF=ZF0*ZC0+ZF1*ZC1
                 ENDIF
              ELSE
                 CALL GAMMA2(A,GA)
                 CALL GAMMA2(B,GB)
                 CALL GAMMA2(C,GC)
                 CALL GAMMA2(C-A,GCA)
                 CALL GAMMA2(C-B,GCB)
                 CALL GAMMA2(C-A-B,GCAB)
                 CALL GAMMA2(A+B-C,GABC)
                 ZC0=GC*GCAB/(GCA*GCB)
                 ZC1=GC*GABC/(GA*GB)*(1.0D0-Z)**(C-A-B)
                 ZHF=(0.0D0,0.0D0)
                 ZR0=ZC0
                 ZR1=ZC1
                 DO 90 K=1,500
                    ZR0=ZR0*(A+K-1.D0)*(B+K-1.D0)/(K*(A+B-C+K))*(1.D0-Z)
                    ZR1=ZR1*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C-A-B+K))
     &                  *(1.0D0-Z)
                    ZHF=ZHF+ZR0+ZR1
                    IF (CDABS(ZHF-ZW).LT.CDABS(ZHF)*EPS) GO TO 95
90                  ZW=ZHF
95               ZHF=ZHF+ZC0+ZC1
              ENDIF
           ELSE
              Z00=(1.0D0,0.0D0)
              IF (C-A.LT.A.AND.C-B.LT.B) THEN
                  Z00=(1.0D0-Z)**(C-A-B)
                  A=C-A
                  B=C-B
              ENDIF
              ZHF=(1.0D0,0.D0)
              ZR=(1.0D0,0.0D0)
              DO 100 K=1,1500
                 ZR=ZR*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*Z
                 ZHF=ZHF+ZR
                 IF (CDABS(ZHF-ZW).LE.CDABS(ZHF)*EPS) GO TO 105
100              ZW=ZHF
105           ZHF=Z00*ZHF
           ENDIF
        ELSE IF (A0.GT.1.0D0) THEN
           MAB=INT(A-B+EPS*DSIGN(1.0D0,A-B))
           IF (DABS(A-B-MAB).LT.EPS.AND.A0.LE.1.1D0) B=B+EPS
           IF (DABS(A-B-MAB).GT.EPS) THEN
              CALL GAMMA2(A,GA)
              CALL GAMMA2(B,GB)
              CALL GAMMA2(C,GC)
              CALL GAMMA2(A-B,GAB)
              CALL GAMMA2(B-A,GBA)
              CALL GAMMA2(C-A,GCA)
              CALL GAMMA2(C-B,GCB)
              ZC0=GC*GBA/(GCA*GB*(-Z)**A)
              ZC1=GC*GAB/(GCB*GA*(-Z)**B)
              ZR0=ZC0
              ZR1=ZC1
              ZHF=(0.0D0,0.0D0)
              DO 110 K=1,500
                 ZR0=ZR0*(A+K-1.0D0)*(A-C+K)/((A-B+K)*K*Z)
                 ZR1=ZR1*(B+K-1.0D0)*(B-C+K)/((B-A+K)*K*Z)
                 ZHF=ZHF+ZR0+ZR1
                 IF (CDABS((ZHF-ZW)/ZHF).LE.EPS) GO TO 115
110              ZW=ZHF
115           ZHF=ZHF+ZC0+ZC1
           ELSE
              IF (A-B.LT.0.0D0) THEN
                 A=BB
                 B=AA
              ENDIF
              CA=C-A
              CB=C-B
              NCA=INT(CA+EPS*DSIGN(1.0D0,CA))
              NCB=INT(CB+EPS*DSIGN(1.0D0,CB))
              IF (DABS(CA-NCA).LT.EPS.OR.DABS(CB-NCB).LT.EPS) C=C+EPS
              CALL GAMMA2(A,GA)
              CALL GAMMA2(C,GC)
              CALL GAMMA2(C-B,GCB)
              CALL PSI_SPEC(A,PA)
              CALL PSI_SPEC(C-A,PCA)
              CALL PSI_SPEC(A-C,PAC)
              MAB=INT(A-B+EPS)
              ZC0=GC/(GA*(-Z)**B)
              CALL GAMMA2(A-B,GM)
              ZF0=GM/GCB*ZC0
              ZR=ZC0
              DO 120 K=1,MAB-1
                 ZR=ZR*(B+K-1.0D0)/(K*Z)
                 T0=A-B-K
                 CALL GAMMA2(T0,G0)
                 CALL GAMMA2(C-B-K,GCBK)
120              ZF0=ZF0+ZR*G0/GCBK
              IF (MAB.EQ.0) ZF0=(0.0D0,0.0D0)
              ZC1=GC/(GA*GCB*(-Z)**A)
              SP=-2.0D0*EL-PA-PCA
              DO 125 J=1,MAB
125              SP=SP+1.0D0/J
              ZP0=SP+CDLOG(-Z)
              SQ=1.0D0
              DO 130 J=1,MAB
130              SQ=SQ*(B+J-1.0D0)*(B-C+J)/J
              ZF1=(SQ*ZP0)*ZC1
              ZR=ZC1
              RK1=1.0D0
              SJ1=0.0D0
              W0=0.0D0
              DO 145 K=1,10000
                 ZR=ZR/Z
                 RK1=RK1*(B+K-1.0D0)*(B-C+K)/(K*K)
                 RK2=RK1
                 DO 135 J=K+1,K+MAB
135                 RK2=RK2*(B+J-1.0D0)*(B-C+J)/J
                 SJ1=SJ1+(A-1.0D0)/(K*(A+K-1.0D0))+(A-C-1.0D0)/
     &               (K*(A-C+K-1.0D0))
                 SJ2=SJ1
                 DO 140 J=K+1,K+MAB
140                 SJ2=SJ2+1.0D0/J
                 ZP=-2.0D0*EL-PA-PAC+SJ2-1.0D0/(K+A-C)
     &              -PI/DTAN(PI*(K+A-C))+CDLOG(-Z)
                 ZF1=ZF1+RK2*ZR*ZP
                 WS=CDABS(ZF1)
                 IF (DABS((WS-W0)/WS).LT.EPS) GO TO 150
145              W0=WS
150           ZHF=ZF0+ZF1
           ENDIF
        ENDIF
        A=AA
        B=BB
        IF (K.GT.150) WRITE(*,160)
160     FORMAT(1X,'Warning! You should check the accuracy')
        RETURN
        END



C       **********************************

        SUBROUTINE ITAIRY(X,APT,BPT,ANT,BNT)
C
C       ======================================================
C       Purpose: Compute the integrals of Airy fnctions with
C                respect to t from 0 and x ( x ≥ 0 )
C       Input  : x   --- Upper limit of the integral
C       Output : APT --- Integration of Ai(t) from 0 and x
C                BPT --- Integration of Bi(t) from 0 and x
C                ANT --- Integration of Ai(-t) from 0 and x
C                BNT --- Integration of Bi(-t) from 0 and x
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(16)
        EPS=1.0D-15
        PI=3.141592653589793D0
        C1=.355028053887817D0
        C2=.258819403792807D0
        SR3=1.732050807568877D0
        IF (X.EQ.0.0D0) THEN
           APT=0.0D0
           BPT=0.0D0
           ANT=0.0D0
           BNT=0.0D0
        ELSE
           IF (DABS(X).LE.9.25D0) THEN
              DO 30 L=0,1
                 X=(-1)**L*X
                 FX=X
                 R=X
                 DO 10 K=1,40
                    R=R*(3.0*K-2.0D0)/(3.0*K+1.0D0)*X/(3.0*K)
     &                *X/(3.0*K-1.0D0)*X
                    FX=FX+R
                    IF (DABS(R).LT.DABS(FX)*EPS) GO TO 15
10               CONTINUE
15               GX=.5D0*X*X
                 R=GX
                 DO 20 K=1,40
                    R=R*(3.0*K-1.0D0)/(3.0*K+2.0D0)*X/(3.0*K)
     &                *X/(3.0*K+1.0D0)*X
                    GX=GX+R
                    IF (DABS(R).LT.DABS(GX)*EPS) GO TO 25
20               CONTINUE
25               ANT=C1*FX-C2*GX
                 BNT=SR3*(C1*FX+C2*GX)
                 IF (L.EQ.0) THEN
                    APT=ANT
                    BPT=BNT
                 ELSE
                    ANT=-ANT
                    BNT=-BNT
                    X=-X
                 ENDIF
30            CONTINUE
           ELSE
              DATA A/.569444444444444D0,.891300154320988D0,
     &             .226624344493027D+01,.798950124766861D+01,
     &             .360688546785343D+02,.198670292131169D+03,
     &             .129223456582211D+04,.969483869669600D+04,
     &             .824184704952483D+05,.783031092490225D+06,
     &             .822210493622814D+07,.945557399360556D+08,
     &             .118195595640730D+10,.159564653040121D+11,
     &             .231369166433050D+12,.358622522796969D+13/
              Q2=1.414213562373095D0
              Q0=.3333333333333333D0
              Q1=.6666666666666667D0
              XE=X*DSQRT(X)/1.5D0
              XP6=1.0D0/DSQRT(6.0D0*PI*XE)
              SU1=1.0D0
              R=1.0D0
              XR1=1.0D0/XE
              DO 35 K=1,16
                 R=-R*XR1
35               SU1=SU1+A(K)*R
              SU2=1.0D0
              R=1.0D0
              DO 40 K=1,16
                 R=R*XR1
40               SU2=SU2+A(K)*R
              APT=Q0-DEXP(-XE)*XP6*SU1
              BPT=2.0D0*DEXP(XE)*XP6*SU2
              SU3=1.0D0
              R=1.0D0
              XR2=1.0D0/(XE*XE)
              DO 45 K=1,8
                 R=-R*XR2
45               SU3=SU3+A(2*K)*R
              SU4=A(1)*XR1
              R=XR1
              DO 50 K=1,7
                 R=-R*XR2
50               SU4=SU4+A(2*K+1)*R
              SU5=SU3+SU4
              SU6=SU3-SU4
              ANT=Q1-Q2*XP6*(SU5*DCOS(XE)-SU6*DSIN(XE))
              BNT=Q2*XP6*(SU5*DSIN(XE)+SU6*DCOS(XE))
           ENDIF
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE IKNA(N,X,NM,BI,DI,BK,DK)
C
C       ========================================================
C       Purpose: Compute modified Bessel functions In(x) and
C                Kn(x), and their derivatives
C       Input:   x --- Argument of In(x) and Kn(x) ( x ≥ 0 )
C                n --- Order of In(x) and Kn(x)
C       Output:  BI(n) --- In(x)
C                DI(n) --- In'(x)
C                BK(n) --- Kn(x)
C                DK(n) --- Kn'(x)
C                NM --- Highest order computed
C       Routines called:
C            (1) IK01A for computing I0(x),I1(x),K0(x) & K1(x)
C            (2) MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BI(0:N),DI(0:N),BK(0:N),DK(0:N)
        NM=N
        IF (X.LE.1.0D-100) THEN
           DO 10 K=0,N
              BI(K)=0.0D0
              DI(K)=0.0D0
              BK(K)=1.0D+300
10            DK(K)=-1.0D+300
           BI(0)=1.0D0
           DI(1)=0.5D0
           RETURN
        ENDIF
        CALL IK01A(X,BI0,DI0,BI1,DI1,BK0,DK0,BK1,DK1)
        BI(0)=BI0
        BI(1)=BI1
        BK(0)=BK0
        BK(1)=BK1
        DI(0)=DI0
        DI(1)=DI1
        DK(0)=DK0
        DK(1)=DK1
        IF (N.LE.1) RETURN
        IF (X.GT.40.0.AND.N.LT.INT(0.25*X)) THEN
           H0=BI0
           H1=BI1
           DO 15 K=2,N
           H=-2.0D0*(K-1.0D0)/X*H1+H0
           BI(K)=H
           H0=H1
15         H1=H
        ELSE
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F0=0.0D0
           F1=1.0D-100
           F=0.0D0
           DO 20 K=M,0,-1
              F=2.0D0*(K+1.0D0)*F1/X+F0
              IF (K.LE.NM) BI(K)=F
              F0=F1
20            F1=F
           S0=BI0/F
           DO 25 K=0,NM
25            BI(K)=S0*BI(K)
        ENDIF
        G0=BK0
        G1=BK1
        DO 30 K=2,NM
           G=2.0D0*(K-1.0D0)/X*G1+G0
           BK(K)=G
           G0=G1
30         G1=G
        DO 40 K=2,NM
           DI(K)=BI(K-1)-K/X*BI(K)
40         DK(K)=-BK(K-1)-K/X*BK(K)
        RETURN
        END



C       **********************************

        SUBROUTINE CJYNA(N,Z,NM,CBJ,CDJ,CBY,CDY)
C
C       =======================================================
C       Purpose: Compute Bessel functions Jn(z), Yn(z) and
C                their derivatives for a complex argument
C       Input :  z --- Complex argument of Jn(z) and Yn(z)
C                n --- Order of Jn(z) and Yn(z)
C       Output:  CBJ(n) --- Jn(z)
C                CDJ(n) --- Jn'(z)
C                CBY(n) --- Yn(z)
C                CDY(n) --- Yn'(z)
C                NM --- Highest order computed
C       Rouitines called:
C            (1) CJY01 to calculate J0(z), J1(z), Y0(z), Y1(z)
C            (2) MSTA1 and MSTA2 to calculate the starting
C                point for backward recurrence
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,E,P,R,W,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:N),CDJ(0:N),CBY(0:N),CDY(0:N)
        PI=3.141592653589793D0
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-100) THEN
           DO 5 K=0,N
              CBJ(K)=(0.0D0,0.0D0)
              CDJ(K)=(0.0D0,0.0D0)
              CBY(K)=-(1.0D+300,0.0D0)
5             CDY(K)=(1.0D+300,0.0D0)
           CBJ(0)=(1.0D0,0.0D0)
           CDJ(1)=(0.5D0,0.0D0)
           RETURN
        ENDIF
        CALL CJY01(Z,CBJ0,CDJ0,CBJ1,CDJ1,CBY0,CDY0,CBY1,CDY1)
        CBJ(0)=CBJ0
        CBJ(1)=CBJ1
        CBY(0)=CBY0
        CBY(1)=CBY1
        CDJ(0)=CDJ0
        CDJ(1)=CDJ1
        CDY(0)=CDY0
        CDY(1)=CDY1
        IF (N.LE.1) RETURN
        IF (N.LT.INT(0.25*A0)) THEN
           CJ0=CBJ0
           CJ1=CBJ1
           DO 70 K=2,N
              CJK=2.0D0*(K-1.0D0)/Z*CJ1-CJ0
              CBJ(K)=CJK
              CJ0=CJ1
70            CJ1=CJK
        ELSE
           M=MSTA1(A0,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(A0,N,15)
           ENDIF
           CF2=(0.0D0,0.0D0)
           CF1=(1.0D-100,0.0D0)
           DO 75 K=M,0,-1
              CF=2.0D0*(K+1.0D0)/Z*CF1-CF2
              IF (K.LE.NM) CBJ(K)=CF
              CF2=CF1
75            CF1=CF
           IF (CDABS(CBJ0).GT.CDABS(CBJ1)) THEN
              CS=CBJ0/CF
           ELSE
              CS=CBJ1/CF2
           ENDIF
           DO 80 K=0,NM
80            CBJ(K)=CS*CBJ(K)
        ENDIF
        DO 85 K=2,NM
85         CDJ(K)=CBJ(K-1)-K/Z*CBJ(K)
        YA0=CDABS(CBY0)
        LB=0
        LB0=0
        CG0=CBY0
        CG1=CBY1
        DO 90 K=2,NM
           CYK=2.0D0*(K-1.0D0)/Z*CG1-CG0
           IF (CDABS(CYK).GT.1.0D+290) GO TO 90
           YAK=CDABS(CYK)
           YA1=CDABS(CG0)
           IF (YAK.LT.YA0.AND.YAK.LT.YA1) LB=K
           CBY(K)=CYK
           CG0=CG1
           CG1=CYK
90      CONTINUE
        IF (LB.LE.4.OR.DIMAG(Z).EQ.0.0D0) GO TO 125
95      IF (LB.EQ.LB0) GO TO 125
        CH2=(1.0D0,0.0D0)
        CH1=(0.0D0,0.0D0)
        LB0=LB
        DO 100 K=LB,1,-1
           CH0=2.0D0*K/Z*CH1-CH2
           CH2=CH1
100        CH1=CH0
        CP12=CH0
        CP22=CH2
        CH2=(0.0D0,0.0D0)
        CH1=(1.0D0,0.0D0)
        DO 105 K=LB,1,-1
           CH0=2.0D0*K/Z*CH1-CH2
           CH2=CH1
105        CH1=CH0
        CP11=CH0
        CP21=CH2
        IF (LB.EQ.NM) CBJ(LB+1)=2.0D0*LB/Z*CBJ(LB)-CBJ(LB-1)
        IF (CDABS(CBJ(0)).GT.CDABS(CBJ(1))) THEN
           CBY(LB+1)=(CBJ(LB+1)*CBY0-2.0D0*CP11/(PI*Z))/CBJ(0)
           CBY(LB)=(CBJ(LB)*CBY0+2.0D0*CP12/(PI*Z))/CBJ(0)
        ELSE
           CBY(LB+1)=(CBJ(LB+1)*CBY1-2.0D0*CP21/(PI*Z))/CBJ(1)
           CBY(LB)=(CBJ(LB)*CBY1+2.0D0*CP22/(PI*Z))/CBJ(1)
        ENDIF
        CYL2=CBY(LB+1)
        CYL1=CBY(LB)
        DO 110 K=LB-1,0,-1
           CYLK=2.0D0*(K+1.0D0)/Z*CYL1-CYL2
           CBY(K)=CYLK
           CYL2=CYL1
110        CYL1=CYLK
        CYL1=CBY(LB)
        CYL2=CBY(LB+1)
        DO 115 K=LB+1,NM-1
           CYLK=2.0D0*K/Z*CYL2-CYL1
           CBY(K+1)=CYLK
           CYL1=CYL2
115        CYL2=CYLK
        DO 120 K=2,NM
           WA=CDABS(CBY(K))
           IF (WA.LT.CDABS(CBY(K-1))) LB=K
120     CONTINUE
        GO TO 95
125     CONTINUE
        DO 130 K=2,NM
130        CDY(K)=CBY(K-1)-K/Z*CBY(K)
        RETURN
        END



C       **********************************

        SUBROUTINE CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)
C
C       =======================================================
C       Purpose: Compute Bessel functions Jn(z), Yn(z) and
C                their derivatives for a complex argument
C       Input :  z --- Complex argument of Jn(z) and Yn(z)
C                n --- Order of Jn(z) and Yn(z)
C       Output:  CBJ(n) --- Jn(z)
C                CDJ(n) --- Jn'(z)
C                CBY(n) --- Yn(z)
C                CDY(n) --- Yn'(z)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 to calculate the starting
C                point for backward recurrence
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:N),CDJ(0:N),CBY(0:N),CDY(0:N),
     &            A(4),B(4),A1(4),B1(4)
        EL=0.5772156649015329D0
        PI=3.141592653589793D0
        R2P=.63661977236758D0
        Y0=DABS(DIMAG(Z))
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBJ(K)=(0.0D0,0.0D0)
              CDJ(K)=(0.0D0,0.0D0)
              CBY(K)=-(1.0D+300,0.0D0)
10            CDY(K)=(1.0D+300,0.0D0)
           CBJ(0)=(1.0D0,0.0D0)
           CDJ(1)=(0.5D0,0.0D0)
           RETURN
        ENDIF
        IF (A0.LE.300.D0.OR.N.GT.80) THEN
           IF (N.EQ.0) NM=1
           M=MSTA1(A0,200)
           IF (M.LT.NM) THEN
              NM=M
           ELSE
              M=MSTA2(A0,NM,15)
           ENDIF
           CBS=(0.0D0,0.0D0)
           CSU=(0.0D0,0.0D0)
           CSV=(0.0D0,0.0D0)
           CF2=(0.0D0,0.0D0)
           CF1=(1.0D-100,0.0D0)
           DO 15 K=M,0,-1
              CF=2.0D0*(K+1.0D0)/Z*CF1-CF2
              IF (K.LE.NM) CBJ(K)=CF
              IF (K.EQ.2*INT(K/2).AND.K.NE.0) THEN
                 IF (Y0.LE.1.0D0) THEN
                    CBS=CBS+2.0D0*CF
                 ELSE
                    CBS=CBS+(-1)**(K/2)*2.0D0*CF
                 ENDIF
                 CSU=CSU+(-1)**(K/2)*CF/K
              ELSE IF (K.GT.1) THEN
                 CSV=CSV+(-1)**(K/2)*K/(K*K-1.0D0)*CF
              ENDIF
              CF2=CF1
15            CF1=CF
           IF (Y0.LE.1.0D0) THEN
              CS0=CBS+CF
           ELSE
              CS0=(CBS+CF)/CDCOS(Z)
           ENDIF
           DO 20 K=0,NM
20            CBJ(K)=CBJ(K)/CS0
           CE=CDLOG(Z/2.0D0)+EL
           CBY(0)=R2P*(CE*CBJ(0)-4.0D0*CSU/CS0)
           CBY(1)=R2P*(-CBJ(0)/Z+(CE-1.0D0)*CBJ(1)-4.0D0*CSV/CS0)
        ELSE
           DATA A/-.7031250000000000D-01,.1121520996093750D+00,
     &            -.5725014209747314D+00,.6074042001273483D+01/
           DATA B/ .7324218750000000D-01,-.2271080017089844D+00,
     &             .1727727502584457D+01,-.2438052969955606D+02/
           DATA A1/.1171875000000000D+00,-.1441955566406250D+00,
     &             .6765925884246826D+00,-.6883914268109947D+01/
           DATA B1/-.1025390625000000D+00,.2775764465332031D+00,
     &             -.1993531733751297D+01,.2724882731126854D+02/
           CT1=Z-0.25D0*PI
           CP0=(1.0D0,0.0D0)
           DO 25 K=1,4
25            CP0=CP0+A(K)*Z**(-2*K)
           CQ0=-0.125D0/Z
           DO 30 K=1,4
30            CQ0=CQ0+B(K)*Z**(-2*K-1)
           CU=CDSQRT(R2P/Z)
           CBJ0=CU*(CP0*CDCOS(CT1)-CQ0*CDSIN(CT1))
           CBY0=CU*(CP0*CDSIN(CT1)+CQ0*CDCOS(CT1))
           CBJ(0)=CBJ0
           CBY(0)=CBY0
           CT2=Z-0.75D0*PI
           CP1=(1.0D0,0.0D0)
           DO 35 K=1,4
35            CP1=CP1+A1(K)*Z**(-2*K)
           CQ1=0.375D0/Z
           DO 40 K=1,4
40            CQ1=CQ1+B1(K)*Z**(-2*K-1)
           CBJ1=CU*(CP1*CDCOS(CT2)-CQ1*CDSIN(CT2))
           CBY1=CU*(CP1*CDSIN(CT2)+CQ1*CDCOS(CT2))
           CBJ(1)=CBJ1
           CBY(1)=CBY1
           DO 45 K=2,NM
              CBJK=2.0D0*(K-1.0D0)/Z*CBJ1-CBJ0
              CBJ(K)=CBJK
              CBJ0=CBJ1
45            CBJ1=CBJK
        ENDIF
        CDJ(0)=-CBJ(1)
        DO 50 K=1,NM
50         CDJ(K)=CBJ(K-1)-K/Z*CBJ(K)
        IF (CDABS(CBJ(0)).GT.1.0D0) THEN
           CBY(1)=(CBJ(1)*CBY(0)-2.0D0/(PI*Z))/CBJ(0)
        ENDIF
        DO 55 K=2,NM
           IF (CDABS(CBJ(K-1)).GE.CDABS(CBJ(K-2))) THEN
              CYY=(CBJ(K)*CBY(K-1)-2.0D0/(PI*Z))/CBJ(K-1)
           ELSE
              CYY=(CBJ(K)*CBY(K-2)-4.0D0*(K-1.0D0)/(PI*Z*Z))/CBJ(K-2)
           ENDIF
           CBY(K)=CYY
55      CONTINUE
        CDY(0)=-CBY(1)
        DO 60 K=1,NM
60         CDY(K)=CBY(K-1)-K/Z*CBY(K)
        RETURN
        END



C       **********************************

        SUBROUTINE IKNB(N,X,NM,BI,DI,BK,DK)
C
C       ============================================================
C       Purpose: Compute modified Bessel functions In(x) and Kn(x),
C                and their derivatives
C       Input:   x --- Argument of In(x) and Kn(x) ( 0 ≤ x ≤ 700 )
C                n --- Order of In(x) and Kn(x)
C       Output:  BI(n) --- In(x)
C                DI(n) --- In'(x)
C                BK(n) --- Kn(x)
C                DK(n) --- Kn'(x)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting point
C                for backward recurrence
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BI(0:N),DI(0:N),BK(0:N),DK(0:N)
        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        NM=N
        IF (X.LE.1.0D-100) THEN
           DO 10 K=0,N
              BI(K)=0.0D0
              DI(K)=0.0D0
              BK(K)=1.0D+300
10            DK(K)=-1.0D+300
           BI(0)=1.0D0
           DI(1)=0.5D0
           RETURN
        ENDIF
        IF (N.EQ.0) NM=1
        M=MSTA1(X,200)
        IF (M.LT.NM) THEN
           NM=M
        ELSE
           M=MSTA2(X,NM,15)
        ENDIF
        BS=0.0D0
        SK0=0.0D0
        F=0.0D0
        F0=0.0D0
        F1=1.0D-100
        DO 15 K=M,0,-1
           F=2.0D0*(K+1.0D0)/X*F1+F0
           IF (K.LE.NM) BI(K)=F
           IF (K.NE.0.AND.K.EQ.2*INT(K/2)) SK0=SK0+4.0D0*F/K
           BS=BS+2.0D0*F
           F0=F1
15         F1=F
        S0=DEXP(X)/(BS-F)
        DO 20 K=0,NM
20         BI(K)=S0*BI(K)
        IF (X.LE.8.0D0) THEN
           BK(0)=-(DLOG(0.5D0*X)+EL)*BI(0)+S0*SK0
           BK(1)=(1.0D0/X-BI(1)*BK(0))/BI(0)
        ELSE
           A0=DSQRT(PI/(2.0D0*X))*DEXP(-X)
           K0=16
           IF (X.GE.25.0) K0=10
           IF (X.GE.80.0) K0=8
           IF (X.GE.200.0) K0=6
           DO 30 L=0,1
              BKL=1.0D0
              VT=4.0D0*L
              R=1.0D0
              DO 25 K=1,K0
                 R=0.125D0*R*(VT-(2.0*K-1.0)**2)/(K*X)
25               BKL=BKL+R
              BK(L)=A0*BKL
30         CONTINUE
        ENDIF
        G0=BK(0)
        G1=BK(1)
        DO 35 K=2,NM
           G=2.0D0*(K-1.0D0)/X*G1+G0
           BK(K)=G
           G0=G1
35         G1=G
        DI(0)=BI(1)
        DK(0)=-BK(1)
        DO 40 K=1,NM
           DI(K)=BI(K-1)-K/X*BI(K)
40         DK(K)=-BK(K-1)-K/X*BK(K)
        RETURN
        END



C       **********************************

        SUBROUTINE LPMN(MM,M,N,X,PM,PD)
C
C       =====================================================
C       Purpose: Compute the associated Legendre functions
C                Pmn(x) and their derivatives Pmn'(x) for
C                real argument
C       Input :  x  --- Argument of Pmn(x)
C                m  --- Order of Pmn(x),  m = 0,1,2,...,n
C                n  --- Degree of Pmn(x), n = 0,1,2,...,N
C                mm --- Physical dimension of PM and PD
C       Output:  PM(m,n) --- Pmn(x)
C                PD(m,n) --- Pmn'(x)
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (D,P,X)
        DIMENSION PM(0:MM,0:N),PD(0:MM,0:N)
        INTRINSIC MIN
        DO 10 I=0,N
        DO 10 J=0,M
           PM(J,I)=0.0D0
10         PD(J,I)=0.0D0
        PM(0,0)=1.0D0
        IF (N.EQ.0) RETURN
        IF (DABS(X).EQ.1.0D0) THEN
           DO 15 I=1,N
              PM(0,I)=X**I
15            PD(0,I)=0.5D0*I*(I+1.0D0)*X**(I+1)
           DO 20 J=1,N
           DO 20 I=1,M
              IF (I.EQ.1) THEN
                 PD(I,J)=DINF()
              ELSE IF (I.EQ.2) THEN
                 PD(I,J)=-0.25D0*(J+2)*(J+1)*J*(J-1)*X**(J+1)
              ENDIF
20         CONTINUE
           RETURN
        ENDIF
        LS=1
        IF (DABS(X).GT.1.0D0) LS=-1
        XQ=DSQRT(LS*(1.0D0-X*X))
C       Ensure connection to the complex-valued function for |x| > 1
        IF (X.LT.-1D0) XQ=-XQ
        XS=LS*(1.0D0-X*X)
        DO 30 I=1,M
30         PM(I,I)=-LS*(2.0D0*I-1.0D0)*XQ*PM(I-1,I-1)
        DO 35 I=0,MIN(M,N-1)
35         PM(I,I+1)=(2.0D0*I+1.0D0)*X*PM(I,I)
        DO 40 I=0,M
        DO 40 J=I+2,N
           PM(I,J)=((2.0D0*J-1.0D0)*X*PM(I,J-1)-
     &             (I+J-1.0D0)*PM(I,J-2))/(J-I)
40      CONTINUE
        PD(0,0)=0.0D0
        DO 45 J=1,N
45         PD(0,J)=LS*J*(PM(0,J-1)-X*PM(0,J))/XS
        DO 50 I=1,M
        DO 50 J=I,N
           PD(I,J)=LS*I*X*PM(I,J)/XS+(J+I)
     &             *(J-I+1.0D0)/XQ*PM(I-1,J)
50      CONTINUE
        RETURN
        END

C       **********************************

        SUBROUTINE MTU0(KF,M,Q,X,CSF,CSD)
C
C       ===============================================================
C       Purpose: Compute Mathieu functions cem(x,q) and sem(x,q)
C                and their derivatives ( q ≥ 0 )
C       Input :  KF  --- Function code
C                        KF=1 for computing cem(x,q) and cem'(x,q)
C                        KF=2 for computing sem(x,q) and sem'(x,q)
C                m   --- Order of Mathieu functions
C                q   --- Parameter of Mathieu functions
C                x   --- Argument of Mathieu functions (in degrees)
C       Output:  CSF --- cem(x,q) or sem(x,q)
C                CSD --- cem'x,q) or sem'x,q)
C       Routines called:
C            (1) CVA2 for computing the characteristic values
C            (2) FCOEF for computing the expansion coefficients
C       ===============================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION FG(251)
        EPS=1.0D-14
        IF (KF.EQ.1.AND.M.EQ.2*INT(M/2)) KD=1
        IF (KF.EQ.1.AND.M.NE.2*INT(M/2)) KD=2
        IF (KF.EQ.2.AND.M.NE.2*INT(M/2)) KD=3
        IF (KF.EQ.2.AND.M.EQ.2*INT(M/2)) KD=4
        CALL CVA2(KD,M,Q,A)
        IF (Q.LE.1.0D0) THEN
           QM=7.5+56.1*SQRT(Q)-134.7*Q+90.7*SQRT(Q)*Q
        ELSE
           QM=17.0+3.1*SQRT(Q)-.126*Q+.0037*SQRT(Q)*Q
        ENDIF
        KM=INT(QM+0.5*M)
        IF(KM.GT.251) THEN
           CSF=DNAN()
           CSD=DNAN()
           RETURN
        END IF
        CALL FCOEF(KD,M,Q,A,FG)
        IC=INT(M/2)+1
        RD=1.74532925199433D-2
        XR=X*RD
        CSF=0.0D0
        DO 10 K=1,KM
           IF (KD.EQ.1) THEN
              CSF=CSF+FG(K)*DCOS((2*K-2)*XR)
           ELSE IF (KD.EQ.2) THEN
              CSF=CSF+FG(K)*DCOS((2*K-1)*XR)
           ELSE IF (KD.EQ.3) THEN
              CSF=CSF+FG(K)*DSIN((2*K-1)*XR)
           ELSE IF (KD.EQ.4) THEN
              CSF=CSF+FG(K)*DSIN(2*K*XR)
           ENDIF
           IF (K.GE.IC.AND.DABS(FG(K)).LT.DABS(CSF)*EPS) GO TO 15
10         CONTINUE
15      CSD=0.0D0
        DO 20 K=1,KM
           IF (KD.EQ.1) THEN
              CSD=CSD-(2*K-2)*FG(K)*DSIN((2*K-2)*XR)
           ELSE IF (KD.EQ.2) THEN
              CSD=CSD-(2*K-1)*FG(K)*DSIN((2*K-1)*XR)
           ELSE IF (KD.EQ.3) THEN
              CSD=CSD+(2*K-1)*FG(K)*DCOS((2*K-1)*XR)
           ELSE IF (KD.EQ.4) THEN
              CSD=CSD+2.0D0*K*FG(K)*DCOS(2*K*XR)
           ENDIF
           IF (K.GE.IC.AND.DABS(FG(K)).LT.DABS(CSD)*EPS) GO TO 25
20         CONTINUE
25      RETURN
        END



C       **********************************

        SUBROUTINE CY01(KF,Z,ZF,ZD)
C
C       ===========================================================
C       Purpose: Compute complex Bessel functions Y0(z), Y1(z)
C                and their derivatives
C       Input :  z  --- Complex argument of Yn(z) ( n=0,1 )
C                KF --- Function choice code
C                    KF=0 for ZF=Y0(z) and ZD=Y0'(z)
C                    KF=1 for ZF=Y1(z) and ZD=Y1'(z)
C                    KF=2 for ZF=Y1'(z) and ZD=Y1''(z)
C       Output:  ZF --- Y0(z) or Y1(z) or Y1'(z)
C                ZD --- Y0'(z) or Y1'(z) or Y1''(z)
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,E,P,R,W)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION A(12),B(12),A1(12),B1(12)
        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        RP2=2.0D0/PI
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z2=Z*Z
        Z1=Z
        IF (A0.EQ.0.0D0) THEN
           CBJ0=(1.0D0,0.0D0)
           CBJ1=(0.0D0,0.0D0)
           CBY0=-(1.0D300,0.0D0)
           CBY1=-(1.0D300,0.0D0)
           CDY0=(1.0D300,0.0D0)
           CDY1=(1.0D300,0.0D0)
           GO TO 70
        ENDIF
        IF (DBLE(Z).LT.0.0) Z1=-Z
        IF (A0.LE.12.0) THEN
           CBJ0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 10 K=1,40
              CR=-0.25D0*CR*Z2/(K*K)
              CBJ0=CBJ0+CR
              IF (CDABS(CR).LT.CDABS(CBJ0)*1.0D-15) GO TO 15
10         CONTINUE
15         CBJ1=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 20 K=1,40
              CR=-0.25D0*CR*Z2/(K*(K+1.0D0))
              CBJ1=CBJ1+CR
              IF (CDABS(CR).LT.CDABS(CBJ1)*1.0D-15) GO TO 25
20         CONTINUE
25         CBJ1=0.5D0*Z1*CBJ1
           W0=0.0D0
           CR=(1.0D0,0.0D0)
           CS=(0.0D0,0.0D0)
           DO 30 K=1,40
              W0=W0+1.0D0/K
              CR=-0.25D0*CR/(K*K)*Z2
              CP=CR*W0
              CS=CS+CP
              IF (CDABS(CP).LT.CDABS(CS)*1.0D-15) GO TO 35
30         CONTINUE
35         CBY0=RP2*(CDLOG(Z1/2.0D0)+EL)*CBJ0-RP2*CS
           W1=0.0D0
           CR=(1.0D0,0.0D0)
           CS=(1.0D0,0.0D0)
           DO 40 K=1,40
              W1=W1+1.0D0/K
              CR=-0.25D0*CR/(K*(K+1))*Z2
              CP=CR*(2.0D0*W1+1.0D0/(K+1.0D0))
              CS=CS+CP
              IF (CDABS(CP).LT.CDABS(CS)*1.0D-15) GO TO 45
40         CONTINUE
45         CBY1=RP2*((CDLOG(Z1/2.0D0)+EL)*CBJ1-1.0D0/Z1-.25D0*Z1*CS)
        ELSE
           DATA A/-.703125D-01,.112152099609375D+00,
     &            -.5725014209747314D+00,.6074042001273483D+01,
     &            -.1100171402692467D+03,.3038090510922384D+04,
     &            -.1188384262567832D+06,.6252951493434797D+07,
     &            -.4259392165047669D+09,.3646840080706556D+11,
     &            -.3833534661393944D+13,.4854014686852901D+15/
           DATA B/ .732421875D-01,-.2271080017089844D+00,
     &             .1727727502584457D+01,-.2438052969955606D+02,
     &             .5513358961220206D+03,-.1825775547429318D+05,
     &             .8328593040162893D+06,-.5006958953198893D+08,
     &             .3836255180230433D+10,-.3649010818849833D+12,
     &             .4218971570284096D+14,-.5827244631566907D+16/
           DATA A1/.1171875D+00,-.144195556640625D+00,
     &             .6765925884246826D+00,-.6883914268109947D+01,
     &             .1215978918765359D+03,-.3302272294480852D+04,
     &             .1276412726461746D+06,-.6656367718817688D+07,
     &             .4502786003050393D+09,-.3833857520742790D+11,
     &             .4011838599133198D+13,-.5060568503314727D+15/
           DATA B1/-.1025390625D+00,.2775764465332031D+00,
     &             -.1993531733751297D+01,.2724882731126854D+02,
     &             -.6038440767050702D+03,.1971837591223663D+05,
     &             -.8902978767070678D+06,.5310411010968522D+08,
     &             -.4043620325107754D+10,.3827011346598605D+12,
     &             -.4406481417852278D+14,.6065091351222699D+16/
           K0=12
           IF (A0.GE.35.0) K0=10
           IF (A0.GE.50.0) K0=8
           CT1=Z1-.25D0*PI
           CP0=(1.0D0,0.0D0)
           DO 50 K=1,K0
50            CP0=CP0+A(K)*Z1**(-2*K)
           CQ0=-0.125D0/Z1
           DO 55 K=1,K0
55            CQ0=CQ0+B(K)*Z1**(-2*K-1)
           CU=CDSQRT(RP2/Z1)
           CBJ0=CU*(CP0*CDCOS(CT1)-CQ0*CDSIN(CT1))
           CBY0=CU*(CP0*CDSIN(CT1)+CQ0*CDCOS(CT1))
           CT2=Z1-.75D0*PI
           CP1=(1.0D0,0.0D0)
           DO 60 K=1,K0
60            CP1=CP1+A1(K)*Z1**(-2*K)
           CQ1=0.375D0/Z1
           DO 65 K=1,K0
65            CQ1=CQ1+B1(K)*Z1**(-2*K-1)
           CBJ1=CU*(CP1*CDCOS(CT2)-CQ1*CDSIN(CT2))
           CBY1=CU*(CP1*CDSIN(CT2)+CQ1*CDCOS(CT2))
        ENDIF
        IF (DBLE(Z).LT.0.0) THEN
           IF (DIMAG(Z).LT.0.0) CBY0=CBY0-2.0D0*CI*CBJ0
           IF (DIMAG(Z).GT.0.0) CBY0=CBY0+2.0D0*CI*CBJ0
           IF (DIMAG(Z).LT.0.0) CBY1=-(CBY1-2.0D0*CI*CBJ1)
           IF (DIMAG(Z).GT.0.0) CBY1=-(CBY1+2.0D0*CI*CBJ1)
           CBJ1=-CBJ1
        ENDIF
        CDY0=-CBY1
        CDY1=CBY0-1.0D0/Z*CBY1
70      IF (KF.EQ.0) THEN
           ZF=CBY0
           ZD=CDY0
        ELSE IF (KF.EQ.1) THEN
           ZF=CBY1
           ZD=CDY1
        ELSE IF (KF.EQ.2) THEN
           ZF=CDY1
           ZD=-CDY1/Z-(1.0D0-1.0D0/(Z*Z))*CBY1
        ENDIF
        RETURN
        END


C       **********************************

        SUBROUTINE FFK(KS,X,FR,FI,FM,FA,GR,GI,GM,GA)
C
C       =======================================================
C       Purpose: Compute modified Fresnel integrals F±(x)
C                and K±(x)
C       Input :  x   --- Argument of F±(x) and K±(x)
C                KS  --- Sign code
C                        KS=0 for calculating F+(x) and K+(x)
C                        KS=1 for calculating F_(x) and K_(x)
C       Output:  FR  --- Re[F±(x)]
C                FI  --- Im[F±(x)]
C                FM  --- |F±(x)|
C                FA  --- Arg[F±(x)]  (Degs.)
C                GR  --- Re[K±(x)]
C                GI  --- Im[K±(x)]
C                GM  --- |K±(x)|
C                GA  --- Arg[K±(x)]  (Degs.)
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        SRD= 57.29577951308233D0
        EPS=1.0D-15
        PI=3.141592653589793D0
        PP2=1.2533141373155D0
        P2P=.7978845608028654D0
        XA=DABS(X)
        X2=X*X
        X4=X2*X2
        IF (X.EQ.0.0D0) THEN
           FR=.5D0*DSQRT(0.5D0*PI)
           FI=(-1)**KS*FR
           FM=DSQRT(0.25D0*PI)
           FA=(-1)**KS*45.0D0
           GR=.5D0
           GI=0.0D0
           GM=.5D0
           GA=0.0D0
        ELSE
           IF (XA.LE.2.5D0) THEN
              XR=P2P*XA
              C1=XR
              DO 10 K=1,50
                 XR=-.5D0*XR*(4.0D0*K-3.0D0)/K/(2.0D0*K-1.0D0)
     &              /(4.0D0*K+1.0D0)*X4
                 C1=C1+XR
                 IF (DABS(XR/C1).LT.EPS) GO TO 15
10            CONTINUE
15            S1=P2P*XA*XA*XA/3.0D0
              XR=S1
              DO 20 K=1,50
                 XR=-.5D0*XR*(4.0D0*K-1.0D0)/K/(2.0D0*K+1.0D0)
     &              /(4.0D0*K+3.0D0)*X4
                 S1=S1+XR
                 IF (DABS(XR/S1).LT.EPS) GO TO 40
20            CONTINUE
           ELSE IF (XA.LT.5.5D0) THEN
              M=INT(42+1.75*X2)
              XSU=0.0D0
              XC=0.0D0
              XS=0.0D0
              XF1=0.0D0
              XF0=1D-100
              DO 25 K=M,0,-1
                 XF=(2.0D0*K+3.0D0)*XF0/X2-XF1
                 IF (K.EQ.2*INT(K/2))  THEN
                    XC=XC+XF
                 ELSE
                    XS=XS+XF
                 ENDIF
                 XSU=XSU+(2.0D0*K+1.0D0)*XF*XF
                 XF1=XF0
25               XF0=XF
              XQ=DSQRT(XSU)
              XW=P2P*XA/XQ
              C1=XC*XW
              S1=XS*XW
           ELSE
              XR=1.0D0
              XF=1.0D0
              DO 30 K=1,12
                 XR=-.25D0*XR*(4.0D0*K-1.0D0)*(4.0D0*K-3.0D0)/X4
30               XF=XF+XR
              XR=1.0D0/(2.0D0*XA*XA)
              XG=XR
              DO 35 K=1,12
                 XR=-.25D0*XR*(4.0D0*K+1.0D0)*(4.0D0*K-1.0D0)/X4
35               XG=XG+XR
              C1=.5D0+(XF*DSIN(X2)-XG*DCOS(X2))/DSQRT(2.0D0*PI)/XA
              S1=.5D0-(XF*DCOS(X2)+XG*DSIN(X2))/DSQRT(2.0D0*PI)/XA
           ENDIF
40         FR=PP2*(.5D0-C1)
           FI0=PP2*(.5D0-S1)
           FI=(-1)**KS*FI0
           FM=DSQRT(FR*FR+FI*FI)
           IF (FR.GE.0.0) THEN
              FA=SRD*DATAN(FI/FR)
           ELSE IF (FI.GT.0.0) THEN
              FA=SRD*(DATAN(FI/FR)+PI)
           ELSE IF (FI.LT.0.0) THEN
              FA=SRD*(DATAN(FI/FR)-PI)
           ENDIF
           XP=X*X+PI/4.0D0
           CS=DCOS(XP)
           SS=DSIN(XP)
           XQ2=1.0D0/DSQRT(PI)
           GR=XQ2*(FR*CS+FI0*SS)
           GI=(-1)**KS*XQ2*(FI0*CS-FR*SS)
           GM=DSQRT(GR*GR+GI*GI)
           IF (GR.GE.0.0) THEN
              GA=SRD*DATAN(GI/GR)
           ELSE IF (GI.GT.0.0) THEN
              GA=SRD*(DATAN(GI/GR)+PI)
           ELSE IF (GI.LT.0.0) THEN
              GA=SRD*(DATAN(GI/GR)-PI)
           ENDIF
           IF (X.LT.0.0D0) THEN
              FR=PP2-FR
              FI=(-1)**KS*PP2-FI
              FM=DSQRT(FR*FR+FI*FI)
              FA=SRD*DATAN(FI/FR)
              GR=DCOS(X*X)-GR
              GI=-(-1)**KS*DSIN(X*X)-GI
              GM=DSQRT(GR*GR+GI*GI)
              GA=SRD*DATAN(GI/GR)
           ENDIF
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE AIRYA(X,AI,BI,AD,BD)
C
C       ======================================================
C       Purpose: Compute Airy functions and their derivatives
C       Input:   x  --- Argument of Airy function
C       Output:  AI --- Ai(x)
C                BI --- Bi(x)
C                AD --- Ai'(x)
C                BD --- Bi'(x)
C       Routine called:
C                AJYIK for computing Jv(x), Yv(x), Iv(x) and
C                Kv(x) with v=1/3 and 2/3
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PIR=0.318309886183891D0
        C1=0.355028053887817D0
        C2=0.258819403792807D0
        SR3=1.732050807568877D0
        Z=XA**1.5/1.5D0
        XQ=DSQRT(XA)
        CALL AJYIK(Z,VJ1,VJ2,VY1,VY2,VI1,VI2,VK1,VK2)
        IF (X.EQ.0.0D0) THEN
           AI=C1
           BI=SR3*C1
           AD=-C2
           BD=SR3*C2
        ELSE IF (X.GT.0.0D0) THEN
           AI=PIR*XQ/SR3*VK1
           BI=XQ*(PIR*VK1+2.0D0/SR3*VI1)
           AD=-XA/SR3*PIR*VK2
           BD=XA*(PIR*VK2+2.0D0/SR3*VI2)
        ELSE
           AI=0.5D0*XQ*(VJ1-VY1/SR3)
           BI=-0.5D0*XQ*(VJ1/SR3+VY1)
           AD=0.5D0*XA*(VJ2+VY2/SR3)
           BD=0.5D0*XA*(VJ2/SR3-VY2)
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE AIRYB(X,AI,BI,AD,BD)
C
C       =======================================================
C       Purpose: Compute Airy functions and their derivatives
C       Input:   x  --- Argument of Airy function
C       Output:  AI --- Ai(x)
C                BI --- Bi(x)
C                AD --- Ai'(x)
C                BD --- Bi'(x)
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION CK(41),DK(41)
        EPS=1.0D-15
        PI=3.141592653589793D0
        C1=0.355028053887817D0
        C2=0.258819403792807D0
        SR3=1.732050807568877D0
        XA=DABS(X)
        XQ=DSQRT(XA)
        XM=8.0D0
        IF (X.GT.0.0D0) XM=5.0D0
        IF (X.EQ.0.0D0) THEN
           AI=C1
           BI=SR3*C1
           AD=-C2
           BD=SR3*C2
           RETURN
        ENDIF
        IF (XA.LE.XM) THEN
           FX=1.0D0
           R=1.0D0
           DO 10 K=1,40
              R=R*X/(3.0D0*K)*X/(3.0D0*K-1.0D0)*X
              FX=FX+R
              IF (DABS(R).LT.DABS(FX)*EPS) GO TO 15
10         CONTINUE
15         GX=X
           R=X
           DO 20 K=1,40
              R=R*X/(3.0D0*K)*X/(3.0D0*K+1.0D0)*X
              GX=GX+R
              IF (DABS(R).LT.DABS(GX)*EPS) GO TO 25
20         CONTINUE
25         AI=C1*FX-C2*GX
           BI=SR3*(C1*FX+C2*GX)
           DF=0.5D0*X*X
           R=DF
           DO 30 K=1,40
              R=R*X/(3.0D0*K)*X/(3.0D0*K+2.0D0)*X
              DF=DF+R
              IF (DABS(R).LT.DABS(DF)*EPS) GO TO 35
30         CONTINUE
35         DG=1.0D0
           R=1.0D0
           DO 40 K=1,40
              R=R*X/(3.0D0*K)*X/(3.0D0*K-2.0D0)*X
              DG=DG+R
              IF (DABS(R).LT.DABS(DG)*EPS) GO TO 45
40         CONTINUE
45         AD=C1*DF-C2*DG
           BD=SR3*(C1*DF+C2*DG)
        ELSE
           XE=XA*XQ/1.5D0
           XR1=1.0D0/XE
           XAR=1.0D0/XQ
           XF=DSQRT(XAR)
           RP=0.5641895835477563D0
           R=1.0D0
           DO 50 K=1,40
              R=R*(6.0D0*K-1.0D0)/216.0D0*(6.0D0*K-3.0D0)
     &          /K*(6.0D0*K-5.0D0)/(2.0D0*K-1.0D0)
              CK(K)=R
50            DK(K)=-(6.0D0*K+1.0D0)/(6.0D0*K-1.0D0)*CK(K)
           KM=INT(24.5-XA)
           IF (XA.LT.6.0) KM=14
           IF (XA.GT.15.0) KM=10
           IF (X.GT.0.0D0) THEN
              SAI=1.0D0
              SAD=1.0D0
              R=1.0D0
              DO 55 K=1,KM
                 R=-R*XR1
                 SAI=SAI+CK(K)*R
55               SAD=SAD+DK(K)*R
              SBI=1.0D0
              SBD=1.0D0
              R=1.0D0
              DO 60 K=1,KM
                 R=R*XR1
                 SBI=SBI+CK(K)*R
60               SBD=SBD+DK(K)*R
              XP1=DEXP(-XE)
              AI=0.5D0*RP*XF*XP1*SAI
              BI=RP*XF/XP1*SBI
              AD=-.5D0*RP/XF*XP1*SAD
              BD=RP/XF/XP1*SBD
           ELSE
              XCS=DCOS(XE+PI/4.0D0)
              XSS=DSIN(XE+PI/4.0D0)
              SSA=1.0D0
              SDA=1.0D0
              R=1.0D0
              XR2=1.0D0/(XE*XE)
              DO 65 K=1,KM
                 R=-R*XR2
                 SSA=SSA+CK(2*K)*R
65               SDA=SDA+DK(2*K)*R
              SSB=CK(1)*XR1
              SDB=DK(1)*XR1
              R=XR1
              DO 70 K=1,KM
                 R=-R*XR2
                 SSB=SSB+CK(2*K+1)*R
70               SDB=SDB+DK(2*K+1)*R
              AI=RP*XF*(XSS*SSA-XCS*SSB)
              BI=RP*XF*(XCS*SSA+XSS*SSB)
              AD=-RP/XF*(XCS*SDA+XSS*SDB)
              BD=RP/XF*(XSS*SDA-XCS*SDB)
           ENDIF
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE SCKA(M,N,C,CV,KD,CK)
C
C       ======================================================
C       Purpose: Compute the expansion coefficients of the
C                prolate and oblate spheroidal functions, c2k
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                cv --- Characteristic value
C                KD --- Function code
C                       KD=1 for prolate; KD=-1 for oblate
C       Output:  CK(k) --- Expansion coefficients ck;
C                          CK(1), CK(2),... correspond to
C                          c0, c2,...
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION CK(200)
        IF (C.LE.1.0D-10) C=1.0D-10
        NM=25+INT((N-M)/2+C)
        CS=C*C*KD
        IP=1
        IF (N-M.EQ.2*INT((N-M)/2)) IP=0
        FS=1.0D0
        F1=0.0D0
        F0=1.0D-100
        KB=0
        CK(NM+1)=0.0D0
        FL=0.0D0
        DO 15 K=NM,1,-1
           F=(((2.0D0*K+M+IP)*(2.0D0*K+M+1.0D0+IP)-CV+CS)*F0
     &       -4.0D0*(K+1.0D0)*(K+M+1.0D0)*F1)/CS
           IF (DABS(F).GT.DABS(CK(K+1))) THEN
              CK(K)=F
              F1=F0
              F0=F
              IF (DABS(F).GT.1.0D+100) THEN
                 DO 5 K1=NM,K,-1
5                   CK(K1)=CK(K1)*1.0D-100
                 F1=F1*1.0D-100
                 F0=F0*1.0D-100
              ENDIF
           ELSE
              KB=K
              FL=CK(K+1)
              F1=1.0D0
              F2=0.25D0*((M+IP)*(M+IP+1.0)-CV+CS)/(M+1.0)*F1
              CK(1)=F1
              IF (KB.EQ.1) THEN
                 FS=F2
              ELSE IF (KB.EQ.2) THEN
                 CK(2)=F2
                 FS=0.125D0*(((M+IP+2.0)*(M+IP+3.0)-CV+CS)*F2
     &              -CS*F1)/(M+2.0)
              ELSE
                 CK(2)=F2
                 DO 10 J=3,KB+1
                    F=0.25D0*(((2.0*J+M+IP-4.0)*(2.0*J+M+IP-
     &                3.0)-CV+CS)*F2-CS*F1)/((J-1.0)*(J+M-1.0))
                    IF (J.LE.KB) CK(J)=F
                    F1=F2
10                  F2=F
                 FS=F
              ENDIF
              GO TO 20
           ENDIF
15      CONTINUE
20      SU1=0.0D0
        DO 25 K=1,KB
25         SU1=SU1+CK(K)
        SU2=0.0D0
        DO 30 K=KB+1,NM
30         SU2=SU2+CK(K)
        R1=1.0D0
        DO 35 J=1,(N+M+IP)/2
35         R1=R1*(J+0.5D0*(N+M+IP))
        R2=1.0D0
        DO 40 J=1,(N-M-IP)/2
40         R2=-R2*J
        IF (KB.EQ.0) THEN
            S0=R1/(2.0D0**N*R2*SU2)
        ELSE
            S0=R1/(2.0D0**N*R2*(FL/FS*SU1+SU2))
        ENDIF
        DO 45 K=1,KB
45         CK(K)=FL/FS*S0*CK(K)
        DO 50 K=KB+1,NM
50         CK(K)=S0*CK(K)
        RETURN
        END



C       **********************************

        SUBROUTINE SCKB(M,N,C,DF,CK)
C
C       ======================================================
C       Purpose: Compute the expansion coefficients of the
C                prolate and oblate spheroidal functions
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                DF(k) --- Expansion coefficients dk
C       Output:  CK(k) --- Expansion coefficients ck;
C                          CK(1), CK(2), ... correspond to
C                          c0, c2, ...
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION DF(200),CK(200)
        IF (C.LE.1.0D-10) C=1.0D-10
        NM=25+INT(0.5*(N-M)+C)
        IP=1
        IF (N-M.EQ.2*INT((N-M)/2)) IP=0
        REG=1.0D0
        IF (M+NM.GT.80) REG=1.0D-200
        FAC=-0.5D0**M
        SW=0.0D0
        DO 35 K=0,NM-1
           FAC=-FAC
           I1=2*K+IP+1
           R=REG
           DO 10 I=I1,I1+2*M-1
10            R=R*I
           I2=K+M+IP
           DO 15 I=I2,I2+K-1
15            R=R*(I+0.5D0)
           SUM=R*DF(K+1)
           DO 20 I=K+1,NM
              D1=2.0D0*I+IP
              D2=2.0D0*M+D1
              D3=I+M+IP-0.5D0
              R=R*D2*(D2-1.0D0)*I*(D3+K)/(D1*(D1-1.0D0)*(I-K)*D3)
              SUM=SUM+R*DF(I+1)
              IF (DABS(SW-SUM).LT.DABS(SUM)*1.0D-14) GOTO 25
20            SW=SUM
25         R1=REG
           DO 30 I=2,M+K
30            R1=R1*I
35         CK(K+1)=FAC*SUM/R1
        RETURN
        END



C       **********************************

        SUBROUTINE CPDLA(N,Z,CDN)
C
C       ===========================================================
C       Purpose: Compute complex parabolic cylinder function Dn(z)
C                for large argument
C       Input:   z   --- Complex argument of Dn(z)
C                n   --- Order of Dn(z) (n = 0,±1,±2,…)
C       Output:  CDN --- Dn(z)
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        CB0=Z**N*CDEXP(-.25D0*Z*Z)
        CR=(1.0D0,0.0D0)
        CDN=(1.0D0,0.0D0)
        DO 10 K=1,16
           CR=-0.5D0*CR*(2.0*K-N-1.0)*(2.0*K-N-2.0)/(K*Z*Z)
           CDN=CDN+CR
           IF (CDABS(CR).LT.CDABS(CDN)*1.0D-12) GO TO 15
10      CONTINUE
15      CDN=CB0*CDN
        RETURN
        END



C       **********************************

        SUBROUTINE FCSZO(KF,NT,ZO)
C
C       ===============================================================
C       Purpose: Compute the complex zeros of Fresnel integral C(z)
C                or S(z) using modified Newton's iteration method
C       Input :  KF  --- Function code
C                        KF=1 for C(z) or KF=2 for S(z)
C                NT  --- Total number of zeros
C       Output:  ZO(L) --- L-th zero of C(z) or S(z)
C       Routines called:
C            (1) CFC for computing Fresnel integral C(z)
C            (2) CFS for computing Fresnel integral S(z)
C       ==============================================================
C
        IMPLICIT DOUBLE PRECISION (E,P,W)
        IMPLICIT COMPLEX *16 (C,Z)
        DIMENSION ZO(NT)
        PI=3.141592653589793D0
        PSQ=0.0D0
        W=0.0D0
        DO 35 NR=1,NT
           IF (KF.EQ.1) PSQ=DSQRT(4.0D0*NR-1.0D0)
           IF (KF.EQ.2) PSQ=2.0D0*NR**(0.5)
           PX=PSQ-DLOG(PI*PSQ)/(PI*PI*PSQ**3.0)
           PY=DLOG(PI*PSQ)/(PI*PSQ)
           Z = DCMPLX(PX, PY)
           IF (KF.EQ.2) THEN
              IF (NR.EQ.2) Z=(2.8334,0.2443)
              IF (NR.EQ.3) Z=(3.4674,0.2185)
              IF (NR.EQ.4) Z=(4.0025,0.2008)
           ENDIF
           IT=0
15         IT=IT+1
           IF (KF.EQ.1) CALL CFC(Z,ZF,ZD)
           IF (KF.EQ.2) CALL CFS(Z,ZF,ZD)
           ZP=(1.0D0,0.0D0)
           DO 20 I=1,NR-1
20            ZP=ZP*(Z-ZO(I))
           ZFD=ZF/ZP
           ZQ=(0.0D0,0.0D0)
           DO 30 I=1,NR-1
              ZW=(1.0D0,0.0D0)
              DO 25 J=1,NR-1
                 IF (J.EQ.I) GO TO 25
                 ZW=ZW*(Z-ZO(J))
25            CONTINUE
30            ZQ=ZQ+ZW
           ZGD=(ZD-ZQ*ZFD)/ZP
           Z=Z-ZFD/ZGD
           W0=W
           W=CDABS(Z)
           IF (IT.LE.50.AND.DABS((W-W0)/W).GT.1.0D-12) GO TO 15
35         ZO(NR)=Z
        RETURN
        END



C       **********************************

        SUBROUTINE E1XA(X,E1)
C
C       ============================================
C       Purpose: Compute exponential integral E1(x)
C       Input :  x  --- Argument of E1(x)
C       Output:  E1 --- E1(x) ( x > 0 )
C       ============================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (X.EQ.0.0) THEN
           E1=1.0D+300
        ELSE IF (X.LE.1.0) THEN
           E1=-DLOG(X)+((((1.07857D-3*X-9.76004D-3)*X+5.519968D-2)*X
     &        -0.24991055D0)*X+0.99999193D0)*X-0.57721566D0
        ELSE
           ES1=(((X+8.5733287401D0)*X+18.059016973D0)*X
     &         +8.6347608925D0)*X+0.2677737343D0
           ES2=(((X+9.5733223454D0)*X+25.6329561486D0)*X
     &         +21.0996530827D0)*X+3.9584969228D0
           E1=DEXP(-X)/X*ES1/ES2
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE LPMV0(V,M,X,PMV)
C
C       =======================================================
C       Purpose: Compute the associated Legendre function
C                Pmv(x) with an integer order and an arbitrary
C                nonnegative degree v
C       Input :  x   --- Argument of Pm(x)  ( -1 ≤ x ≤ 1 )
C                m   --- Order of Pmv(x)
C                v   --- Degree of Pmv(x)
C       Output:  PMV --- Pmv(x)
C       Routine called:  PSI_SPEC for computing Psi function
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        EPS=1.0D-14
        NV=INT(V)
        V0=V-NV
        IF (X.EQ.-1.0D0.AND.V.NE.NV) THEN
           IF (M.EQ.0) PMV=-1.0D+300
           IF (M.NE.0) PMV=1.0D+300
           RETURN
        ENDIF
        C0=1.0D0
        IF (M.NE.0) THEN
           RG=V*(V+M)
           DO 10 J=1,M-1
10            RG=RG*(V*V-J*J)
           XQ=DSQRT(1.0D0-X*X)
           R0=1.0D0
           DO 15 J=1,M
15            R0=.5D0*R0*XQ/J
           C0=R0*RG
        ENDIF
        IF (V0.EQ.0.0D0) THEN
           PMV=1.0D0
           R=1.0D0
           DO 20 K=1,NV-M
              R=0.5D0*R*(-NV+M+K-1.0D0)*(NV+M+K)/(K*(K+M))
     &          *(1.0D0+X)
20            PMV=PMV+R
           PMV=(-1)**NV*C0*PMV
        ELSE
           IF (X.GE.-0.35D0) THEN
              PMV=1.0D0
              R=1.0D0
              DO 25 K=1,100
                 R=0.5D0*R*(-V+M+K-1.0D0)*(V+M+K)/(K*(M+K))*(1.0D0-X)
                 PMV=PMV+R
                 IF (K.GT.12.AND.DABS(R/PMV).LT.EPS) GO TO 30
25            CONTINUE
30            PMV=(-1)**M*C0*PMV
           ELSE
              VS=DSIN(V*PI)/PI
              PV0=0.0D0
              IF (M.NE.0) THEN
                 QR=DSQRT((1.0D0-X)/(1.0D0+X))
                 R2=1.0D0
                 DO 35 J=1,M
35                  R2=R2*QR*J
                 S0=1.0D0
                 R1=1.0D0
                 DO 40 K=1,M-1
                    R1=0.5D0*R1*(-V+K-1)*(V+K)/(K*(K-M))*(1.0D0+X)
40                  S0=S0+R1
                 PV0=-VS*R2/M*S0
              ENDIF
              CALL PSI_SPEC(V,PSV)
              PA=2.0D0*(PSV+EL)+PI/DTAN(PI*V)+1.0D0/V
              S1=0.0D0
              DO 45 J=1,M
45               S1=S1+(J*J+V*V)/(J*(J*J-V*V))
              PMV=PA+S1-1.0D0/(M-V)+DLOG(0.5D0*(1.0D0+X))
              R=1.0D0
              DO 60 K=1,100
                 R=0.5D0*R*(-V+M+K-1.0D0)*(V+M+K)/(K*(K+M))*(1.0D0+X)
                 S=0.0D0
                 DO 50 J=1,M
50                  S=S+((K+J)**2+V*V)/((K+J)*((K+J)**2-V*V))
                 S2=0.0D0
                 DO 55 J=1,K
55                  S2=S2+1.0D0/(J*(J*J-V*V))
                 PSS=PA+S+2.0D0*V*V*S2-1.0D0/(M+K-V)
     &               +DLOG(0.5D0*(1.0D0+X))
                 R2=PSS*R
                 PMV=PMV+R2
                 IF (DABS(R2/PMV).LT.EPS) GO TO 65
60            CONTINUE
65            PMV=PV0+PMV*VS*C0
           ENDIF
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE LPMV(V,M,X,PMV)
C
C       =======================================================
C       Purpose: Compute the associated Legendre function
C                Pmv(x) with an integer order and an arbitrary
C                degree v, using recursion for large degrees
C       Input :  x   --- Argument of Pm(x)  ( -1 ≤ x ≤ 1 )
C                m   --- Order of Pmv(x)
C                v   --- Degree of Pmv(x)
C       Output:  PMV --- Pmv(x)
C       Routine called:  LPMV0
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (X.EQ.-1.0D0.AND.V.NE.INT(V)) THEN
           IF (M.EQ.0) PMV=-DINF()
           IF (M.NE.0) PMV=DINF()
           RETURN
        ENDIF
        VX=V
        MX=M
        IF (V.LT.0) THEN
           VX=-VX-1
        ENDIF
        NEG_M=0
        IF (M.LT.0) THEN
           IF ((VX+M+1).GT.0D0.OR.VX.NE.INT(VX)) THEN
              NEG_M=1
              MX=-M
           ELSE
C             We don't handle cases where DLMF 14.9.3 doesn't help
              PMV=DNAN()
              RETURN
           END IF
        ENDIF
        NV=INT(VX)
        V0=VX-NV
        IF (NV.GT.2.AND.NV.GT.MX) THEN
C          Up-recursion on degree, AMS 8.5.3
           CALL LPMV0(V0+MX, MX, X, P0)
           CALL LPMV0(V0+MX+1, MX, X, P1)
           PMV = P1
           DO 10 J=MX+2,NV
              PMV = ((2*(V0+J)-1)*X*P1 - (V0+J-1+MX)*P0) / (V0+J-MX)
              P0 = P1
              P1 = PMV
10         CONTINUE
        ELSE
           CALL LPMV0(VX, MX, X, PMV)
        ENDIF
        IF (NEG_M.NE.0.AND.ABS(PMV).LT.1.0D+300) THEN
C          DLMF 14.9.3
           CALL GAMMA2(VX-MX+1, G1)
           CALL GAMMA2(VX+MX+1, G2)
           PMV = PMV*G1/G2 * (-1)**MX
        ENDIF
        END


C       **********************************

        SUBROUTINE CGAMA(X,Y,KF,GR,GI)
C
C       =========================================================
C       Purpose: Compute the gamma function Г(z) or ln[Г(z)]
C                for a complex argument
C       Input :  x  --- Real part of z
C                y  --- Imaginary part of z
C                KF --- Function code
C                       KF=0 for ln[Г(z)]
C                       KF=1 for Г(z)
C       Output:  GR --- Real part of ln[Г(z)] or Г(z)
C                GI --- Imaginary part of ln[Г(z)] or Г(z)
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(10)
        PI=3.141592653589793D0
        DATA A/8.333333333333333D-02,-2.777777777777778D-03,
     &         7.936507936507937D-04,-5.952380952380952D-04,
     &         8.417508417508418D-04,-1.917526917526918D-03,
     &         6.410256410256410D-03,-2.955065359477124D-02,
     &         1.796443723688307D-01,-1.39243221690590D+00/
        IF (Y.EQ.0.0D0.AND.X.EQ.INT(X).AND.X.LE.0.0D0) THEN
           GR=1.0D+300
           GI=0.0D0
           RETURN
        ELSE IF (X.LT.0.0D0) THEN
           X1=X
           Y1=Y
           X=-X
           Y=-Y
        ELSE
           Y1=0.0D0
           X1=X
        ENDIF
        X0=X
        NA=0
        IF (X.LE.7.0) THEN
           NA=INT(7-X)
           X0=X+NA
        ENDIF
        Z1=DSQRT(X0*X0+Y*Y)
        TH=DATAN(Y/X0)
        GR=(X0-.5D0)*DLOG(Z1)-TH*Y-X0+0.5D0*DLOG(2.0D0*PI)
        GI=TH*(X0-0.5D0)+Y*DLOG(Z1)-Y
        DO 10 K=1,10
           T=Z1**(1-2*K)
           GR=GR+A(K)*T*DCOS((2.0D0*K-1.0D0)*TH)
10         GI=GI-A(K)*T*DSIN((2.0D0*K-1.0D0)*TH)
        IF (X.LE.7.0) THEN
           GR1=0.0D0
           GI1=0.0D0
           DO 15 J=0,NA-1
              GR1=GR1+.5D0*DLOG((X+J)**2+Y*Y)
15            GI1=GI1+DATAN(Y/(X+J))
           GR=GR-GR1
           GI=GI-GI1
        ENDIF
        IF (X1.LT.0.0D0) THEN
           Z1=DSQRT(X*X+Y*Y)
           TH1=DATAN(Y/X)
           SR=-DSIN(PI*X)*DCOSH(PI*Y)
           SI=-DCOS(PI*X)*DSINH(PI*Y)
           Z2=DSQRT(SR*SR+SI*SI)
           TH2=DATAN(SI/SR)
           IF (SR.LT.0.0D0) TH2=PI+TH2
           GR=DLOG(PI/(Z1*Z2))-GR
           GI=-TH1-TH2-GI
           X=X1
           Y=Y1
        ENDIF
        IF (KF.EQ.1) THEN
           G0=DEXP(GR)
           GR=G0*DCOS(GI)
           GI=G0*DSIN(GI)
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE ASWFB(M,N,C,X,KD,CV,S1F,S1D)
C
C       ===========================================================
C       Purpose: Compute the prolate and oblate spheroidal angular
C                functions of the first kind and their derivatives
C       Input :  m  --- Mode parameter,  m = 0,1,2,...
C                n  --- Mode parameter,  n = m,m+1,...
C                c  --- Spheroidal parameter
C                x  --- Argument of angular function, |x| < 1.0
C                KD --- Function code
C                       KD=1 for prolate;  KD=-1 for oblate
C                cv --- Characteristic value
C       Output:  S1F --- Angular function of the first kind
C                S1D --- Derivative of the angular function of
C                        the first kind
C       Routines called:
C            (1) SDMN for computing expansion coefficients dk
C            (2) LPMNS for computing associated Legendre function
C                of the first kind Pmn(x)
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION DF(200),PM(0:251),PD(0:251)
        EPS=1.0D-14
        IP=1
        IF (N-M.EQ.2*INT((N-M)/2)) IP=0
        NM=25+INT((N-M)/2+C)
        NM2=2*NM+M
        CALL SDMN(M,N,C,CV,KD,DF)
        CALL LPMNS(M,NM2,X,PM,PD)
        SW=0.0D0
        SU1=0.0D0
        DO 10 K=1,NM
           MK=M+2*(K-1)+IP
           SU1=SU1+DF(K)*PM(MK)
           IF (DABS(SW-SU1).LT.DABS(SU1)*EPS) GOTO 15
10         SW=SU1
15      S1F=(-1)**M*SU1
        SU1=0.0D0
        DO 20 K=1,NM
           MK=M+2*(K-1)+IP
           SU1=SU1+DF(K)*PD(MK)
           IF (DABS(SW-SU1).LT.DABS(SU1)*EPS) GOTO 25
20         SW=SU1
25      S1D=(-1)**M*SU1
        RETURN
        END



C       **********************************

        SUBROUTINE CHGUS(A,B,X,HU,ID)
C
C       ======================================================
C       Purpose: Compute confluent hypergeometric function
C                U(a,b,x) for small argument x
C       Input  : a  --- Parameter
C                b  --- Parameter ( b <> 0,-1,-2,...)
C                x  --- Argument
C       Output:  HU --- U(a,b,x)
C                ID --- Estimated number of significant digits
C       Routine called: GAMMA2 for computing gamma function
C       ======================================================
C
C       DLMF 13.2.42 with prefactors rewritten according to
C       DLMF 5.5.3, M(a, b, x) with DLMF 13.2.2
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        ID=-100
        PI=3.141592653589793D0
        CALL GAMMA2(A,GA)
        CALL GAMMA2(B,GB)
        XG1=1.0D0+A-B
        CALL GAMMA2(XG1,GAB)
        XG2=2.0D0-B
        CALL GAMMA2(XG2,GB2)
        HU0=PI/DSIN(PI*B)
        R1=HU0/(GAB*GB)
        R2=HU0*X**(1.0D0-B)/(GA*GB2)
        HU=R1-R2
        HMAX=0.0D0
        HMIN=1.0D+300
        H0=0.0D0
        DO 10 J=1,150
           R1=R1*(A+J-1.0D0)/(J*(B+J-1.0D0))*X
           R2=R2*(A-B+J)/(J*(1.0D0-B+J))*X
           HU=HU+R1-R2
           HUA=DABS(HU)
           IF (HUA.GT.HMAX) HMAX=HUA
           IF (HUA.LT.HMIN) HMIN=HUA
           IF (DABS(HU-H0).LT.DABS(HU)*1.0D-15) GO TO 15
10         H0=HU
15      D1=LOG10(HMAX)
        D2=0.0D0
        IF (HMIN.NE.0.0) D2=LOG10(HMIN)
        ID=15-ABS(D1-D2)
        RETURN
        END



C       **********************************

        SUBROUTINE ITTH0(X,TTH)
C
C       ===========================================================
C       Purpose: Evaluate the integral H0(t)/t with respect to t
C                from x to infinity
C       Input :  x   --- Lower limit  ( x ≥ 0 )
C       Output:  TTH --- Integration of H0(t)/t from x to infinity
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        S=1.0D0
        R=1.0D0
        IF (X.LT.24.5D0) THEN
           DO 10 K=1,60
              R=-R*X*X*(2.0*K-1.0D0)/(2.0*K+1.0D0)**3
              S=S+R
              IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 15
10         CONTINUE
15         TTH=PI/2.0D0-2.0D0/PI*X*S
        ELSE
           DO 20 K=1,10
              R=-R*(2.0*K-1.0D0)**3/((2.0*K+1.0D0)*X*X)
              S=S+R
              IF (DABS(R).LT.DABS(S)*1.0D-12) GO TO 25
20            CONTINUE
25         TTH=2.0D0/(PI*X)*S
           T=8.0D0/X
           XT=X+.25D0*PI
           F0=(((((.18118D-2*T-.91909D-2)*T+.017033D0)*T
     &        -.9394D-3)*T-.051445D0)*T-.11D-5)*T+.7978846D0
           G0=(((((-.23731D-2*T+.59842D-2)*T+.24437D-2)*T
     &        -.0233178D0)*T+.595D-4)*T+.1620695D0)*T
           TTY=(F0*DSIN(XT)-G0*DCOS(XT))/(DSQRT(X)*X)
           TTH=TTH+TTY
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE LGAMA(KF,X,GL)
C
C       ==================================================
C       Purpose: Compute gamma function Г(x) or ln[Г(x)]
C       Input:   x  --- Argument of Г(x) ( x > 0 )
C                KF --- Function code
C                       KF=1 for Г(x); KF=0 for ln[Г(x)]
C       Output:  GL --- Г(x) or ln[Г(x)]
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(10)
        DATA A/8.333333333333333D-02,-2.777777777777778D-03,
     &         7.936507936507937D-04,-5.952380952380952D-04,
     &         8.417508417508418D-04,-1.917526917526918D-03,
     &         6.410256410256410D-03,-2.955065359477124D-02,
     &         1.796443723688307D-01,-1.39243221690590D+00/
        X0=X
        N=0
        IF (X.EQ.1.0.OR.X.EQ.2.0) THEN
           GL=0.0D0
           GO TO 20
        ELSE IF (X.LE.7.0) THEN
           N=INT(7-X)
           X0=X+N
        ENDIF
        X2=1.0D0/(X0*X0)
        XP=6.283185307179586477D0
        GL0=A(10)
        DO 10 K=9,1,-1
10         GL0=GL0*X2+A(K)
        GL=GL0/X0+0.5D0*DLOG(XP)+(X0-.5D0)*DLOG(X0)-X0
        IF (X.LE.7.0) THEN
           DO 15 K=1,N
              GL=GL-DLOG(X0-1.0D0)
15            X0=X0-1.0D0
        ENDIF
20      IF (KF.EQ.1) GL=DEXP(GL)
        RETURN
        END

C       **********************************

        SUBROUTINE LQNA(N,X,QN,QD)
C
C       =====================================================
C       Purpose: Compute Legendre functions Qn(x) and Qn'(x)
C       Input :  x  --- Argument of Qn(x) ( -1 ≤ x ≤ 1 )
C                n  --- Degree of Qn(x) ( n = 0,1,2,… )
C       Output:  QN(n) --- Qn(x)
C                QD(n) --- Qn'(x)
C                ( 1.0D+300 stands for infinity )
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (Q,X)
        DIMENSION QN(0:N),QD(0:N)
        IF (DABS(X).EQ.1.0D0) THEN
           DO 10 K=0,N
              QN(K)=1.0D+300
              QD(K)=-1.0D+300
10         CONTINUE
        ELSE IF (DABS(X).LT.1.0D0) THEN
           Q0=0.5D0*DLOG((1.0D0+X)/(1.0D0-X))
           Q1=X*Q0-1.0D0
           QN(0)=Q0
           QN(1)=Q1
           QD(0)=1.0D0/(1.0D0-X*X)
           QD(1)=QN(0)+X*QD(0)
           DO 15 K=2,N
              QF=((2*K-1)*X*Q1-(K-1)*Q0)/K
              QN(K)=QF
              QD(K)=(QN(K-1)-X*QF)*K/(1.0D0-X*X)
              Q0=Q1
15            Q1=QF
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE DVLA(VA,X,PD)
C
C       ====================================================
C       Purpose: Compute parabolic cylinder functions Dv(x)
C                for large argument
C       Input:   x  --- Argument
C                va --- Order
C       Output:  PD --- Dv(x)
C       Routines called:
C             (1) VVLA for computing Vv(x) for large |x|
C             (2) GAMMA2 for computing Г(x)
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EPS=1.0D-12
        EP=DEXP(-.25*X*X)
        A0=DABS(X)**VA*EP
        R=1.0D0
        PD=1.0D0
        DO 10 K=1,16
           R=-0.5D0*R*(2.0*K-VA-1.0)*(2.0*K-VA-2.0)/(K*X*X)
           PD=PD+R
           IF (DABS(R/PD).LT.EPS) GO TO 15
10      CONTINUE
15      PD=A0*PD
        IF (X.LT.0.0D0) THEN
            X1=-X
            CALL VVLA(VA,X1,VL)
            CALL GAMMA2(-VA,GL)
            PD=PI*VL/GL+DCOS(PI*VA)*PD
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE IK01A(X,BI0,DI0,BI1,DI1,BK0,DK0,BK1,DK1)
C
C       =========================================================
C       Purpose: Compute modified Bessel functions I0(x), I1(1),
C                K0(x) and K1(x), and their derivatives
C       Input :  x   --- Argument ( x ≥ 0 )
C       Output:  BI0 --- I0(x)
C                DI0 --- I0'(x)
C                BI1 --- I1(x)
C                DI1 --- I1'(x)
C                BK0 --- K0(x)
C                DK0 --- K0'(x)
C                BK1 --- K1(x)
C                DK1 --- K1'(x)
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(12),B(12),A1(8)
        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        X2=X*X
        IF (X.EQ.0.0D0) THEN
           BI0=1.0D0
           BI1=0.0D0
           BK0=1.0D+300
           BK1=1.0D+300
           DI0=0.0D0
           DI1=0.5D0
           DK0=-1.0D+300
           DK1=-1.0D+300
           RETURN
        ELSE IF (X.LE.18.0D0) THEN
           BI0=1.0D0
           R=1.0D0
           DO 15 K=1,50
              R=0.25D0*R*X2/(K*K)
              BI0=BI0+R
              IF (DABS(R/BI0).LT.1.0D-15) GO TO 20
15         CONTINUE
20         BI1=1.0D0
           R=1.0D0
           DO 25 K=1,50
              R=0.25D0*R*X2/(K*(K+1))
              BI1=BI1+R
              IF (DABS(R/BI1).LT.1.0D-15) GO TO 30
25         CONTINUE
30         BI1=0.5D0*X*BI1
        ELSE
           DATA A/0.125D0,7.03125D-2,
     &            7.32421875D-2,1.1215209960938D-1,
     &            2.2710800170898D-1,5.7250142097473D-1,
     &            1.7277275025845D0,6.0740420012735D0,
     &            2.4380529699556D01,1.1001714026925D02,
     &            5.5133589612202D02,3.0380905109224D03/
           DATA B/-0.375D0,-1.171875D-1,
     &            -1.025390625D-1,-1.4419555664063D-1,
     &            -2.7757644653320D-1,-6.7659258842468D-1,
     &            -1.9935317337513D0,-6.8839142681099D0,
     &            -2.7248827311269D01,-1.2159789187654D02,
     &            -6.0384407670507D02,-3.3022722944809D03/
           K0=12
           IF (X.GE.35.0) K0=9
           IF (X.GE.50.0) K0=7
           CA=DEXP(X)/DSQRT(2.0D0*PI*X)
           BI0=1.0D0
           XR=1.0D0/X
           DO 35 K=1,K0
35            BI0=BI0+A(K)*XR**K
           BI0=CA*BI0
           BI1=1.0D0
           DO 40 K=1,K0
40            BI1=BI1+B(K)*XR**K
           BI1=CA*BI1
        ENDIF
        WW=0.0D0
        IF (X.LE.9.0D0) THEN
           CT=-(DLOG(X/2.0D0)+EL)
           BK0=0.0D0
           W0=0.0D0
           R=1.0D0
           DO 65 K=1,50
              W0=W0+1.0D0/K
              R=0.25D0*R/(K*K)*X2
              BK0=BK0+R*(W0+CT)
              IF (DABS((BK0-WW)/BK0).LT.1.0D-15) GO TO 70
65            WW=BK0
70         BK0=BK0+CT
        ELSE
           DATA A1/0.125D0,0.2109375D0,
     &             1.0986328125D0,1.1775970458984D01,
     &             2.1461706161499D02,5.9511522710323D03,
     &             2.3347645606175D05,1.2312234987631D07/
           CB=0.5D0/X
           XR2=1.0D0/X2
           BK0=1.0D0
           DO 75 K=1,8
75            BK0=BK0+A1(K)*XR2**K
           BK0=CB*BK0/BI0
        ENDIF
        BK1=(1.0D0/X-BI1*BK0)/BI0
        DI0=BI1
        DI1=BI0-BI1/X
        DK0=-BK1
        DK1=-BK0-BK1/X
        RETURN
        END

C       **********************************

        SUBROUTINE CPBDN(N,Z,CPB,CPD)
C
C       ==================================================
C       Purpose: Compute the parabolic cylinder functions
C                 Dn(z) and Dn'(z) for a complex argument
C       Input:   z --- Complex argument of Dn(z)
C                n --- Order of Dn(z)  ( n=0,±1,±2,… )
C       Output:  CPB(|n|) --- Dn(z)
C                CPD(|n|) --- Dn'(z)
C       Routines called:
C            (1) CPDSA for computing Dn(z) for a small |z|
C            (2) CPDLA for computing Dn(z) for a large |z|
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CPB(0:*),CPD(0:*)
        PI=3.141592653589793D0
        X=DBLE(Z)
        A0=CDABS(Z)
        C0=(0.0D0,0.0D0)
        CA0=CDEXP(-0.25D0*Z*Z)
        N0=0
        IF (N.GE.0) THEN
           CF0=CA0
           CF1=Z*CA0
           CPB(0)=CF0
           CPB(1)=CF1
           DO 10 K=2,N
              CF=Z*CF1-(K-1.0D0)*CF0
              CPB(K)=CF
              CF0=CF1
10            CF1=CF
        ELSE
           N0=-N
           IF (X.LE.0.0.OR.CDABS(Z).EQ.0.0) THEN
              CF0=CA0
              CPB(0)=CF0
              Z1=-Z
              IF (A0.LE.7.0) THEN
                 CALL CPDSA(-1,Z1,CF1)
              ELSE
                 CALL CPDLA(-1,Z1,CF1)
              ENDIF
              CF1=DSQRT(2.0D0*PI)/CA0-CF1
              CPB(1)=CF1
              DO 15 K=2,N0
                 CF=(-Z*CF1+CF0)/(K-1.0D0)
                 CPB(K)=CF
                 CF0=CF1
15               CF1=CF
           ELSE
              IF (A0.LE.3.0) THEN
                 CALL CPDSA(-N0,Z,CFA)
                 CPB(N0)=CFA
                 N1=N0+1
                 CALL CPDSA(-N1,Z,CFB)
                 CPB(N1)=CFB
                 NM1=N0-1
                 DO 20 K=NM1,0,-1
                    CF=Z*CFA+(K+1.0D0)*CFB
                    CPB(K)=CF
                    CFB=CFA
20                  CFA=CF
              ELSE
                 M=100+ABS(N)
                 CFA=C0
                 CFB=(1.0D-30,0.0D0)
                 DO 25 K=M,0,-1
                    CF=Z*CFB+(K+1.0D0)*CFA
                    IF (K.LE.N0) CPB(K)=CF
                    CFA=CFB
25                  CFB=CF
                 CS0=CA0/CF
                 DO 30 K=0,N0
30                  CPB(K)=CS0*CPB(K)
              ENDIF
           ENDIF
        ENDIF
        CPD(0)=-0.5D0*Z*CPB(0)
        IF (N.GE.0) THEN
           DO 35 K=1,N
35            CPD(K)=-0.5D0*Z*CPB(K)+K*CPB(K-1)
        ELSE
           DO 40 K=1,N0
40            CPD(K)=0.5D0*Z*CPB(K)-CPB(K-1)
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE IK01B(X,BI0,DI0,BI1,DI1,BK0,DK0,BK1,DK1)
C
C       =========================================================
C       Purpose: Compute modified Bessel functions I0(x), I1(1),
C                K0(x) and K1(x), and their derivatives
C       Input :  x   --- Argument ( x ≥ 0 )
C       Output:  BI0 --- I0(x)
C                DI0 --- I0'(x)
C                BI1 --- I1(x)
C                DI1 --- I1'(x)
C                BK0 --- K0(x)
C                DK0 --- K0'(x)
C                BK1 --- K1(x)
C                DK1 --- K1'(x)
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (X.EQ.0.0D0) THEN
           BI0=1.0D0
           BI1=0.0D0
           BK0=1.0D+300
           BK1=1.0D+300
           DI0=0.0D0
           DI1=0.5D0
           DK0=-1.0D+300
           DK1=-1.0D+300
           RETURN
        ELSE IF (X.LE.3.75D0) THEN
           T=X/3.75D0
           T2=T*T
           BI0=(((((.0045813D0*T2+.0360768D0)*T2+.2659732D0)
     &         *T2+1.2067492D0)*T2+3.0899424D0)*T2
     &         +3.5156229D0)*T2+1.0D0
           BI1=X*((((((.00032411D0*T2+.00301532D0)*T2
     &         +.02658733D0)*T2+.15084934D0)*T2+.51498869D0)
     &         *T2+.87890594D0)*T2+.5D0)
        ELSE
           T=3.75D0/X
           BI0=((((((((.00392377D0*T-.01647633D0)*T
     &         +.02635537D0)*T-.02057706D0)*T+.916281D-2)*T
     &         -.157565D-2)*T+.225319D-2)*T+.01328592D0)*T
     &         +.39894228D0)*DEXP(X)/DSQRT(X)
           BI1=((((((((-.420059D-2*T+.01787654D0)*T
     &         -.02895312D0)*T+.02282967D0)*T-.01031555D0)*T
     &         +.163801D-2)*T-.00362018D0)*T-.03988024D0)*T
     &         +.39894228D0)*DEXP(X)/DSQRT(X)
        ENDIF
        IF (X.LE.2.0D0) THEN
           T=X/2.0D0
           T2=T*T
           BK0=(((((.0000074D0*T2+.0001075D0)*T2+.00262698D0)
     &         *T2+.0348859D0)*T2+.23069756D0)*T2+.4227842D0)
     &         *T2-.57721566D0-BI0*DLOG(T)
           BK1=((((((-.00004686D0*T2-.00110404D0)*T2
     &         -.01919402D0)*T2-.18156897D0)*T2-.67278579D0)
     &         *T2+.15443144D0)*T2+1.0D0)/X+BI1*DLOG(T)
        ELSE
           T=2.0D0/X
           T2=T*T
           BK0=((((((.00053208D0*T-.0025154D0)*T+.00587872D0)
     &         *T-.01062446D0)*T+.02189568D0)*T-.07832358D0)
     &         *T+1.25331414D0)*DEXP(-X)/DSQRT(X)
           BK1=((((((-.00068245D0*T+.00325614D0)*T
     &         -.00780353D0)*T+.01504268D0)*T-.0365562D0)*T+
     &         .23498619D0)*T+1.25331414D0)*DEXP(-X)/DSQRT(X)
        ENDIF
        DI0=BI1
        DI1=BI0-BI1/X
        DK0=-BK1
        DK1=-BK0-BK1/X
        RETURN
        END

C       **********************************

        SUBROUTINE BETA(P,Q,BT)
C
C       ==========================================
C       Purpose: Compute the beta function B(p,q)
C       Input :  p  --- Parameter  ( p > 0 )
C                q  --- Parameter  ( q > 0 )
C       Output:  BT --- B(p,q)
C       Routine called: GAMMA2 for computing Г(x)
C       ==========================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        CALL GAMMA2(P,GP)
        CALL GAMMA2(Q,GQ)
        PPQ=P+Q
        CALL GAMMA2(PPQ,GPQ)
        BT=GP*GQ/GPQ
        RETURN
        END



C       **********************************

        SUBROUTINE LPN(N,X,PN,PD)
C
C       ===============================================
C       Purpose: Compute Legendre polynomials Pn(x)
C                and their derivatives Pn'(x)
C       Input :  x --- Argument of Pn(x)
C                n --- Degree of Pn(x) ( n = 0,1,...)
C       Output:  PN(n) --- Pn(x)
C                PD(n) --- Pn'(x)
C       ===============================================
C
        IMPLICIT DOUBLE PRECISION (P,X)
        DIMENSION PN(0:N),PD(0:N)
        PN(0)=1.0D0
        PN(1)=X
        PD(0)=0.0D0
        PD(1)=1.0D0
        P0=1.0D0
        P1=X
        DO 10 K=2,N
           PF=(2.0D0*K-1.0D0)/K*X*P1-(K-1.0D0)/K*P0
           PN(K)=PF
           IF (DABS(X).EQ.1.0D0) THEN
              PD(K)=0.5D0*X**(K+1)*K*(K+1.0D0)
           ELSE
              PD(K)=K*(P1-X*PF)/(1.0D0-X*X)
           ENDIF
           P0=P1
10         P1=PF
        RETURN
        END

C       **********************************

        SUBROUTINE FCOEF(KD,M,Q,A,FC)
C
C       =====================================================
C       Purpose: Compute expansion coefficients for Mathieu
C                functions and modified Mathieu functions
C       Input :  m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions
C                KD --- Case code
C                       KD=1 for cem(x,q)  ( m = 0,2,4,...)
C                       KD=2 for cem(x,q)  ( m = 1,3,5,...)
C                       KD=3 for sem(x,q)  ( m = 1,3,5,...)
C                       KD=4 for sem(x,q)  ( m = 2,4,6,...)
C                A  --- Characteristic value of Mathieu
C                       functions for given m and q
C       Output:  FC(k) --- Expansion coefficients of Mathieu
C                       functions ( k= 1,2,...,KM )
C                       FC(1),FC(2),FC(3),... correspond to
C                       A0,A2,A4,... for KD=1 case, A1,A3,
C                       A5,... for KD=2 case, B1,B3,B5,...
C                       for KD=3 case and B2,B4,B6,... for
C                       KD=4 case
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION FC(251)
        DO 5 I=1,251
5          FC(I)=0.0D0
        IF (DABS(Q).LE.1.0D-7) THEN
C          Expansion up to order Q^1 (Abramowitz & Stegun 20.2.27-28)
           IF (KD.EQ.1) THEN
              JM=M/2 + 1
           ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
              JM=(M-1)/2+1
           ELSE IF (KD.EQ.4) THEN
              JM=M/2
           END IF
C          Check for overflow
           IF (JM+1.GT.251) GOTO 6
C          Proceed using the simplest expansion
           IF (KD.EQ.1.OR.KD.EQ.2) THEN
              IF (M.EQ.0) THEN
                 FC(1) = 1/SQRT(2.0D0)
                 FC(2) = -Q/2.0D0/SQRT(2.0D0)
              ELSE IF (M.EQ.1) THEN
                 FC(1) = 1.0D0
                 FC(2) = -Q/8.0D0
              ELSE IF (M.EQ.2) THEN
                 FC(1) = Q/4.0D0
                 FC(2) = 1.0D0
                 FC(3) = -Q/12.0D0
              ELSE
                 FC(JM) = 1.0D0
                 FC(JM+1) = -Q/(4.0D0 * (M + 1))
                 FC(JM-1) =  Q/(4.0D0 * (M - 1))
              END IF
           ELSE IF (KD.EQ.3.OR.KD.EQ.4) THEN
              IF (M.EQ.1) THEN
                 FC(1) = 1.0D0
                 FC(2) = -Q/8.0D0
              ELSE IF (M.EQ.2) THEN
                 FC(1) = 1.0D0
                 FC(2) = -Q/12.0D0
              ELSE
                 FC(JM) = 1.0D0
                 FC(JM+1) = -Q/(4.0D0 * (M + 1))
                 FC(JM-1) =  Q/(4.0D0 * (M - 1))
              END IF
           ENDIF
           RETURN
        ELSE IF (Q.LE.1.0D0) THEN
           QM=7.5+56.1*SQRT(Q)-134.7*Q+90.7*SQRT(Q)*Q
        ELSE
           QM=17.0+3.1*SQRT(Q)-.126*Q+.0037*SQRT(Q)*Q
        ENDIF
        KM=INT(QM+0.5*M)
        IF (KM.GT.251) THEN
C          Overflow, generate NaNs
 6         FNAN=DNAN()
           DO 7 I=1,251
 7            FC(I)=FNAN
           RETURN
        ENDIF
        KB=0
        S=0.0D0
        F=1.0D-100
        U=0.0D0
        FC(KM)=0.0D0
        F2=0.0D0
        IF (KD.EQ.1) THEN
           DO 25 K=KM,3,-1
              V=U
              U=F
              F=(A-4.0D0*K*K)*U/Q-V
              IF (DABS(F).LT.DABS(FC(K+1))) THEN
                 KB=K
                 FC(1)=1.0D-100
                 SP=0.0D0
                 F3=FC(K+1)
                 FC(2)=A/Q*FC(1)
                 FC(3)=(A-4.0D0)*FC(2)/Q-2.0D0*FC(1)
                 U=FC(2)
                 F1=FC(3)
                 DO 15 I=3,KB
                    V=U
                    U=F1
                    F1=(A-4.0D0*(I-1.0D0)**2)*U/Q-V
                    FC(I+1)=F1
                    IF (I.EQ.KB) F2=F1
                    IF (I.NE.KB) SP=SP+F1*F1
15               CONTINUE
                 SP=SP+2.0D0*FC(1)**2+FC(2)**2+FC(3)**2
                 SS=S+SP*(F3/F2)**2
                 S0=DSQRT(1.0D0/SS)
                 DO 20 J=1,KM
                    IF (J.LE.KB+1) THEN
                       FC(J)=S0*FC(J)*F3/F2
                    ELSE
                       FC(J)=S0*FC(J)
                    ENDIF
20               CONTINUE
                 GO TO 85
              ELSE
                 FC(K)=F
                 S=S+F*F
              ENDIF
25         CONTINUE
           FC(2)=Q*FC(3)/(A-4.0D0-2.0D0*Q*Q/A)
           FC(1)=Q/A*FC(2)
           S=S+2.0D0*FC(1)**2+FC(2)**2
           S0=DSQRT(1.0D0/S)
           DO 30 K=1,KM
30            FC(K)=S0*FC(K)
        ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
           DO 35 K=KM,3,-1
              V=U
              U=F
              F=(A-(2.0D0*K-1)**2)*U/Q-V
              IF (DABS(F).GE.DABS(FC(K))) THEN
                 FC(K-1)=F
                 S=S+F*F
              ELSE
                 KB=K
                 F3=FC(K)
                 GO TO 45
              ENDIF
35         CONTINUE
           FC(1)=Q/(A-1.0D0-(-1)**KD*Q)*FC(2)
           S=S+FC(1)*FC(1)
           S0=DSQRT(1.0D0/S)
           DO 40 K=1,KM
40            FC(K)=S0*FC(K)
           GO TO 85
45         FC(1)=1.0D-100
           FC(2)=(A-1.0D0-(-1)**KD*Q)/Q*FC(1)
           SP=0.0D0
           U=FC(1)
           F1=FC(2)
           DO 50 I=2,KB-1
              V=U
              U=F1
              F1=(A-(2.0D0*I-1.0D0)**2)*U/Q-V
              IF (I.NE.KB-1) THEN
                 FC(I+1)=F1
                 SP=SP+F1*F1
              ELSE
                 F2=F1
              ENDIF
50         CONTINUE
           SP=SP+FC(1)**2+FC(2)**2
           SS=S+SP*(F3/F2)**2
           S0=1.0D0/DSQRT(SS)
           DO 55 J=1,KM
              IF (J.LT.KB) FC(J)=S0*FC(J)*F3/F2
              IF (J.GE.KB) FC(J)=S0*FC(J)
55         CONTINUE
        ELSE IF (KD.EQ.4) THEN
           DO 60 K=KM,3,-1
              V=U
              U=F
              F=(A-4.0D0*K*K)*U/Q-V
              IF (DABS(F).GE.DABS(FC(K))) THEN
                 FC(K-1)=F
                 S=S+F*F
              ELSE
                 KB=K
                 F3=FC(K)
                 GO TO 70
              ENDIF
60         CONTINUE
           FC(1)=Q/(A-4.0D0)*FC(2)
           S=S+FC(1)*FC(1)
           S0=DSQRT(1.0D0/S)
           DO 65 K=1,KM
65            FC(K)=S0*FC(K)
           GO TO 85
70         FC(1)=1.0D-100
           FC(2)=(A-4.0D0)/Q*FC(1)
           SP=0.0D0
           U=FC(1)
           F1=FC(2)
           DO 75 I=2,KB-1
              V=U
              U=F1
              F1=(A-4.0D0*I*I)*U/Q-V
              IF (I.NE.KB-1) THEN
                 FC(I+1)=F1
                 SP=SP+F1*F1
              ELSE
                 F2=F1
              ENDIF
75         CONTINUE
           SP=SP+FC(1)**2+FC(2)**2
           SS=S+SP*(F3/F2)**2
           S0=1.0D0/DSQRT(SS)
           DO 80 J=1,KM
              IF (J.LT.KB) FC(J)=S0*FC(J)*F3/F2
              IF (J.GE.KB) FC(J)=S0*FC(J)
80         CONTINUE
        ENDIF
85      IF (FC(1).LT.0.0D0) THEN
           DO 90 J=1,KM
90            FC(J)=-FC(J)
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE SPHI(N,X,NM,SI,DI)
C
C       ========================================================
C       Purpose: Compute modified spherical Bessel functions
C                of the first kind, in(x) and in'(x)
C       Input :  x --- Argument of in(x)
C                n --- Order of in(x) ( n = 0,1,2,... )
C       Output:  SI(n) --- in(x)
C                DI(n) --- in'(x)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SI(0:N),DI(0:N)
        NM=N
        IF (DABS(X).LT.1.0D-100) THEN
           DO 10 K=0,N
              SI(K)=0.0D0
10            DI(K)=0.0D0
           SI(0)=1.0D0
           DI(1)=0.333333333333333D0
           RETURN
        ENDIF
        SI(0)=DSINH(X)/X
        SI(1)=-(DSINH(X)/X-DCOSH(X))/X
        SI0=SI(0)
        IF (N.GE.2) THEN
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F=0.0D0
           F0=0.0D0
           F1=1.0D0-100
           DO 15 K=M,0,-1
              F=(2.0D0*K+3.0D0)*F1/X+F0
              IF (K.LE.NM) SI(K)=F
              F0=F1
15            F1=F
           CS=SI0/F
           DO 20 K=0,NM
20            SI(K)=CS*SI(K)
        ENDIF
        DI(0)=SI(1)
        DO 25 K=1,NM
25         DI(K)=SI(K-1)-(K+1.0D0)/X*SI(K)
        RETURN
        END



C       **********************************

        SUBROUTINE PBWA(A,X,W1F,W1D,W2F,W2D)
C
C       ======================================================
C       Purpose: Compute parabolic cylinder functions W(a,±x)
C                and their derivatives
C       Input  : a --- Parameter  ( 0 ≤ |a| ≤ 5 )
C                x --- Argument of W(a,±x)  ( 0 ≤ |x| ≤ 5 )
C       Output : W1F --- W(a,x)
C                W1D --- W'(a,x)
C                W2F --- W(a,-x)
C                W2D --- W'(a,-x)
C       Routine called:
C               CGAMA for computing complex gamma function
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX *16 (C,Z)
        DIMENSION H(100),D(100)
        EPS=1.0D-15
        P0=0.59460355750136D0
        IF (A.EQ.0.0D0) THEN
           G1=3.625609908222D0
           G2=1.225416702465D0
        ELSE
           X1=0.25D0
           Y1=0.5D0*A
           CALL CGAMA(X1,Y1,1,UGR,UGI)
           G1=DSQRT(UGR*UGR+UGI*UGI)
           X2=0.75D0
           CALL CGAMA(X2,Y1,1,VGR,VGI)
           G2=DSQRT(VGR*VGR+VGI*VGI)
        ENDIF
        F1=DSQRT(G1/G2)
        F2=DSQRT(2.0D0*G2/G1)
        H0=1.0D0
        H1=A
        H(1)=A
        DO 10 L1=4,200,2
           M=L1/2
           HL=A*H1-0.25D0*(L1-2.0D0)*(L1-3.0D0)*H0
           H(M)=HL
           H0=H1
10         H1=HL
        Y1F=1.0D0
        R=1.0D0
        DO 15 K=1,100
           R=0.5D0*R*X*X/(K*(2.0D0*K-1.0D0))
           R1=H(K)*R
           Y1F=Y1F+R1
           IF (DABS(R1/Y1F).LE.EPS.AND.K.GT.30) GO TO 20
15      CONTINUE
20      Y1D=A
        R=1.0D0
        DO 25 K=1,100
           R=0.5D0*R*X*X/(K*(2.0D0*K+1.0D0))
           R1=H(K+1)*R
           Y1D=Y1D+R1
           IF (DABS(R1/Y1D).LE.EPS.AND.K.GT.30) GO TO 30
25      CONTINUE
30      Y1D=X*Y1D
        D1=1.0D0
        D2=A
        D(1)=1.0D0
        D(2)=A
        DO 40 L2=5,160,2
           M=(L2+1)/2
           DL=A*D2-0.25D0*(L2-2.0D0)*(L2-3.0D0)*D1
           D(M)=DL
           D1=D2
40         D2=DL
        Y2F=1.0D0
        R=1.0D0
        DO 45 K=1,100
           R=0.5D0*R*X*X/(K*(2.0D0*K+1.0D0))
           R1=D(K+1)*R
           Y2F=Y2F+R1
           IF (DABS(R1/Y2F).LE.EPS.AND.K.GT.30) GO TO 50
45      CONTINUE
50      Y2F=X*Y2F
        Y2D=1.0D0
        R=1.0D0
        DO 55 K=1,100
           R=0.5D0*R*X*X/(K*(2.0D0*K-1.0D0))
           R1=D(K+1)*R
           Y2D=Y2D+R1
           IF (DABS(R1/Y2D).LE.EPS.AND.K.GT.30) GO TO 60
55      CONTINUE
60      W1F=P0*(F1*Y1F-F2*Y2F)
        W2F=P0*(F1*Y1F+F2*Y2F)
        W1D=P0*(F1*Y1D-F2*Y2D)
        W2D=P0*(F1*Y1D+F2*Y2D)
        RETURN
        END



C       **********************************

        SUBROUTINE RMN1(M,N,C,X,DF,KD,R1F,R1D)
C
C       =======================================================
C       Purpose: Compute prolate and oblate spheroidal radial
C                functions of the first kind for given m, n,
C                c and x
C       Routines called:
C            (1) SCKB for computing expansion coefficients c2k
C            (2) SPHJ for computing the spherical Bessel
C                functions of the first kind
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION CK(200),DF(200),SJ(0:251),DJ(0:251)
        EPS=1.0D-14
        IP=1
        NM1=INT((N-M)/2)
        IF (N-M.EQ.2*NM1) IP=0
        NM=25+NM1+INT(C)
        REG=1.0D0
        IF (M+NM.GT.80) REG=1.0D-200
        R0=REG
        DO 10 J=1,2*M+IP
10         R0=R0*J
        R=R0
        SUC=R*DF(1)
        SW=0.0D0
        DO 15 K=2,NM
           R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
           SUC=SUC+R*DF(K)
           IF (K.GT.NM1.AND.DABS(SUC-SW).LT.DABS(SUC)*EPS) GO TO 20
15         SW=SUC
20      CONTINUE
        IF (X.EQ.0.0) THEN
           CALL SCKB(M,N,C,DF,CK)
           SUM=0.0D0
           SW1=0.0D0
           DO 25 J=1,NM
              SUM=SUM+CK(J)
              IF (DABS(SUM-SW1).LT.DABS(SUM)*EPS) GO TO 30
25            SW1=SUM
30         R1=1.0D0
           DO 35 J=1,(N+M+IP)/2
35            R1=R1*(J+0.5D0*(N+M+IP))
           R2=1.0D0
           DO 40 J=1,M
40            R2=2.0D0*C*R2*J
           R3=1.0D0
           DO 45 J=1,(N-M-IP)/2
45            R3=R3*J
           SA0=(2.0*(M+IP)+1.0)*R1/(2.0**N*C**IP*R2*R3)
           IF (IP.EQ.0) THEN
              R1F=SUM/(SA0*SUC)*DF(1)*REG
              R1D=0.0D0
           ELSE IF (IP.EQ.1) THEN
              R1F=0.0D0
              R1D=SUM/(SA0*SUC)*DF(1)*REG
           ENDIF
           RETURN
        ENDIF
        CX=C*X
        NM2=2*NM+M
        CALL SPHJ(NM2,CX,NM2,SJ,DJ)
        A0=(1.0D0-KD/(X*X))**(0.5D0*M)/SUC
        R1F=0.0D0
        SW=0.0D0
        LG=0
        DO 50 K=1,NM
           L=2*K+M-N-2+IP
           IF (L.EQ.4*INT(L/4)) LG=1
           IF (L.NE.4*INT(L/4)) LG=-1
           IF (K.EQ.1) THEN
              R=R0
           ELSE
              R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
           ENDIF
           NP=M+2*K-2+IP
           R1F=R1F+LG*R*DF(K)*SJ(NP)
           IF (K.GT.NM1.AND.DABS(R1F-SW).LT.DABS(R1F)*EPS) GO TO 55
50         SW=R1F
55      R1F=R1F*A0
        B0=KD*M/X**3.0D0/(1.0-KD/(X*X))*R1F
        SUD=0.0D0
        SW=0.0D0
        DO 60 K=1,NM
           L=2*K+M-N-2+IP
           IF (L.EQ.4*INT(L/4)) LG=1
           IF (L.NE.4*INT(L/4)) LG=-1
           IF (K.EQ.1) THEN
              R=R0
           ELSE
              R=R*(M+K-1.0)*(M+K+IP-1.5D0)/(K-1.0D0)/(K+IP-1.5D0)
           ENDIF
           NP=M+2*K-2+IP
           SUD=SUD+LG*R*DF(K)*DJ(NP)
           IF (K.GT.NM1.AND.DABS(SUD-SW).LT.DABS(SUD)*EPS) GO TO 65
60         SW=SUD
65      R1D=B0+A0*C*SUD
        RETURN
        END



C       **********************************

        SUBROUTINE DVSA(VA,X,PD)
C
C       ===================================================
C       Purpose: Compute parabolic cylinder function Dv(x)
C                for small argument
C       Input:   x  --- Argument
C                va --- Order
C       Output:  PD --- Dv(x)
C       Routine called: GAMMA2 for computing Г(x)
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EPS=1.0D-15
        PI=3.141592653589793D0
        SQ2=DSQRT(2.0D0)
        EP=DEXP(-.25D0*X*X)
        VA0=0.5D0*(1.0D0-VA)
        IF (VA.EQ.0.0) THEN
           PD=EP
        ELSE
           IF (X.EQ.0.0) THEN
              IF (VA0.LE.0.0.AND.VA0.EQ.INT(VA0)) THEN
                 PD=0.0D0
              ELSE
                 CALL GAMMA2(VA0,GA0)
                 PD=DSQRT(PI)/(2.0D0**(-.5D0*VA)*GA0)
              ENDIF
           ELSE
              CALL GAMMA2(-VA,G1)
              A0=2.0D0**(-0.5D0*VA-1.0D0)*EP/G1
              VT=-.5D0*VA
              CALL GAMMA2(VT,G0)
              PD=G0
              R=1.0D0
              DO 10 M=1,250
                 VM=.5D0*(M-VA)
                 CALL GAMMA2(VM,GM)
                 R=-R*SQ2*X/M
                 R1=GM*R
                 PD=PD+R1
                 IF (DABS(R1).LT.DABS(PD)*EPS) GO TO 15
10            CONTINUE
15            PD=A0*PD
           ENDIF
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE E1Z(Z,CE1)
C
C       ====================================================
C       Purpose: Compute complex exponential integral E1(z)
C       Input :  z   --- Argument of E1(z)
C       Output:  CE1 --- E1(z)
C       ====================================================
C
        IMPLICIT COMPLEX*16 (C,Z)
        IMPLICIT DOUBLE PRECISION (A,D-H,O-Y)
        PI=3.141592653589793D0
        EL=0.5772156649015328D0
        X=DBLE(Z)
        A0=CDABS(Z)
C       Continued fraction converges slowly near negative real axis,
C       so use power series in a wedge around it until radius 40.0
        XT=-2*DABS(DIMAG(Z))
        IF (A0.EQ.0.0D0) THEN
           CE1=(1.0D+300,0.0D0)
        ELSE IF (A0.LE.5.0.OR.X.LT.XT.AND.A0.LT.40.0) THEN
C          Power series
           CE1=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 10 K=1,500
              CR=-CR*K*Z/(K+1.0D0)**2
              CE1=CE1+CR
              IF (CDABS(CR).LE.CDABS(CE1)*1.0D-15) GO TO 15
10         CONTINUE
15         CONTINUE
           IF (X.LE.0.0.AND.DIMAG(Z).EQ.0.0) THEN
C             Careful on the branch cut -- avoid signed zeros
              CE1=-EL-CDLOG(-Z)+Z*CE1-PI*(0.0D0,1.0D0)
           ELSE
              CE1=-EL-CDLOG(Z)+Z*CE1
           ENDIF
        ELSE
C          Continued fraction http://dlmf.nist.gov/6.9
C
C                           1     1     1     2     2     3     3
C          E1 = exp(-z) * ----- ----- ----- ----- ----- ----- ----- ...
C                         Z +   1 +   Z +   1 +   Z +   1 +   Z +
           ZC=0D0
           ZD=1/Z
           ZDC=1*ZD
           ZC=ZC + ZDC
           DO 20 K=1,500
              ZD=1/(ZD*K + 1)
              ZDC=(1*ZD - 1)*ZDC
              ZC=ZC + ZDC

              ZD=1/(ZD*K + Z)
              ZDC=(Z*ZD - 1)*ZDC
              ZC=ZC + ZDC

              IF (CDABS(ZDC).LE.CDABS(ZC)*1.0D-15.AND.K.GT.20) GO TO 25
20         CONTINUE
25         CE1=CDEXP(-Z)*ZC
           IF (X.LE.0.0.AND.DIMAG(Z).EQ.0.0) CE1=CE1-PI*(0.0D0,1.0D0)
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE ITJYB(X,TJ,TY)
C
C       =======================================================
C       Purpose: Integrate Bessel functions J0(t) and Y0(t)
C                with respect to t from 0 to x ( x ≥ 0 )
C       Input :  x  --- Upper limit of the integral
C       Output:  TJ --- Integration of J0(t) from 0 to x
C                TY --- Integration of Y0(t) from 0 to x
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           TJ=0.0D0
           TY=0.0D0
        ELSE IF (X.LE.4.0D0) THEN
           X1=X/4.0D0
           T=X1*X1
           TJ=(((((((-.133718D-3*T+.2362211D-2)*T
     &        -.025791036D0)*T+.197492634D0)*T-1.015860606D0)
     &        *T+3.199997842D0)*T-5.333333161D0)*T+4.0D0)*X1
           TY=((((((((.13351D-4*T-.235002D-3)*T+.3034322D-2)*
     &        T-.029600855D0)*T+.203380298D0)*T-.904755062D0)
     &        *T+2.287317974D0)*T-2.567250468D0)*T
     &        +1.076611469D0)*X1
           TY=2.0D0/PI*DLOG(X/2.0D0)*TJ-TY
        ELSE IF (X.LE.8.0D0) THEN
           XT=X-.25D0*PI
           T=16.0D0/(X*X)
           F0=((((((.1496119D-2*T-.739083D-2)*T+.016236617D0)
     &        *T-.022007499D0)*T+.023644978D0)
     &        *T-.031280848D0)*T+.124611058D0)*4.0D0/X
           G0=(((((.1076103D-2*T-.5434851D-2)*T+.01242264D0)
     &        *T-.018255209)*T+.023664841D0)*T-.049635633D0)
     &        *T+.79784879D0
           TJ=1.0D0-(F0*DCOS(XT)-G0*DSIN(XT))/DSQRT(X)
           TY=-(F0*DSIN(XT)+G0*DCOS(XT))/DSQRT(X)
        ELSE
           T=64.0D0/(X*X)
           XT=X-.25D0*PI
           F0=(((((((-.268482D-4*T+.1270039D-3)*T
     &        -.2755037D-3)*T+.3992825D-3)*T-.5366169D-3)*T
     &        +.10089872D-2)*T-.40403539D-2)*T+.0623347304D0)
     &        *8.0D0/X
           G0=((((((-.226238D-4*T+.1107299D-3)*T-.2543955D-3)
     &        *T+.4100676D-3)*T-.6740148D-3)*T+.17870944D-2)
     &        *T-.01256424405D0)*T+.79788456D0
           TJ=1.0D0-(F0*DCOS(XT)-G0*DSIN(XT))/DSQRT(X)
           TY=-(F0*DSIN(XT)+G0*DCOS(XT))/DSQRT(X)
        ENDIF
        RETURN
        END


C       **********************************

        SUBROUTINE CHGUL(A,B,X,HU,ID)
C
C       =======================================================
C       Purpose: Compute the confluent hypergeometric function
C                U(a,b,x) for large argument x
C       Input  : a  --- Parameter
C                b  --- Parameter
C                x  --- Argument
C       Output:  HU --- U(a,b,x)
C                ID --- Estimated number of significant digits
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL IL1,IL2
        ID=-100
        AA=A-B+1.0D0
        IL1=A.EQ.INT(A).AND.A.LE.0.0
        IL2=AA.EQ.INT(AA).AND.AA.LE.0.0
        NM=0
        IF (IL1) NM=ABS(A)
        IF (IL2) NM=ABS(AA)
C       IL1: DLMF 13.2.7 with k=-s-a
C       IL2: DLMF 13.2.8
        IF (IL1.OR.IL2) THEN
           HU=1.0D0
           R=1.0D0
           DO 10 K=1,NM
              R=-R*(A+K-1.0D0)*(A-B+K)/(K*X)
              HU=HU+R
10         CONTINUE
           HU=X**(-A)*HU
           ID=10
        ELSE
C       DLMF 13.7.3
           HU=1.0D0
           R=1.0D0
           DO 15 K=1,25
              R=-R*(A+K-1.0D0)*(A-B+K)/(K*X)
              RA=DABS(R)
              IF (K.GT.5.AND.RA.GE.R0.OR.RA.LT.1.0D-15) GO TO 20
              R0=RA
15            HU=HU+R
20         ID=ABS(LOG10(RA))
           HU=X**(-A)*HU
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE GMN(M,N,C,X,BK,GF,GD)
C
C       ===========================================================
C       Purpose: Compute gmn(-ic,ix) and its derivative for oblate
C                radial functions with a small argument
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BK(200)
        EPS=1.0D-14
        IP=1
        IF (N-M.EQ.2*INT((N-M)/2)) IP=0
        NM=25+INT(0.5*(N-M)+C)
        XM=(1.0D0+X*X)**(-0.5D0*M)
        GF0=0.0D0
        GW=0.0D0
        DO 10 K=1,NM
           GF0=GF0+BK(K)*X**(2.0*K-2.0)
           IF (DABS((GF0-GW)/GF0).LT.EPS.AND.K.GE.10) GO TO 15
10         GW=GF0
15      GF=XM*GF0*X**(1-IP)
        GD1=-M*X/(1.0D0+X*X)*GF
        GD0=0.0D0
        DO 20 K=1,NM
           IF (IP.EQ.0) THEN
              GD0=GD0+(2.0D0*K-1.0)*BK(K)*X**(2.0*K-2.0)
           ELSE
              GD0=GD0+2.0D0*K*BK(K+1)*X**(2.0*K-1.0)
           ENDIF
           IF (DABS((GD0-GW)/GD0).LT.EPS.AND.K.GE.10) GO TO 25
20         GW=GD0
25      GD=GD1+XM*GD0
        RETURN
        END



C       **********************************

        SUBROUTINE ITJYA(X,TJ,TY)
C
C       ==========================================================
C       Purpose: Integrate Bessel functions J0(t) & Y0(t) with
C                respect to t from 0 to x
C       Input :  x  --- Upper limit of the integral ( x >= 0 )
C       Output:  TJ --- Integration of J0(t) from 0 to x
C                TY --- Integration of Y0(t) from 0 to x
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(18)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        EPS=1.0D-12
        IF (X.EQ.0.0D0) THEN
           TJ=0.0D0
           TY=0.0D0
        ELSE IF (X.LE.20.0D0) THEN
           X2=X*X
           TJ=X
           R=X
           DO 10 K=1,60
              R=-.25D0*R*(2*K-1.0D0)/(2*K+1.0D0)/(K*K)*X2
              TJ=TJ+R
              IF (DABS(R).LT.DABS(TJ)*EPS) GO TO 15
10         CONTINUE
15         TY1=(EL+DLOG(X/2.0D0))*TJ
           RS=0.0D0
           TY2=1.0D0
           R=1.0D0
           DO 20 K=1,60
              R=-.25D0*R*(2*K-1.0D0)/(2*K+1.0D0)/(K*K)*X2
              RS=RS+1.0D0/K
              R2=R*(RS+1.0D0/(2.0D0*K+1.0D0))
              TY2=TY2+R2
              IF (DABS(R2).LT.DABS(TY2)*EPS) GO TO 25
20         CONTINUE
25         TY=(TY1-X*TY2)*2.0D0/PI
        ELSE
           A0=1.0D0
           A1=5.0D0/8.0D0
           A(1)=A1
           DO 30 K=1,16
              AF=((1.5D0*(K+.5D0)*(K+5.0D0/6.0D0)*A1-.5D0
     &           *(K+.5D0)*(K+.5D0)*(K-.5D0)*A0))/(K+1.0D0)
              A(K+1)=AF
              A0=A1
30            A1=AF
           BF=1.0D0
           R=1.0D0
           DO 35 K=1,8
              R=-R/(X*X)
35            BF=BF+A(2*K)*R
           BG=A(1)/X
           R=1.0D0/X
           DO 40 K=1,8
              R=-R/(X*X)
40            BG=BG+A(2*K+1)*R
           XP=X+.25D0*PI
           RC=DSQRT(2.0D0/(PI*X))
           TJ=1.0D0-RC*(BF*DCOS(XP)+BG*DSIN(XP))
           TY=RC*(BG*DCOS(XP)-BF*DSIN(XP))
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE STVLV(V,X,SLV)
C
C       ======================================================
C       Purpose:  Compute modified Struve function Lv(x) with
C                 an arbitrary order v
C       Input :   v   --- Order of Lv(x)  ( |v| ≤ 20 )
C                 x   --- Argument of Lv(x) ( x ≥ 0 )
C       Output:   SLV --- Lv(x)
C       Routine called: GAMMA2 to compute the gamma function
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           IF (V.GT.-1.0.OR.INT(V)-V.EQ.0.5D0) THEN
              SLV=0.0D0
           ELSE IF (V.LT.-1.0D0) THEN
              SLV=(-1)**(INT(0.5D0-V)-1)*1.0D+300
           ELSE IF (V.EQ.-1.0D0) THEN
              SLV=2.0D0/PI
           ENDIF
           RETURN
        ENDIF
        IF (X.LE.40.0D0) THEN
           V0=V+1.5D0
           CALL GAMMA2(V0,GA)
           S=2.0D0/(DSQRT(PI)*GA)
           R1=1.0D0
           DO 10 K=1,100
              VA=K+1.5D0
              CALL GAMMA2(VA,GA)
              VB=V+K+1.5D0
              CALL GAMMA2(VB,GB)
              R1=R1*(0.5D0*X)**2
              R2=R1/(GA*GB)
              S=S+R2
              IF (DABS(R2/S).LT.1.0D-12) GO TO 15
10         CONTINUE
15         SLV=(0.5D0*X)**(V+1.0D0)*S
        ELSE
           SA=-1.0D0/PI*(0.5D0*X)**(V-1.0)
           V0=V+0.5D0
           CALL GAMMA2(V0,GA)
           S=-DSQRT(PI)/GA
           R1=-1.0D0
           DO 20 K=1,12
              VA=K+0.5D0
              CALL GAMMA2(VA,GA)
              VB=-K+V+0.5D0
              CALL GAMMA2(VB,GB)
              R1=-R1/(0.5D0*X)**2
              S=S+R1*GA/GB
20         CONTINUE
           S0=SA*S
           U=DABS(V)
           N=INT(U)
           U0=U-N
           BIV0=0.0D0
           DO 35 L=0,1
              VT=U0+L
              R=1.0D0
              BIV=1.0D0
              DO 25 K=1,16
                 R=-0.125*R*(4.0*VT*VT-(2.0*K-1.0D0)**2)/(K*X)
                 BIV=BIV+R
                 IF (DABS(R/BIV).LT.1.0D-12) GO TO 30
25            CONTINUE
30            IF (L.EQ.0) BIV0=BIV
35         CONTINUE
           BF=0.0D0
           BF0=BIV0
           BF1=BIV
           DO 40 K=2,N
              BF=-2.0D0*(K-1.0+U0)/X*BF1+BF0
              BF0=BF1
40            BF1=BF
           IF (N.EQ.0) BIV=BIV0
           IF (N.GT.1) BIV=BF
           SLV=DEXP(X)/DSQRT(2.0D0*PI*X)*BIV+S0
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE RCTY(N,X,NM,RY,DY)
C
C       ========================================================
C       Purpose: Compute Riccati-Bessel functions of the second
C                kind and their derivatives
C       Input:   x --- Argument of Riccati-Bessel function
C                n --- Order of yn(x)
C       Output:  RY(n) --- x·yn(x)
C                DY(n) --- [x·yn(x)]'
C                NM --- Highest order computed
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION RY(0:N),DY(0:N)
        NM=N
        IF (X.LT.1.0D-60) THEN
           DO 10 K=0,N
              RY(K)=-1.0D+300
10            DY(K)=1.0D+300
           RY(0)=-1.0D0
           DY(0)=0.0D0
           RETURN
        ENDIF
        RY(0)=-DCOS(X)
        RY(1)=RY(0)/X-DSIN(X)
        RF0=RY(0)
        RF1=RY(1)
        DO 15 K=2,N
           RF2=(2.0D0*K-1.0D0)*RF1/X-RF0
           IF (DABS(RF2).GT.1.0D+300) GO TO 20
           RY(K)=RF2
           RF0=RF1
15         RF1=RF2
20      NM=K-1
        DY(0)=DSIN(X)
        DO 25 K=1,NM
25         DY(K)=-K*RY(K)/X+RY(K-1)
        RETURN
        END

C       **********************************

        SUBROUTINE LPNI(N,X,PN,PD,PL)
C
C       =====================================================
C       Purpose: Compute Legendre polynomials Pn(x), Pn'(x)
C                and the integral of Pn(t) from 0 to x
C       Input :  x --- Argument of Pn(x)
C                n --- Degree of Pn(x) ( n = 0,1,... )
C       Output:  PN(n) --- Pn(x)
C                PD(n) --- Pn'(x)
C                PL(n) --- Integral of Pn(t) from 0 to x
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (P,R,X)
        DIMENSION PN(0:N),PD(0:N),PL(0:N)
        PN(0)=1.0D0
        PN(1)=X
        PD(0)=0.0D0
        PD(1)=1.0D0
        PL(0)=X
        PL(1)=0.5D0*X*X
        P0=1.0D0
        P1=X
        DO 15 K=2,N
           PF=(2.0D0*K-1.0D0)/K*X*P1-(K-1.0D0)/K*P0
           PN(K)=PF
           IF (DABS(X).EQ.1.0D0) THEN
              PD(K)=0.5D0*X**(K+1)*K*(K+1.0D0)
           ELSE
              PD(K)=K*(P1-X*PF)/(1.0D0-X*X)
           ENDIF
           PL(K)=(X*PN(K)-PN(K-1))/(K+1.0D0)
           P0=P1
           P1=PF
           IF (K.EQ.2*INT(K/2)) GO TO 15
           R=1.0D0/(K+1.0D0)
           N1=(K-1)/2
           DO 10 J=1,N1
10            R=(0.5D0/J-1.0D0)*R
           PL(K)=PL(K)+R
15      CONTINUE
        RETURN
        END

C       **********************************

        SUBROUTINE KLVNA(X,BER,BEI,GER,GEI,DER,DEI,HER,HEI)
C
C       ======================================================
C       Purpose: Compute Kelvin functions ber x, bei x, ker x
C                and kei x, and their derivatives  ( x > 0 )
C       Input :  x   --- Argument of Kelvin functions
C       Output:  BER --- ber x
C                BEI --- bei x
C                GER --- ker x
C                GEI --- kei x
C                DER --- ber'x
C                DEI --- bei'x
C                HER --- ker'x
C                HEI --- kei'x
C       ================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        EPS=1.0D-15
        IF (X.EQ.0.0D0) THEN
           BER=1.0D0
           BEI=0.0D0
           GER=1.0D+300
           GEI=-0.25D0*PI
           DER=0.0D0
           DEI=0.0D0
           HER=-1.0D+300
           HEI=0.0D0
           RETURN
        ENDIF
        X2=0.25D0*X*X
        X4=X2*X2
        IF (DABS(X).LT.10.0D0) THEN
           BER=1.0D0
           R=1.0D0
           DO 10 M=1,60
              R=-0.25D0*R/(M*M)/(2.0D0*M-1.0D0)**2*X4
              BER=BER+R
              IF (DABS(R).LT.DABS(BER)*EPS) GO TO 15
10         CONTINUE
15         BEI=X2
           R=X2
           DO 20 M=1,60
              R=-0.25D0*R/(M*M)/(2.0D0*M+1.0D0)**2*X4
              BEI=BEI+R
              IF (DABS(R).LT.DABS(BEI)*EPS) GO TO 25
20         CONTINUE
25         GER=-(DLOG(X/2.0D0)+EL)*BER+0.25D0*PI*BEI
           R=1.0D0
           GS=0.0D0
           DO 30 M=1,60
              R=-0.25D0*R/(M*M)/(2.0D0*M-1.0D0)**2*X4
              GS=GS+1.0D0/(2.0D0*M-1.0D0)+1.0D0/(2.0D0*M)
              GER=GER+R*GS
              IF (DABS(R*GS).LT.DABS(GER)*EPS) GO TO 35
30         CONTINUE
35         GEI=X2-(DLOG(X/2.0D0)+EL)*BEI-0.25D0*PI*BER
           R=X2
           GS=1.0D0
           DO 40 M=1,60
              R=-0.25D0*R/(M*M)/(2.0D0*M+1.0D0)**2*X4
              GS=GS+1.0D0/(2.0D0*M)+1.0D0/(2.0D0*M+1.0D0)
              GEI=GEI+R*GS
              IF (DABS(R*GS).LT.DABS(GEI)*EPS) GO TO 45
40         CONTINUE
45         DER=-0.25D0*X*X2
           R=DER
           DO 50 M=1,60
              R=-0.25D0*R/M/(M+1.0D0)/(2.0D0*M+1.0D0)**2*X4
              DER=DER+R
              IF (DABS(R).LT.DABS(DER)*EPS) GO TO 55
50         CONTINUE
55         DEI=0.5D0*X
           R=DEI
           DO 60 M=1,60
              R=-0.25D0*R/(M*M)/(2.D0*M-1.D0)/(2.D0*M+1.D0)*X4
              DEI=DEI+R
              IF (DABS(R).LT.DABS(DEI)*EPS) GO TO 65
60            CONTINUE
65         R=-0.25D0*X*X2
           GS=1.5D0
           HER=1.5D0*R-BER/X-(DLOG(X/2.D0)+EL)*DER+0.25*PI*DEI
           DO 70 M=1,60
              R=-0.25D0*R/M/(M+1.0D0)/(2.0D0*M+1.0D0)**2*X4
              GS=GS+1.0D0/(2*M+1.0D0)+1.0D0/(2*M+2.0D0)
              HER=HER+R*GS
              IF (DABS(R*GS).LT.DABS(HER)*EPS) GO TO 75
70         CONTINUE
75         R=0.5D0*X
           GS=1.0D0
           HEI=0.5D0*X-BEI/X-(DLOG(X/2.D0)+EL)*DEI-0.25*PI*DER
           DO 80 M=1,60
              R=-0.25D0*R/(M*M)/(2*M-1.0D0)/(2*M+1.0D0)*X4
              GS=GS+1.0D0/(2.0D0*M)+1.0D0/(2*M+1.0D0)
              HEI=HEI+R*GS
              IF (DABS(R*GS).LT.DABS(HEI)*EPS) RETURN
80         CONTINUE
        ELSE
           PP0=1.0D0
           PN0=1.0D0
           QP0=0.0D0
           QN0=0.0D0
           R0=1.0D0
           KM=18
           IF (DABS(X).GE.40.0) KM=10
           FAC=1.0D0
           DO 85 K=1,KM
              FAC=-FAC
              XT=0.25D0*K*PI-INT(0.125D0*K)*2.0D0*PI
              CS=COS(XT)
              SS=SIN(XT)
              R0=0.125D0*R0*(2.0D0*K-1.0D0)**2/K/X
              RC=R0*CS
              RS=R0*SS
              PP0=PP0+RC
              PN0=PN0+FAC*RC
              QP0=QP0+RS
85            QN0=QN0+FAC*RS
           XD=X/DSQRT(2.0D0)
           XE1=DEXP(XD)
           XE2=DEXP(-XD)
           XC1=1.D0/DSQRT(2.0D0*PI*X)
           XC2=DSQRT(.5D0*PI/X)
           CP0=DCOS(XD+0.125D0*PI)
           CN0=DCOS(XD-0.125D0*PI)
           SP0=DSIN(XD+0.125D0*PI)
           SN0=DSIN(XD-0.125D0*PI)
           GER=XC2*XE2*(PN0*CP0-QN0*SP0)
           GEI=XC2*XE2*(-PN0*SP0-QN0*CP0)
           BER=XC1*XE1*(PP0*CN0+QP0*SN0)-GEI/PI
           BEI=XC1*XE1*(PP0*SN0-QP0*CN0)+GER/PI
           PP1=1.0D0
           PN1=1.0D0
           QP1=0.0D0
           QN1=0.0D0
           R1=1.0D0
           FAC=1.0D0
           DO 90 K=1,KM
              FAC=-FAC
              XT=0.25D0*K*PI-INT(0.125D0*K)*2.0D0*PI
              CS=DCOS(XT)
              SS=DSIN(XT)
              R1=0.125D0*R1*(4.D0-(2.0D0*K-1.0D0)**2)/K/X
              RC=R1*CS
              RS=R1*SS
              PP1=PP1+FAC*RC
              PN1=PN1+RC
              QP1=QP1+FAC*RS
              QN1=QN1+RS
90         CONTINUE
           HER=XC2*XE2*(-PN1*CN0+QN1*SN0)
           HEI=XC2*XE2*(PN1*SN0+QN1*CN0)
           DER=XC1*XE1*(PP1*CP0+QP1*SP0)-HEI/PI
           DEI=XC1*XE1*(PP1*SP0-QP1*CP0)+HER/PI
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE CHGUBI(A,B,X,HU,ID)
C
C       ======================================================
C       Purpose: Compute confluent hypergeometric function
C                U(a,b,x) with integer b ( b = ±1,±2,... )
C       Input  : a  --- Parameter
C                b  --- Parameter
C                x  --- Argument
C       Output:  HU --- U(a,b,x)
C                ID --- Estimated number of significant digits
C       Routines called:
C            (1) GAMMA2 for computing gamma function Г(x)
C            (2) PSI_SPEC for computing psi function
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        ID=-100
        EL=0.5772156649015329D0
        N=ABS(B-1)
        RN1=1.0D0
        RN=1.0D0
        DO 10 J=1,N
           RN=RN*J
           IF (J.EQ.N-1) RN1=RN
10      CONTINUE
        CALL PSI_SPEC(A,PS)
        CALL GAMMA2(A,GA)
        IF (B.GT.0.0) THEN
           A0=A
           A1=A-N
           A2=A1
           CALL GAMMA2(A1,GA1)
           UA=(-1)**(N-1)/(RN*GA1)
           UB=RN1/GA*X**(-N)
        ELSE
           A0=A+N
           A1=A0
           A2=A
           CALL GAMMA2(A1,GA1)
           UA=(-1)**(N-1)/(RN*GA)*X**N
           UB=RN1/GA1
        ENDIF
        HM1=1.0D0
        R=1.0D0
        HMAX=0.0D0
        HMIN=1.0D+300
        H0=0D0
        DO 15 K=1,150
           R=R*(A0+K-1.0D0)*X/((N+K)*K)
           HM1=HM1+R
           HU1=DABS(HM1)
           IF (HU1.GT.HMAX) HMAX=HU1
           IF (HU1.LT.HMIN) HMIN=HU1
           IF (DABS(HM1-H0).LT.DABS(HM1)*1.0D-15) GO TO 20
15         H0=HM1
20      DA1=LOG10(HMAX)
        DA2=0.0D0
        IF (HMIN.NE.0.0) DA2=LOG10(HMIN)
        ID=15-ABS(DA1-DA2)
        HM1=HM1*DLOG(X)
        S0=0.0D0
        DO 25 M=1,N
           IF (B.GE.0.0) S0=S0-1.0D0/M
25         IF (B.LT.0.0) S0=S0+(1.0D0-A)/(M*(A+M-1.0D0))
        HM2=PS+2.0D0*EL+S0
        R=1.0D0
        HMAX=0.0D0
        HMIN=1.0D+300
        DO 50 K=1,150
           S1=0.0D0
           S2=0.0D0
           IF (B.GT.0.0) THEN
              DO 30 M=1,K
30               S1=S1-(M+2.0D0*A-2.0D0)/(M*(M+A-1.0D0))
              DO 35 M=1,N
35               S2=S2+1.0D0/(K+M)
           ELSE
              DO 40 M=1,K+N
40               S1=S1+(1.0D0-A)/(M*(M+A-1.0D0))
              DO 45 M=1,K
45               S2=S2+1.0D0/M
           ENDIF
           HW=2.0D0*EL+PS+S1-S2
           R=R*(A0+K-1.0D0)*X/((N+K)*K)
           HM2=HM2+R*HW
           HU2=DABS(HM2)
           IF (HU2.GT.HMAX) HMAX=HU2
           IF (HU2.LT.HMIN) HMIN=HU2
           IF (DABS((HM2-H0)/HM2).LT.1.0D-15) GO TO 55
50         H0=HM2
55      DB1=LOG10(HMAX)
        DB2=0.0D0
        IF (HMIN.NE.0.0) DB2=LOG10(HMIN)
        ID1=15-ABS(DB1-DB2)
        IF (ID1.LT.ID) ID=ID1
        HM3=1.0D0
        IF (N.EQ.0) HM3=0.0D0
        R=1.0D0
        DO 60 K=1,N-1
           R=R*(A2+K-1.0D0)/((K-N)*K)*X
60         HM3=HM3+R
        SA=UA*(HM1+HM2)
        SB=UB*HM3
        HU=SA+SB
        ID2=0.0D0
        IF (SA.NE.0.0) ID1=INT(LOG10(ABS(SA)))
        IF (HU.NE.0.0) ID2=INT(LOG10(ABS(HU)))
        IF (SA*SB.LT.0.0) ID=ID-ABS(ID1-ID2)
        RETURN
        END



C       **********************************

        SUBROUTINE CYZO(NT,KF,KC,ZO,ZV)
C
C       ===========================================================
C       Purpose : Compute the complex zeros of Y0(z), Y1(z) and
C                 Y1'(z), and their associated values at the zeros
C                 using the modified Newton's iteration method
C       Input:    NT --- Total number of zeros/roots
C                 KF --- Function choice code
C                        KF=0 for  Y0(z) & Y1(z0)
C                        KF=1 for  Y1(z) & Y0(z1)
C                        KF=2 for  Y1'(z) & Y1(z1')
C                 KC --- Choice code
C                        KC=0 for complex roots
C                        KC=1 for real roots
C       Output:   ZO(L) --- L-th zero of Y0(z) or Y1(z) or Y1'(z)
C                 ZV(L) --- Value of Y0'(z) or Y1'(z) or Y1(z)
C                           at the L-th zero
C       Routine called: CY01 for computing Y0(z) and Y1(z), and
C                       their derivatives
C       ===========================================================
        IMPLICIT DOUBLE PRECISION (H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION ZO(NT),ZV(NT)
        X=0.0D0
        Y=0.0D0
        H=0.0D0
        IF (KC.EQ.0) THEN
           X=-2.4D0
           Y=0.54D0
           H=3.14D0
        ELSE IF (KC.EQ.1) THEN
           X=0.89
           Y=0.0
           H=-3.14
        ENDIF
        IF (KF.EQ.1) X=-0.503
        IF (KF.EQ.2) X=0.577
        ZERO = DCMPLX(X, Y)
        Z=ZERO
        W=0.0D0
        DO 35 NR=1,NT
           IF (NR.NE.1) Z=ZO(NR-1)-H
           IT=0
15         IT=IT+1
           CALL CY01(KF,Z,ZF,ZD)
           ZP=(1.0D0,0.0D0)
           DO 20 I=1,NR-1
20            ZP=ZP*(Z-ZO(I))
           ZFD=ZF/ZP
           ZQ=(0.0D0,0.0D0)
           DO 30 I=1,NR-1
              ZW=(1.0D0,0.0D0)
              DO 25 J=1,NR-1
                 IF (J.EQ.I) GO TO 25
                 ZW=ZW*(Z-ZO(J))
25            CONTINUE
              ZQ=ZQ+ZW
30         CONTINUE
           ZGD=(ZD-ZQ*ZFD)/ZP
           Z=Z-ZFD/ZGD
           W0=W
           W=CDABS(Z)
           IF (IT.LE.50.AND.DABS((W-W0)/W).GT.1.0D-12) GO TO 15
           ZO(NR)=Z
35      CONTINUE
        DO 40 I=1,NT
           Z=ZO(I)
           IF (KF.EQ.0.OR.KF.EQ.2) THEN
              CALL CY01(1,Z,ZF,ZD)
              ZV(I)=ZF
           ELSE IF (KF.EQ.1) THEN
              CALL CY01(0,Z,ZF,ZD)
              ZV(I)=ZF
           ENDIF
40      CONTINUE
        RETURN
        END



C       **********************************

        SUBROUTINE KLVNB(X,BER,BEI,GER,GEI,DER,DEI,HER,HEI)
C
C       ======================================================
C       Purpose: Compute Kelvin functions ber x, bei x, ker x
C                and kei x, and their derivatives  ( x > 0 )
C       Input :  x   --- Argument of Kelvin functions
C       Output:  BER --- ber x
C                BEI --- bei x
C                GER --- ker x
C                GEI --- kei x
C                DER --- ber'x
C                DEI --- bei'x
C                HER --- ker'x
C                HEI --- kei'x
C       ================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           BER=1.0D0
           BEI=0.0D0
           GER=1.0D+300
           GEI=-.25D0*PI
           DER=0.0D0
           DEI=0.0D0
           HER=-1.0D+300
           HEI=0.0D0
        ELSE IF (X.LT.8.0D0) THEN
           T=X/8.0D0
           T2=T*T
           U=T2*T2
           BER=((((((-.901D-5*U+.122552D-2)*U-.08349609D0)*U
     &         +2.64191397D0)*U-32.36345652D0)*U
     &         +113.77777774D0)*U-64.0D0)*U+1.0D0
           BEI=T*T*((((((.11346D-3*U-.01103667D0)*U
     &         +.52185615D0)*U-10.56765779D0)*U
     &         +72.81777742D0)*U-113.77777774D0)*U+16.0D0)
           GER=((((((-.2458D-4*U+.309699D-2)*U-.19636347D0)
     &         *U+5.65539121D0)*U-60.60977451D0)*U+
     &         171.36272133D0)*U-59.05819744D0)*U-.57721566D0
           GER=GER-DLOG(.5D0*X)*BER+.25D0*PI*BEI
           GEI=T2*((((((.29532D-3*U-.02695875D0)*U
     &         +1.17509064D0)*U-21.30060904D0)*U
     &         +124.2356965D0)*U-142.91827687D0)*U
     &         +6.76454936D0)
           GEI=GEI-DLOG(.5D0*X)*BEI-.25D0*PI*BER
           DER=X*T2*((((((-.394D-5*U+.45957D-3)*U
     &         -.02609253D0)*U+.66047849D0)*U-6.0681481D0)*U
     &         +14.22222222D0)*U-4.0D0)
           DEI=X*((((((.4609D-4*U-.379386D-2)*U+.14677204D0)
     &         *U-2.31167514D0)*U+11.37777772D0)*U
     &         -10.66666666D0)*U+.5D0)
           HER=X*T2*((((((-.1075D-4*U+.116137D-2)*U
     &         -.06136358D0)*U+1.4138478D0)*U-11.36433272D0)
     &         *U+21.42034017D0)*U-3.69113734D0)
           HER=HER-DLOG(.5D0*X)*DER-BER/X+.25D0*PI*DEI
           HEI=X*((((((.11997D-3*U-.926707D-2)*U
     &         +.33049424D0)*U-4.65950823D0)*U+19.41182758D0)
     &         *U-13.39858846D0)*U+.21139217D0)
           HEI=HEI-DLOG(.5D0*X)*DEI-BEI/X-.25D0*PI*DER
        ELSE
           T=8.0D0/X
           TNR=0.0D0
           TNI=0.0D0
           DO 10 L=1,2
              V=(-1)**L*T
              TPR=((((.6D-6*V-.34D-5)*V-.252D-4)*V-.906D-4)
     &            *V*V+.0110486D0)*V
              TPI=((((.19D-5*V+.51D-5)*V*V-.901D-4)*V
     &            -.9765D-3)*V-.0110485D0)*V-.3926991D0
              IF (L.EQ.1) THEN
                 TNR=TPR
                 TNI=TPI
              ENDIF
10         CONTINUE
           YD=X/DSQRT(2.0D0)
           YE1=DEXP(YD+TPR)
           YE2=DEXP(-YD+TNR)
           YC1=1.0D0/DSQRT(2.0D0*PI*X)
           YC2=DSQRT(PI/(2.0D0*X))
           CSP=DCOS(YD+TPI)
           SSP=DSIN(YD+TPI)
           CSN=DCOS(-YD+TNI)
           SSN=DSIN(-YD+TNI)
           GER=YC2*YE2*CSN
           GEI=YC2*YE2*SSN
           FXR=YC1*YE1*CSP
           FXI=YC1*YE1*SSP
           BER=FXR-GEI/PI
           BEI=FXI+GER/PI
           PNR=0.0D0
           PNI=0.0D0
           DO 15 L=1,2
              V=(-1)**L*T
              PPR=(((((.16D-5*V+.117D-4)*V+.346D-4)*V+.5D-6)
     &            *V-.13813D-2)*V-.0625001D0)*V+.7071068D0
              PPI=(((((-.32D-5*V-.24D-5)*V+.338D-4)*V+
     &           .2452D-3)*V+.13811D-2)*V-.1D-6)*V+.7071068D0
              IF (L.EQ.1) THEN
                 PNR=PPR
                 PNI=PPI
              ENDIF
15         CONTINUE
           HER=GEI*PNI-GER*PNR
           HEI=-(GEI*PNR+GER*PNI)
           DER=FXR*PPR-FXI*PPI-HEI/PI
           DEI=FXI*PPR+FXR*PPI+HER/PI
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE RMN2SO(M,N,C,X,CV,DF,KD,R2F,R2D)
C
C       =============================================================
C       Purpose: Compute oblate radial functions of the second kind
C                with a small argument, Rmn(-ic,ix) & Rmn'(-ic,ix)
C       Routines called:
C            (1) SCKB for computing the expansion coefficients c2k
C            (2) KMN for computing the joining factors
C            (3) QSTAR for computing the factor defined in (15.7.3)
C            (4) CBK for computing the the expansion coefficient
C                defined in (15.7.6)
C            (5) GMN for computing the function defined in (15.7.4)
C            (6) RMN1 for computing the radial function of the first
C                kind
C       =============================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BK(200),CK(200),DF(200),DN(200)
        IF (DABS(DF(1)).LE.1.0D-280) THEN
           R2F=1.0D+300
           R2D=1.0D+300
           RETURN
        ENDIF
        EPS=1.0D-14
        PI=3.141592653589793D0
        NM=25+INT((N-M)/2+C)
        IP=1
        IF (N-M.EQ.2*INT((N-M)/2)) IP=0
        CALL SCKB(M,N,C,DF,CK)
        CALL KMN(M,N,C,CV,KD,DF,DN,CK1,CK2)
        CALL QSTAR(M,N,C,CK,CK1,QS,QT)
        CALL CBK(M,N,C,CV,QT,CK,BK)
        IF (X.EQ.0.0D0) THEN
           SUM=0.0D0
           SW=0.0D0
           DO 10 J=1,NM
              SUM=SUM+CK(J)
              IF (DABS(SUM-SW).LT.DABS(SUM)*EPS) GO TO 15
10            SW=SUM
15         IF (IP.EQ.0) THEN
              R1F=SUM/CK1
              R2F=-0.5D0*PI*QS*R1F
              R2D=QS*R1F+BK(1)
           ELSE IF (IP.EQ.1) THEN
              R1D=SUM/CK1
              R2F=BK(1)
              R2D=-0.5D0*PI*QS*R1D
           ENDIF
           RETURN
        ELSE
           CALL GMN(M,N,C,X,BK,GF,GD)
           CALL RMN1(M,N,C,X,DF,KD,R1F,R1D)
           H0=DATAN(X)-0.5D0*PI
           R2F=QS*R1F*H0+GF
           R2D=QS*(R1D*H0+R1F/(1.0D0+X*X))+GD
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE CSPHIK(N,Z,NM,CSI,CDI,CSK,CDK)
C
C       =======================================================
C       Purpose: Compute modified spherical Bessel functions
C                and their derivatives for a complex argument
C       Input :  z --- Complex argument
C                n --- Order of in(z) & kn(z) ( n = 0,1,2,... )
C       Output:  CSI(n) --- in(z)
C                CDI(n) --- in'(z)
C                CSK(n) --- kn(z)
C                CDK(n) --- kn'(z)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       =======================================================
C
        IMPLICIT COMPLEX*16 (C,Z)
        DOUBLE PRECISION A0,PI
        DIMENSION CSI(0:N),CDI(0:N),CSK(0:N),CDK(0:N)
        PI=3.141592653589793D0
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-60) THEN
           DO 10 K=0,N
              CSI(K)=0.0D0
              CDI(K)=0.0D0
              CSK(K)=1.0D+300
10            CDK(K)=-1.0D+300
           CSI(0)=1.0D0
           CDI(1)=0.3333333333333333D0
           RETURN
        ENDIF
        CI = DCMPLX(0.0D0, 1.0D0)
        CSINH=CDSIN(CI*Z)/CI
        CCOSH=CDCOS(CI*Z)
        CSI0=CSINH/Z
        CSI1=(-CSINH/Z+CCOSH)/Z
        CSI(0)=CSI0
        CSI(1)=CSI1
        IF (N.GE.2) THEN
           M=MSTA1(A0,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(A0,N,15)
           ENDIF
           CF0=0.0D0
           CF1=1.0D0-100
           DO 15 K=M,0,-1
              CF=(2.0D0*K+3.0D0)*CF1/Z+CF0
              IF (K.LE.NM) CSI(K)=CF
              CF0=CF1
15            CF1=CF
           IF (CDABS(CSI0).GT.CDABS(CSI1)) CS=CSI0/CF
           IF (CDABS(CSI0).LE.CDABS(CSI1)) CS=CSI1/CF0
           DO 20 K=0,NM
20            CSI(K)=CS*CSI(K)
        ENDIF
        CDI(0)=CSI(1)
        DO 25 K=1,NM
25         CDI(K)=CSI(K-1)-(K+1.0D0)*CSI(K)/Z
        CSK(0)=0.5D0*PI/Z*CDEXP(-Z)
        CSK(1)=CSK(0)*(1.0D0+1.0D0/Z)
        DO 30 K=2,NM
           IF (CDABS(CSI(K-1)).GT.CDABS(CSI(K-2))) THEN
              CSK(K)=(0.5D0*PI/(Z*Z)-CSI(K)*CSK(K-1))/CSI(K-1)
           ELSE
              CSK(K)=(CSI(K)*CSK(K-2)+(K-0.5D0)*PI/Z**3)/CSI(K-2)
           ENDIF
30      CONTINUE
        CDK(0)=-CSK(1)
        DO 35 K=1,NM
35         CDK(K)=-CSK(K-1)-(K+1.0D0)*CSK(K)/Z
        RETURN
        END



C       **********************************

        SUBROUTINE BJNDD(N,X,BJ,DJ,FJ)
C
C       =====================================================
C       Purpose: Compute Bessel functions Jn(x) and their
C                first and second derivatives ( n= 0,1,… )
C       Input:   x ---  Argument of Jn(x)  ( x ≥ 0 )
C                n ---  Order of Jn(x)
C       Output:  BJ(n+1) ---  Jn(x)
C                DJ(n+1) ---  Jn'(x)
C                FJ(n+1) ---  Jn"(x)
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BJ(101),DJ(101),FJ(101)
        DO 10 NT=1,900
           MT=INT(0.5*LOG10(6.28*NT)-NT*LOG10(1.36*DABS(X)/NT))
           IF (MT.GT.20) GO TO 15
10      CONTINUE
15      M=NT
        BS=0.0D0
        F=0.0D0
        F0=0.0D0
        F1=1.0D-35
        DO 20 K=M,0,-1
           F=2.0D0*(K+1.0D0)*F1/X-F0
           IF (K.LE.N) BJ(K+1)=F
           IF (K.EQ.2*INT(K/2)) BS=BS+2.0D0*F
           F0=F1
20         F1=F
        DO 25 K=0,N
25         BJ(K+1)=BJ(K+1)/(BS-F)
        DJ(1)=-BJ(2)
        FJ(1)=-1.0D0*BJ(1)-DJ(1)/X
        DO 30 K=1,N
           DJ(K+1)=BJ(K)-K*BJ(K+1)/X
30         FJ(K+1)=(K*K/(X*X)-1.0D0)*BJ(K+1)-DJ(K+1)/X
        RETURN
        END

C       **********************************


        SUBROUTINE SPHJ(N,X,NM,SJ,DJ)
C       MODIFIED to ALLOW N=0 CASE (ALSO IN CSPHJY, SPHY)
C
C       =======================================================
C       Purpose: Compute spherical Bessel functions jn(x) and
C                their derivatives
C       Input :  x --- Argument of jn(x)
C                n --- Order of jn(x)  ( n = 0,1,… )
C       Output:  SJ(n) --- jn(x)
C                DJ(n) --- jn'(x)
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SJ(0:N),DJ(0:N)
        NM=N
        IF (DABS(X).LT.1.0D-100) THEN
           DO 10 K=0,N
              SJ(K)=0.0D0
10            DJ(K)=0.0D0
           SJ(0)=1.0D0
           IF (N.GT.0) THEN
              DJ(1)=.3333333333333333D0
           ENDIF
           RETURN
        ENDIF
        SJ(0)=DSIN(X)/X
        DJ(0)=(DCOS(X)-DSIN(X)/X)/X
        IF (N.LT.1) THEN
           RETURN
        ENDIF
        SJ(1)=(SJ(0)-DCOS(X))/X
        IF (N.GE.2) THEN
           SA=SJ(0)
           SB=SJ(1)
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F=0.0D0
           F0=0.0D0
           F1=1.0D0-100
           DO 15 K=M,0,-1
              F=(2.0D0*K+3.0D0)*F1/X-F0
              IF (K.LE.NM) SJ(K)=F
              F0=F1
15            F1=F
           CS=0.0D0
           IF (DABS(SA).GT.DABS(SB)) CS=SA/F
           IF (DABS(SA).LE.DABS(SB)) CS=SB/F0
           DO 20 K=0,NM
20            SJ(K)=CS*SJ(K)
        ENDIF
        DO 25 K=1,NM
25         DJ(K)=SJ(K-1)-(K+1.0D0)*SJ(K)/X
        RETURN
        END



C       **********************************

        SUBROUTINE OTHPL(KF,N,X,PL,DPL)
C
C       ==========================================================
C       Purpose: Compute orthogonal polynomials: Tn(x) or Un(x),
C                or Ln(x) or Hn(x), and their derivatives
C       Input :  KF --- Function code
C                       KF=1 for Chebyshev polynomial Tn(x)
C                       KF=2 for Chebyshev polynomial Un(x)
C                       KF=3 for Laguerre polynomial Ln(x)
C                       KF=4 for Hermite polynomial Hn(x)
C                n ---  Order of orthogonal polynomials
C                x ---  Argument of orthogonal polynomials
C       Output:  PL(n) --- Tn(x) or Un(x) or Ln(x) or Hn(x)
C                DPL(n)--- Tn'(x) or Un'(x) or Ln'(x) or Hn'(x)
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION PL(0:N),DPL(0:N)
        A=2.0D0
        B=0.0D0
        C=1.0D0
        Y0=1.0D0
        Y1=2.0D0*X
        DY0=0.0D0
        DY1=2.0D0
        PL(0)=1.0D0
        PL(1)=2.0D0*X
        DPL(0)=0.0D0
        DPL(1)=2.0D0
        IF (KF.EQ.1) THEN
           Y1=X
           DY1=1.0D0
           PL(1)=X
           DPL(1)=1.0D0
        ELSE IF (KF.EQ.3) THEN
           Y1=1.0D0-X
           DY1=-1.0D0
           PL(1)=1.0D0-X
           DPL(1)=-1.0D0
        ENDIF
        DO 10 K=2,N
           IF (KF.EQ.3) THEN
              A=-1.0D0/K
              B=2.0D0+A
              C=1.0D0+A
           ELSE IF (KF.EQ.4) THEN
              C=2.0D0*(K-1.0D0)
           ENDIF
           YN=(A*X+B)*Y1-C*Y0
           DYN=A*Y1+(A*X+B)*DY1-C*DY0
           PL(K)=YN
           DPL(K)=DYN
           Y0=Y1
           Y1=YN
           DY0=DY1
10         DY1=DYN
        RETURN
        END

C       **********************************

        SUBROUTINE KLVNZO(NT,KD,ZO)
C
C       ====================================================
C       Purpose: Compute the zeros of Kelvin functions
C       Input :  NT  --- Total number of zeros
C                KD  --- Function code
C                KD=1 to 8 for ber x, bei x, ker x, kei x,
C                          ber'x, bei'x, ker'x and kei'x,
C                          respectively.
C       Output:  ZO(M) --- the M-th zero of Kelvin function
C                          for code KD
C       Routine called:
C                KLVNA for computing Kelvin functions and
C                their derivatives
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION ZO(NT),RT0(8)
        RT0(1)=2.84891
        RT0(2)=5.02622
        RT0(3)=1.71854
        RT0(4)=3.91467
        RT0(5)=6.03871
        RT0(6)=3.77268
        RT0(7)=2.66584
        RT0(8)=4.93181
        RT=RT0(KD)
        DO 15 M=1,NT
10         CALL KLVNA(RT,BER,BEI,GER,GEI,DER,DEI,HER,HEI)
           IF (KD.EQ.1) THEN
              RT=RT-BER/DER
           ELSE IF (KD.EQ.2) THEN
              RT=RT-BEI/DEI
           ELSE IF (KD.EQ.3) THEN
              RT=RT-GER/HER
           ELSE IF (KD.EQ.4) THEN
              RT=RT-GEI/HEI
           ELSE IF (KD.EQ.5) THEN
              DDR=-BEI-DER/RT
              RT=RT-DER/DDR
           ELSE IF (KD.EQ.6) THEN
              DDI=BER-DEI/RT
              RT=RT-DEI/DDI
           ELSE IF (KD.EQ.7) THEN
              GDR=-GEI-HER/RT
              RT=RT-HER/GDR
           ELSE
              GDI=GER-HEI/RT
              RT=RT-HEI/GDI
           ENDIF
           IF (DABS(RT-RT0(KD)).GT.5.0D-10) THEN
              RT0(KD)=RT
              GO TO 10
           ENDIF
           ZO(M)=RT
15         RT=RT+4.44D0
        RETURN
        END



C       **********************************

        SUBROUTINE RSWFO(M,N,C,X,CV,KF,R1F,R1D,R2F,R2D)
C
C       ==========================================================
C       Purpose: Compute oblate radial functions of the first
C                and second kinds, and their derivatives
C       Input :  m  --- Mode parameter,  m = 0,1,2,...
C                n  --- Mode parameter,  n = m,m+1,m+2,...
C                c  --- Spheroidal parameter
C                x  --- Argument (x ≥ 0)
C                cv --- Characteristic value
C                KF --- Function code
C                       KF=1 for the first kind
C                       KF=2 for the second kind
C                       KF=3 for both the first and second kinds
C       Output:  R1F --- Radial function of the first kind
C                R1D --- Derivative of the radial function of
C                        the first kind
C                R2F --- Radial function of the second kind
C                R2D --- Derivative of the radial function of
C                        the second kind
C       Routines called:
C            (1) SDMN for computing expansion coefficients dk
C            (2) RMN1 for computing prolate or oblate radial
C                function of the first kind
C            (3) RMN2L for computing prolate or oblate radial
C                function of the second kind for a large argument
C            (4) RMN2SO for computing oblate radial functions of
C                the second kind for a small argument
C       ==========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION DF(200)
        KD=-1
        CALL SDMN(M,N,C,CV,KD,DF)
        IF (KF.NE.2) THEN
           CALL RMN1(M,N,C,X,DF,KD,R1F,R1D)
        ENDIF
        IF (KF.GT.1) THEN
           ID=10
           IF (X.GT.1.0D-8) THEN
              CALL RMN2L(M,N,C,X,DF,KD,R2F,R2D,ID)
           ENDIF
           IF (ID.GT.-1) THEN
              CALL RMN2SO(M,N,C,X,CV,DF,KD,R2F,R2D)
           ENDIF
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE CH12N(N,Z,NM,CHF1,CHD1,CHF2,CHD2)
C
C       ====================================================
C       Purpose: Compute Hankel functions of the first and
C                second kinds and their derivatives for a
C                complex argument
C       Input :  z --- Complex argument
C                n --- Order of Hn(1)(z) and Hn(2)(z)
C       Output:  CHF1(n) --- Hn(1)(z)
C                CHD1(n) --- Hn(1)'(z)
C                CHF2(n) --- Hn(2)(z)
C                CHD2(n) --- Hn(2)'(z)
C                NM --- Highest order computed
C       Routines called:
C             (1) CJYNB for computing Jn(z) and Yn(z)
C             (2) CIKNB for computing In(z) and Kn(z)
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBJ(0:250),CDJ(0:250),CBY(0:250),CDY(0:250),
     &            CBI(0:250),CDI(0:250),CBK(0:250),CDK(0:250)
        DIMENSION CHF1(0:N),CHD1(0:N),CHF2(0:N),CHD2(0:N)
        CI=(0.0D0,1.0D0)
        PI=3.141592653589793D0
        IF (DIMAG(Z).LT.0.0D0) THEN
           CALL CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)
           DO 10 K=0,NM
              CHF1(K)=CBJ(K)+CI*CBY(K)
10            CHD1(K)=CDJ(K)+CI*CDY(K)
           ZI=CI*Z
           CALL CIKNB(N,ZI,NM,CBI,CDI,CBK,CDK)
           CFAC=-2.0D0/(PI*CI)
           DO 15 K=0,NM
              CHF2(K)=CFAC*CBK(K)
              CHD2(K)=CFAC*CI*CDK(K)
15            CFAC=CFAC*CI
        ELSE IF (DIMAG(Z).GT.0.0D0) THEN
           ZI=-CI*Z
           CALL CIKNB(N,ZI,NM,CBI,CDI,CBK,CDK)
           CF1=-CI
           CFAC=2.0D0/(PI*CI)
           DO 20 K=0,NM
              CHF1(K)=CFAC*CBK(K)
              CHD1(K)=-CFAC*CI*CDK(K)
20            CFAC=CFAC*CF1
           CALL CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)
           DO 25 K=0,NM
              CHF2(K)=CBJ(K)-CI*CBY(K)
25            CHD2(K)=CDJ(K)-CI*CDY(K)
        ELSE
           CALL CJYNB(N,Z,NM,CBJ,CDJ,CBY,CDY)
           DO 30 K=0,NM
              CHF1(K)=CBJ(K)+CI*CBY(K)
              CHD1(K)=CDJ(K)+CI*CDY(K)
              CHF2(K)=CBJ(K)-CI*CBY(K)
30            CHD2(K)=CDJ(K)-CI*CDY(K)
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE JYZO(N,NT,RJ0,RJ1,RY0,RY1)
C
C       ======================================================
C       Purpose: Compute the zeros of Bessel functions Jn(x),
C                Yn(x), and their derivatives
C       Input :  n  --- Order of Bessel functions  (n >= 0)
C                NT --- Number of zeros (roots)
C       Output:  RJ0(L) --- L-th zero of Jn(x),  L=1,2,...,NT
C                RJ1(L) --- L-th zero of Jn'(x), L=1,2,...,NT
C                RY0(L) --- L-th zero of Yn(x),  L=1,2,...,NT
C                RY1(L) --- L-th zero of Yn'(x), L=1,2,...,NT
C       Routine called: JYNDD for computing Jn(x), Yn(x), and
C                       their first and second derivatives
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION RJ0(NT),RJ1(NT),RY0(NT),RY1(NT)
        PI=3.141592653589793D0
C       -- Newton method for j_{N,L}
C       1) initial guess for j_{N,1}
        IF (N.LE.20) THEN
           X=2.82141+1.15859*N
        ELSE
C          Abr & Stg (9.5.14)
           X=N+1.85576*N**0.33333+1.03315/N**0.33333
        ENDIF
        L=0
C       2) iterate
        XGUESS=X
10      X0=X
        CALL JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
        X=X-BJN/DJN
        IF (X-X0.LT.-1) X=X0-1
        IF (X-X0.GT.1) X=X0+1
        IF (DABS(X-X0).GT.1.0D-11) GO TO 10
C       3) initial guess for j_{N,L+1}
        IF (L.GE.1 .AND. X.LE.RJ0(L)+0.5) THEN
           X=XGUESS+PI
           XGUESS=X
           GO TO 10
        END IF
        L=L+1
        RJ0(L)=X
C       XXX: should have a better initial guess for large N ~> 100 here
        X=X+PI+MAX((0.0972d0+0.0679*N-0.000354*N**2)/L, 0d0)
        IF (L.LT.NT) GO TO 10
C       -- Newton method for j_{N,L}'
        IF (N.LE.20) THEN
           X=0.961587+1.07703*N
        ELSE
           X=N+0.80861*N**0.33333+0.07249/N**0.33333
        ENDIF
        IF (N.EQ.0) X=3.8317
        L=0
        XGUESS=X
15      X0=X
        CALL JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
        X=X-DJN/FJN
        IF (X-X0.LT.-1) X=X0-1
        IF (X-X0.GT.1) X=X0+1
        IF (DABS(X-X0).GT.1.0D-11) GO TO 15
        IF (L.GE.1 .AND. X.LE.RJ1(L)+0.5) THEN
           X=XGUESS+PI
           XGUESS=X
           GO TO 15
        END IF
        L=L+1
        RJ1(L)=X
C       XXX: should have a better initial guess for large N ~> 100 here
        X=X+PI+MAX((0.4955d0+0.0915*N-0.000435*N**2)/L, 0d0)
        IF (L.LT.NT) GO TO 15
C       -- Newton method for y_{N,L}
        IF (N.LE.20) THEN
           X=1.19477+1.08933*N
        ELSE
           X=N+0.93158*N**0.33333+0.26035/N**0.33333
        ENDIF
        L=0
        XGUESS=X
20      X0=X
        CALL JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
        X=X-BYN/DYN
        IF (X-X0.LT.-1) X=X0-1
        IF (X-X0.GT.1) X=X0+1
        IF (DABS(X-X0).GT.1.0D-11) GO TO 20
        IF (L.GE.1 .AND. X.LE.RY0(L)+0.5) THEN
           X=XGUESS+PI
           XGUESS=X
           GO TO 20
        END IF
        L=L+1
        RY0(L)=X
C       XXX: should have a better initial guess for large N ~> 100 here
        X=X+PI+MAX((0.312d0+0.0852*N-0.000403*N**2)/L,0d0)
        IF (L.LT.NT) GO TO 20
C       -- Newton method for y_{N,L}'
        IF (N.LE.20) THEN
           X=2.67257+1.16099*N
        ELSE
           X=N+1.8211*N**0.33333+0.94001/N**0.33333
        ENDIF
        L=0
        XGUESS=X
25      X0=X
        CALL JYNDD(N,X,BJN,DJN,FJN,BYN,DYN,FYN)
        X=X-DYN/FYN
        IF (DABS(X-X0).GT.1.0D-11) GO TO 25
        IF (L.GE.1 .AND. X.LE.RY1(L)+0.5) THEN
           X=XGUESS+PI
           XGUESS=X
           GO TO 25
        END IF
        L=L+1
        RY1(L)=X
C       XXX: should have a better initial guess for large N ~> 100 here
        X=X+PI+MAX((0.197d0+0.0643*N-0.000286*N**2)/L,0d0)
        IF (L.LT.NT) GO TO 25
        RETURN
        END



C       **********************************

        SUBROUTINE IKV(V,X,VM,BI,DI,BK,DK)
C
C       =======================================================
C       Purpose: Compute modified Bessel functions Iv(x) and
C                Kv(x), and their derivatives
C       Input :  x --- Argument ( x ≥ 0 )
C                v --- Order of Iv(x) and Kv(x)
C                      ( v = n+v0, n = 0,1,2,..., 0 ≤ v0 < 1 )
C       Output:  BI(n) --- In+v0(x)
C                DI(n) --- In+v0'(x)
C                BK(n) --- Kn+v0(x)
C                DK(n) --- Kn+v0'(x)
C                VM --- Highest order computed
C       Routines called:
C            (1) GAMMA2 for computing the gamma function
C            (2) MSTA1 and MSTA2 to compute the starting
C                point for backward recurrence
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION BI(0:*),DI(0:*),BK(0:*),DK(0:*)
        PI=3.141592653589793D0
        X2=X*X
        N=INT(V)
        V0=V-N
        IF (N.EQ.0) N=1
        IF (X.LT.1.0D-100) THEN
           DO 10 K=0,N
              BI(K)=0.0D0
              DI(K)=0.0D0
              BK(K)=-1.0D+300
10            DK(K)=1.0D+300
           IF (V.EQ.0.0) THEN
              BI(0)=1.0D0
              DI(1)=0.5D0
           ENDIF
           VM=V
           RETURN
        ENDIF
        PIV=PI*V0
        VT=4.0D0*V0*V0
        IF (V0.EQ.0.0D0) THEN
           A1=1.0D0
        ELSE
           V0P=1.0D0+V0
           CALL GAMMA2(V0P,GAP)
           A1=(0.5D0*X)**V0/GAP
        ENDIF
        K0=14
        IF (X.GE.35.0) K0=10
        IF (X.GE.50.0) K0=8
        IF (X.LE.18.0) THEN
           BI0=1.0D0
           R=1.0D0
           DO 15 K=1,30
              R=0.25D0*R*X2/(K*(K+V0))
              BI0=BI0+R
              IF (DABS(R/BI0).LT.1.0D-15) GO TO 20
15         CONTINUE
20         BI0=BI0*A1
        ELSE
           CA=DEXP(X)/DSQRT(2.0D0*PI*X)
           SUM=1.0D0
           R=1.0D0
           DO 25 K=1,K0
              R=-0.125D0*R*(VT-(2.0D0*K-1.0D0)**2.0)/(K*X)
25            SUM=SUM+R
           BI0=CA*SUM
        ENDIF
        M=MSTA1(X,200)
        IF (M.LT.N) THEN
           N=M
        ELSE
           M=MSTA2(X,N,15)
        ENDIF
        F=0.0D0
        F2=0.0D0
        F1=1.0D-100
        WW=0.0D0
        DO 30 K=M,0,-1
           F=2.0D0*(V0+K+1.0D0)/X*F1+F2
           IF (K.LE.N) BI(K)=F
           F2=F1
30         F1=F
        CS=BI0/F
        DO 35 K=0,N
35         BI(K)=CS*BI(K)
        DI(0)=V0/X*BI(0)+BI(1)
        DO 40 K=1,N
40         DI(K)=-(K+V0)/X*BI(K)+BI(K-1)
        IF (X.LE.9.0D0) THEN
           IF (V0.EQ.0.0D0) THEN
              CT=-DLOG(0.5D0*X)-0.5772156649015329D0
              CS=0.0D0
              W0=0.0D0
              R=1.0D0
              DO 45 K=1,50
                 W0=W0+1.0D0/K
                 R=0.25D0*R/(K*K)*X2
                 CS=CS+R*(W0+CT)
                 WA=DABS(CS)
                 IF (DABS((WA-WW)/WA).LT.1.0D-15) GO TO 50
45               WW=WA
50            BK0=CT+CS
           ELSE
              V0N=1.0D0-V0
              CALL GAMMA2(V0N,GAN)
              A2=1.0D0/(GAN*(0.5D0*X)**V0)
              A1=(0.5D0*X)**V0/GAP
              SUM=A2-A1
              R1=1.0D0
              R2=1.0D0
              DO 55 K=1,120
                 R1=0.25D0*R1*X2/(K*(K-V0))
                 R2=0.25D0*R2*X2/(K*(K+V0))
                 SUM=SUM+A2*R1-A1*R2
                 WA=DABS(SUM)
                 IF (DABS((WA-WW)/WA).LT.1.0D-15) GO TO 60
55               WW=WA
60            BK0=0.5D0*PI*SUM/DSIN(PIV)
           ENDIF
        ELSE
           CB=DEXP(-X)*DSQRT(0.5D0*PI/X)
           SUM=1.0D0
           R=1.0D0
           DO 65 K=1,K0
              R=0.125D0*R*(VT-(2.0*K-1.0)**2.0)/(K*X)
65            SUM=SUM+R
           BK0=CB*SUM
        ENDIF
        BK1=(1.0D0/X-BI(1)*BK0)/BI(0)
        BK(0)=BK0
        BK(1)=BK1
        DO 70 K=2,N
           BK2=2.0D0*(V0+K-1.0D0)/X*BK1+BK0
           BK(K)=BK2
           BK0=BK1
70         BK1=BK2
        DK(0)=V0/X*BK(0)-BK(1)
        DO 80 K=1,N
80         DK(K)=-(K+V0)/X*BK(K)-BK(K-1)
        VM=N+V0
        RETURN
        END



C       **********************************

        SUBROUTINE SDMN(M,N,C,CV,KD,DF)
C
C       =====================================================
C       Purpose: Compute the expansion coefficients of the
C                prolate and oblate spheroidal functions, dk
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                cv --- Characteristic value
C                KD --- Function code
C                       KD=1 for prolate; KD=-1 for oblate
C       Output:  DF(k) --- Expansion coefficients dk;
C                          DF(1), DF(2), ... correspond to
C                          d0, d2, ... for even n-m and d1,
C                          d3, ... for odd n-m
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(200),D(200),G(200),DF(200)
        NM=25+INT(0.5*(N-M)+C)
        IF (C.LT.1.0D-10) THEN
           DO 5 I=1,NM
5             DF(I)=0D0
           DF((N-M)/2+1)=1.0D0
           RETURN
        ENDIF
        CS=C*C*KD
        IP=1
        K=0
        IF (N-M.EQ.2*INT((N-M)/2)) IP=0
        DO 10 I=1,NM+2
           IF (IP.EQ.0) K=2*(I-1)
           IF (IP.EQ.1) K=2*I-1
           DK0=M+K
           DK1=M+K+1
           DK2=2*(M+K)
           D2K=2*M+K
           A(I)=(D2K+2.0)*(D2K+1.0)/((DK2+3.0)*(DK2+5.0))*CS
           D(I)=DK0*DK1+(2.0*DK0*DK1-2.0*M*M-1.0)/((DK2-1.0)
     &          *(DK2+3.0))*CS
           G(I)=K*(K-1.0)/((DK2-3.0)*(DK2-1.0))*CS
10      CONTINUE
        FS=1.0D0
        F1=0.0D0
        F0=1.0D-100
        KB=0
        DF(NM+1)=0.0D0
        FL=0.0D0
        DO 30 K=NM,1,-1
           F=-((D(K+1)-CV)*F0+A(K+1)*F1)/G(K+1)
           IF (DABS(F).GT.DABS(DF(K+1))) THEN
              DF(K)=F
              F1=F0
              F0=F
              IF (DABS(F).GT.1.0D+100) THEN
                 DO 12 K1=K,NM
12                  DF(K1)=DF(K1)*1.0D-100
                 F1=F1*1.0D-100
                 F0=F0*1.0D-100
              ENDIF
           ELSE
              KB=K
              FL=DF(K+1)
              F1=1.0D-100
              F2=-(D(1)-CV)/A(1)*F1
              DF(1)=F1
              IF (KB.EQ.1) THEN
                 FS=F2
              ELSE IF (KB.EQ.2) THEN
                 DF(2)=F2
                 FS=-((D(2)-CV)*F2+G(2)*F1)/A(2)
              ELSE
                 DF(2)=F2
                 DO 20 J=3,KB+1
                    F=-((D(J-1)-CV)*F2+G(J-1)*F1)/A(J-1)
                    IF (J.LE.KB) DF(J)=F
                    IF (DABS(F).GT.1.0D+100) THEN
                       DO 15 K1=1,J
15                        DF(K1)=DF(K1)*1.0D-100
                       F=F*1.0D-100
                       F2=F2*1.0D-100
                    ENDIF
                    F1=F2
20                  F2=F
                 FS=F
              ENDIF
              GO TO 35
           ENDIF
30      CONTINUE
35      SU1=0.0D0
        R1=1.0D0
        DO 40 J=M+IP+1,2*(M+IP)
40         R1=R1*J
        SU1=DF(1)*R1
        DO 45 K=2,KB
           R1=-R1*(K+M+IP-1.5D0)/(K-1.0D0)
45           SU1=SU1+R1*DF(K)
        SU2=0.0D0
        SW=0.0D0
        DO 50 K=KB+1,NM
           IF (K.NE.1) R1=-R1*(K+M+IP-1.5D0)/(K-1.0D0)
           SU2=SU2+R1*DF(K)
           IF (DABS(SW-SU2).LT.DABS(SU2)*1.0D-14) GOTO 55
50         SW=SU2
55      R3=1.0D0
        DO 60 J=1,(M+N+IP)/2
60         R3=R3*(J+0.5D0*(N+M+IP))
        R4=1.0D0
        DO 65 J=1,(N-M-IP)/2
65         R4=-4.0D0*R4*J
        S0=R3/(FL*(SU1/FS)+SU2)/R4
        DO 70 K=1,KB
70         DF(K)=FL/FS*S0*DF(K)
        DO 75 K=KB+1,NM
75         DF(K)=S0*DF(K)
        RETURN
        END




C       **********************************

        SUBROUTINE AJYIK(X,VJ1,VJ2,VY1,VY2,VI1,VI2,VK1,VK2)
C
C       =======================================================
C       Purpose: Compute Bessel functions Jv(x) and Yv(x),
C                and modified Bessel functions Iv(x) and
C                Kv(x), and their derivatives with v=1/3,2/3
C       Input :  x --- Argument of Jv(x),Yv(x),Iv(x) and
C                      Kv(x) ( x ≥ 0 )
C       Output:  VJ1 --- J1/3(x)
C                VJ2 --- J2/3(x)
C                VY1 --- Y1/3(x)
C                VY2 --- Y2/3(x)
C                VI1 --- I1/3(x)
C                VI2 --- I2/3(x)
C                VK1 --- K1/3(x)
C                VK2 --- K2/3(x)
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        IF (X.EQ.0.0D0) THEN
           VJ1=0.0D0
           VJ2=0.0D0
           VY1=-1.0D+300
           VY2=1.0D+300
           VI1=0.0D0
           VI2=0.0D0
           VK1=-1.0D+300
           VK2=-1.0D+300
           RETURN
        ENDIF
        PI=3.141592653589793D0
        RP2=.63661977236758D0
        GP1=.892979511569249D0
        GP2=.902745292950934D0
        GN1=1.3541179394264D0
        GN2=2.678938534707747D0
        VV0=0.444444444444444D0
        UU0=1.1547005383793D0
        X2=X*X
        K0=12
        IF (X.GE.35.0) K0=10
        IF (X.GE.50.0) K0=8
        IF (X.LE.12.0) THEN
           DO 25 L=1,2
              VL=L/3.0D0
              VJL=1.0D0
              R=1.0D0
              DO 15 K=1,40
                 R=-0.25D0*R*X2/(K*(K+VL))
                 VJL=VJL+R
                 IF (DABS(R).LT.1.0D-15) GO TO 20
15            CONTINUE
20            A0=(0.5D0*X)**VL
              IF (L.EQ.1) VJ1=A0/GP1*VJL
              IF (L.EQ.2) VJ2=A0/GP2*VJL
25         CONTINUE
        ELSE
           DO 40 L=1,2
              VV=VV0*L*L
              PX=1.0D0
              RP=1.0D0
              DO 30 K=1,K0
                 RP=-0.78125D-2*RP*(VV-(4.0*K-3.0)**2.0)*(VV-
     &              (4.0*K-1.0)**2.0)/(K*(2.0*K-1.0)*X2)
30               PX=PX+RP
              QX=1.0D0
              RQ=1.0D0
              DO 35 K=1,K0
                 RQ=-0.78125D-2*RQ*(VV-(4.0*K-1.0)**2.0)*(VV-
     &              (4.0*K+1.0)**2.0)/(K*(2.0*K+1.0)*X2)
35               QX=QX+RQ
              QX=0.125D0*(VV-1.0)*QX/X
              XK=X-(0.5D0*L/3.0D0+0.25D0)*PI
              A0=DSQRT(RP2/X)
              CK=DCOS(XK)
              SK=DSIN(XK)
              IF (L.EQ.1) THEN
                 VJ1=A0*(PX*CK-QX*SK)
                 VY1=A0*(PX*SK+QX*CK)
              ELSE IF (L.EQ.2) THEN
                 VJ2=A0*(PX*CK-QX*SK)
                 VY2=A0*(PX*SK+QX*CK)
              ENDIF
40         CONTINUE
        ENDIF
        IF (X.LE.12.0D0) THEN
           UJ1=0.0D0
           UJ2=0.0D0
           DO 55 L=1,2
              VL=L/3.0D0
              VJL=1.0D0
              R=1.0D0
              DO 45 K=1,40
                 R=-0.25D0*R*X2/(K*(K-VL))
                 VJL=VJL+R
                 IF (DABS(R).LT.1.0D-15) GO TO 50
45            CONTINUE
50            B0=(2.0D0/X)**VL
              IF (L.EQ.1) UJ1=B0*VJL/GN1
              IF (L.EQ.2) UJ2=B0*VJL/GN2
55         CONTINUE
           PV1=PI/3.0D0
           PV2=PI/1.5D0
           VY1=UU0*(VJ1*DCOS(PV1)-UJ1)
           VY2=UU0*(VJ2*DCOS(PV2)-UJ2)
        ENDIF
        IF (X.LE.18.0) THEN
           DO 70 L=1,2
              VL=L/3.0D0
              VIL=1.0D0
              R=1.0D0
              DO 60 K=1,40
                 R=0.25D0*R*X2/(K*(K+VL))
                 VIL=VIL+R
                 IF (DABS(R).LT.1.0D-15) GO TO 65
60            CONTINUE
65            A0=(0.5D0*X)**VL
              IF (L.EQ.1) VI1=A0/GP1*VIL
              IF (L.EQ.2) VI2=A0/GP2*VIL
70         CONTINUE
        ELSE
           C0=DEXP(X)/DSQRT(2.0D0*PI*X)
           DO 80 L=1,2
              VV=VV0*L*L
              VSL=1.0D0
              R=1.0D0
              DO 75 K=1,K0
                 R=-0.125D0*R*(VV-(2.0D0*K-1.0D0)**2.0)/(K*X)
75               VSL=VSL+R
              IF (L.EQ.1) VI1=C0*VSL
              IF (L.EQ.2) VI2=C0*VSL
80         CONTINUE
        ENDIF
        IF (X.LE.9.0D0) THEN
           GN=0.0D0
           DO 95 L=1,2
              VL=L/3.0D0
               IF (L.EQ.1) GN=GN1
               IF (L.EQ.2) GN=GN2
               A0=(2.0D0/X)**VL/GN
               SUM=1.0D0
               R=1.0D0
               DO 85 K=1,60
                  R=0.25D0*R*X2/(K*(K-VL))
                  SUM=SUM+R
                  IF (DABS(R).LT.1.0D-15) GO TO 90
85             CONTINUE
90            IF (L.EQ.1) VK1=0.5D0*UU0*PI*(SUM*A0-VI1)
              IF (L.EQ.2) VK2=0.5D0*UU0*PI*(SUM*A0-VI2)
95         CONTINUE
        ELSE
           C0=DEXP(-X)*DSQRT(0.5D0*PI/X)
           DO 105 L=1,2
              VV=VV0*L*L
              SUM=1.0D0
              R=1.0D0
              DO 100 K=1,K0
                 R=0.125D0*R*(VV-(2.0*K-1.0)**2.0)/(K*X)
100              SUM=SUM+R
              IF (L.EQ.1) VK1=C0*SUM
              IF (L.EQ.2) VK2=C0*SUM
105        CONTINUE
        ENDIF
        RETURN
        END



C       **********************************

        SUBROUTINE CIKVB(V,Z,VM,CBI,CDI,CBK,CDK)
C
C       ===========================================================
C       Purpose: Compute the modified Bessel functions Iv(z), Kv(z)
C                and their derivatives for an arbitrary order and
C                complex argument
C       Input :  z --- Complex argument z
C                v --- Real order of Iv(z) and Kv(z)
C                      ( v =n+v0, n = 0,1,2,..., 0 ≤ v0 < 1 )
C       Output:  CBI(n) --- In+v0(z)
C                CDI(n) --- In+v0'(z)
C                CBK(n) --- Kn+v0(z)
C                CDK(n) --- Kn+v0'(z)
C                VM --- Highest order computed
C       Routines called:
C            (1) GAMMA2 for computing the gamma function
C            (2) MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBI(0:*),CDI(0:*),CBK(0:*),CDK(0:*)
        Z1=Z
        Z2=Z*Z
        A0=CDABS(Z)
        PI=3.141592653589793D0
        CI=(0.0D0,1.0D0)
        N=INT(V)
        V0=V-N
        PIV=PI*V0
        VT=4.0D0*V0*V0
        IF (N.EQ.0) N=1
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBI(K)=0.0D0
              CDI(K)=0.0D0
              CBK(K)=-1.0D+300
10            CDK(K)=1.0D+300
           IF (V0.EQ.0.0) THEN
              CBI(0)=(1.0D0,0.0D0)
              CDI(1)=(0.5D0,0.0D0)
           ENDIF
           VM=V
           RETURN
        ENDIF
        K0=14
        IF (A0.GE.35.0) K0=10
        IF (A0.GE.50.0) K0=8
        IF (DBLE(Z).LT.0.0) Z1=-Z
        IF (A0.LT.18.0) THEN
           IF (V0.EQ.0.0) THEN
              CA1=(1.0D0,0.0D0)
           ELSE
              V0P=1.0D0+V0
              CALL GAMMA2(V0P,GAP)
              CA1=(0.5D0*Z1)**V0/GAP
           ENDIF
           CI0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 15 K=1,50
              CR=0.25D0*CR*Z2/(K*(K+V0))
              CI0=CI0+CR
              IF (CDABS(CR/CI0).LT.1.0D-15) GO TO 20
15         CONTINUE
20         CBI0=CI0*CA1
        ELSE
           CA=CDEXP(Z1)/CDSQRT(2.0D0*PI*Z1)
           CS=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 25 K=1,K0
              CR=-0.125D0*CR*(VT-(2.0D0*K-1.0D0)**2.0)/(K*Z1)
25            CS=CS+CR
           CBI0=CA*CS
        ENDIF
        M=MSTA1(A0,200)
        IF (M.LT.N) THEN
           N=M
        ELSE
           M=MSTA2(A0,N,15)
        ENDIF
        CF2=(0.0D0,0.0D0)
        CF1=(1.0D-100,0.0D0)
        DO 30 K=M,0,-1
           CF=2.0D0*(V0+K+1.0D0)/Z1*CF1+CF2
           IF (K.LE.N) CBI(K)=CF
           CF2=CF1
30         CF1=CF
        CS=CBI0/CF
        DO 35 K=0,N
35         CBI(K)=CS*CBI(K)
        IF (A0.LE.9.0) THEN
           IF (V0.EQ.0.0) THEN
              CT=-CDLOG(0.5D0*Z1)-0.5772156649015329D0
              CS=(0.0D0,0.0D0)
              W0=0.0D0
              CR=(1.0D0,0.0D0)
              DO 40 K=1,50
                 W0=W0+1.0D0/K
                 CR=0.25D0*CR/(K*K)*Z2
                 CP=CR*(W0+CT)
                 CS=CS+CP
                 IF (K.GE.10.AND.CDABS(CP/CS).LT.1.0D-15) GO TO 45
40            CONTINUE
45            CBK0=CT+CS
           ELSE
              V0N=1.0D0-V0
              CALL GAMMA2(V0N,GAN)
              CA2=1.0D0/(GAN*(0.5D0*Z1)**V0)
              CA1=(0.5D0*Z1)**V0/GAP
              CSU=CA2-CA1
              CR1=(1.0D0,0.0D0)
              CR2=(1.0D0,0.0D0)
              DO 50 K=1,50
                 CR1=0.25D0*CR1*Z2/(K*(K-V0))
                 CR2=0.25D0*CR2*Z2/(K*(K+V0))
                 CP=CA2*CR1-CA1*CR2
                 CSU=CSU+CP
                 IF (K.GE.10.AND.CDABS(CP/CSU).LT.1.0D-15) GO TO 55
50            CONTINUE
55            CBK0=0.5D0*PI*CSU/DSIN(PIV)
           ENDIF
        ELSE
           CB=CDEXP(-Z1)*CDSQRT(0.5D0*PI/Z1)
           CS=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 60 K=1,K0
              CR=0.125D0*CR*(VT-(2.0D0*K-1.0D0)**2.0)/(K*Z1)
60            CS=CS+CR
           CBK0=CB*CS
        ENDIF
        CBK(0)=CBK0
        IF (DBLE(Z).LT.0.0) THEN
           DO 65 K=0,N
              CVK=CDEXP((K+V0)*PI*CI)
              IF (DIMAG(Z).LT.0.0D0) THEN
                 CBK(K)=CVK*CBK(K)+PI*CI*CBI(K)
                 CBI(K)=CBI(K)/CVK
              ELSE IF (DIMAG(Z).GT.0.0) THEN
                 CBK(K)=CBK(K)/CVK-PI*CI*CBI(K)
                 CBI(K)=CVK*CBI(K)
              ENDIF
65         CONTINUE
        ENDIF
        DO 70 K=1,N
           CKK=(1.0D0/Z-CBI(K)*CBK(K-1))/CBI(K-1)
           CBK(K)=CKK
70      CONTINUE
        CDI(0)=V0/Z*CBI(0)+CBI(1)
        CDK(0)=V0/Z*CBK(0)-CBK(1)
        DO 80 K=1,N
           CDI(K)=-(K+V0)/Z*CBI(K)+CBI(K-1)
80         CDK(K)=-(K+V0)/Z*CBK(K)-CBK(K-1)
        VM=N+V0
        RETURN
        END



C       **********************************

        SUBROUTINE CIKVA(V,Z,VM,CBI,CDI,CBK,CDK)
C
C       ============================================================
C       Purpose: Compute the modified Bessel functions Iv(z), Kv(z)
C                and their derivatives for an arbitrary order and
C                complex argument
C       Input :  z --- Complex argument
C                v --- Real order of Iv(z) and Kv(z)
C                      ( v = n+v0, n = 0,1,2,…, 0 ≤ v0 < 1 )
C       Output:  CBI(n) --- In+v0(z)
C                CDI(n) --- In+v0'(z)
C                CBK(n) --- Kn+v0(z)
C                CDK(n) --- Kn+v0'(z)
C                VM --- Highest order computed
C       Routines called:
C            (1) GAMMA2 for computing the gamma function
C            (2) MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       ============================================================
C
        IMPLICIT DOUBLE PRECISION (A,G,P,R,V,W)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBI(0:*),CDI(0:*),CBK(0:*),CDK(0:*)
        PI=3.141592653589793D0
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z1=Z
        Z2=Z*Z
        N=INT(V)
        V0=V-N
        PIV=PI*V0
        VT=4.0D0*V0*V0
        IF (N.EQ.0) N=1
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBI(K)=0.0D0
              CDI(K)=0.0D0
              CBK(K)=-1.0D+300
10            CDK(K)=1.0D+300
           IF (V0.EQ.0.0) THEN
              CBI(0)=(1.0D0,0.0D0)
              CDI(1)=(0.5D0,0.0D0)
           ENDIF
           VM=V
           RETURN
        ENDIF
        K0=14
        IF (A0.GE.35.0) K0=10
        IF (A0.GE.50.0) K0=8
        IF (DBLE(Z).LT.0.0) Z1=-Z
        IF (A0.LT.18.0) THEN
           IF (V0.EQ.0.0) THEN
              CA1=(1.0D0,0.0D0)
           ELSE
              V0P=1.0D0+V0
              CALL GAMMA2(V0P,GAP)
              CA1=(0.5D0*Z1)**V0/GAP
           ENDIF
           CI0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 15 K=1,50
              CR=0.25D0*CR*Z2/(K*(K+V0))
              CI0=CI0+CR
              IF (CDABS(CR).LT.CDABS(CI0)*1.0D-15) GO TO 20
15         CONTINUE
20         CBI0=CI0*CA1
        ELSE
           CA=CDEXP(Z1)/CDSQRT(2.0D0*PI*Z1)
           CS=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 25 K=1,K0
              CR=-0.125D0*CR*(VT-(2.0D0*K-1.0D0)**2.0)/(K*Z1)
25            CS=CS+CR
           CBI0=CA*CS
        ENDIF
        M=MSTA1(A0,200)
        IF (M.LT.N) THEN
           N=M
        ELSE
           M=MSTA2(A0,N,15)
        ENDIF
        CF2=(0.0D0,0.0D0)
        CF1=(1.0D-100,0.0D0)
        DO 30 K=M,0,-1
           CF=2.0D0*(V0+K+1.0D0)/Z1*CF1+CF2
           IF (K.LE.N) CBI(K)=CF
           CF2=CF1
30         CF1=CF
        CS=CBI0/CF
        DO 35 K=0,N
35         CBI(K)=CS*CBI(K)
        IF (A0.LE.9.0) THEN
           IF (V0.EQ.0.0) THEN
              CT=-CDLOG(0.5D0*Z1)-0.5772156649015329D0
              CS=(0.0D0,0.0D0)
              W0=0.0D0
              CR=(1.0D0,0.0D0)
              DO 40 K=1,50
                 W0=W0+1.0D0/K
                 CR=0.25D0*CR/(K*K)*Z2
                 CP=CR*(W0+CT)
                 CS=CS+CP
                 IF (K.GE.10.AND.CDABS(CP/CS).LT.1.0D-15) GO TO 45
40            CONTINUE
45            CBK0=CT+CS
           ELSE
              V0N=1.0D0-V0
              CALL GAMMA2(V0N,GAN)
              CA2=1.0D0/(GAN*(0.5D0*Z1)**V0)
              CA1=(0.5D0*Z1)**V0/GAP
              CSU=CA2-CA1
              CR1=(1.0D0,0.0D0)
              CR2=(1.0D0,0.0D0)
              WS0=0.0D0
              DO 50 K=1,50
                 CR1=0.25D0*CR1*Z2/(K*(K-V0))
                 CR2=0.25D0*CR2*Z2/(K*(K+V0))
                 CSU=CSU+CA2*CR1-CA1*CR2
                 WS=CDABS(CSU)
                 IF (K.GE.10.AND.DABS(WS-WS0)/WS.LT.1.0D-15) GO TO 55
                 WS0=WS
50            CONTINUE
55            CBK0=0.5D0*PI*CSU/DSIN(PIV)
           ENDIF
        ELSE
           CB=CDEXP(-Z1)*CDSQRT(0.5D0*PI/Z1)
           CS=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 60 K=1,K0
              CR=0.125D0*CR*(VT-(2.0D0*K-1.0D0)**2.0)/(K*Z1)
60            CS=CS+CR
           CBK0=CB*CS
        ENDIF
        CBK1=(1.0D0/Z1-CBI(1)*CBK0)/CBI(0)
        CBK(0)=CBK0
        CBK(1)=CBK1
        CG0=CBK0
        CG1=CBK1
        DO 65 K=2,N
           CGK=2.0D0*(V0+K-1.0D0)/Z1*CG1+CG0
           CBK(K)=CGK
           CG0=CG1
65         CG1=CGK
        IF (DBLE(Z).LT.0.0) THEN
           DO 70 K=0,N
              CVK=CDEXP((K+V0)*PI*CI)
              IF (DIMAG(Z).LT.0.0D0) THEN
                 CBK(K)=CVK*CBK(K)+PI*CI*CBI(K)
                 CBI(K)=CBI(K)/CVK
              ELSE IF (DIMAG(Z).GT.0.0) THEN
                 CBK(K)=CBK(K)/CVK-PI*CI*CBI(K)
                 CBI(K)=CVK*CBI(K)
              ENDIF
70         CONTINUE
        ENDIF
        CDI(0)=V0/Z*CBI(0)+CBI(1)
        CDK(0)=V0/Z*CBK(0)-CBK(1)
        DO 75 K=1,N
           CDI(K)=-(K+V0)/Z*CBI(K)+CBI(K-1)
75         CDK(K)=-(K+V0)/Z*CBK(K)-CBK(K-1)
        VM=N+V0
        RETURN
        END



C       **********************************

        SUBROUTINE CFC(Z,ZF,ZD)
C
C       =========================================================
C       Purpose: Compute complex Fresnel integral C(z) and C'(z)
C       Input :  z --- Argument of C(z)
C       Output:  ZF --- C(z)
C                ZD --- C'(z)
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (E,P,W)
        IMPLICIT COMPLEX *16 (C,S,Z)
        EPS=1.0D-14
        PI=3.141592653589793D0
        W0=CDABS(Z)
        ZP=0.5D0*PI*Z*Z
        ZP2=ZP*ZP
        Z0=(0.0D0,0.0D0)
        IF (Z.EQ.Z0) THEN
           C=Z0
        ELSE IF (W0.LE.2.5) THEN
           CR=Z
           C=CR
           WA0=0.0D0
           DO 10 K=1,80
              CR=-.5D0*CR*(4.0D0*K-3.0D0)/K/(2.0D0*K-1.0D0)
     &          /(4.0D0*K+1.0D0)*ZP2
              C=C+CR
              WA=CDABS(C)
              IF (DABS((WA-WA0)/WA).LT.EPS.AND.K.GT.10) GO TO 30
10            WA0=WA
        ELSE IF (W0.GT.2.5.AND.W0.LT.4.5) THEN
           M=85
           C=Z0
           CF1=Z0
           CF0=(1.0D-100,0.0D0)
           DO 15 K=M,0,-1
              CF=(2.0D0*K+3.0D0)*CF0/ZP-CF1
              IF (K.EQ.INT(K/2)*2) C=C+CF
              CF1=CF0
15            CF0=CF
           C=CDSQRT(2.0D0/(PI*ZP))*CDSIN(ZP)/CF*C
        ELSE
           CR=(1.0D0,0.0D0)
           CF=(1.0D0,0.0D0)
           DO 20 K=1,20
              CR=-.25D0*CR*(4.0D0*K-1.0D0)*(4.0D0*K-3.0D0)/ZP2
20            CF=CF+CR
           CR=1.0D0/(PI*Z*Z)
           CG=CR
           DO 25 K=1,12
              CR=-.25D0*CR*(4.0D0*K+1.0D0)*(4.0D0*K-1.0D0)/ZP2
25            CG=CG+CR
           C=.5D0+(CF*CDSIN(ZP)-CG*CDCOS(ZP))/(PI*Z)
        ENDIF
30      ZF=C
        ZD=CDCOS(0.5*PI*Z*Z)
        RETURN
        END



C       **********************************

        SUBROUTINE FCS(X,C,S)
C
C       =================================================
C       Purpose: Compute Fresnel integrals C(x) and S(x)
C       Input :  x --- Argument of C(x) and S(x)
C       Output:  C --- C(x)
C                S --- S(x)
C       =================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        EPS=1.0D-15
        PI=3.141592653589793D0
        XA=DABS(X)
        PX=PI*XA
        T=.5D0*PX*XA
        T2=T*T
        IF (XA.EQ.0.0) THEN
           C=0.0D0
           S=0.0D0
        ELSE IF (XA.LT.2.5D0) THEN
           R=XA
           C=R
           DO 10 K=1,50
              R=-.5D0*R*(4.0D0*K-3.0D0)/K/(2.0D0*K-1.0D0)
     &          /(4.0D0*K+1.0D0)*T2
              C=C+R
              IF (DABS(R).LT.DABS(C)*EPS) GO TO 15
10         CONTINUE
15         S=XA*T/3.0D0
           R=S
           DO 20 K=1,50
              R=-.5D0*R*(4.0D0*K-1.0D0)/K/(2.0D0*K+1.0D0)
     &          /(4.0D0*K+3.0D0)*T2
              S=S+R
              IF (DABS(R).LT.DABS(S)*EPS) GO TO 40
20         CONTINUE
        ELSE IF (XA.LT.4.5D0) THEN
           M=INT(42.0+1.75*T)
           SU=0.0D0
           C=0.0D0
           S=0.0D0
           F1=0.0D0
           F0=1.0D-100
           DO 25 K=M,0,-1
              F=(2.0D0*K+3.0D0)*F0/T-F1
              IF (K.EQ.INT(K/2)*2) THEN
                 C=C+F
              ELSE
                 S=S+F
              ENDIF
              SU=SU+(2.0D0*K+1.0D0)*F*F
              F1=F0
25            F0=F
           Q=DSQRT(SU)
           C=C*XA/Q
           S=S*XA/Q
        ELSE
           R=1.0D0
           F=1.0D0
           DO 30 K=1,20
              R=-.25D0*R*(4.0D0*K-1.0D0)*(4.0D0*K-3.0D0)/T2
30            F=F+R
           R=1.0D0/(PX*XA)
           G=R
           DO 35 K=1,12
              R=-.25D0*R*(4.0D0*K+1.0D0)*(4.0D0*K-1.0D0)/T2
35            G=G+R
           T0=T-INT(T/(2.0D0*PI))*2.0D0*PI
           C=.5D0+(F*DSIN(T0)-G*DCOS(T0))/PX
           S=.5D0-(F*DCOS(T0)+G*DSIN(T0))/PX
        ENDIF
40      IF (X.LT.0.0D0) THEN
           C=-C
           S=-S
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE RCTJ(N,X,NM,RJ,DJ)
C
C       ========================================================
C       Purpose: Compute Riccati-Bessel functions of the first
C                kind and their derivatives
C       Input:   x --- Argument of Riccati-Bessel function
C                n --- Order of jn(x)  ( n = 0,1,2,... )
C       Output:  RJ(n) --- x·jn(x)
C                DJ(n) --- [x·jn(x)]'
C                NM --- Highest order computed
C       Routines called:
C                MSTA1 and MSTA2 for computing the starting
C                point for backward recurrence
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION RJ(0:N),DJ(0:N)
        NM=N
        IF (DABS(X).LT.1.0D-100) THEN
           DO 10 K=0,N
              RJ(K)=0.0D0
10            DJ(K)=0.0D0
           DJ(0)=1.0D0
           RETURN
        ENDIF
        RJ(0)=DSIN(X)
        RJ(1)=RJ(0)/X-DCOS(X)
        RJ0=RJ(0)
        RJ1=RJ(1)
        CS=0.0D0
        F=0.0D0
        IF (N.GE.2) THEN
           M=MSTA1(X,200)
           IF (M.LT.N) THEN
              NM=M
           ELSE
              M=MSTA2(X,N,15)
           ENDIF
           F0=0.0D0
           F1=1.0D-100
           DO 15 K=M,0,-1
              F=(2.0D0*K+3.0D0)*F1/X-F0
              IF (K.LE.NM) RJ(K)=F
              F0=F1
15            F1=F
           IF (DABS(RJ0).GT.DABS(RJ1)) CS=RJ0/F
           IF (DABS(RJ0).LE.DABS(RJ1)) CS=RJ1/F0
           DO 20 K=0,NM
20            RJ(K)=CS*RJ(K)
        ENDIF
        DJ(0)=DCOS(X)
        DO 25 K=1,NM
25         DJ(K)=-K*RJ(K)/X+RJ(K-1)
        RETURN
        END



C       **********************************

        SUBROUTINE HERZO(N,X,W)
C
C       ========================================================
C       Purpose : Compute the zeros of Hermite polynomial Ln(x)
C                 in the interval [-∞,∞], and the corresponding
C                 weighting coefficients for Gauss-Hermite
C                 integration
C       Input :   n    --- Order of the Hermite polynomial
C                 X(n) --- Zeros of the Hermite polynomial
C                 W(n) --- Corresponding weighting coefficients
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION X(N),W(N)
        HN=1.0D0/N
        ZL=-1.1611D0+1.46D0*N**0.5
        Z=0.0D0
        HF=0.0D0
        HD=0.0D0
        DO 40 NR=1,N/2
           IF (NR.EQ.1) Z=ZL
           IF (NR.NE.1) Z=Z-HN*(N/2+1-NR)
           IT=0
10         IT=IT+1
           Z0=Z
           F0=1.0D0
           F1=2.0D0*Z
           DO 15 K=2,N
              HF=2.0D0*Z*F1-2.0D0*(K-1.0D0)*F0
              HD=2.0D0*K*F1
              F0=F1
15            F1=HF
           P=1.0D0
           DO 20 I=1,NR-1
20            P=P*(Z-X(I))
           FD=HF/P
           Q=0.0D0
           DO 30 I=1,NR-1
              WP=1.0D0
              DO 25 J=1,NR-1
                 IF (J.EQ.I) GO TO 25
                 WP=WP*(Z-X(J))
25            CONTINUE
30            Q=Q+WP
           GD=(HD-Q*FD)/P
           Z=Z-FD/GD
           IF (IT.LE.40.AND.DABS((Z-Z0)/Z).GT.1.0D-15) GO TO 10
           X(NR)=Z
           X(N+1-NR)=-Z
           R=1.0D0
           DO 35 K=1,N
35            R=2.0D0*R*K
           W(NR)=3.544907701811D0*R/(HD*HD)
40         W(N+1-NR)=W(NR)
        IF (N.NE.2*INT(N/2)) THEN
           R1=1.0D0
           R2=1.0D0
           DO 45 J=1,N
              R1=2.0D0*R1*J
              IF (J.GE.(N+1)/2) R2=R2*J
45         CONTINUE
           W(N/2+1)=0.88622692545276D0*R1/(R2*R2)
           X(N/2+1)=0.0D0
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE JY01B(X,BJ0,DJ0,BJ1,DJ1,BY0,DY0,BY1,DY1)
C
C       =======================================================
C       Purpose: Compute Bessel functions J0(x), J1(x), Y0(x),
C                Y1(x), and their derivatives
C       Input :  x   --- Argument of Jn(x) & Yn(x) ( x ≥ 0 )
C       Output:  BJ0 --- J0(x)
C                DJ0 --- J0'(x)
C                BJ1 --- J1(x)
C                DJ1 --- J1'(x)
C                BY0 --- Y0(x)
C                DY0 --- Y0'(x)
C                BY1 --- Y1(x)
C                DY1 --- Y1'(x)
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           BJ0=1.0D0
           BJ1=0.0D0
           DJ0=0.0D0
           DJ1=0.5D0
           BY0=-1.0D+300
           BY1=-1.0D+300
           DY0=1.0D+300
           DY1=1.0D+300
           RETURN
        ELSE IF (X.LE.4.0D0) THEN
           T=X/4.0D0
           T2=T*T
           BJ0=((((((-.5014415D-3*T2+.76771853D-2)*T2
     &         -.0709253492D0)*T2+.4443584263D0)*T2
     &         -1.7777560599D0)*T2+3.9999973021D0)
     &         *T2-3.9999998721D0)*T2+1.0D0
           BJ1=T*(((((((-.1289769D-3*T2+.22069155D-2)
     &         *T2-.0236616773D0)*T2+.1777582922D0)*T2
     &         -.8888839649D0)*T2+2.6666660544D0)*T2
     &         -3.9999999710D0)*T2+1.9999999998D0)
           BY0=(((((((-.567433D-4*T2+.859977D-3)*T2
     &         -.94855882D-2)*T2+.0772975809D0)*T2
     &         -.4261737419D0)*T2+1.4216421221D0)*T2
     &         -2.3498519931D0)*T2+1.0766115157D0)*T2
     &         +.3674669052D0
           BY0=2.0D0/PI*DLOG(X/2.0D0)*BJ0+BY0
           BY1=((((((((.6535773D-3*T2-.0108175626D0)*T2
     &         +.107657606D0)*T2-.7268945577D0)*T2
     &         +3.1261399273D0)*T2-7.3980241381D0)*T2
     &         +6.8529236342D0)*T2+.3932562018D0)*T2
     &         -.6366197726D0)/X
           BY1=2.0D0/PI*DLOG(X/2.0D0)*BJ1+BY1
        ELSE
           T=4.0D0/X
           T2=T*T
           A0=DSQRT(2.0D0/(PI*X))
           P0=((((-.9285D-5*T2+.43506D-4)*T2-.122226D-3)*T2
     &        +.434725D-3)*T2-.4394275D-2)*T2+.999999997D0
           Q0=T*(((((.8099D-5*T2-.35614D-4)*T2+.85844D-4)*T2
     &        -.218024D-3)*T2+.1144106D-2)*T2-.031249995D0)
           TA0=X-.25D0*PI
           BJ0=A0*(P0*DCOS(TA0)-Q0*DSIN(TA0))
           BY0=A0*(P0*DSIN(TA0)+Q0*DCOS(TA0))
           P1=((((.10632D-4*T2-.50363D-4)*T2+.145575D-3)*T2
     &        -.559487D-3)*T2+.7323931D-2)*T2+1.000000004D0
           Q1=T*(((((-.9173D-5*T2+.40658D-4)*T2-.99941D-4)*T2
     &        +.266891D-3)*T2-.1601836D-2)*T2+.093749994D0)
           TA1=X-.75D0*PI
           BJ1=A0*(P1*DCOS(TA1)-Q1*DSIN(TA1))
           BY1=A0*(P1*DSIN(TA1)+Q1*DCOS(TA1))
        ENDIF
        DJ0=-BJ1
        DJ1=BJ0-BJ1/X
        DY0=-BY1
        DY1=BY0-BY1/X
        RETURN
        END

C       **********************************

        SUBROUTINE ENXB(N,X,EN)
C
C       ===============================================
C       Purpose: Compute exponential integral En(x)
C       Input :  x --- Argument of En(x)
C                n --- Order of En(x)  (n = 0,1,2,...)
C       Output:  EN(n) --- En(x)
C       ===============================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION EN(0:N)
        IF (X.EQ.0.0) THEN
           EN(0)=1.0D+300
           EN(1)=1.0D+300
           DO 10 K=2,N
10            EN(K)=1.0D0/(K-1.0)
           RETURN
        ELSE IF (X.LE.1.0) THEN
           EN(0)=DEXP(-X)/X
           S0=0.0D0
           DO 40 L=1,N
              RP=1.0D0
              DO 15 J=1,L-1
15               RP=-RP*X/J
              PS=-0.5772156649015328D0
              DO 20 M=1,L-1
20               PS=PS+1.0D0/M
              ENS=RP*(-DLOG(X)+PS)
              S=0.0D0
              DO 30 M=0,20
                 IF (M.EQ.L-1) GO TO 30
                 R=1.0D0
                 DO 25 J=1,M
25                  R=-R*X/J
                 S=S+R/(M-L+1.0D0)
                 IF (DABS(S-S0).LT.DABS(S)*1.0D-15) GO TO 35
                 S0=S
30            CONTINUE
35            EN(L)=ENS-S
40         CONTINUE
        ELSE
           EN(0)=DEXP(-X)/X
           M=15+INT(100.0/X)
           DO 50 L=1,N
              T0=0.0D0
              DO 45 K=M,1,-1
45               T0=(L+K-1.0D0)/(1.0D0+K/(X+T0))
              T=1.0D0/(X+T0)
50            EN(L)=DEXP(-X)*T
        ENDIF
        END

C       **********************************

        SUBROUTINE SPHK(N,X,NM,SK,DK)
C
C       =====================================================
C       Purpose: Compute modified spherical Bessel functions
C                of the second kind, kn(x) and kn'(x)
C       Input :  x --- Argument of kn(x)  ( x ≥ 0 )
C                n --- Order of kn(x) ( n = 0,1,2,... )
C       Output:  SK(n) --- kn(x)
C                DK(n) --- kn'(x)
C                NM --- Highest order computed
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SK(0:N),DK(0:N)
        PI=3.141592653589793D0
        NM=N
        IF (X.LT.1.0D-60) THEN
           DO 10 K=0,N
              SK(K)=1.0D+300
10            DK(K)=-1.0D+300
           RETURN
        ENDIF
        SK(0)=0.5D0*PI/X*DEXP(-X)
        SK(1)=SK(0)*(1.0D0+1.0D0/X)
        F0=SK(0)
        F1=SK(1)
        DO 15 K=2,N
           F=(2.0D0*K-1.0D0)*F1/X+F0
           SK(K)=F
           IF (DABS(F).GT.1.0D+300) GO TO 20
           F0=F1
15         F1=F
20      NM=K-1
        DK(0)=-SK(1)
        DO 25 K=1,NM
25         DK(K)=-SK(K-1)-(K+1.0D0)/X*SK(K)
        RETURN
        END

C       **********************************

        SUBROUTINE ENXA(N,X,EN)
C
C       ============================================
C       Purpose: Compute exponential integral En(x)
C       Input :  x --- Argument of En(x) ( x ≤ 20 )
C                n --- Order of En(x)
C       Output:  EN(n) --- En(x)
C       Routine called: E1XB for computing E1(x)
C       ============================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION EN(0:N)
        EN(0)=DEXP(-X)/X
        CALL E1XB(X,E1)
        EN(1)=E1
        DO 10 K=2,N
           EK=(DEXP(-X)-X*E1)/(K-1.0D0)
           EN(K)=EK
10         E1=EK
        RETURN
        END



C       **********************************

        SUBROUTINE GAIH(X,GA)
C
C       =====================================================
C       Purpose: Compute gamma function Г(x)
C       Input :  x  --- Argument of Г(x), x = n/2, n=1,2,…
C       Output:  GA --- Г(x)
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X).AND.X.GT.0.0) THEN
           GA=1.0D0
           M1=INT(X-1.0)
           DO 10 K=2,M1
10            GA=GA*K
        ELSE IF (X+.5D0.EQ.INT(X+.5D0).AND.X.GT.0.0) THEN
           M=INT(X)
           GA=DSQRT(PI)
           DO 15 K=1,M
15            GA=0.5D0*GA*(2.0D0*K-1.0D0)
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE PBVV(V,X,VV,VP,PVF,PVD)
C
C       ===================================================
C       Purpose: Compute parabolic cylinder functions Vv(x)
C                and their derivatives
C       Input:   x --- Argument of Vv(x)
C                v --- Order of Vv(x)
C       Output:  VV(na) --- Vv(x)
C                VP(na) --- Vv'(x)
C                ( na = |n|, v = n+v0, |v0| < 1
C                  n = 0,±1,±2,… )
C                PVF --- Vv(x)
C                PVD --- Vv'(x)
C       Routines called:
C             (1) VVSA for computing Vv(x) for small |x|
C             (2) VVLA for computing Vv(x) for large |x|
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION VV(0:*),VP(0:*)
        PI=3.141592653589793D0
        XA=DABS(X)
        VH=V
        V=V+DSIGN(1.0D0,V)
        NV=INT(V)
        V0=V-NV
        NA=ABS(NV)
        QE=DEXP(0.25D0*X*X)
        Q2P=DSQRT(2.0D0/PI)
        JA=0
        IF (NA.GE.1) JA=1
        F=0.0D0
        IF (V.LE.0.0) THEN
           IF (V0.EQ.0.0) THEN
              IF (XA.LE.7.5) CALL VVSA(V0,X,PV0)
              IF (XA.GT.7.5) CALL VVLA(V0,X,PV0)
              F0=Q2P*QE
              F1=X*F0
              VV(0)=PV0
              VV(1)=F0
              VV(2)=F1
           ELSE
              DO 10 L=0,JA
                 V1=V0-L
                 IF (XA.LE.7.5) CALL VVSA(V1,X,F1)
                 IF (XA.GT.7.5) CALL VVLA(V1,X,F1)
                 IF (L.EQ.0) F0=F1
10            CONTINUE
              VV(0)=F0
              VV(1)=F1
           ENDIF
           KV=2
           IF (V0.EQ.0.0) KV=3
           DO 15 K=KV,NA
              F=X*F1+(K-V0-2.0D0)*F0
              VV(K)=F
              F0=F1
15            F1=F
        ELSE
           IF (X.GE.0.0.AND.X.LE.7.5D0) THEN
              V2=V
              IF (V2.LT.1.0) V2=V2+1.0D0
              CALL VVSA(V2,X,F1)
              V1=V2-1.0D0
              KV=INT(V2)
              CALL VVSA(V1,X,F0)
              VV(KV)=F1
              VV(KV-1)=F0
              DO 20 K=KV-2,0,-1
                 F=X*F0-(K+V0+2.0D0)*F1
                 IF (K.LE.NA) VV(K)=F
                 F1=F0
20               F0=F
           ELSE IF (X.GT.7.5D0) THEN
              CALL VVLA(V0,X,PV0)
              M=100+ABS(NA)
              VV(1)=PV0
              F1=0.0D0
              F0=1.0D-40
              DO 25 K=M,0,-1
                 F=X*F0-(K+V0+2.0D0)*F1
                 IF (K.LE.NA) VV(K)=F
                 F1=F0
25               F0=F
              S0=PV0/F
              DO 30 K=0,NA
30               VV(K)=S0*VV(K)
           ELSE
              IF (XA.LE.7.5D0) THEN
                 CALL VVSA(V0,X,F0)
                 V1=V0+1.0
                 CALL VVSA(V1,X,F1)
              ELSE
                 CALL VVLA(V0,X,F0)
                 V1=V0+1.0D0
                 CALL VVLA(V1,X,F1)
              ENDIF
              VV(0)=F0
              VV(1)=F1
              DO 35 K=2,NA
                 F=(X*F1-F0)/(K+V0)
                 VV(K)=F
                 F0=F1
35               F1=F
           ENDIF
        ENDIF
        DO 40 K=0,NA-1
           V1=V0+K
           IF (V.GE.0.0D0) THEN
              VP(K)=0.5D0*X*VV(K)-(V1+1.0D0)*VV(K+1)
           ELSE
              VP(K)=-0.5D0*X*VV(K)+VV(K+1)
           ENDIF
40      CONTINUE
        PVF=VV(NA-1)
        PVD=VP(NA-1)
        V=VH
        RETURN
        END



C       **********************************

        SUBROUTINE CLQMN(MM,M,N,X,Y,NTYPE,CQM,CQD)
C
C       =======================================================
C       Purpose: Compute the associated Legendre functions of
C                the second kind, Qmn(z) and Qmn'(z) for a
C                complex argument
C       Input :  x  --- Real part of z
C                y  --- Imaginary part of z
C                m  --- Order of Qmn(z)  ( m = 0,1,2,… )
C                n  --- Degree of Qmn(z) ( n = 0,1,2,… )
C                mm --- Physical dimension of CQM and CQD
C       Output:  CQM(m,n) --- Qmn(z)
C                CQD(m,n) --- Qmn'(z)
C       =======================================================
C
        IMPLICIT DOUBLE PRECISION (X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CQM(0:MM,0:N),CQD(0:MM,0:N)
        Z = DCMPLX(X, Y)
        IF (DABS(X).EQ.1.0D0.AND.Y.EQ.0.0D0) THEN
           DO 10 I=0,M
           DO 10 J=0,N
              CQM(I,J)=(1.0D+300,0.0D0)
              CQD(I,J)=(1.0D+300,0.0D0)
10         CONTINUE
           RETURN
        ENDIF
        XC=CDABS(Z)
        IF (NTYPE.EQ.2) THEN
           LS=1
        ELSE
           LS=-1
        END IF
        ZQ=CDSQRT(LS*(1.0D0-Z*Z))
        ZS=LS*(1.0D0-Z*Z)
        CQ0=0.5D0*CDLOG(LS*(1.0D0+Z)/(1.0D0-Z))
        IF ((XC.LT.1.0001D0).OR.(NTYPE.EQ.2)) THEN
           CQM(0,0)=CQ0
           CQM(0,1)=Z*CQ0-1.0D0
           CQM(1,0)=-1.0D0/ZQ
           CQM(1,1)=-LS*ZQ*(CQ0+Z/(1.0D0-Z*Z))
C       DLMF 14.10.3 applied to Q
           DO 15 I=0,1
           DO 15 J=2,N
              CQM(I,J)=((2.0D0*J-1.0D0)*Z*CQM(I,J-1)
     &                -(J+I-1.0D0)*CQM(I,J-2))/(J-I)
15         CONTINUE
C       NTYPE=2: DLMF 14.10.1 applied to Q
C       NTYPE=3: DLMF 14.10.6 applied to Q
           DO 20 J=0,N
           DO 20 I=2,M
              CQM(I,J)=-2.0D0*(I-1.0D0)*Z/ZQ*CQM(I-1,J)-LS*
     &                 (J+I-1.0D0)*(J-I+2.0D0)*CQM(I-2,J)
20         CONTINUE
        ELSE
           IF (XC.GT.1.1) THEN
              KM=40+M+N
           ELSE
              KM=(40+M+N)*INT(-1.0-1.8*LOG(XC-1.0))
           ENDIF
C       backward recursion with DLMF 14.10.3 applied to Q with mu=0
           CQF2=(0.0D0,0.0D0)
           CQF1=(1.0D0,0.0D0)
           DO 25 K=KM,0,-1
              CQF0=((2*K+3.0D0)*Z*CQF1-(K+2.0D0)*CQF2)/(K+1.0D0)
              IF (K.LE.N) CQM(0,K)=CQF0
              CQF2=CQF1
25            CQF1=CQF0
           DO 30 K=0,N
30            CQM(0,K)=CQ0*CQM(0,K)/CQF0
           IF (M.GT.0) THEN
C       backward recursion with DLMF 14.10.3 applied to Q with mu=1
              CQF2=0.0D0
              CQF1=1.0D0
              DO 35 K=KM,0,-1
                 CQF0=((2*K+3.0D0)*Z*CQF1-(K+1.0D0)*CQF2)/(K+2.0D0)
                 IF (K.LE.N) CQM(1,K)=CQF0
                 CQF2=CQF1
35               CQF1=CQF0
              CQ10=-1.0D0/ZQ
              DO 40 K=0,N
40               CQM(1,K)=CQ10*CQM(1,K)/CQF0
              DO 45 J=0,N
                 CQ0=CQM(0,J)
                 CQ1=CQM(1,J)
C       NTYPE=2: DLMF 14.10.1 applied to Q
C       NTYPE=3: DLMF 14.10.6 applied to Q
                 DO 45 I=0,M-2
                    CQF=-2.0D0*(I+1)*Z/ZQ*CQ1-LS*(J-I)*(J+I+1.0D0)*CQ0
                    CQM(I+2,J)=CQF
                    CQ0=CQ1
                    CQ1=CQF
45            CONTINUE
           ENDIF
        ENDIF
        CQD(0,0)=LS/ZS
C       DLMF 14.10.5 applied to Q with mu=0
        DO 50 J=1,N
50         CQD(0,J)=LS*J*(CQM(0,J-1)-Z*CQM(0,J))/ZS
C       NTYPE=2: derivative of DLMF 14.7.9
C                and DLMF 14.10.1 applied to Q
C       NTYPE=3: derivative of LDMF 14.7.12
C                and DLMF 14.10.6 applied to Q
        DO 55 J=0,N
        DO 55 I=1,M
           CQD(I,J)=LS*I*Z/ZS*CQM(I,J)+(I+J)*(J-I+1.0D0)
     &              /ZQ*CQM(I-1,J)
55      CONTINUE
        RETURN
        END


C       **********************************

        SUBROUTINE SEGV(M,N,C,KD,CV,EG)
C
C       =========================================================
C       Purpose: Compute the characteristic values of spheroidal
C                wave functions
C       Input :  m  --- Mode parameter
C                n  --- Mode parameter
C                c  --- Spheroidal parameter
C                KD --- Function code
C                       KD=1 for Prolate; KD=-1 for Oblate
C       Output:  CV --- Characteristic value for given m, n and c
C                EG(L) --- Characteristic value for mode m and n'
C                          ( L = n' - m + 1 )
C       =========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION B(100),H(100),D(300),E(300),F(300),CV0(100),
     &            A(300),G(300),EG(200)
        IF (C.LT.1.0D-10) THEN
           DO 5 I=1,N-M+1
5             EG(I)=(I+M)*(I+M-1.0D0)
           GO TO 70
        ENDIF
        ICM=(N-M+2)/2
        NM=10+INT(0.5*(N-M)+C)
        CS=C*C*KD
        K=0
        DO 60 L=0,1
           DO 10 I=1,NM
              IF (L.EQ.0) K=2*(I-1)
              IF (L.EQ.1) K=2*I-1
              DK0=M+K
              DK1=M+K+1
              DK2=2*(M+K)
              D2K=2*M+K
              A(I)=(D2K+2.0)*(D2K+1.0)/((DK2+3.0)*(DK2+5.0))*CS
              D(I)=DK0*DK1+(2.0*DK0*DK1-2.0*M*M-1.0)/((DK2-1.0)
     &             *(DK2+3.0))*CS
10            G(I)=K*(K-1.0)/((DK2-3.0)*(DK2-1.0))*CS
           DO 15 K=2,NM
              E(K)=DSQRT(A(K-1)*G(K))
15            F(K)=E(K)*E(K)
           F(1)=0.0D0
           E(1)=0.0D0
           XA=D(NM)+DABS(E(NM))
           XB=D(NM)-DABS(E(NM))
           NM1=NM-1
           DO 20 I=1,NM1
              T=DABS(E(I))+DABS(E(I+1))
              T1=D(I)+T
              IF (XA.LT.T1) XA=T1
              T1=D(I)-T
              IF (T1.LT.XB) XB=T1
20         CONTINUE
           DO 25 I=1,ICM
              B(I)=XA
25            H(I)=XB
           DO 55 K=1,ICM
              DO 30 K1=K,ICM
                 IF (B(K1).LT.B(K)) THEN
                    B(K)=B(K1)
                    GO TO 35
                 ENDIF
30            CONTINUE
35            IF (K.NE.1.AND.H(K).LT.H(K-1)) H(K)=H(K-1)
40            X1=(B(K)+H(K))/2.0D0
              CV0(K)=X1
              IF (DABS((B(K)-H(K))/X1).LT.1.0D-14) GO TO 50
              J=0
              S=1.0D0
              DO 45 I=1,NM
                 IF (S.EQ.0.0D0) S=S+1.0D-30
                 T=F(I)/S
                 S=D(I)-T-X1
                 IF (S.LT.0.0D0) J=J+1
45            CONTINUE
              IF (J.LT.K) THEN
                 H(K)=X1
              ELSE
                 B(K)=X1
                 IF (J.GE.ICM) THEN
                    B(ICM)=X1
                 ELSE
                    IF (H(J+1).LT.X1) H(J+1)=X1
                    IF (X1.LT.B(J)) B(J)=X1
                 ENDIF
              ENDIF
              GO TO 40
50            CV0(K)=X1
              IF (L.EQ.0) EG(2*K-1)=CV0(K)
              IF (L.EQ.1) EG(2*K)=CV0(K)
55         CONTINUE
60      CONTINUE
70      CV=EG(N-M+1)
        RETURN
        END


C       **********************************

        SUBROUTINE CIKNB(N,Z,NM,CBI,CDI,CBK,CDK)
C
C       ============================================================
C       Purpose: Compute modified Bessel functions In(z) and Kn(z),
C                and their derivatives for a complex argument
C       Input:   z --- Complex argument
C                n --- Order of In(z) and Kn(z)
C       Output:  CBI(n) --- In(z)
C                CDI(n) --- In'(z)
C                CBK(n) --- Kn(z)
C                CDK(n) --- Kn'(z)
C                NM --- Highest order computed
C       Routones called:
C                MSTA1 and MSTA2 to compute the starting point for
C                backward recurrence
C       ===========================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBI(0:N),CDI(0:N),CBK(0:N),CDK(0:N)
        PI=3.141592653589793D0
        EL=0.57721566490153D0
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBI(K)=(0.0D0,0.0D0)
              CBK(K)=(1.0D+300,0.0D0)
              CDI(K)=(0.0D0,0.0D0)
10            CDK(K)=-(1.0D+300,0.0D0)
           CBI(0)=(1.0D0,0.0D0)
           CDI(1)=(0.5D0,0.0D0)
           RETURN
        ENDIF
        Z1=Z
        CI=(0.0D0,1.0D0)
        IF (DBLE(Z).LT.0.0) Z1=-Z
        IF (N.EQ.0) NM=1
        M=MSTA1(A0,200)
        IF (M.LT.NM) THEN
           NM=M
        ELSE
           M=MSTA2(A0,NM,15)
        ENDIF
        CBS=0.0D0
        CSK0=0.0D0
        CF0=0.0D0
        CF1=1.0D-100
        DO 15 K=M,0,-1
           CF=2.0D0*(K+1.0D0)*CF1/Z1+CF0
           IF (K.LE.NM) CBI(K)=CF
           IF (K.NE.0.AND.K.EQ.2*INT(K/2)) CSK0=CSK0+4.0D0*CF/K
           CBS=CBS+2.0D0*CF
           CF0=CF1
15         CF1=CF
        CS0=CDEXP(Z1)/(CBS-CF)
        DO 20 K=0,NM
20         CBI(K)=CS0*CBI(K)
        IF (A0.LE.9.0) THEN
           CBK(0)=-(CDLOG(0.5D0*Z1)+EL)*CBI(0)+CS0*CSK0
           CBK(1)=(1.0D0/Z1-CBI(1)*CBK(0))/CBI(0)
        ELSE
           CA0=CDSQRT(PI/(2.0D0*Z1))*CDEXP(-Z1)
           K0=16
           IF (A0.GE.25.0) K0=10
           IF (A0.GE.80.0) K0=8
           IF (A0.GE.200.0) K0=6
           DO 30 L=0,1
              CBKL=1.0D0
              VT=4.0D0*L
              CR=(1.0D0,0.0D0)
              DO 25 K=1,K0
                 CR=0.125D0*CR*(VT-(2.0*K-1.0)**2)/(K*Z1)
25               CBKL=CBKL+CR
              CBK(L)=CA0*CBKL
30         CONTINUE
        ENDIF
        CG0=CBK(0)
        CG1=CBK(1)
        DO 35 K=2,NM
           CG=2.0D0*(K-1.0D0)/Z1*CG1+CG0
           CBK(K)=CG
           CG0=CG1
35         CG1=CG
        IF (DBLE(Z).LT.0.0) THEN
           FAC=1.0D0
           DO 45 K=0,NM
              IF (DIMAG(Z).LT.0.0) THEN
                 CBK(K)=FAC*CBK(K)+CI*PI*CBI(K)
              ELSE
                 CBK(K)=FAC*CBK(K)-CI*PI*CBI(K)
              ENDIF
              CBI(K)=FAC*CBI(K)
              FAC=-FAC
45         CONTINUE
        ENDIF
        CDI(0)=CBI(1)
        CDK(0)=-CBK(1)
        DO 50 K=1,NM
           CDI(K)=CBI(K-1)-K/Z*CBI(K)
50         CDK(K)=-CBK(K-1)-K/Z*CBK(K)
        RETURN
        END


C       **********************************

        SUBROUTINE CIKNA(N,Z,NM,CBI,CDI,CBK,CDK)
C
C       ========================================================
C       Purpose: Compute modified Bessel functions In(z), Kn(x)
C                and their derivatives for a complex argument
C       Input :  z --- Complex argument of In(z) and Kn(z)
C                n --- Order of In(z) and Kn(z)
C       Output:  CBI(n) --- In(z)
C                CDI(n) --- In'(z)
C                CBK(n) --- Kn(z)
C                CDK(n) --- Kn'(z)
C                NM --- Highest order computed
C       Routines called:
C             (1) CIK01 to compute I0(z), I1(z) K0(z) & K1(z)
C             (2) MSTA1 and MSTA2 to compute the starting
C                 point for backward recurrence
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,P,W,X,Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION CBI(0:N),CDI(0:N),CBK(0:N),CDK(0:N)
        A0=CDABS(Z)
        NM=N
        IF (A0.LT.1.0D-100) THEN
           DO 10 K=0,N
              CBI(K)=(0.0D0,0.0D0)
              CDI(K)=(0.0D0,0.0D0)
              CBK(K)=-(1.0D+300,0.0D0)
10            CDK(K)=(1.0D+300,0.0D0)
           CBI(0)=(1.0D0,0.0D0)
           CDI(1)=(0.5D0,0.0D0)
           RETURN
        ENDIF
        CALL CIK01(Z,CBI0,CDI0,CBI1,CDI1,CBK0,CDK0,CBK1,CDK1)
        CBI(0)=CBI0
        CBI(1)=CBI1
        CBK(0)=CBK0
        CBK(1)=CBK1
        CDI(0)=CDI0
        CDI(1)=CDI1
        CDK(0)=CDK0
        CDK(1)=CDK1
        IF (N.LE.1) RETURN
        M=MSTA1(A0,200)
        IF (M.LT.N) THEN
           NM=M
        ELSE
           M=MSTA2(A0,N,15)
        ENDIF
        CF2=(0.0D0,0.0D0)
        CF1=(1.0D-100,0.0D0)
        DO 45 K=M,0,-1
           CF=2.0D0*(K+1.0D0)/Z*CF1+CF2
           IF (K.LE.NM) CBI(K)=CF
           CF2=CF1
45         CF1=CF
        CS=CBI0/CF
        DO 50 K=0,NM
50         CBI(K)=CS*CBI(K)
        DO 60 K=2,NM
           IF (CDABS(CBI(K-1)).GT.CDABS(CBI(K-2))) THEN
              CKK=(1.0D0/Z-CBI(K)*CBK(K-1))/CBI(K-1)
           ELSE
              CKK=(CBI(K)*CBK(K-2)+2.0D0*(K-1.0D0)/(Z*Z))/CBI(K-2)
           ENDIF
60         CBK(K)=CKK
        DO 70 K=2,NM
           CDI(K)=CBI(K-1)-K/Z*CBI(K)
70         CDK(K)=-CBK(K-1)-K/Z*CBK(K)
        RETURN
        END



C       **********************************

        SUBROUTINE MTU12(KF,KC,M,Q,X,F1R,D1R,F2R,D2R)
C
C       ==============================================================
C       Purpose: Compute modified Mathieu functions of the first and
C                second kinds, Mcm(1)(2)(x,q) and Msm(1)(2)(x,q),
C                and their derivatives
C       Input:   KF --- Function code
C                       KF=1 for computing Mcm(x,q)
C                       KF=2 for computing Msm(x,q)
C                KC --- Function Code
C                       KC=1 for computing the first kind
C                       KC=2 for computing the second kind
C                            or Msm(2)(x,q) and Msm(2)'(x,q)
C                       KC=3 for computing both the first
C                            and second kinds
C                m  --- Order of Mathieu functions
C                q  --- Parameter of Mathieu functions ( q ≥ 0 )
C                x  --- Argument of Mathieu functions
C       Output:  F1R --- Mcm(1)(x,q) or Msm(1)(x,q)
C                D1R --- Derivative of Mcm(1)(x,q) or Msm(1)(x,q)
C                F2R --- Mcm(2)(x,q) or Msm(2)(x,q)
C                D2R --- Derivative of Mcm(2)(x,q) or Msm(2)(x,q)
C       Routines called:
C            (1) CVA2 for computing the characteristic values
C            (2) FCOEF for computing expansion coefficients
C            (3) JYNB for computing Jn(x), Yn(x) and their
C                derivatives
C       ==============================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION FG(251),BJ1(0:251),DJ1(0:251),BJ2(0:251),DJ2(0:251),
     &            BY1(0:251),DY1(0:251),BY2(0:251),DY2(0:251)
        EPS=1.0D-14
        IF (KF.EQ.1.AND.M.EQ.2*INT(M/2)) KD=1
        IF (KF.EQ.1.AND.M.NE.2*INT(M/2)) KD=2
        IF (KF.EQ.2.AND.M.NE.2*INT(M/2)) KD=3
        IF (KF.EQ.2.AND.M.EQ.2*INT(M/2)) KD=4
        CALL CVA2(KD,M,Q,A)
        IF (Q.LE.1.0D0) THEN
           QM=7.5+56.1*SQRT(Q)-134.7*Q+90.7*SQRT(Q)*Q
        ELSE
           QM=17.0+3.1*SQRT(Q)-.126*Q+.0037*SQRT(Q)*Q
        ENDIF
        KM=INT(QM+0.5*M)
        IF(KM.GE.251) THEN
           F1R=DNAN()
           D1R=DNAN()
           F2R=DNAN()
           D2R=DNAN()
           RETURN
        END IF
        CALL FCOEF(KD,M,Q,A,FG)
        IC=INT(M/2)+1
        IF (KD.EQ.4) IC=M/2
        C1=DEXP(-X)
        C2=DEXP(X)
        U1=DSQRT(Q)*C1
        U2=DSQRT(Q)*C2
        CALL JYNB(KM+1,U1,NM,BJ1,DJ1,BY1,DY1)
        CALL JYNB(KM+1,U2,NM,BJ2,DJ2,BY2,DY2)
        W1=0.0D0
        W2=0.0D0
        IF (KC.EQ.2) GO TO 50
        F1R=0.0D0
        DO 30 K=1,KM
           IF (KD.EQ.1) THEN
              F1R=F1R+(-1)**(IC+K)*FG(K)*BJ1(K-1)*BJ2(K-1)
           ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
              F1R=F1R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BJ2(K)
     &            +(-1)**KD*BJ1(K)*BJ2(K-1))
           ELSE
              F1R=F1R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BJ2(K+1)
     &            -BJ1(K+1)*BJ2(K-1))
           ENDIF
           IF (K.GE.5.AND.DABS(F1R-W1).LT.DABS(F1R)*EPS) GO TO 35
30         W1=F1R
35      F1R=F1R/FG(1)
        D1R=0.0D0
        DO 40 K=1,KM
           IF (KD.EQ.1) THEN
              D1R=D1R+(-1)**(IC+K)*FG(K)*(C2*BJ1(K-1)*DJ2(K-1)
     &            -C1*DJ1(K-1)*BJ2(K-1))
           ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
              D1R=D1R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DJ2(K)
     &            +(-1)**KD*BJ1(K)*DJ2(K-1))-C1*(DJ1(K-1)*BJ2(K)
     &            +(-1)**KD*DJ1(K)*BJ2(K-1)))
           ELSE
              D1R=D1R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DJ2(K+1)
     &            -BJ1(K+1)*DJ2(K-1))-C1*(DJ1(K-1)*BJ2(K+1)
     &            -DJ1(K+1)*BJ2(K-1)))
           ENDIF
           IF (K.GE.5.AND.DABS(D1R-W2).LT.DABS(D1R)*EPS) GO TO 45
40         W2=D1R
45      D1R=D1R*DSQRT(Q)/FG(1)
        IF (KC.EQ.1) RETURN
50      F2R=0.0D0
        DO 55 K=1,KM
           IF (KD.EQ.1) THEN
              F2R=F2R+(-1)**(IC+K)*FG(K)*BJ1(K-1)*BY2(K-1)
           ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
              F2R=F2R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BY2(K)
     &            +(-1)**KD*BJ1(K)*BY2(K-1))
           ELSE
              F2R=F2R+(-1)**(IC+K)*FG(K)*(BJ1(K-1)*BY2(K+1)
     &            -BJ1(K+1)*BY2(K-1))
           ENDIF
           IF (K.GE.5.AND.DABS(F2R-W1).LT.DABS(F2R)*EPS) GO TO 60
55         W1=F2R
60      F2R=F2R/FG(1)
        D2R=0.0D0
        DO 65 K=1,KM
           IF (KD.EQ.1) THEN
              D2R=D2R+(-1)**(IC+K)*FG(K)*(C2*BJ1(K-1)*DY2(K-1)
     &            -C1*DJ1(K-1)*BY2(K-1))
           ELSE IF (KD.EQ.2.OR.KD.EQ.3) THEN
              D2R=D2R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DY2(K)
     &            +(-1)**KD*BJ1(K)*DY2(K-1))-C1*(DJ1(K-1)*BY2(K)
     &            +(-1)**KD*DJ1(K)*BY2(K-1)))
           ELSE
              D2R=D2R+(-1)**(IC+K)*FG(K)*(C2*(BJ1(K-1)*DY2(K+1)
     &            -BJ1(K+1)*DY2(K-1))-C1*(DJ1(K-1)*BY2(K+1)
     &            -DJ1(K+1)*BY2(K-1)))
           ENDIF
           IF (K.GE.5.AND.DABS(D2R-W2).LT.DABS(D2R)*EPS) GO TO 70
65         W2=D2R
70         D2R=D2R*DSQRT(Q)/FG(1)
        RETURN
        END



C       **********************************

        SUBROUTINE CIK01(Z,CBI0,CDI0,CBI1,CDI1,CBK0,CDK0,CBK1,CDK1)
C
C       ==========================================================
C       Purpose: Compute modified Bessel functions I0(z), I1(z),
C                K0(z), K1(z), and their derivatives for a
C                complex argument
C       Input :  z --- Complex argument
C       Output:  CBI0 --- I0(z)
C                CDI0 --- I0'(z)
C                CBI1 --- I1(z)
C                CDI1 --- I1'(z)
C                CBK0 --- K0(z)
C                CDK0 --- K0'(z)
C                CBK1 --- K1(z)
C                CDK1 --- K1'(z)
C       ==========================================================
C
        IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
        IMPLICIT COMPLEX*16 (C,Z)
        DIMENSION A(12),B(12),A1(10)
        PI=3.141592653589793D0
        CI=(0.0D0,1.0D0)
        A0=CDABS(Z)
        Z2=Z*Z
        Z1=Z
        IF (A0.EQ.0.0D0) THEN
           CBI0=(1.0D0,0.0D0)
           CBI1=(0.0D0,0.0D0)
           CDI0=(0.0D0,0.0D0)
           CDI1=(0.5D0,0.0D0)
           CBK0=(1.0D+300,0.0D0)
           CBK1=(1.0D+300,0.0D0)
           CDK0=-(1.0D+300,0.0D0)
           CDK1=-(1.0D+300,0.0D0)
           RETURN
        ENDIF
        IF (DBLE(Z).LT.0.0) Z1=-Z
        IF (A0.LE.18.0) THEN
           CBI0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 10 K=1,50
              CR=0.25D0*CR*Z2/(K*K)
              CBI0=CBI0+CR
              IF (CDABS(CR/CBI0).LT.1.0D-15) GO TO 15
10         CONTINUE
15         CBI1=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO 20 K=1,50
              CR=0.25D0*CR*Z2/(K*(K+1))
              CBI1=CBI1+CR
              IF (CDABS(CR/CBI1).LT.1.0D-15) GO TO 25
20         CONTINUE
25         CBI1=0.5D0*Z1*CBI1
        ELSE
           DATA A/0.125D0,7.03125D-2,
     &            7.32421875D-2,1.1215209960938D-1,
     &            2.2710800170898D-1,5.7250142097473D-1,
     &            1.7277275025845D0,6.0740420012735D0,
     &            2.4380529699556D01,1.1001714026925D02,
     &            5.5133589612202D02,3.0380905109224D03/
           DATA B/-0.375D0,-1.171875D-1,
     &            -1.025390625D-1,-1.4419555664063D-1,
     &            -2.7757644653320D-1,-6.7659258842468D-1,
     &            -1.9935317337513D0,-6.8839142681099D0,
     &            -2.7248827311269D01,-1.2159789187654D02,
     &            -6.0384407670507D02,-3.3022722944809D03/
           K0=12
           IF (A0.GE.35.0) K0=9
           IF (A0.GE.50.0) K0=7
           CA=CDEXP(Z1)/CDSQRT(2.0D0*PI*Z1)
           CBI0=(1.0D0,0.0D0)
           ZR=1.0D0/Z1
           DO 30 K=1,K0
30            CBI0=CBI0+A(K)*ZR**K
           CBI0=CA*CBI0
           CBI1=(1.0D0,0.0D0)
           DO 35 K=1,K0
35            CBI1=CBI1+B(K)*ZR**K
           CBI1=CA*CBI1
        ENDIF
        IF (A0.LE.9.0) THEN
           CS=(0.0D0,0.0D0)
           CT=-CDLOG(0.5D0*Z1)-0.5772156649015329D0
           W0=0.0D0
           CR=(1.0D0,0.0D0)
           DO 40 K=1,50
              W0=W0+1.0D0/K
              CR=0.25D0*CR/(K*K)*Z2
              CS=CS+CR*(W0+CT)
              IF (CDABS((CS-CW)/CS).LT.1.0D-15) GO TO 45
40            CW=CS
45         CBK0=CT+CS
        ELSE
           DATA A1/0.125D0,0.2109375D0,
     &             1.0986328125D0,1.1775970458984D01,
     &             2.1461706161499D02,5.9511522710323D03,
     &             2.3347645606175D05,1.2312234987631D07,
     &             8.401390346421D08,7.2031420482627D10/
           CB=0.5D0/Z1
           ZR2=1.0D0/Z2
           CBK0=(1.0D0,0.0D0)
           DO 50 K=1,10
50            CBK0=CBK0+A1(K)*ZR2**K
           CBK0=CB*CBK0/CBI0
        ENDIF
        CBK1=(1.0D0/Z1-CBI1*CBK0)/CBI0
        IF (DBLE(Z).LT.0.0) THEN
           IF (DIMAG(Z).LT.0.0) CBK0=CBK0+CI*PI*CBI0
           IF (DIMAG(Z).GT.0.0) CBK0=CBK0-CI*PI*CBI0
           IF (DIMAG(Z).LT.0.0) CBK1=-CBK1+CI*PI*CBI1
           IF (DIMAG(Z).GT.0.0) CBK1=-CBK1-CI*PI*CBI1
           CBI1=-CBI1
        ENDIF
        CDI0=CBI1
        CDI1=CBI0-1.0D0/Z*CBI1
        CDK0=-CBK1
        CDK1=-CBK0-1.0D0/Z*CBK1
        RETURN
        END

C       **********************************

        SUBROUTINE CPSI(X,Y,PSR,PSI)
C
C       =============================================
C       Purpose: Compute the psi function for a
C                complex argument
C       Input :  x   --- Real part of z
C                y   --- Imaginary part of z
C       Output:  PSR --- Real part of psi(z)
C                PSI --- Imaginary part of psi(z)
C       =============================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(8)
        DATA A/-.8333333333333D-01,.83333333333333333D-02,
     &       -.39682539682539683D-02,.41666666666666667D-02,
     &       -.75757575757575758D-02,.21092796092796093D-01,
     &       -.83333333333333333D-01,.4432598039215686D0/
        PI=3.141592653589793D0
        IF (Y.EQ.0.0D0.AND.X.EQ.INT(X).AND.X.LE.0.0D0) THEN
           PSR=1.0D+300
           PSI=0.0D0
        ELSE
           X1=X
           Y1=Y
           IF (X.LT.0.0D0) THEN
              X=-X
              Y=-Y
           ENDIF
           X0=X
           N=0
           IF (X.LT.8.0D0) THEN
              N=8-INT(X)
              X0=X+N
           ENDIF
           TH=0.0D0
           IF (X0.EQ.0.0D0.AND.Y.NE.0.0D0) TH=0.5D0*PI
           IF (X0.NE.0.0D0) TH=DATAN(Y/X0)
           Z2=X0*X0+Y*Y
           Z0=DSQRT(Z2)
           PSR=DLOG(Z0)-0.5D0*X0/Z2
           PSI=TH+0.5D0*Y/Z2
           DO 10 K=1,8
              PSR=PSR+A(K)*Z2**(-K)*DCOS(2.0D0*K*TH)
10            PSI=PSI-A(K)*Z2**(-K)*DSIN(2.0D0*K*TH)
           IF (X.LT.8.0D0) THEN
              RR=0.0D0
              RI=0.0D0
              DO 20 K=1,N
                 RR=RR+(X0-K)/((X0-K)**2.0D0+Y*Y)
20               RI=RI+Y/((X0-K)**2.0D0+Y*Y)
              PSR=PSR-RR
              PSI=PSI+RI
           ENDIF
           IF (X1.LT.0.0D0) THEN
              TN=DTAN(PI*X)
              TM=DTANH(PI*Y)
              CT2=TN*TN+TM*TM
              PSR=PSR+X/(X*X+Y*Y)+PI*(TN-TN*TM*TM)/CT2
              PSI=PSI-Y/(X*X+Y*Y)-PI*TM*(1.0D0+TN*TN)/CT2
              X=X1
              Y=Y1
           ENDIF
        ENDIF
        RETURN
        END

C       **********************************

        SUBROUTINE SPHY(N,X,NM,SY,DY)
C
C       ======================================================
C       Purpose: Compute spherical Bessel functions yn(x) and
C                their derivatives
C       Input :  x --- Argument of yn(x) ( x ≥ 0 )
C                n --- Order of yn(x) ( n = 0,1,… )
C       Output:  SY(n) --- yn(x)
C                DY(n) --- yn'(x)
C                NM --- Highest order computed
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION SY(0:N),DY(0:N)
        NM=N
        IF (X.LT.1.0D-60) THEN
           DO 10 K=0,N
              SY(K)=-1.0D+300
10            DY(K)=1.0D+300
           RETURN
        ENDIF
        SY(0)=-DCOS(X)/X
        F0=SY(0)
        DY(0)=(DSIN(X)+DCOS(X)/X)/X
        IF (N.LT.1) THEN
           RETURN
        ENDIF
        SY(1)=(SY(0)-DSIN(X))/X
        F1=SY(1)
        DO 15 K=2,N
           F=(2.0D0*K-1.0D0)*F1/X-F0
           SY(K)=F
           IF (DABS(F).GE.1.0D+300) GO TO 20
           F0=F1
15         F1=F
20      NM=K-1
        DO 25 K=1,NM
25         DY(K)=SY(K-1)-(K+1.0D0)*SY(K)/X
        RETURN
        END



C       **********************************

        SUBROUTINE JELP(U,HK,ESN,ECN,EDN,EPH)
C
C       ========================================================
C       Purpose: Compute Jacobian elliptic functions sn u, cn u
C                and dn u
C       Input  : u   --- Argument of Jacobian elliptic fuctions
C                Hk  --- Modulus k ( 0 ≤ k ≤ 1 )
C       Output : ESN --- sn u
C                ECN --- cn u
C                EDN --- dn u
C                EPH --- phi ( in degrees )
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION R(40)
        PI=3.14159265358979D0
        A0=1.0D0
        B0=DSQRT(1.0D0-HK*HK)
        DO 10 N=1,40
           A=(A0+B0)/2.0D0
           B=DSQRT(A0*B0)
           C=(A0-B0)/2.0D0
           R(N)=C/A
           IF (C.LT.1.0D-7) GO TO 15
           A0=A
10         B0=B
15      DN=2.0D0**N*A*U
        D=0.0D0
        DO 20 J=N,1,-1
           T=R(J)*DSIN(DN)
           SA=DATAN(T/DSQRT(DABS(1.0D0-T*T)))
           D=.5D0*(DN+SA)
20         DN=D
        EPH=D*180.0D0/PI
        ESN=DSIN(D)
        ECN=DCOS(D)
        EDN=DSQRT(1.0D0-HK*HK*ESN*ESN)
        RETURN
        END

C       **********************************

        SUBROUTINE STVHV(V,X,HV)
C
C       =====================================================
C       Purpose: Compute Struve function Hv(x) with an
C                arbitrary order v
C       Input :  v  --- Order of Hv(x)  ( -8.0 ≤ v ≤ 12.5 )
C                x  --- Argument of Hv(x) ( x ≥ 0 )
C       Output:  HV --- Hv(x)
C       Note: numerically unstable away from the above range for `v`
C       Routine called: GAMMA2 to compute the gamma function
C       =====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PI=3.141592653589793D0
        IF (X.EQ.0.0D0) THEN
           IF (V.GT.-1.0.OR.INT(V)-V.EQ.0.5D0) THEN
              HV=0.0D0
           ELSE IF (V.LT.-1.0D0) THEN
              HV=(-1)**(INT(0.5D0-V)-1)*1.0D+300
           ELSE IF (V.EQ.-1.0D0) THEN
              HV=2.0D0/PI
           ENDIF
           RETURN
        ENDIF
        BYV=0.0D0
        BF=0.0D0
        QU0=0.0D0
        PU0=0.0D0
        IF (X.LE.20.0D0) THEN
C          Power series for Hv (Abramowitz & Stegun 12.1.3)
           V0=V+1.5D0
           CALL GAMMA2(V0,GA)
           S=2.0D0/(DSQRT(PI)*GA)
           R1=1.0D0
           DO 10 K=1,100
              VA=K+1.5D0
              CALL GAMMA2(VA,GA)
              VB=V+K+1.5D0
              CALL GAMMA2(VB,GB)
              R1=-R1*(0.5D0*X)**2
              R2=R1/(GA*GB)
              S=S+R2
              IF (DABS(R2).LT.DABS(S)*1.0D-12) GO TO 15
10         CONTINUE
15         HV=(0.5D0*X)**(V+1.0D0)*S
        ELSE
C          Asymptotic large |z| expansion for Hv - Yv  (Abm & Stg 12.1.29)
           SA=(0.5D0*X)**(V-1.0)/PI
           V0=V+0.5D0
           CALL GAMMA2(V0,GA)
           S=DSQRT(PI)/GA
           R1=1.0D0
           DO 20 K=1,12
              VA=K+0.5D0
              CALL GAMMA2(VA,GA)
              VB=-K+V+0.5D0
              CALL GAMMA2(VB,GB)
              R1=R1/(0.5D0*X)**2
              S=S+R1*GA/GB
20         CONTINUE
           S0=SA*S

C          Compute Y_(|v|-N)   (Abm & Stg 9.2.6)
           U=DABS(V)
           N=INT(U)
           U0=U-N
           DO 35 L=0,1
              VT=4.0D0*(U0+L)**2
              R1=1.0D0
              PU1=1.0D0
              DO 25 K=1,12
                 R1=-0.0078125D0*R1*(VT-(4.0*K-3.0D0)**2)*
     &             (VT-(4.0D0*K-1.0)**2)/((2.0D0*K-1.0)*K*X*X)
                 PU1=PU1+R1
25            CONTINUE
              QU1=1.0D0
              R2=1.0D0
              DO 30 K=1,12
                 R2=-0.0078125D0*R2*(VT-(4.0D0*K-1.0)**2)*
     &             (VT-(4.0D0*K+1.0)**2)/((2.0D0*K+1.0)*K*X*X)
                 QU1=QU1+R2
30            CONTINUE
              QU1=0.125D0*(VT-1.0D0)/X*QU1
              IF (L.EQ.0) THEN
                 PU0=PU1
                 QU0=QU1
              ENDIF
35         CONTINUE
           T0=X-(0.5*U0+0.25D0)*PI
           T1=X-(0.5*U0+0.75D0)*PI
           SR=DSQRT(2.0D0/(PI*X))
           BY0=SR*(PU0*DSIN(T0)+QU0*DCOS(T0))
           BY1=SR*(PU1*DSIN(T1)+QU1*DCOS(T1))

C          Compute Y_|v|   (Abm & Stg 9.1.27)
           BF0=BY0
           BF1=BY1
           DO 40 K=2,N
              BF=2.0D0*(K-1.0+U0)/X*BF1-BF0
              BF0=BF1
40            BF1=BF
           IF (N.EQ.0) BYV=BY0
           IF (N.EQ.1) BYV=BY1
           IF (N.GT.1) BYV=BF

C          Compute Y_v  (handle the case v < 0 appropriately)
           IF (V .LT. 0) THEN
              IF (U0 .EQ. 0) THEN
C                Use symmetry (Abm & Stg 9.1.5)
                 BYV=(-1)**N*BYV
              ELSE
C                Use relation between Yv & Jv (Abm & Stg 9.1.6)

C                Compute J_(|v|-N) (Abm & Stg 9.2.5)
                 BJ0=SR*(PU0*DCOS(T0)-QU0*DSIN(T0))
                 BJ1=SR*(PU1*DCOS(T1)-QU1*DSIN(T1))
C                Forward recurrence for J_|v| (Abm & Stg 9.1.27)
C                It's OK for the limited range -8.0 ≤ v ≤ 12.5,
C                since x >= 20 here; but would be unstable for v <~ -20
                 BF0=BJ0
                 BF1=BJ1
                 DO 50 K=2,N
                    BF=2.0D0*(K-1.0+U0)/X*BF1-BF0
                    BF0=BF1
50                  BF1=BF
                 IF (N.EQ.0) BJV=BJ0
                 IF (N.EQ.1) BJV=BJ1
                 IF (N.GT.1) BJV=BF

C                Compute Y_v    (Abm & Stg 9.1.6)
                 BYV = DCOS(V*PI)*BYV + DSIN(-V*PI)*BJV
              END IF
           END IF

C          Compute H_v
           HV=BYV+S0
        ENDIF
        RETURN
        END



C       **********************************
