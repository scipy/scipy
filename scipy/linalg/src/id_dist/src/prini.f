C 
C 
C 
C 
        SUBROUTINE PRINI(IP1,IQ1)
        save
        CHARACTER *1 MES(1), AA(1)
        REAL *4 A(1)
        REAL *8 A2(1)
        REAL *8 A4(1)
        INTEGER *4 IA(1)
        INTEGER *2 IA2(1)
        IP=IP1
        IQ=IQ1
  
        RETURN
  
C 
C 
C 
C 
C 
        ENTRY PRIN(MES,A,N)
        CALL  MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1200)(A(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1200)(A(J),J=1,N)
 1200 FORMAT(6(2X,E11.5))
         RETURN
C 
C 
C 
C 
        ENTRY PRIN2(MES,A2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1400)(A2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1400)(A2(J),J=1,N)
 1400 FORMAT(6(2X,E11.5))
        RETURN
C 
C 
C 
C 
        ENTRY PRIN2_long(MES,A2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1450)(A2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1450)(A2(J),J=1,N)
 1450 FORMAT(2(2X,E22.16))
        RETURN
C 
C 
C 
C 
        ENTRY PRINQ(MES,A4,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1500)(A4(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1500)(A4(J),J=1,N)
 1500 FORMAT(6(2X,e11.5))
        RETURN
C 
C 
C 
C 
        ENTRY PRINF(MES,IA,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA(J),J=1,N)
 1600 FORMAT(10(1X,I7))
        RETURN
C 
C 
C 
C 
        ENTRY PRINF2(MES,IA2,N)
        CALL MESSPR(MES,IP,IQ)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,1600)(IA2(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,1600)(IA2(J),J=1,N)
        RETURN
C 
C 
C 
C 
        ENTRY PRINA(MES,AA,N)
        CALL MESSPR(MES,IP,IQ)
 2000 FORMAT(1X,80A1)
        IF(IP.NE.0 .AND. N.NE.0) WRITE(IP,2000)(AA(J),J=1,N)
        IF(IQ.NE.0 .AND. N.NE.0) WRITE(IQ,2000)(AA(J),J=1,N)
        RETURN
        END
c 
c 
c 
c 
c 
        SUBROUTINE MESSPR(MES,IP,IQ)
        save
        CHARACTER *1 MES(1),AST
        DATA AST/'*'/
C 
C         DETERMINE THE LENGTH OF THE MESSAGE
C 
        I1=0
        DO 1400 I=1,10000
        IF(MES(I).EQ.AST) GOTO 1600
        I1=I
 1400 CONTINUE
 1600 CONTINUE
         IF ( (I1.NE.0) .AND. (IP.NE.0) )
     1     WRITE(IP,1800) (MES(I),I=1,I1)
         IF ( (I1.NE.0) .AND. (IQ.NE.0) )
     1     WRITE(IQ,1800) (MES(I),I=1,I1)
 1800 FORMAT(1X,80A1)
         RETURN
         END
C 
C 
C 
C 
C 
        SUBROUTINE ZTIME(I)
        save
        J=1
        J=7-I+J
CCCC    I=MRUN(J)
        RETURN
        END
c 
c 
c 
c 
c 
        subroutine msgmerge(a,b,c)
        save
        character *1 a(1),b(1),c(1),ast
        data ast/'*'/
c 
        do 1200 i=1,1000
c 
        if(a(i) .eq. ast) goto 1400
        c(i)=a(i)
        iadd=i
 1200 continue
c 
 1400 continue
c 
        do 1800 i=1,1000
c 
        c(iadd+i)=b(i)
        if(b(i) .eq. ast) return
 1800 continue
        return
        end
c 
c 
c 
c 
c 
  
        subroutine fileflush(iw)
        implicit real *8 (a-h,o-z)
c 
        save
        close(iw)
        open(iw,status='old')
        do 1400 i=1,1000000
c 
        read(iw,1200,end=1600)
 1200 format(1a1)
 1400 continue
 1600 continue
c 
        return
        end
  
  
c 
c 
c 
c 
c 
        subroutine mach_zero(zero_mach)
        implicit real *8 (a-h,o-z)
        save
c
        zero_mach=100       
c
        d1=1.1
        d3=1.1
        d=1.11
        do 1200 i=1,1000
c

        d=d/2
        d2=d1+d
        call mach_zero0(d2,d3,d4)
c
        if(d4 .eq. 0) goto 1400
c
 1200 continue
 1400 continue
c
        zero_mach=d
        return
        end

c 
c 
c 
c 
c 
        subroutine mach_zero0(a,b,c)
        implicit real *8 (a-h,o-z)
        save
c
        c=b-a

        return
        end
