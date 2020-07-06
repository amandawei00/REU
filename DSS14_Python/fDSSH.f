c Example 
c      PROGRAM DSS14EXAMPLE
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c      INTEGER FINI14
c
c      COMMON / FRAGINI14 / FINI14
c
c      Q2=10.d0
c
c      write(6,*)  '    z        Q2      zD_utot    zD_g    relerr-u  rel  
c     1err-g'
c
c
c      do i=0,100
c         z=0.05+0.95d0*i/101.d0
c         uerr=0.d0
c         gerr=0.d0
c
cc   set 
c        IS=0      ! DSS14 BEST FIT  
cc        IS=1     ! POSITIVE SHIFT HESSIAN EIGENVECTOR #1 
cc        IS=-1    ! NEGATIVE SHIFT HESSIAN EIGENVECTOR #1  
c
cc   order 
c        IO=1 ! NLO       
cc   charge
c        IC=1 ! charge +        
c
cc      calls best fit 
c        fINI14=0
c        call fDSSH(IS,IL,IC,IO,Z,Q2,U1,UB1,D1,DB1,S1,SB1,C1,B1,GL1)
c 
c 
cc    % CL
c        IL=0    ! 90 % CL
cc         IL=1    ! 68 % CL    
c 
c 
c 
c
cc      loops over eigenvector sets
c        do j=1,28
c
c           fINI14=0               !resets hadron-type/order/set        
c           IS=j   
c           call fDSSH(IS,IL,IC,IO,Z,Q2,U2,UB2,D2,DB2,S2,SB2,C2,B2,GL2)
c           fINI14=0               !resets hadron-type/order/set        
c           IS=-j           
c           call fDSSH(IS,IL,IC,IO,Z,Q2,U3,UB3,D3,DB3,S3,SB3,C3,B3,GL3)       
c           uerr=uerr+(U2+UB2-U3-UB3)**2
c           gerr=gerr+(GL2-GL3)**2        
c        enddo      
c        uerr=0.5d0*DSQRT(uerr)
c        gerr=0.5d0*DSQRT(gerr)
c       
cc        write(6,90) z,q2,u1+ub1,gl1,uerr/(u1+ub1),gerr/gl1 
c        write(6,90) z, gerr/gl1, -gerr/gl1
cc        write(6,90) z, gl1, gl1+gerr,gl1-gerr
c 90     FORMAT (10(1PE12.3))
c
c       enddo
c       end



*************************************************************************
*                                                                                     
*        DSS14  PION  HESSIAN SETS   (NLO only)            
*                                                                                                                                                                           
*  D.de Florian, M. Epele, R. Hernadez Pinto R.Sassot, M.Stratmann         
*                                                  arXiv:1410.6027               
*                                                                                                                 
*     CALL fDSSH (IS,IL,IC,IO, X, Q2, U, UB, D, DB, S, SB, C, B, GL)               
*                                                                                                                 
*  INPUT:    
*  IS = SET VERSION    0: DSS14  BEST FIT
*                      +/- 1, +/- 2, ... , +/- 28:  DSS14 VARIATIONS
*
*  IL = CL    0: 90% CL
*             1: 68% CL 
*                                                                                                    
*  IC = Hadron Charge  0: 0 (as average of + and -)                    
*                      1: +                                        
*                     -1: -                                        
*                                                                  
*  IO= Order       1: NLO                          
*                                                              
*                                                                  
*             X                    (between  0.05   and  1.0)       
*            Q2 = scale in GeV**2 (between  1.0    and  1.D5)      
*             (for values outside the allowed range the program    
*              writes a warning and extrapolates to the x and      
*              Q2 values requested)                                
*                                                                  
*   OUTPUT: U, UB, D, DB, S, SB,   C,           B,       GL        
*           U Ubar D Dbar S Sbar Charm=Cbar Bottom=Bbar Gluon      
*           Always X times the distribution is returned            
*                                                                  
*                                                                  
*   COMMON:  The main program or the calling routine has to have   
*            a common block  COMMON / FRAGINI14 / FINI14 , and  FINI14   
*            has always to be zero when DSS is called for the      
*            first time or when the SET has been changed.          
*                                                                  
********************************************************************

      SUBROUTINE fDSSH(IS,IL,IC,IO, X, Q2, U, UB, D, DB, S, SB,C,B,GL)
Cf2py intent(in) is
Cf2py intent(in) il
Cf2py intent(in) ic
Cf2py intent(in) io
Cf2py intent(in) x
Cf2py intent(in) q2
Cf2py intent(out) gl
Cf2py intent(out) u
Cf2py intent(out) ub
Cf2py intent(out) d
Cf2py intent(out) db
Cf2py intent(out) s
Cf2py intent(out) sb
Cf2py intent(out) c
Cf2py intent(out) b
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IS
      INTEGER FINI14,IL
      CHARACTER(len=255) root
      PARAMETER (NPART=9, NX=47, NQ=24, NARG=2)
      DIMENSION XUTOTF(NX,NQ), XDTOTF(NX,NQ), XSTOTF(NX,NQ)
      DIMENSION XUVALF(NX,NQ), XDVALF(NX,NQ), XSVALF(NX,NQ)
      DIMENSION XCTOTF(NX,NQ), XBTOTF(NX,NQ)
      DIMENSION XGF(NX,NQ), PARTON (NPART,NQ,NX-1)
      DIMENSION QS(NQ), XB(NX), XT(NARG), NA(NARG), ARRF(NX+NQ) 
      COMMON /root/root
      
      character *10 fname1(4)
      character *2  fname2(40)
      character *5  fname3
      character *20 fname,file

      DATA fname1 / '90pi-plus-','90pi-minu-','68pi-plus-','68pi-minu-'/
      DATA fname2 / '01','02','03','04','05','06','07','08','09','10',
     1     '11','12','13','14','15','16','17','18','19','20','21','22',
     1     '23','24','25','26','27','28','29','30','31','32','33','34',
     1     '35','36','37','38','39','40'/
      
      COMMON / FRAGINI14 / FINI14
      SAVE XUTOTF, XDTOTF, XSTOTF, XCTOTF, XBTOTF, XGF, NA, ARRF
      SAVE XUVALF, XDVALF, XSVALF
*...BJORKEN-X AND Q**2 VALUES OF THE GRID :
       DATA QS / 1.d0, 1.25D0, 1.5D0, 2.5D0, 
     1           4.0D0, 6.4D0, 1.0D1, 1.5D1, 2.5D1, 4.0D1, 6.4D1,
     2           1.0D2, 1.8D2, 3.2D2, 5.8D2, 1.0D3, 1.8D3,
     3           3.2D3, 5.8D3, 1.0D4, 1.8D4, 3.2D4, 5.8D4, 1.0D5/
       DATA XB /0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
     4        0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     5        0.3, 0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475,  0.5, 
     6        0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675, 0.7,
     7        0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 
     8        0.925, 0.95, 0.975, 1.0/     
     
*...CHECK OF X AND Q2 VALUES : 
       IF ( (X.LT.0.05D0) .OR. (X.GT.1.0D0) ) THEN
           WRITE(6,91) 
  91       FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE')
C          STOP
       ENDIF
       IF ( (Q2.LT.1.D0) .OR. (Q2.GT.1.D5) ) THEN
           WRITE(6,92) 
  92       FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE')
       ENDIF
*...INITIALIZATION :
*    SELECTION AND READING OF THE GRID :
       IF (FINI14.NE.0) GOTO 16

      fname3='.grid'

       IF (IS.EQ.0) THEN
       IIREAD=10
       OPEN(IIREAD,FILE='grids/PI_DSS14.GRID')   
       ELSE
       IIREAD=11          
       IF(IS.GT.0) THEN
       IBAND=1+2*IL
       KBAND=IS
       ELSE
       IBAND=2+2*IL
       KBAND=-IS
      
       ENDIF
            call strcat(fname1(IBAND),fname2(KBAND),fname)
            call strcat(fname,fname3,file)
            
            OPEN(IIREAD,FILE=FILE)  !PARAMETER-FILE            
      END IF


C
       DO 15 M = 1, NX-1 
       DO 15 N = 1, NQ
       READ(IIREAD,90) PARTON(1,N,M), PARTON(2,N,M), PARTON(3,N,M), 
     1                 PARTON(4,N,M), PARTON(5,N,M), PARTON(6,N,M),
     2                 PARTON(7,N,M), PARTON(8,N,M), PARTON(9,N,M)
  90   FORMAT (9(1PE10.3))
  15   CONTINUE
       CLOSE(IIREAD)
C
      FINI14 = 1
*....ARRAYS FOR THE INTERPOLATION SUBROUTINE :
      DO 10 IQ = 1, NQ
      DO 20 IX = 1, NX-1
        XB0 = XB(IX) 
        XB1 = 1.D0-XB(IX)
        XUTOTF(IX,IQ) = PARTON(1,IQ,IX) / (XB1**4 * XB0**0.5)
        XDTOTF(IX,IQ) = PARTON(2,IQ,IX) / (XB1**4 * XB0**0.5)
        XSTOTF(IX,IQ) = PARTON(3,IQ,IX) / (XB1**4 * XB0**0.5) 
        XCTOTF(IX,IQ) = PARTON(4,IQ,IX) / (XB1**7 * XB0**0.3) 
        XBTOTF(IX,IQ) = PARTON(5,IQ,IX) / (XB1**7 * XB0**0.3)
        XGF(IX,IQ)    = PARTON(6,IQ,IX) / (XB1**4 * XB0**0.3)
        XUVALF(IX,IQ) = PARTON(7,IQ,IX) / (XB1**4 * XB0**0.5)
        XDVALF(IX,IQ) = PARTON(8,IQ,IX) / (XB1**4 * XB0**0.5)
        XSVALF(IX,IQ) = PARTON(9,IQ,IX) / (XB1**4 * XB0**0.5)
  20  CONTINUE
        XUTOTF(NX,IQ) = 0.D0
        XDTOTF(NX,IQ) = 0.D0
        XSTOTF(NX,IQ) = 0.D0
        XCTOTF(NX,IQ) = 0.D0
        XBTOTF(NX,IQ) = 0.D0
        XGF(NX,IQ)    = 0.D0
        XUVALF(NX,IQ) = 0.D0
        XDVALF(NX,IQ) = 0.D0
        XSVALF(NX,IQ) = 0.D0
  10  CONTINUE  
      NA(1) = NX
      NA(2) = NQ
      DO 30 IX = 1, NX
        ARRF(IX) = DLOG(XB(IX))
  30  CONTINUE
      DO 40 IQ = 1, NQ
        ARRF(NX+IQ) = DLOG(QS(IQ))
  40  CONTINUE
  16  CONTINUE
*...INTERPOLATION :
      XT(1) = DLOG(X)
      XT(2) = DLOG(Q2)
      UTOT = FINT14(NARG,XT,NA,ARRF,XUTOTF) * (1.D0-X)**4 * X**0.5
      DTOT = FINT14(NARG,XT,NA,ARRF,XDTOTF) * (1.D0-X)**4 * X**0.5 
      STOT = FINT14(NARG,XT,NA,ARRF,XSTOTF) * (1.D0-X)**4 * X**0.5
      CTOT = FINT14(NARG,XT,NA,ARRF,XCTOTF) * (1.D0-X)**7 * X**0.3
      BTOT = FINT14(NARG,XT,NA,ARRF,XBTOTF) * (1.D0-X)**7 * X**0.3
      GL   = FINT14(NARG,XT,NA,ARRF,XGF)    * (1.D0-X)**4 * X**0.3
      UVAL = FINT14(NARG,XT,NA,ARRF,XUVALF) * (1.D0-X)**4 * X**0.5
      DVAL = FINT14(NARG,XT,NA,ARRF,XDVALF) * (1.D0-X)**4 * X**0.5 
      SVAL = FINT14(NARG,XT,NA,ARRF,XSVALF) * (1.D0-X)**4 * X**0.5
       
       Up  = (UTOT+UVAL)/2.
       UBp = (UTOT-UVAL)/2.
       Dp  = (DTOT+DVAL)/2.
       DBp = (DTOT-DVAL)/2.
       Sp  = (STOT+SVAL)/2.
       SBp = (STOT-SVAL)/2.
       Cp  =  CTOT/2.
       Bp  =  BTOT/2.
              
       IF (IC.EQ.1) THEN
       U  = Up
       UB = UBp
       D  = Dp 
       DB = DBp
       S  = Sp
       SB = SBp
       C  = Cp
       B  = Bp
       ELSEIF (IC.EQ.-1) THEN
       U  = UBp
       UB = Up
       D  = DBp 
       DB = Dp
       S  = SBp
       SB = Sp
       C  = Cp
       B  = Bp
       ELSEIF (IC.EQ.0) THEN
       U  = (UBp+Up)/2.
       UB =  U
       D  = (DBp+Dp)/2. 
       DB =  D
       S  = (SBp+Sp)/2.
       SB =  S
       C  =  Cp
       B  =  Bp 
       ELSE
         WRITE(6,94)
 94      FORMAT (2X,' WRONG CHARGE')
         STOP
       END IF 
 60   RETURN
       END
*
*...CERN LIBRARY ROUTINE E104 (INTERPOLATION) :
*
      FUNCTION FINT14(NARG,ARG,NENT,ENT,TABLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARG(2),NENT(2),ENT(59),TABLE(840)
      DIMENSION D(5),NCOMB(5),IENT(5)
      KD=1
      M=1
      JA=1
         DO 5 I=1,NARG
      NCOMB(I)=1
      JB=JA-1+NENT(I)
         DO 2 J=JA,JB
      IF (ARG(I).LE.ENT(J)) GO TO 3
    2 CONTINUE
      J=JB
    3 IF (J.NE.JA) GO TO 4
      J=J+1
    4 JR=J-1
      D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR))
      IENT(I)=J-JA
      KD=KD+IENT(I)*M
      M=M*NENT(I)
    5 JA=JB+1
      FINT14=0.D0
   10 FAC=1.D0
      IADR=KD
      IFADR=1
         DO 15 I=1,NARG
      IF (NCOMB(I).EQ.0) GO TO 12
      FAC=FAC*(1.D0-D(I))
      GO TO 15
   12 FAC=FAC*D(I)
      IADR=IADR-IFADR
   15 IFADR=IFADR*NENT(I)
      FINT14=FINT14+FAC*TABLE(IADR)
      IL=NARG
   40 IF (NCOMB(IL).EQ.0) GO TO 80
      NCOMB(IL)=0
      IF (IL.EQ.NARG) GO TO 10
      IL=IL+1
         DO 50  K=IL,NARG
   50 NCOMB(K)=1
      GO TO 10
   80 IL=IL-1
      IF(IL.NE.0) GO TO 40
      RETURN
      END
      
c----------------------------------------------------------------------------------------
      subroutine strcat(str1,str2,str)
Cf2py intent(in) str1
Cf2py intent(in) str2
Cf2py intent(out) str
c concatenates str1 and str2 into str. Ignores trailing blanks of str1,str2
      character *(*) str1,str2,str
      l1=istrl(str1)
      l2=istrl(str2)
      l =len(str)
      if(l.lt.l1+l2) then
          write(*,*) 'error: l1+l2>l in strcat'
          write(*,*) 'l1=',l1,' str1=',str1
          write(*,*) 'l2=',l2,' str2=',str2
          write(*,*) 'l=',l
          stop
      endif
      if(l1.ne.0) str(1:l1)=str1(1:l1)
      if(l2.ne.0) str(l1+1:l1+l2)=str2(1:l2)
      if(l1+l2+1.le.l) str(l1+l2+1:l)= ' '
      end
c
      function istrl(string)
c returns the position of the last non-blank character in string
      character * (*) string
      i = len(string)
      dowhile(i.gt.0.and.string(i:i).eq.' ')
         i=i-1
      enddo
      istrl = i
      end
c----------------------------------------------------------------------------------------
