      PROGRAM MASTER
      IMPLICIT NONE
      CHARACTER*72 DSS07_PI_DIR, DSS07_PI_GRIDS
      CHARACTER*72 DSS07_KA_DIR, DSS07_KA_GRIDS
      CHARACTER*72 DSS07_PR_DIR, DSS07_PR_GRIDS
      CHARACTER*72 DSS07_HA_DIR, DSS07_HA_GRIDS
      LOGICAL DIRECTORY_EXISTS
      INTEGER I_PI, I_KA, I_PR, I_HA
      INTEGER NUM_XVALUES
      INTEGER NUM_QVALUES
      INTEGER I, J
      CHARACTER*72 MKDIR
      REAL*8 X(35)
      REAL*8 QS(24),Q(24)
      INTEGER PARTICLE_ID(-5:21)
      INTEGER IH,IC,IO
      REAL*8 xx,QQ
      REAL*8 GL,U,UB,D,DB,S,SB,C,CB,B,BB

C-----CHOOSE GRIDS TO EXPORT
      I_PI = 1
      I_KA = 1
      I_PR = 1
      I_HA = 1

C      DATA QS / 0.100000E+01, 0.111803E+01, 0.122474E+01, 0.158114E+01,
C     1 0.200000E+01, 0.252982E+01, 0.316228E+01, 0.387298E+01, 0.500000E+01,
C     2 0.632456E+01, 0.800000E+01, 0.141421E+02, 0.134164E+02, 0.178885E+02,
C     3 0.240832E+02, 0.316228E+02, 0.424264E+02, 0.979796E+02, 0.761577E+02,
C     4 0.100000E+03, 0.134164E+03, 0.178885E+03, 0.240832E+03, 0.316228E+03 /
                       
      DATA QS / 1.D0, 1.11803, 1.22474, 1.58114, 2., 2.52982, 3.16228,
     1  3.87298, 5., 6.32456, 8., 14.1421, 13.4164, 17.8885, 24.0832,
     2  31.6228, 42.4264, 97.9796, 76.1577, 100., 134.164, 178.885,
     3  240.832, 316.228/  

        DO I = 1, 24
        Q(I) = SQRT(QS(I))
       END DO

C       DATA X / 0.100000E-01, 0.200000E-01, 0.300000E-01, 0.400000E-01,
C     1 0.500000E-01, 0.600000E-01, 0.700000E-01, 0.800000E-01, 0.900000E-01,
C     2 0.950000E-01, 0.100000E+00, 0.125000E+00, 0.150000E+00, 0.175000E+00,
C     3 0.200000E+00, 0.225000E+00, 0.250000E+00, 0.275000E+00, 0.300000E+00,
C     4 0.325000E+00, 0.350000E+00, 0.375000E+00, 0.400000E+00, 0.450000E+00,
C     5 0.500000E+00, 0.550000E+00, 0.600000E+00, 0.650000E+00, 0.700000E+00,
C     6 0.750000E+00, 0.800000E+00, 0.850000E+00, 0.900000E+00, 0.930000E+00,
C     7 0.100000E+01 /

      DATA X /0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
     1  0.095, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275,
     2  0.3, 0.325, 0.35, 0.375, 0.4, 0.45, 0.5, 0.55,
     3  0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.93, 1.0/

C-----PARTICLE ID NUMBERS: Based on Standard PDG ID Numbers:
C------https://pdg.lbl.gov/2006/reviews/pdf-files/montecarlo-web.pdf
      PARTICLE_ID( 1) =  1  ! d
      PARTICLE_ID(-1) = -1  ! dbar
      PARTICLE_ID( 2) =  2  ! u
      PARTICLE_ID(-2) = -2  ! ubar
      PARTICLE_ID( 3) =  3  ! s
      PARTICLE_ID(-3) = -3  ! sbar
      PARTICLE_ID( 4) =  4  ! c
      PARTICLE_ID(-4) = -4  ! cbar
      PARTICLE_ID( 5) =  5  ! b
      PARTICLE_ID(-5) = -5  ! bbar
      PARTICLE_ID(21) = 21  ! g

*-----Set Grid Directory Names
      DSS07_PI_DIR = 'DSS07PI/'
      DSS07_KA_DIR = 'DSS07PI/'
      DSS07_PR_DIR = 'DSS07PI/'
      DSS07_HA_DIR = 'DSS07PI/'

*-----GRID FILE NAMES:
      DSS07_PI_GRIDS = 'PILO.GRID'
      DSS07_KA_GRIDS = 'KALO.GRID'
      DSS07_PR_GRIDS = 'PROLO.GRID'
      DSS07_HA_GRIDS = 'HALO.GRID'


*-----Create Directories if they do not exist
      MKDIR = 'mkdir'

C------PION
      IF (I_PI.eq.1) THEN
        INQUIRE(FILE=DSS07_PI_DIR, EXIST=DIRECTORY_EXISTS)
        IF (.not.DIRECTORY_EXISTS) THEN
            CALL SYSTEM(MKDIR//DSS07_PI_DIR)
        ENDIF
*-------Generate LHAPDF Grids
        OPEN(UNIT = 4, FILE = trim(DSS07_PI_DIR)//DSS07_PI_GRIDS)
        WRITE(4,*) 'PdfType: central'
        WRITE(4,*) 'Format: lhagrid1'
        WRITE(4,*) '---'
        WRITE(4,100) X
        WRITE(4,101) Q
        WRITE(4,102) PARTICLE_ID(-5:-1),PARTICLE_ID(1:5),PARTICLE_ID(21)


C-------PION
        IH = 1
        IC = 1
        IO = 1

        DO I = 1, 35
          xx = x(I)
          DO J = 1, 24
          QQ = QS(J)
            CALL fDSS (IH,IC,IO, xx, QQ, GL,U,UB,D,DB,S,SB,C,CB,B,BB)
            WRITE(4,103) BB,CB,SB,UB,DB,D,U,
     +                   S ,C,B,GL
          ENDDO
        ENDDO


        WRITE(4,*) '---'
        CLOSE(4)
      ENDIF


*------FORMATTING
  100 FORMAT(81(E12.6E2,' '))
  101 FORMAT(31(E12.6E2,' '))
  102 FORMAT(' ',5(I2,' '),5(I1,' '),I2)
  103 FORMAT(11(E12.6E2,'  '))
  202 FORMAT(' ',5(I2,' '),5(I1,' '),I2)



      END PROGRAM MASTER
