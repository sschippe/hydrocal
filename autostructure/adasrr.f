C N. R. BADNELL     PROGRAM ADASRR       NRB v3.12              21/11/22
C
C***********************************************************************
C
C         POST-PROCESSOR FOR  ** AUTOSTRUCTURE ** (ADAS ONLY)
C         ***************************************************
C
C CALCULATES CA/LS/IC (+ HYBRID) PARTIAL RR RATE COEFFS IN ADF48 FORMAT
C
C***********************************************************************
cparc                                                               !par
cparc                    + Parallel +                               !par
cparc                                                               !par
cparc   Only useful for totals since gives separate adf48 per proc  !par
cparc                                                               !par
cpar!***************************************************************!par
cparc                                                               !par
cpar      module comm_interface                                     !par
cparc                                                               !par
cpar      use mpi                                                   !par
cparc                                                               !par
cpar      implicit none                                             !par
cparc                                                               !par
cpar      public comm_init          ! Initialize MPI                !par
cpar      public comm_barrier       ! MPI barrier                   !par
cpar      public comm_finalize      ! Terminate MPI                 !par
cpar      integer*4, public  :: iam                                 !par
cpar      integer*4, public  :: nproc                               !par
cparc                                                               !par
cpar      SAVE                                                      !par
cparc                                                               !par
cpar      private                                                   !par
cpar      integer*4 :: mpicom                                       !par
cparc                                                               !par
cpar      CONTAINS                                                  !par
cparc                                                               !par
cpar!---------------------------------------------------------------!par
cpar      subroutine comm_init()                                    !par
cparc                                                               !par
cpar      implicit none                                             !par
cparc                                                               !par
cpar      integer*4 :: ier                                          !par
cparc                                                               !par
cpar      mpicom = MPI_COMM_WORLD                                   !par
cparc                                                               !par
cpar      call mpi_init(ier)                                        !par
cpar      call mpi_comm_rank(mpicom, iam, ier)                      !par
cpar      call mpi_comm_size(mpicom, nproc, ier)                    !par
cparc                                                               !par
cpar      return                                                    !par
cparc                                                               !par
cpar      end subroutine comm_init                                  !par
cparc                                                               !par
cpar!---------------------------------------------------------------!par
cpar      subroutine comm_barrier()                                 !par
cparc                                                               !par
cpar      implicit none                                             !par
cparc                                                               !par
cpar      integer*4 :: ier                                          !par
cparc                                                               !par
cpar      call mpi_barrier(mpicom, ier)                             !par
cparc                                                               !par
cpar      return                                                    !par
cparc                                                               !par
cpar      end subroutine comm_barrier                               !par
cpar!---------------------------------------------------------------!par
cparc                                                               !par
cpar      subroutine comm_finalize()                                !par
cparc                                                               !par
cpar      implicit none                                             !par
cparc                                                               !par
cpar      integer*4 :: ier                                          !par
cparc                                                               !par
cpar      call mpi_finalize(ier)                                    !par
cparc                                                               !par
cpar      return                                                    !par
cparc                                                               !par
cpar      end subroutine comm_finalize                              !par
cpar!---------------------------------------------------------------!par
cparc                                                               !par
cpar      end module comm_interface                                 !par
cparc                                                               !par
cpar!***************************************************************!par
C
      PROGRAM MAIN
cparc                                                               !par
cpar      use comm_interface, only : iam,nproc,comm_init,           !par
cpar     A                           comm_barrier,comm_finalize     !par
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      CHARACTER*2 NAM0
C
C SUN TIME
      REAL*4 TARRY(2),TIME
C
      NAM0=''
cparc                                                               !par
cpar      call comm_init()                                          !par
cpar      write(0,*)'Starting proc', iam                            !par
cparc                                                               !par
cpar      ic1=iam/10                                                !par
cpar      ic2=iam-10*ic1                                            !par
cpar      ich0=ichar('0')                                           !par
cpar      ic1=ic1+ich0                                              !par
cpar      ic2=ic2+ich0                                              !par
cpar      nam0=char(ic1)//char(ic2)                                 !par
cpar      OPEN(5,FILE='adasin')                                     !par
C
C      OPEN(5,FILE='adasin')                                      !STDIN
      OPEN(6,FILE='adasout'//nam0)                               !STDOUT
C      OPEN(9,FILE='CAVES/TERMS/LEVELS')         !OPT TARGET SYMM/ENERGY
C      OPEN(10,FILE='adf48'//nam0)                   !ADAS FORMAT OUTPUT
C      OPEN(11,FILE='adf37')                     !ADAS NUMERIC E-DISTRIB
      OPEN(14,FILE='XRRTOT')                       !CROSS SECTION TOTALS
C
C on, onu FILES CONTAIN LEVEL INFO.
C      THEIR NV,LV SET MUST MATCH THOSE OF opn, opnu.
C      THIS IS THE CASE WHEN OUTPUT FROM THE SAME RUN.
C
C      OPEN(70,FILE='on')             !AUTOS DATA FILE (FORMATTED)
C OR
C      OPEN(70,FILE='onu',FORM='UNFORMATTED') !AUTOS DATA FILE (UNFORM)
C
C opn, opnu CONTAIN PHOTOIONIZATION DATA.
C      OPEN(80,FILE='opn')             !AUTOS DATA FILE (FORMATTED)
C OR
C      OPEN(80,FILE='opnu',FORM='UNFORMATTED') !AUTOS DATA FILE (UNFORM)
C
      CALL POSTP(NAM0)
C
C SUN TIME
      DUM=DTIME(TARRY)
      TIME=TARRY(1)
C
C CRAY TIME
CCRAY CALL SECOND(TIME)
C
      TIME=TIME/60.0
      WRITE(6,999)TIME
 999  FORMAT(//1X,'CPU TIME=',F9.3,' MIN')
C
      WRITE(6,*) 'PROGRAM ADASRR: NORMAL END'
C
      CLOSE(6)
C      CLOSE(9)
C      CLOSE(10)
C      CLOSE(11)
      CLOSE(14)
C      CLOSE(70)
C      CLOSE(80)
cparc                                                               !par
cpar      write(0,*)'Ending proc', iam                              !par
cpar      call comm_barrier()                                       !par
cpar      call comm_finalize()                                      !par
C
      END
C
C***********************************************************************
C
      SUBROUTINE POSTP(NAMX)
cparc                                                               !par
cpar      use comm_interface, only : iam,nproc,comm_init,           !par
cpar     A                           comm_barrier,comm_finalize     !par
cpar      use mpi                                                   !par
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (NDIM2=12)
      PARAMETER (NDIM5=150)
      PARAMETER (NDIM27=150)
      PARAMETER (NDIM25=NDIM27)
      PARAMETER (NDIM31=10)
      PARAMETER (NDIM37=19)
C
      PARAMETER (DZERO=0.0D0)
      PARAMETER (DONE=1.0D0)
c      PARAMETER (DTWO=2.0D0)
C
      PARAMETER (CONRYK=1.5789D5) !1.578885D5)      !RYDBERGS TO KELVIN
      PARAMETER (CONRYEV=13.6058)                   !RYDBERGS TO EV
C
      CHARACTER NAME*30,DATE*30,DATE8*8,NAMX*2,NAM0*2
C
      CHARACTER*4 NAMEJ,COD(20),MERGE
      CHARACTER*3 PIG
C
      REAL*8 KAPPA
C
      LOGICAL BCA,BLS,BIC,EX
C
      DIMENSION IWT(NDIM5),IWS(NDIM5),IWL(NDIM5)
     X         ,EI(NDIM5),IWJ(NDIM5),LCP(NDIM5),TEMP(NDIM37)
C
      ALLOCATABLE :: TEMPE37(:),E37(:,:),F37(:,:)
C
      COMMON /ECOR/E1C(NDIM5),TOLB,TOLB0
      COMMON /NRBTST/BSP0,INAG,NCP,NLAG,IMESH
      COMMON /COMMS/COD,NAME,DATE,nam0
C
      NAMELIST/ONE/NTAR1,NTAR2,IPRINT,NCUT,LCUT,NMIN,LMIN,NMAX,LMAX
     X            ,JCFJ,NRSLMX,NLMAX,LLMAX,NCORE,NZ,MR5
C
      NAMELIST/TWO/EBDMIN,EBDMAX,TOLR,TOLB,DELTAF,PFACT
     X            ,ITYPE,JTEMP,IRDT,JTHETA
     X            ,ICON,KAPPA,XDRY                !DEFINE NON-MAXWELLIAN
     X            ,IREL,NR1,NENG,INAG,NCP,NLAG,BSP,IMESH
     X            ,PIG,PMIN,TEAPOT,IOLDW
C
      nam0=namx
C
      NAME='NAME: ADASRR'
      CALL DATE_AND_TIME(DATE8)
      DATE='DATE: '//DATE8(7:7)//DATE8(8:8)//'/'//DATE8(5:5)//
     X               DATE8(6:6)//'/'//DATE8(3:3)//DATE8(4:4)
C
C READ HEADER TO DETERMINE IF CA, LS OR IC RUN (/CA/, /LS/ OR /IC/),
C THE REST OF THE LINE IS FOR COMMENT, IT IS ADDED TO THE END OF ADF48.
C
      READ(5,1000)COD
      WRITE(6,1001)COD
C
      BCA=COD(1).EQ.'/CA/'
      BLS=COD(1).EQ.'/LS/'
      BIC=COD(1).EQ.'/IC/'
      IF(.NOT.BCA.AND..NOT.BLS.AND..NOT.BIC)THEN
        WRITE(6,1002)COD(1)
        STOP 'ERROR: INPUT ERROR ON UNIT5'
      ENDIF
C
C
C**********************NAMELIST-ONE*************************************
C
C           ONLY THE STARRED '***' INPUT IS COMPULSORY
C-----------------------------------------------------------------------
C *** NTAR1 = NO OF INITIAL STATES POPULATED.
C *** NTAR2 = NO OF FINAL (PARENT) STATES
C             NEW DEFAULT: ALL FOUND IN FIXED FORMAT WEIGHT LIST
C           > 0 FIRST NTAR1 OF THEM=INITIAL
C           < 0 THESE ARE CONFIGURATION RESOLVED (REQUIRES TARGET o_str)
C           = 0 TOTAL RATE COEFFICIENTS ONLY, NO ADF48 FILE
C
C     FOR HEADER /CA/, /LS/ OR /IC/ THEN, *AFTER* NAMELIST-TWO, READ
C     NTAR2 CFG/TERM/LEVEL INFO AS OUTPUT BY AUTOS CAVES/TERMS/LEVELS
C     FILES. VIZ.
C   ***            IW  IPAR       CA
C       OR
C   ***      IWS  IWL  IPAR       LS
C       OR
C   ***      IWJ  IPAR IWS  IWL   IC  (HERE IWJ=2J).
C     THIS INFO IS WRITTEN TO THE ADF48 FILE.
C
C ***    END OF COMPULSORY INPUT    ***
C
C *   HOWEVER, IT MAY BE NECESSARY TO TELL ADASRR THE MAX CORE N IF IT
C     IS UNABLE TO DETERMINE IT ITSELF: SET NCORE =  "    "    " .
C
C-----------------------------------------------------------------------
C
C     NZ SET .GT.0 FOR BARE ION OF CHARGE NZ I.E. PURE HYDROGENIC.
C     NRSLMX -- IS THE MAX N OF RESOLVED DATA, DEFAULT=8.
C     NLMAX  -- IS THE MAX N OF BUNDLED-NL DATA (NOT LS) DEFAULT=10 (IC)
C                   ALSO USED BY CA, DEFAULT=100.
C     LLMAX  -- IS THE MAX L OF BUNDLED-NL DATA (NOT LS) DEFAULT=9.
C
C     JCFJ .GT. 0   NEGLECTS CAPTURE INTO CONFIGS .GT. JCFJ.
C                            DEFAULT: INCLUDE ALL.
C
C     IPRINT=UNIT6 PRINT LEVEL
C           .GE. 0, DETAILED PRINTOUT OF EACH PARTIAL CROSS SECTION
C           .EQ.-1, NL CROSS SECTIONS
C           .EQ.-2,  L CROSS SECTIONS
C           .LE.-3, TOTAL CROSS SECTION ONLY.
C
C     NCUT(OR MAX) .GT. 0 IGNORES CONTRIBUTIONS FROM N .GT. NCUT(OR MAX)
C     LCUT(OR MAX) .GE. 0         "         "        L .GT. LCUT(OR MAX)
C     NMIN .GT. 0, IGNORES CONTRIBUTIONS FROM N .LT. NMIN
C     LMIN .GE. 0,         "         "        L .LT. LMIN
C                            DEFAULT: INCLUDE ALL.
C     NOTE: IF NMIN,NMAX DO NOT MATCH A REPRESENTATIVE-N THEN THE TOTAL
C           IS TRUNCATED AT THE FIRST/LAST REP-N WHICH SATISFIES IT,
C           EVEN IF THERE ARE ADDITIONAL REP-N ON FILE. THUS, IT DIFFERS
C           FROM MDRCS13 IN THIS RESPECT BECAUSE PARTIALS ARE IN CONTROL
C
C
C****************************END-ONE***********************************
C
C
      MR5=0
      NTAR1=1
      NTAR2=NDIM5-1
      IPRINT=-1
      NCUT=-66
      LCUT=-77
      JCFJ=0
      NRSLMX=8
      NLMAX=-9999
      LLMAX=-9999
      NMIN=-10
      LMIN=-10
      NMAX=-10
      LMAX=-10
      NCORE=0
      NZ=0
C
C-----------------------------------------------------------------------
C
      READ(5,ONE)
C
C-----------------------------------------------------------------------
C
      IF(NZ.GT.0)THEN            !BARE ION ELSE MUST SPECIFY WNP(2) ETC.
        NTAR1=1
        NTAR2=1
      ENDIF
C
      IF(NTAR1.LT.1) THEN
        WRITE(6,849)NTAR1
 849    FORMAT(/' NTAR1 MUST BE .GT. 0, BUT INPUT NTAR1=',I3)
        STOP 'ERROR: NTAR1 MUST BE .GT. 0'
      ENDIF
      IF(NTAR1.GT.NDIM2)THEN
        WRITE(6,847)NTAR1
 847    FORMAT(/' INCREASE NDIM2 TO AT LEAST',I4)
        STOP 'ERROR: INCREASE NDIM2'
      ENDIF
      NBINI=NTAR1+1                                   !HISTORIC INTERNAL
      NBINM=NTAR1
C
      NTAR2O=NTAR2
      IF(NTAR2.LT.0)THEN                              !HYBRID PARENTS
        NBINR=NTAR2-1
      ELSE
        IF(NTAR2*NTAR2.LT.NTAR1*NTAR2)NTAR2=NTAR1    !.LT.ALLOWS NTAR2=0
COLD        IF(NTAR2.EQ.0)NTAR2=NDIM5-1         !=0 NOW FLAGS TOTAL ONLY
        NBINR=NTAR2+1
      ENDIF
      NBINRM=NTAR2
      IF(IABS(NBINR).GT.NDIM5)THEN
        WRITE(6,848)IABS(NBINR)
 848    FORMAT(/' INCREASE NDIM5 TO AT LEAST',I5)
        STOP 'ERROR: INCREASE NDIM5'
      ENDIF
C
      IF(JCFJ.LE.0)JCFJ=999
      JCFX=0
      NAMEJ='JCF*'
      IF(JCFJ.NE.999)THEN
        JCFX=JCFJ
        NAMEJ='JCFJ'
      ENDIF
C
      IF(NMAX.GT.0)NCUT=NMAX
      IF(LMAX.GT.-1)LCUT=LMAX
      NMAX=NCUT
      LMAX=LCUT
      IF(NCUT.LT.1)NCUT=1000
      IF(LCUT.LT.0)LCUT=150
C
      WRITE(6,11) NTAR1,NTAR2,NMIN,LMIN,NCUT,LCUT,NAMEJ,JCFX,NCORE
C
      IF(NZ.GT.0)THEN
        WRITE(6,846)NZ
 846    FORMAT(/'*** RR OF BARE NUCLEUS OF CHARGE:',I3/)
      ENDIF
C
C*************************NAMELIST-TWO**********************************
C
C  NOTHING IS COMPULSORY, BUT SEE JTHETA REGARDS TYPES 3 & 5 ADF48
C
C-----------------------------------------------------------------------
C
C     JTHETA.NE.0 RESTRICT NO. OF ADAS DEFAULT TEMPERATURES TO |JTHETA|
C                 WILL ALSO RESTRICT ANY ADF37 TEMPS.
C           .GT.0 WRITE A TYPE-3 adf48 FILE OF TEMP vs RATE COEFFICIENT.
C           .EQ.0 WRITE A TYPE-5 adf48 FILE OF E vs E*SIGMA.
C           .LT.0 WRITE A TYPE-5 adf48_5 AND A TYPE-3 adf48_3 FILE.
C         DEFAULT: 999 ALL POSSIBLE TEMPERATURES.
C                 CAN BE OVERRIDDEN BY ITYPE:
C
C     ITYPE .EQ.5 SETS JTHETA=0, I.E. TYPE-5 (FOR SCRIPT WORK)
C              -5 SETS JTHETA=-99, I.E. TYPE-5 AND TYPE-3.
C           .EQ.3 SETS JTHETA POSITIVE, I.E. TYPE-3 (ONLY), BUT ALLOWS
C                 THE NUMBER OF TEMPS TO BE RESTRICTED BY JTHETA.
C            ELSE DOES NOTHING, SO JTHETA PREVAILS (DEFAULT IS TYPE-3
C                 AT ALL ADAS TEMPERATURES.)
C
C     ANY ELECTRON DISTRIBUTION IS MAXWELLIAN BY DEFAULT (ICON=0), BUT
C
C     KAPPA   -- IF .GT. 1.5 THEN USES A KAPPA-DISTRIBUTION (ICON=1).
C
C     XDRY    -- IF .GE. 1.0 USES A DRUYVESTEYN DISTRIBUTION (ICON=2).
C
C     ICON    -- IF .EQ. 3 THEN READ AN adf37 NUMERICAL DISTRIBUTION.
C                DEFAULT: -1 IS RESET TO MAXWELLIAN, KAPPA OR XDRY,
C                I.E. USER SHOULD ONLY SET FOR .EQ. 3.
C                OTHER VALUES ARE RESERVED FOR POSSIBLE FUTURE USE.
C
C     TEAPOT       SET .GT. 0 TO OBSERVED GROUND STATE IONIZATION
C                  POTENTIAL (RYDBERGS). CAN BE USED TO CORRECT
C                  CONVERSION FROM PI TO PR AT LOW-E FOR LOW-CHARGE.
C                  ONLY APPLIED FOR LV=-1 AND ASSUMES GROUND STATE
C                  OF ATOM AND ION PRESENT. ***BETTER to USE ECORLS,
C                  ECORIC in AUTOSTRUCTURE RUN. DEFAULT -1.
C
C  ***NR1          SET IT EQUAL TO THE LARGEST CORE PRINCIPAL QUANTUM
C                  NUMBER (NCORE) +1. THIS IS OFTEN THE LOWEST NV.
C                  PARENTAGE IS DETERMINED BY BUILDING STATES UPON
C                  THE NTAR2 PARENTS FROM NR1 UPWARDS. IF NR1.EQ.0 THEN
C                  IT IS DETERMINED ON ENERGY GROUNDS ALONE (NOT GOOD).
C                  ALSO NO PP HIGH-L RECOMBINATION THEN, SO CODE STOPS.
C
C-----------------------------------------------------------------------
C
C  ***VARIABLES HERE-ON ARE NORMALLY FOR TESTING ETC. ONLY***
C
C     EBDMIN, EBDMAX -- ONLY EVALUATE RR RATE COEFFICIENT TO
C         FINAL ATOM (BINDING) ENERGIES (RELATIVE TO LOWEST
C         CONTINUUM) IN THIS ENERGY (RYD) RANGE.
C         DEFAULT: (-1.0D30, 1.0D30) I.E. NO RESTRICTION.
C
C     TOLR  CONTROLS RR STABILIZATION, NORMALLY SET INTERNALLY, TO
C           FINAL STATES UP TO TOLR RYD ABOVE THE IONIZATION LIMIT.
C           THIS INCLUDES METASTABLE AUTOIONIZING FINAL STATES,
C                         C-R PARTIAL DATA ONLY, TOTALS STILL CORONAL.
C                         SET .EQ. 0.0 TO FORCE NON-AUTOIONIZING ONLY.
C
C     TOLB=MAX(1.5D-7,5.0D-9*DZ*NZ),  DEFAULT.
C         SET TOLB COARSER TO HANDLE USER SUPPLIED IMBALANCED CONTINUUM
C         EXPANSIONS, I.E. IF NOT ALL PARTIAL WAVES HAVE SAME TARGET CI.
C
C     DELTAF IS THE MINIMUM VALUE OF THE OSCILLATOR STRENGTH RETAINED
C            FOR OUTPUT TO ADF48. DEFAULT ALL.
C
C     PFACT: IF PFACT*PCS(IE).GT.PCS(IE-1) THEN DROP ENERGY IE FROM
C         PCS TABULATION. TEST IS DONE ON FIRST PCS ONLY.
C         DEFAULT: -1., ALL ENERGIES INCLUDED.
C     PFACT=1.2 IS A REASONABLE TEST.
C
C     JTEMP.EQ.0 USE ADAS DEFAULT TEMPERATURES, IF REQUIRED.
C          .GT.0 READ-IN JTEMP TEMPERATURES IN KELVIN.
C          .LT.0 READ-IN -JTEMP TEMPERATURES IN LOG(K).
C           N.B. JTEMP.NE.0 IGNORES ANY JTHETA.NE.0, BUT SIGN IS USED.
C
C     IRDT.EQ.0 HISTORIC TEMP(J=1,JTEMP) READ AFTER ALL TARGET INFO ETC.
C         .NE.0 STRAIGHT AFTER NAMELIST, FOR EASE OF USE WITH SCRIPT.
C
C     IREL.NE.0 THEN APPLY RELATIVISTIC (JUTTNER) CORRECTION
C               TO MAXWELL RATE COEFFICIENTS.
C         .GT.0 ALSO ASSUMES RELATIVISTIC ORBITALS, >DIPOLE ETC AND
C               SO SWITCHES-OFF COULOMB REPLACEMENT/TOP-UP BY ENERGY,
C               TOP-UP IN ANGULAR MOMENTUM IS STILL DONE, SINCE SMALL
C               AT RELATIVISTIC ENERGIES. DROPS LAST TWO ADAS TEMPS.
C         .LT.0 ASSUMES NON-RELATIVISTIC DIPOLE, I.E. ASYMPTOTIC TO
C               COULOMB IN ENERGY, NO COULOMB/ENERGY/TEMPERATURE
C               RESTRICTIONS, ONLY JUTTNER APPLIED.
C         .EQ.0 DEFAULT: DOES NOTHING. WE LEAVE JUTTNER TO ADAS TO APPLY
C
C     IOLDW=0 DEFAULT: FIXED FORMAT READ OF STAT. WEIGHTS; AND ANY CALC.
C          =1 USE OLD FREE-FORMATTED READ OF STAT. WEIGHTS; AND NO CALC.
C
C-----------------------------------------------------------------------
C DO NOT CHANGE NEXT 8 VARIABLES UNLESS YOU KNOW WHAT YOU ARE DOING:
C THEY CONTROL ACCURACY OF MAXWELLIAN INTEGRATION OF CROSS SECTIONS.
C-----------------------------------------------------------------------
C
C     NENG.NE.0 RESTRICT PI ENERGIES USED IN MAXWELLIAN TO FIRST |NENG|.
C         IF NENG.GT.0 HYDROGENIC TOP-UP IN ENERGY IS THEN USED THERE-ON
C                      AS DW NOT ACCURATE TO SUFFICIENTLY HIGH ENERGY.
C         DEFAULT .EQ. 0 IS RESET TO 999 I.E. HYDROGENIC TOP-UP ON.
C
C     PMIN: SWITCH TO HYDROGENIC TOP-UP IN ENERGY WHEN COULOMB PI XSCTN
C           FALLS BELOW PMIN (Z-SCALED CM**2), VALUE GAUGE DEPENDENT VIZ
C           1.D-24 (LEN), 1.D-30 (VEL), 1.D-36 (ACC) DEFAULTS FOR PIG.
C
C     PIG = 'LEN', 'VEL', 'ACC' SPECIFIES THE PI GAUGE USED BY AUTOS RUN
C            DEFAULT: 'ACC', I.E. ASSUMES AUTOS RUN WITH VARIABLE GAUGE.
C
C     IMESH .LE. 0 USES THE DEFAULT ENERGY MESH OF AUTOS FOR HYDROGENIC
C                  PI, I.E. 3-POINTS PER DECADE.
C           .GT. 0 USES 5-POINTS PER DECADE, SHOULD CHECK DW MESH THEN.
C
C     INAG SWITCHES BETWEEN QUADRATURE RULES! DEFAULT INAG=0 (NCP=5).
C
C     .EQ. 0 NCP-PT NEWTON-COTES, NCP=3 OR 5: 3 USES DOUBLE STEP OF 5,
C                                 NEGATIVE HALVES STEP-LENGTH.
C     .LT. 0 MID-POINT RULE, DOES NOT CONVERGE HIGH-N, HIGH-T PARTIALS
C     .GT. 0 GAUSS-LAGUERRE, DOES NOT CONVERGE PARTIALS AT HIGH-T
C
C     BSP .NE. 0 USE SPLINE TO FIT RR CROSS SECTIONS, SCALING THEM BY
C                A FACTOR OF (BSP+E**L+.5).
C         .LT. 0 DETERMINES BSP INTERNALLY .ELSE.
C         .GT. 0 USES SAID VALUE FOR BSP.
C     BSP .EQ. 0 USE NLAG POINT LAGRANGE INTERPOLATION (DEFAULT).
C
C     NLAG .EQ. NO. OF LAGRANGE INTERP. PTS OF CROSS SECTIONS. DEFAULT 4
C          FOR AUTOS 3-POINTS PER DECADE. *** DO NOT INCREASE *** UNLESS
C          USING A FINER AUTOS MESH ELSE SPANS TOO WIDE AN ENERGY RANGE.
C
C***************************END-TWO*************************************
C
      TOLR=1.D-10
      TOLB=-DONE
      EBDMIN=-1.0D30
      EBDMAX=1.0D30
      PFACT=-DONE
      DELTAF=-1.D-5
      NR1=999                 !NEW DEFAULT
      ITYPE=0
      JTHETA=999
      JTEMP=0
      IRDT=0
      ICON=-1
      KAPPA=-DONE
      XDRY=-DONE
      IREL=0
      TEAPOT=-1
      IOLDW=0
C
C DO NOT CHANGE NEXT 8 VARIABLES UNLESS YOU KNOW WHAT YOU ARE DOING
      NENG=0
      INAG=0
      NCP=0
      NLAG=0
      BSP=DZERO
      IMESH=0
      PMIN=-DONE
      PIG='ACC'
      IPIG=-1                                 !NOT USED
C
C-----------------------------------------------------------------------
C
      READ(5,TWO)
C
C-----------------------------------------------------------------------
C
      IF(ABS(ITYPE).EQ.5)THEN                 !MAP TO JTHETA TYPE-5 FLAG
        IF(ITYPE.GT.0)THEN                    !TYPE-5 ONLY
          JTHETA=0
        ELSE                                  !TYPE-3 AS WELL
          JTHETA=-99
        ENDIF
      ELSEIF(ITYPE.EQ.3)THEN                  !SWITCH-OFF CALC OF BOTH
        IF(JTHETA.EQ.0)THEN                   !SWITCH TYPE-5 TO TYPE-3
          JTHETA=999                          !ALL TEMPS
        ELSE
          JTHETA=ABS(JTHETA)                  !BUT ALLOW RESTICT NUMBER
        ENDIF
      ENDIF
C
      IF(IRDT*JTEMP.NE.0)THEN            !ALTERNATE TEMP READ FOR SCRIPT
        JJTEMP=ABS(JTEMP)
        IF(JJTEMP.GT.NDIM37)THEN
          WRITE(6,*)'TOO MANY TEMPS; INCREASE NDIM37 TO:',JJTEMP
          STOP 'ERROR: TOO MANY TEMPS; INCREASE NDIM37 TO JTEMP'
        ENDIF
        READ(5,*)(TEMP(K),K=1,JJTEMP)
      ENDIF
C
C CHECK FOR SENSIBLE NON-MAXWELLIAN INPUT
C
      DKMIN=1.5D0+1.D-10                                      !MIN KAPPA
      DKMAX=1.0D10                                            !MAX KAPPA
      DXMIN=0.1D0                                              !MIN XDRY
      DXMAX=1.0D3                                              !MAX XDRY
C
      ICON0=ICON
      IF(ICON0.LT.0)THEN
        IF(KAPPA.GE.DKMIN)THEN
          ICON=1
        ELSEIF(XDRY.GE.DXMIN)THEN
          ICON=2
        ELSE
          ICON=0
        ENDIF
      ELSEIF(ICON0.EQ.1)THEN
        IF(KAPPA.LT.DKMIN)
     X  STOP 'MUST SET A VALID KAPPA VALUE FOR THE K-DISTRIBUTION'
      ELSEIF(ICON0.EQ.2)THEN
        IF(XDRY.LT.DXMIN) STOP
     X  'MUST SET A VALID XDRY VALUE FOR THE DRUYVESTEYN DISTRIBUTION'
      ELSEIF(ICON0.EQ.3)THEN
C
        INQUIRE (FILE='adf37',EXIST=EX)
        IF(.NOT.EX)STOP '*** USER ICON=3 BUT adf37 FILE NOT FOUND...'
        ICON=3
C
        OPEN(11,FILE='adf37')
        READ(11,*)                !SKIP HEADER
C
        READ(11,*)IFORM37         !=1 SUPERPOSITION =2 NUMERICAL
        IF(IFORM37.EQ.1)STOP 'NOT CODED FOR SUPERPOSITION'
        READ(11,*)IUNITS37        !=1 KELVIN, =2 EV
        IF(IUNITS37.EQ.1)THEN
          CONE=CONRYK
        ELSE
          CONE=CONRYEV
        ENDIF
C
        READ(11,*)NENG37          !NO. ENERGIES FOR ADF37 F(E)
C
        READ(11,*)JTEMP37         !NO. TEMPERATURES
        IF(JTHETA.NE.0)JTEMP37=MIN(JTEMP37,ABS(JTHETA))   !USER RESTRICT
        IF(JTEMP37.GT.NDIM37)THEN
          WRITE(0,*)'*** FOR ADF37, INCREASE NDIM37 TO',JTEMP37
          STOP'*** INCREASE NDIM37 FOR ADF37'
        ENDIF
C
        ALLOCATE (TEMPE37(JTEMP37))
        ALLOCATE (E37(NENG37,JTEMP37),F37(NENG37,JTEMP37))
C
        READ(11,*)IDUM            !CUT-OFF THRESHOLD (UNUSED)
        READ(11,*)IDUM,DUM        !HIGH ENERGY BEHAVIOUR (UNUSED)
C
        DO J=1,JTEMP37
          TEMPE37(J)=J                              !SIMPLE INDEX LABEL
          READ(11,*)(E37(N,J),N=1,NENG37)
          READ(11,*)(F37(N,J),N=1,NENG37)
          DO N=1,NENG37                             !CONVERT TO RYD
            IF(E37(N,J).EQ.DZERO)THEN
              F37(N,J)=DZERO
            ELSE
              E37(N,J)=E37(N,J)/CONE
              F37(N,J)=F37(N,J)*CONE/SQRT(E37(N,J)) !SO FBAR NOW
            ENDIF
          ENDDO
          CALL ADF37(IPRINT,J,TEMPE37(J),ICON,NENG37,E37(1,J),F37(1,J))
        ENDDO
        CLOSE(11)
        JTEMP=0                 !SWITCH-OFF ANY USER READS
      ELSEIF(ICON0.GT.3)THEN
        STOP 'ICON.GT.3 IS RESERVED FOR POSSIBLE FUTURE USE'
      ENDIF
C
      IF(XDRY.GT.DXMAX)STOP
     X 'MUST SET A SENSIBLE XDRY VALUE FOR THE DRUYVESTEYN DISTRIBUTION'
      IF(KAPPA.GT.DKMAX)THEN                       !RE-SET TO MAXWELLIAN
        KAPPA=-1
        ICON=0
      ENDIF
C
      IF(JTHETA.NE.0)THEN
        JT=MIN(ABS(JTHETA),NDIM37)
        JTHETA=SIGN(JT,JTHETA)
        IF(ICON.EQ.0)WRITE(6,520)
        IF(ICON.EQ.1)WRITE(6,521)KAPPA
        IF(ICON.EQ.2)WRITE(6,522)XDRY
        IF(ICON.EQ.3)WRITE(6,523)
      ENDIF
C
      IF(ICON*IREL.NE.0)WRITE(6,776)
C
      IF(ICON*INAG.NE.0)THEN
        WRITE(6,*)'*** ATTENTION: ONLY NEWTON-COTES QUADRATURE CODED'
     X           ,' FOR NON-MAXWELLIANS, RE-SETTING YOUR INAG=0'
        INAG=0
      ENDIF
C
C THESE ARE SET IN NAMELIST ONE, BUT NOW USE IREL FROM TWO...
      IF(NLMAX.EQ.-9999)THEN                          !UNSET BY USER, SO
        IF(BCA.OR.NTAR2.LT.0)THEN               !CA OR FINAL CF RESOLVED
          IF(IREL.GT.0)THEN              !TRUE BOUND ONLY - SYNC WITH AS
            NLMAX=20
          ELSE
            NLMAX=100
          ENDIF
        ELSE                                             !LS/IC RESOLVED
          NLMAX=10
        ENDIF
      ENDIF
      IF(NTAR2.EQ.0)NLMAX=0                              !TOTAL ONLY
      NLMAX=MIN(NLMAX,NCUT)
C
      LLMAX0=LLMAX
      IF(LLMAX.EQ.-9999)THEN
C        LLMAX=9
        LLMAX=NLMAX-1
      ENDIF
      IF(LCUT.LT.NCUT)LLMAX=MIN(LCUT+0,ABS(LLMAX))  !+0 AS LCUT IS LOWER
      IF(LLMAX+1.GT.NDIM31)THEN
        WRITE(6,*)
        WRITE(6,*)'NOTE: REPRESENTATIVE NL LIMITED TO L=',NDIM31-1
     X   ,', INCREASE NDIM31 TO OBTAIN MORE'
C        WRITE(0,*)'NOTE: REPRESENTATIVE NL LIMITED TO L=',NDIM31-1
C     X   ,', INCREASE NDIM31 TO OBTAIN MORE'
        LLMAX=NDIM31-1
      ENDIF
      IF(LLMAX0.EQ.-9999)LLMAX=-LLMAX
C
      IF(NRSLMX.LT.0)NRSLMX=8
      NRSLMX=MIN(NRSLMX,NCUT)
      IF(NRSLMX.GT.NDIM25)STOP 'ERROR: INCREASE NDIM25 TO NRSLMX'
      IF(NRSLMX.GT.NDIM31)STOP 'ERROR: INCREASE NDIM31 TO NRSLMX'
C
      IF(PMIN.LE.DZERO)THEN            !PMIN=1.D-5 FOR DIP INT, ELSE PCS
        IF(PIG.EQ.'LEN')PMIN=1.D-24           !IPIG=1
        IF(PIG.EQ.'VEL')PMIN=1.D-30           !IPIG=0
        IF(PIG.EQ.'ACC')PMIN=1.D-36           !IPIG=-1
      ENDIF
C
      BSP0=BSP
      IF(NZ.LE.0.AND.NCP.EQ.0)NCP=5
      IF(NZ.GT.0)THEN
        IF(NCP.EQ.0)NCP=-25
        IF(IMESH.EQ.0)IMESH=1
        IF(NLAG.LE.0)THEN
          NLAG=4
          IF(IMESH.GT.0)NLAG=6
        ENDIF
      ENDIF
C
      NLAG=(NLAG/2)*2
      IF(NLAG.LE.0)NLAG=4
      IF(IMESH.LE.0)NLAG=MIN(NLAG,4)
      IF(IMESH.GT.0)NLAG=MIN(NLAG,6)
C
      IF(INAG.GT.0)INAG=MAX(INAG,6)
      IF(INAG.LT.0)INAG=MIN(INAG,-10)
      IF(INAG.NE.0)THEN
        WRITE(*,*) '*** WARNING *** WARNING *** WARNING *** WARNING ***'
        WRITE(*,*) 'INAG.NE.0: HIGH- N/T PARTIALS MAY NOT BE CONVERGED'
        WRITE(*,*) '*** WARNING *** WARNING *** WARNING *** WARNING ***'
        WRITE(6,*) '*** WARNING *** WARNING *** WARNING *** WARNING ***'
        WRITE(6,*) 'INAG.NE.0: HIGH- N/T PARTIALS MAY NOT BE CONVERGED'
        WRITE(6,*) '*** WARNING *** WARNING *** WARNING *** WARNING ***'
      ENDIF
C

      NP0=NENG
      IF(NP0.EQ.0)THEN
        NP0=999
        IF(IREL.GT.0)NP0=-NP0                    !SWITCH-OFF TOP-UP IN E
      ENDIF
C
      IF(NTAR2O.EQ.0.AND.IOLDW.EQ.1)THEN         !RE-INSTATE OLD DEFAULT
        NTAR2=NTAR1
        NBINRM=NBINM
        NBINR=NBINI
      ENDIF
      TOLB0=TOLB
      IF(NTAR2.LT.0)TOLR=DZERO
      IF(TOLR.LT.0.0D0)TOLR=0.0D0
C
      WRITE(6,13)EBDMIN,EBDMAX,NR1,TOLB,TOLR,DELTAF
C
      IF(NR1.NE.999.AND.NCORE.GT.0)THEN
        IF(IABS(NR1).NE.NCORE+1)THEN
          WRITE(6,*)'*** ERROR, NCORE AND NR1 INCONSISTENT:',NCORE,NR1
          STOP '*** ERROR, NCORE AND NR1 INCONSISTENT'
        ENDIF
      ENDIF
C
C READ ELECTRON TARGET INFO
C
      IF(NTAR2.GT.0)THEN
        NBINP=NBINRM
      ELSE
        NBINP=NBINM
      ENDIF
C
      IF(MR5.LE.0)THEN
        MR5=5
      ELSE
        IF(MR5.NE.5)MR5=9
      ENDIF
C
      EI(1)=0.0D0
      IF(BLS)THEN
        IF(IOLDW.EQ.1)WRITE(6,816)
        IF(IOLDW.EQ.0)WRITE(6,814)
        IF(MR5.NE.5)THEN
          INQUIRE(FILE='TERMS',EXIST=EX)
          IF(.NOT.EX)GO TO 822
          OPEN(MR5,FILE='TERMS')
          READ(MR5,*,END=822)
        ENDIF
        DO I=1,NBINP
          LCP(I)=0
          IF(IOLDW.EQ.1)READ(MR5,*)IWS(I),IWL(I),IPAR
          IF(IOLDW.EQ.0)READ(MR5,992)IWS(I),IWL(I),IPAR
     X                              ,LCP(I),NI,E1C(I),MERGE
          IF(IWS(I).EQ.0)THEN
            IF(NBINP+1.LT.NDIM5)THEN
              WRITE(6,*)'*** PREMATURE END OF STAT. WEIGHT INPUT:'
              WRITE(6,*)'*** REQUESTED=',NBINP,' BUT FOUND=',I-1
              STOP '*** PREMATURE END OF STAT. WEIGHT INPUT'
            ELSEIF(NTAR2.GT.0)THEN      !ASSUME INPUT UNSET (NTAR2=0)
              NBINR=I
              NBINRM=NBINR-1
              NBINP=NBINRM
              GO TO 994                 !TERMINATOR
            ENDIF
          ENDIF
          IF(IOLDW.EQ.1)WRITE(6,817)IWS(I),IWL(I),IPAR
          IF(IOLDW.EQ.0)WRITE(6,817)IWS(I),IWL(I),IPAR,E1C(I)
        ENDDO
C        IF(E1C(2).EQ.DZERO)IOLDW=1
        IF(IOLDW.EQ.0)THEN
          E1C(NBINP+1)=DZERO
          READ(MR5,992,END=994)ITEST,IDUM,IDUM,IDUM,IDUM,ETARG1
          IF(ITEST.EQ.0)GO TO 994       !TERMINATOR
          E1C(NBINP+1)=ETARG1
          IF(MR5.EQ.5)THEN              !SKIP ANY EXTRA TARGET INFO
            DO I=1,9999
              READ(5,992,END=994)ITEST,IDUM,IDUM,IDUM,IDUM,ETARG1
              IF(ITEST.EQ.0)GO TO 994   !TERMINATOR
            ENDDO
          ENDIF
        ENDIF
      ENDIF
C
      IF(BCA)THEN
        WRITE(6,811)
        IF(MR5.NE.5)THEN
          INQUIRE(FILE='CAVES',EXIST=EX)
          IF(.NOT.EX)GO TO 822
          OPEN(MR5,FILE='CAVES')
          READ(MR5,*,END=822)
        ENDIF
        DO I=1,NBINP
          IWL(I)=0                               !TO USE BLS
          READ(MR5,*)IWS(I),IPAR,LCP(I),E1C(I)   !5,* SINCE NO OLD STYLE
          IF(IWS(I).EQ.0)THEN
            IF(NBINP+1.LT.NDIM5)THEN
              WRITE(6,*)'*** PREMATURE END OF STAT. WEIGHT INPUT:'
              WRITE(6,*)'*** REQUESTED=',NBINP,' BUT FOUND=',I-1
              STOP '*** PREMATURE END OF STAT. WEIGHT INPUT'
            ELSEIF(NTAR2.GT.0)THEN      !ASSUME INPUT UNSET (NTAR2=0)
              NBINR=I
              NBINRM=NBINR-1
              NBINP=NBINRM
              GO TO 994                 !TERMINATOR
            ENDIF
          ENDIF
          WRITE(6,837)IWS(I),IPAR,E1C(I)
        ENDDO
        E1C(NBINP+1)=DZERO
        READ(MR5,*,END=994)ITEST,IDUM,IDUM,ETARG1
        IF(ITEST.EQ.0)GO TO 994         !TERMINATOR
        E1C(NBINP+1)=ETARG1
        IF(MR5.EQ.5)THEN
          DO I=1,9999                   !SKIP ANY EXTRA TARGET INFO
            READ(5,*,END=994)ITEST,IDUM,IDUM,ETARG1
            IF(ITEST.EQ.0)GO TO 994     !TERMINATOR
          ENDDO
        ENDIF
      ENDIF
C
      IF(BIC)THEN
        IF(IOLDW.EQ.1)WRITE(6,818)
        IF(IOLDW.EQ.0)WRITE(6,813)
        IF(MR5.NE.5)THEN
          INQUIRE(FILE='LEVELS',EXIST=EX)
          IF(.NOT.EX)GO TO 822
          OPEN(MR5,FILE='LEVELS')
          READ(MR5,*,END=822)
        ENDIF
        DO I=1,NBINP
          LCP(I)=0
          IF(IOLDW.EQ.1)READ(5,*)IWJ(I),IPAR,IWS(I),IWL(I)
          IF(IOLDW.EQ.0)READ(5,993)IWJ(I),IPAR,IWS(I),IWL(I)
     X                            ,LCP(I),NI,E1C(I),MERGE
          IF(IWS(I).EQ.0)THEN
            IF(NBINP+1.LT.NDIM5)THEN
              WRITE(6,*)'*** PREMATURE END OF STAT. WEIGHT INPUT:'
              WRITE(6,*)'*** REQUESTED=',NBINP,' BUT FOUND=',I-1
              STOP '*** PREMATURE END OF STAT. WEIGHT INPUT'
            ELSE                        !ASSUME INPUT UNSET (NTAR2=0)
              NBINR=I
              NBINRM=NBINR-1
              NBINP=NBINRM
              GO TO 994                 !TERMINATOR
            ENDIF
          ENDIF
        ENDDO
        DO I=1,NBINP
          IPAR=0
          IF(IWS(I).LT.0)IPAR=1
          IWS(I)=IABS(IWS(I))
          IF(IOLDW.EQ.1)WRITE(6,819)IWJ(I),IPAR,IWS(I),IWL(I)
          IF(IOLDW.EQ.0)WRITE(6,819)IWJ(I),IPAR,IWS(I),IWL(I),E1C(I)
        ENDDO
C        IF(E1C(2).EQ.DZERO)IOLDW=1
        IF(IOLDW.EQ.0)THEN
          E1C(NBINP+1)=DZERO
          READ(MR5,993,END=994)IDUM,IDUM,ITEST,IDUM,IDUM,IDUM,ETARG1
          IF(ITEST.EQ.0)GO TO 994       !TERMINATOR
          E1C(NBINP+1)=ETARG1
          IF(MR5.EQ.5)THEN              !SKIP ANY EXTRA TARGET INFO
            DO I=1,9999
              READ(5,993,END=994)IDUM,IDUM,ITEST,IDUM,IDUM,IDUM,ETARG1
              IF(ITEST.EQ.0)GO TO 994   !TERMINATOR
            ENDDO
          ENDIF
        ENDIF
      ENDIF
C
 994  IF(MR5.NE.5)CLOSE(MR5)                               !=9
c
      if(nbinp.gt.9999)then
        write(6,*)
     x       'Target indexing needs adjusting for ntar.gt.9999'
        stop 'Target indexing needs adjusting for ntar.gt.9999'
      endif
C
C CHECK/SET TOLB
C
      IF(IOLDW.EQ.0)THEN
        TOLB=1.D10
        DO I=2,NBINP+1
          T=E1C(I)-E1C(I-1)
          IF(T.GT.0.0D0.AND.T.LT.TOLB)TOLB=T
        ENDDO
        IF(TOLB.LT.TOLB0)THEN   !WRITE WARNING, BUT ALLOW USER SET VALUE
          WRITE(6,829)TOLB0,TOLB
          TOLB=TOLB0                              !RE-INSTATE USER VALUE
        ELSE
          IF(TOLB0.LE.0.0D0.AND.TOLB.GT.0.0D0)THEN
            TOLB=TOLB/2
          ELSE                             !ORIGINAL INPUT (MAYBE UNSET)
            TOLB=TOLB0
          ENDIF
        ENDIF
      ENDIF
C
C INITIALIZE TARGET WEIGHTS
C
      DO I=1,NBINP
        IWS(I)=IABS(IWS(I))
        IF(IWS(I).EQ.0.AND.IWL(I).NE.0)THEN
          IWT(I)=IWL(I)                      !2J+1
          IWJ(I)=IWL(I)-1                    !2J
          IWL(I)=0
        ELSE
          IF(BLS.OR.BCA)THEN                 !IWS=W & IWL=0 FOR CA
            IWT(I)=IWS(I)*(2*IWL(I)+1)
          ENDIF
          IF(BIC)THEN
            IWT(I)=IWJ(I)+1                  !2J+1
          ENDIF
        ENDIF
      ENDDO
C
C SET TEMPERATURES (CURRENTLY, ASSUME USER INPUT ALWAYS IN KELVIN)
C
      IF(JTEMP.NE.0)THEN
        JJTEMP=ABS(JTEMP)
        IF(JJTEMP.GT.NDIM37)THEN
          WRITE(6,*)'TOO MANY TEMPS; INCREASE NDIM37 TO:',JJTEMP
          STOP 'ERROR: TOO MANY TEMPS; INCREASE NDIM37 TO JTEMP'
        ENDIF
        IF(IRDT.EQ.0)READ(5,*)(TEMP(K),K=1,JJTEMP)  !IF NOT ALREADY READ
        IF(JTEMP.LT.0)THEN
          DO K=1,JJTEMP
            TEMP(K)=10**TEMP(K)
          ENDDO
          JTEMP=-JTEMP
        ENDIF
        IF(JTHETA.LT.0)JTEMP=-JTEMP
        IF(INAG.EQ.0)THEN
          TMIN=1.D20
          TMAX=-1.D20
          DO K=1,JJTEMP
            TMIN=MIN(TMIN,TEMP(K))
            TMAX=MAX(TMAX,TEMP(K))
          ENDDO
          WRITE(6,833)TMIN
 833      FORMAT('***ATTENTION: IF TMIN=',1PD10.2,' .LT. 10*ZA**2 K,'
     X              ,' THEN NEWTON-COATES QUADRATURE MAYBE INACCURATE.')
          WRITE(*,*)'***ATTENTION: IF TMIN .LT. 10*ZA**2 K, THEN',
     X              ' NEWTON-COATES QUADRATURE MAYBE INACCURATE.'
          WRITE(6,834)TMAX
 834      FORMAT('***WARNING: IF TMAX=',1PD10.2,' .GT. 10**7*ZA**2 K,'
     X     ,' THEN NEWTON-COATES QUADRATURE IN ERROR.'/' TRY INAG.NE.0,'
     X          ,' SLOWER AND INACCURATE AT HIGH-T, ESP. FOR PARTIALS.')
          WRITE(*,*)'***WARNING: IF TMAX .GT. 10**7*ZA**2 K, THEN',
     X              ' NEWTON-COATES QUADRATURE IN ERROR.'
        ENDIF
        DO K=1,JJTEMP
          TEMP(K)=TEMP(K)/CONRYK         !SINCE USER INPUT ASSSUMED IN K
        ENDDO
        JTEMP37=JJTEMP
        NENG37=1
      ELSEIF(ICON.EQ.3)THEN                      !NOW TRANSFER TEMP INFO
        IF(JTHETA.NE.0)THEN
          JTEMP=SIGN(JTEMP37,JTHETA)
          DO J=1,JTEMP37
            TEMP(J)=TEMPE37(J)
          ENDDO
        ELSE
          JTEMP=0
        ENDIF
      ELSE
        TEMP(1)=-DONE
        JTEMP=JTHETA
        IF(JTEMP.EQ.0)THEN                       !NONE (JUST SIGMA)
          JTEMP37=1
        ELSE
          JTEMP37=ABS(JTEMP)
        ENDIF
        NENG37=1
      ENDIF
      IF(ICON0.NE.3)ALLOCATE (E37(NENG37,JTEMP37),F37(NENG37,JTEMP37))
C
C SUM OVER CROSS SECTIONS
C
      CALL CROSSJ(NBINM,NBINR,NMIN,LMIN,NCUT,LCUT,NRSLMX,DELTAF,PMIN
     X     ,IPIG,EI,IWT,NR1,IPRINT,TOLR,EBDMIN,EBDMAX,PFACT,NMAX,LMAX
     X     ,NCORE,IWS,IWL,IWJ,LCP,BCA,BLS,BIC,NLMAX,LLMAX,JTEMP,TEMP
     X     ,IREL,JCFJ,NZ,NP0,TEAPOT,ICON,KAPPA,XDRY,NENG37,E37,F37)
C
C COMMENTS
C
      WRITE(6,770)
      IF(IREL.NE.0.AND.JTEMP.NE.0)WRITE(6,775)
      WRITE(6,1005)(COD(I),I=2,20)
      WRITE(6,790)NAME,DATE
C
      RETURN
C
  822 WRITE(6,735)
      STOP '*** TARGET SYMM INFO NOT FOUND ON FILE!!!'
C
  11  FORMAT(/' NTAR1=',I3,3X,'NTAR2=',I3,3X,'NMIN=',I4,3X,'LMIN=',I3
     X,3X,'NMAX=',I4,3X,'LMAX=',I3,3X,A4,'=',I4,3X,'NCORE=',I3)
  13  FORMAT(/1X,'EBDMIN=',1PE10.3,3X,'EBDMAX=',E10.3,3X,'NR1=',I3
     X,3X,'TOLB=',F12.8,3X,'TOLR=',F10.4,3X,'DELTAF=',F10.7)
  520 FORMAT(/' CONVOLUTION ELECTRON DISTRIBUTION IS MAXWELLIAN')
  521 FORMAT(/' CONVOLUTION ELECTRON DISTRIBUTION IS KAPPA, WITH K='
     X,F7.2)
  522 FORMAT(/' CONVOLUTION ELECTRON DISTRIBUTION IS DRUYVESTEYN,'
     X,' WITH X=',F6.1)
  523 FORMAT(/' CONVOLUTION ELECTRON DISTRIBUTION IS NUMERICAL (adf37)')
  735 FORMAT(/' *** ERROR: EMAX.LT.0 FLAGS READ OF TARGET SYMMS FROM'
     X,' FILES CAVES/TERMS/LEVELS, BUT NONE FOUND/EMPTY')
  770 FORMAT('C',110('-')/'C')
  775 FORMAT('C     JUTTNER RELATIVISTIC CORRECTION APPLIED TO THE',
     X ' DISTRIBUTION'/'C')
  776 FORMAT(/'C *** ATTENTION: YOU ARE APPLYING A JUTTNER RELATIVISTIC'
     X      ,' CORRECTION TO A NON-MAXWELLIAN DISTRIBUTION...'/'C')
  790 FORMAT('C'/'C',1X,A30/'C',1X,A30/'C',110('-'))
  811 FORMAT(/1X,'    W  P        E1C(RYD)')
  813 FORMAT(/1X,'2J P',3X,'(2S+1) L          E1C(RYD)')
  814 FORMAT(/1X,'(2S+1) L  P        E1C(RYD)')
  816 FORMAT(/1X,'(2S+1) L  P')
  817 FORMAT(I6,2I3,3X,F13.6)
  818 FORMAT(/1X,'2J P',3X,'(2S+1) L')
  819 FORMAT(I3,I2,3X,I5,I3,3X,F15.8)
  829 FORMAT(/' *** WARNING: YOUR INPUT TOLB IS LARGER THAN THE',
     X ' MINIMUM TARGET SPLITTING:',1P2E10.3/' *** RECOMMEND',
     X ' UNSETTING TOLB AND LET CODE DETERMINE IT!'/)
  837 FORMAT(I6,I3,3X,F13.6)
  992 FORMAT(3I2,I5,I5,F18.6,3X,A4)
  993 FORMAT(2I2,1X,I3,I2,2I5,3X,F15.8,3X,A4)
 1000 FORMAT(20A4)
 1001 FORMAT(1X,50('-'),'ADASRR',50('-')//1X,20A4//1X,50('-'),'(V3.12)'
     X,50('-')//)
 1002 FORMAT(' ****INPUT CODE ERROR: ONLY /CA/, /LS/ OR /IC/ ARE'
     X  ,' ALLOWED, WHILE YOUR INPUT IS "',A4,'"')
 1005 FORMAT('C',19A4)
C
      END
C
C***********************************************************************
C
      SUBROUTINE CROSSJ(NBINM,NBINR,NMN,LMN,NCUT,LCUT,NRSLMX,DELTAF,PMN0
     X       ,IPIG,EI,IWT,NR1,IPRINT,TOLR,EBDMIN,EBDMAX,PFACT,NMAX,LMAX
     X       ,NCORE,IWS,IWL,IWJ,LCP,BCA,BLS,BIC,NLMAX,LLMAX,JTHETA,TEMP
     X       ,IREL,JCFJ,NZ,NP0,TEAPOT,ICON,TKAPPA,XDRY,NENG37,E37,F37)
cparc                                                               !par
cpar      use comm_interface, only : iam                            !par
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(NDIM0=1500       !INITIAL METASTABLES * FINAL PARENTS
     X         ,NDIM2=12         !INITIAL TARGET METASTABLES
     X         ,NDIM5=150        !FINAL PARENT METASTABLES
     X         ,NDIM10=25000     !N+1 TERMS/LEVELS PER NL, INC CORR ****
     X         ,NDIM13=150000    !N+1 ENERGIES PER NL, EXC CORR ETC ****
     X         ,NDIM14=750       !CONFIGS
     X         ,NDIM15=1001      !HYDROGENIC POST-PROC RAD L-VALUES
     X         ,NDIM17=16000     !TERM/LEVEL RESOLVED FINAL STATES
     X         ,NDIM18=50000     !DUMMY L-RESOLVED FINAL STATES
     X         ,NDIM20=92        !SEQUENCE LABEL
     X         ,NDIM26=150       !NRSLMX NL-VALUES
     X         ,NDIM27=150       !REPRESENTAIVE N-VALUES
     X         ,NDIM25=NDIM27
     X         ,NDIM31=10        !MAX L+1 FOR BUNDLED-NL
     X         ,NDIM33=46        !PI X-SCTNS
     X         ,NDIM37=19        !ADAS/USER/ADF37 TEMPS (MIN 19)
     X         )
C
      PARAMETER (JTHTIC=19)      !NO. NEW ADAS STANDARD TEMPS
C
      PARAMETER (LMAXZ=31)
      PARAMETER (NLIT=60)
      PARAMETER (NLAB=20)
      PARAMETER (MXORB0=50)            !NO. NL READ FROM AS
C
      PARAMETER (IZERO=0)
C
      PARAMETER (DZERO=0.0D0)
      PARAMETER (DONE=1.0D0)
      PARAMETER (DTWO=2.0D0)
      PARAMETER (DFOUR=4.0D0)
      PARAMETER (DTWELV=12.0D0)
      PARAMETER (DKCM=109737.4D0)
      PARAMETER (D1M4=1.D-4)
      PARAMETER (D1M18=1.D-18)
C
c      PARAMETER (TINY=1.D-5)
C      PARAMETER (PI=3.14159265359D0)
C
      PARAMETER (DFSC=DONE/137.03599976D0)
      PARAMETER (DALF=DFSC*DFSC)
      PARAMETER (ALFP=2.568D-18)       !=4*pi*a_0^2*alpha
      PARAMETER (ALF4=1.D18*DALF/4)    != 1.33128D13
      PARAMETER (CONRYK=1.5789D5) !1.578885D5)      !RYDBERGS TO KELVIN
C
      PARAMETER (MXTABC=25)            !SET FOR EDUM IN DATA BELOW
      PARAMETER (MXTABF=41)            !SET FOR EDUM IN DATA BELOW
C
C NO. OF VALUES IN N-MESH (FOR BARE IONS)
C      PARAMETER (NTOP0=59)      !AUTO-GENERATED N-MESH, TO SYNC WITH DR
      PARAMETER (NTOP0=66)                      !2005 RR DATABASE N-MESH
C
      REAL*4 ALF,ALFN,AN,BN,BNL,RNE,RNLE
C
      integer*4 MTEST4,MBLNK4              !keep I*4 for backward compat
C
      INTEGER SS,SSR,QS0,QL0,QSB,QLB,QL,QN,QSP,QLP,QSR,QLR,QNV,QLV
     X,SSZ,QLZ,QSZ,QND,QLD,QSH,QLH
C
      CHARACTER
     X LAB2*2,LSQ*2,INAM*2,JNAM*2
     X,O*3,OP*4
     X,LAB4*4,CMBLNK*4,CMSTAR*4
     X,LAB5*5
     X,O1*8,INDX*8
     X,O1U*9
     X,FILNAM*9,OP1*9
     X,FILNAMP*10,OP1U*10
     X,F101*30,F732*23
      CHARACTER*1 LAB1,CLABL(NLAB),CLIT(-1:NLIT)
C
      CHARACTER NAME*30,DATE*30,COD(20)*4,nam0*2
C
      LOGICAL BPRNTP1,BPRNTP2,BPRNT0,BPRNT1,BPRNT2,BINT,BFAST,BDROP
     X,BBIN,BPASS1,BFORM,BFORMP,BRAD,BCA,BLS,BIC,BEX,BEXP
     X,BHYBRD,BCAH,BRSLF,BTEST,BDOWN,BTWO
     x,bare,bflagp,bnmesh,bskp
C
      DIMENSION
     X IWT(NDIM5),IWS(NDIM5),IWL(NDIM5),IWJ(NDIM5),LCP(NDIM5)
     X,EI(NDIM5)
C
      DIMENSION
     X LIT(0:NLIT),LABL(NLAB),LSQ(NDIM20),NTOP(NTOP0)
     X,EDUMC(MXTABC),EDUMF(MXTABF)
C
      DIMENSION THTIC(JTHTIC),TEMP(*)
C
      DIMENSION E37(NENG37,*),F37(NENG37,*)
C
      ALLOCATABLE
     X ENERGP(:),PCS(:,:)
     X,ENERGP0(:),PCS0(:)
     X,MENGP(:),EDUM(:)
C
      ALLOCATABLE
     X IK(:),SS(:),LL(:),JJ(:),LCF(:)
     X,IWP(:),ENERG(:),IRSOL(:),IRSOL0(:)
     X,JK(:),ITAR(:)
C
      ALLOCATABLE
     X TC(:,:),TCP(:,:),IBN(:)
     X,QN(:),QL(:),QND(:),QLD(:)
C
      ALLOCATABLE
     X LMP(:),QSP(:,:),QLP(:,:)
     X,EII(:),WNP(:)
C
      ALLOCATABLE
     X QS0(:),QL0(:),QSB(:,:),QLB(:,:),LMX(:)
     X,QSH(:,:),QLH(:,:),LMH(:),JJH(:)
     X,WNH(:),JKH(:),ITARH(:)
C
      ALLOCATABLE
     X ALF(:,:),ALFN(:,:,:)
     X,CP(:),CM(:),JDUM(:)
     X,NLSWCH(:,:)
     X,RNE(:,:),RNLE(:,:,:)
     X,FREL(:)
C
      ALLOCATABLE                             !BUNDLED ADF48 ONLY
     X BNL(:,:,:,:,:)
     X,BN(:,:,:,:)
C
      ALLOCATABLE                             !RESOLVED ADF48 ONLY
     X AN(:,:,:)
     X,WNR(:),LMR(:),QSR(:,:),QLR(:,:)
     X,QNV(:),QLV(:),SSR(:)
     X,LLR(:),JJR(:),JVR(:),ITARR(:)
C
      ALLOCATABLE                             !RESOLVED ADF48 ONLY
     X ITARZ(:),SSZ(:),LLZ(:),JJZ(:)
     X,LMZ(:),QSZ(:,:),QLZ(:,:)
C
      COMMON /ECOR/E1C(NDIM5),TOLB,TOLB0
      COMMON /NRBTST/BSP0,INAG,NCP,NLAG,IMESH
      COMMON /COMMS/COD,NAME,DATE,nam0
C
      DATA CLABL /'S','P','D','F','G','H','I','J','K','L','M','N','O'
     X,'P','Q','R','S','T','U','*'/, CMBLNK/'    '/, CMSTAR/'****'/
      DATA CLIT
     X/'*','0','1','2','3','4','5','6','7','8','9','A','B','C','D','E',
     X 'F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U',
     X 'V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k',
     X 'l','m','n','o','p','q','r','s','t','u','v','w','x','y'/
      DATA LSQ
     X/'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE','NA','MG','AL'
     X,'SI','P ','S ','CL','AR','K ','CA','SC','TI','V ','CR','MN','FE'
     X,'CO','NI','CU','ZN','GA','GE','AS','SE','BR','KR','RB','SR','Y '
     X,'ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN','SN','SB','TE'
     X,'I ','XE','CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB'
     X,'DY','HO','ER','TM','YB','LU','HF','TA',' W','RE','OS','IR','PT'
     X,'AU','HG','TL','PB','BI','PO','AT','RN','FR','RA','AC','TH','PA'
     X,'U '/
C
      DATA
     XTHTIC/1.0D1,2.0D1,5.0D1,1.0D2,2.0D2,5.0D2,1.0D3,2.0D3,5.0D3,1.0D4
     X     ,2.0D4,5.0D4,1.0D5,2.0D5,5.0D5,1.0D6,2.0D6,5.0D6,1.0D7/
C
C FOR BARE IONS (& NO AS RUN)
C THIS IS THE AS DEFAULT AUTOMATICALLY GENERATED N-MESH, TO SYNC WITH DR
C NTOP0=59
C      DATA NTOP /1,2,3,4,5,6,7,8,9,10,11,13,16,19,23,28,34,41,49,59,70,
C     X          85,105,125,150,180,220,270,320,390,470,560,680,820,999,
C THIS IS THE HISTORIC 2005 RR DATABASE N-MESH
C NTOP0=66
      DATA NTOP /1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     X           16,18,20,22,25,29,35,39,45,49,55,61,70,81,
     X           100,115,140,162,200,235,300,353,450,535,700,811,999,
     X           1150,1400,1620,2000,2350,3000,3530,4500,5350,7000,8110,
     X           9999,11500,14000,16200,20000,23500,30000,35300,45000,
     X          53500,70000,81100,99999/
C
C DEFAULT MESH, MXTABC=25
      DATA EDUMC/0.0D+0,2.1D-6,4.6D-6,1.0D-5,2.1D-5,4.6D-5,
     X           1.0D-4,2.1D-4,4.6D-4,1.0D-3,2.1D-3,4.6D-3,
     X           1.0D-2,2.1D-2,4.6D-2,1.0D-1,2.1D-1,4.6D-1,
     X           1.0D+0,2.1D+0,4.6D+0,1.0D+1,2.1D+1,4.6D+1,1.0D+2/
C TEST FINER MESH, MXTABF=41
      DATA EDUMF/0.0D+0,1.6D-6,2.5D-6,4.0D-6,6.3D-6,
     X           1.0D-5,1.6D-5,2.5D-5,4.0D-5,6.3D-5,
     X           1.0D-4,1.6D-4,2.5D-4,4.0D-4,6.3D-4,
     X           1.0D-3,1.6D-3,2.5D-3,4.0D-3,6.3D-3,
     X           1.0D-2,1.6D-2,2.5D-2,4.0D-2,6.3D-2,
     X           1.0D-1,1.6D-1,2.5D-1,4.0D-1,6.3D-1,
     X           1.0D+0,1.6D+0,2.5D+0,4.0D+0,6.3D+0,
     X           1.0D+1,1.6D+1,2.5D+1,4.0D+1,6.3D+1,1.0D+2/
C TEST COARSER MESH, MXTABC=10 OR 13
c      DATA EDUMC/0.D0,1.D-6,1.D-5,1.D-4,1.D-3,1.D-2,1.D-1,1.D0,1.D1/
c     X,1.D2,1.D3,1.D4,1.D5/
C
      PI=ACOS(-DONE)
C
C FIX FOR FORTRAN 90 COMPILERS THAT DON'T ALLOW ASSIGNMENT OF CHARACTERS
C TO INTEGER VARIABLES, REQUIRED FOR HISTORIC BACKWARDS COMPATIBILITY
C
      OPEN(90,STATUS='SCRATCH',FORM='FORMATTED')
      WRITE(90,1111)CMSTAR,(CLIT(I),I=0,NLIT)
 1111 FORMAT(A4,80A1)
      BACKSPACE(90)
      READ(90,1111)MSTAR,(LIT(I),I=0,NLIT)
      WRITE(90,1111)CMBLNK,(CLABL(I),I=1,NLAB)
      BACKSPACE(90)
      READ(90,1111)MBLNK,(LABL(I),I=1,NLAB)
      BACKSPACE(90)
      READ(90,1111)MBLNK4
      CLOSE(90)
C
C THESE SWITCHES CONTROL RESOLUTION
C
      bare=.false.
      bflagp=.false.             !flag if master parent not "found".
      bskp=.false.
C
      BHYBRD=NBINR.LE.1         !FINAL PARENT BY CONFIG (=1 FOR NTAR2=0)
      NBINR=IABS(NBINR)
      NBINRM=NBINR-1
      IF(BHYBRD)THEN
        NRX=0
      ELSE
        NRX=NBINM
      ENDIF
      BRSLF=.NOT.BCA.AND..NOT.BHYBRD
C
      LLMAX0=LLMAX                                       !HOLD FOR PRINT
      IF(LLMAX.LT.0)THEN
        LLMAX0=-9999
        LLMAX=-LLMAX
      ENDIF
C
      MXLL=-1
C
      NSYM=1
      IF(BLS)THEN
        NSYM=2
      ELSE
        BLS=BCA
      ENDIF
C
      JTEMP=ABS(JTHETA)
C
      IF(TEMP(1).GT.DZERO)THEN
C ALREADY IN RYD NOW
      ELSE
        IF(JTEMP.GT.0)THEN
          JTEMP=MIN(JTEMP,JTHTIC)
          IF(IREL.GT.0.AND.JTEMP.EQ.JTHTIC)JTEMP=JTEMP-2
          JTHETA=SIGN(JTEMP,JTHETA)
        ENDIF
        IF(JTEMP.GT.NDIM37)THEN                       !FOR NON-ALLOCATED
          WRITE(6,*)'TOO MANY TEMPS; INCREASE NDIM37 TO:',JTEMP
          STOP 'ERROR: TOO MANY TEMPS; INCREASE NDIM37 TO JTEMP'
        ENDIF
      ENDIF
C
C WOULD NEED TO PRE-READ TO GET NENGP SO CAN USE EXACT ALLOCATION AND
C THERE IS THE POSSIBLE INCREASE IN ENGP DUE TO TOP-UP IN E (SEE BFIRST)
      JTMAX=NDIM33                      !NOT JTEMP AS WE STORE PCS FIRST
C
      ALLOCATE(
     X ENERGP(NDIM33),PCS(NDIM33,NDIM13)
     X,ENERGP0(NDIM33),PCS0(NDIM33)
     X,MENGP(NDIM33),EDUM(NDIM33)
     X,STAT=IERR)
      IF(IERR.NE.0)STOP 'ALLOCATION 1 FAILS'
C
      ALLOCATE(
     X IK(NDIM13),SS(NDIM13),LL(NDIM13),JJ(NDIM13),LCF(NDIM13)
     X,IWP(NDIM13),ENERG(NDIM13),IRSOL(NDIM13),IRSOL0(NDIM13)
     X,JK(NDIM10),ITAR(NDIM10)
     X,STAT=IERR)
      IF(IERR.NE.0)STOP 'ALLOCATION 2 FAILS'
C
      ALLOCATE(
     X TC(JTMAX,NDIM2),TCP(JTMAX,NDIM2),IBN(0:NDIM27)
     X,QN(NDIM26),QL(NDIM26),QND(NDIM26),QLD(NDIM26)
     X,STAT=IERR)
      IF(IERR.NE.0)STOP 'ALLOCATION 3 FAILS'
C
      ALLOCATE(
     X LMP(NDIM5),QSP(NDIM5,10),QLP(NDIM5,10)
     X,EII(NDIM5),WNP(NDIM5)
     X,STAT=IERR)
      IF(IERR.NE.0)STOP 'ALLOCATION 4 FAILS'
C
      ALLOCATE(
     X QS0(10),QL0(10),QSB(NDIM14,10),QLB(NDIM14,10),LMX(NDIM14)
     X,QSH(NDIM14,10),QLH(NDIM14,10),LMH(NDIM14),JJH(NDIM14)
     X,WNH(NDIM14),JKH(NDIM14),ITARH(NDIM14)
     X,STAT=IERR)
      IF(IERR.NE.0)STOP 'ALLOCATION 5 FAILS'
C
      ALLOCATE(
     X ALF(JTMAX,NDIM2),ALFN(JTMAX,NDIM27,NDIM2),FREL(JTMAX)
     X,CP(NDIM15),CM(NDIM15),JDUM(NDIM15)
     X,NLSWCH(NDIM31,NDIM25)
     X,RNE(JTMAX,NDIM25),RNLE(JTMAX,NDIM31,NDIM25)
     X,STAT=IERR)
      IF(IERR.NE.0)STOP 'ALLOCATION 6 FAILS'
C
      IF(NBINRM.GT.0)THEN
C
      ALLOCATE(                               !BUNDLED ADF48 ONLY
     X BNL(JTMAX,NDIM31,NDIM27,NSYM,NDIM0)
     X,BN(JTMAX,NDIM27,NSYM,NDIM0)
     X,STAT=IERR)
      IF(IERR.NE.0)STOP 'ALLOCATION 7 FAILS'
C
      ENDIF
C
      IF(BRSLF)THEN
C
      ALLOCATE(                               !RESOLVED ADF48 ONLY
     X AN(JTMAX,NDIM17,NDIM2)
     X,WNR(NDIM17),LMR(NDIM17),QSR(NDIM17,10),QLR(NDIM17,10)
     X,QNV(NDIM17),QLV(NDIM17),SSR(NDIM17)
     X,LLR(NDIM17),JJR(NDIM17),JVR(NDIM17),ITARR(NDIM17)
     X,STAT=IERR)
      IF(IERR.NE.0)STOP 'ALLOCATION 8 FAILS'
C
      ALLOCATE(                               !RESOLVED ADF48 ONLY
     X ITARZ(NDIM18),SSZ(NDIM18),LLZ(NDIM18),JJZ(NDIM18)
     X,LMZ(NDIM18),QSZ(NDIM18,10),QLZ(NDIM18,10)
     X,STAT=IERR)
      IF(IERR.NE.0)STOP 'ALLOCATION 9 FAILS'
C
      ENDIF
c
      iam0=0
cpar      iam0=iam                                                  !par
C
C
C***********
C INITIALIZE
C***********
C
      INR1=IABS(NR1)
      BRAD=NR1.GT.0
      NCMX0=MAX(0,NCORE)
      IFLAGE=0
      IFLAGB=0
      IFLGL1=0
      IFLGL2=0
      IFLAGR=0
      NFLAG2=NDIM5+1
      NVINT=100
      IF(NCUT.Le.999.AND.NVINT.LT.NCUT)NVINT=NCUT
C      NMN0=NMN
      IB00=0
C
C INITIALIZE NL FOR OUTER ELECTRON STABILIZATION TO STANDARD ORDER+MXORB
C
      J=MXORB0
      DO N=1,NRSLMX
        DO L=1,N
          J=J+1
          IF(J.GT.NDIM26)THEN
            NRSLMX=N-1
            WRITE(6,*)'***WARNING, NRSLMX REDUCED TO:',NRSLMX
            WRITE(6,*)'   INCREASE NDIM26 IF NEED BE'
            GO TO 424
          ENDIF
          QN(J)=N
          QL(J)=L-1
        ENDDO
      ENDDO
C
 424  IF(MXORB0.GT.NLIT)WRITE(6,*)'***WARNING: MIGHT NOT BE ABLE '
     X                  ,' TO DECODE ORBITAL, INCREASE LIT SPEC.'
C
C DUPLICATE OF POSTP BUT GUESS CATCHES NDIM2/5 SET DIFFERENT
C
      IF(NBINM.GT.NDIM2)THEN
        WRITE(6,847)NBINM
 847    FORMAT(/' INCREASE NDIM2 TO AT LEAST',I3)
        STOP 'ERROR: INCREASE NDIM2'
      ENDIF
      IF(NBINR.GT.NDIM5)THEN
        WRITE(6,848)NBINR
 848    FORMAT(/' INCREASE NDIM5 TO AT LEAST',I3)
        STOP 'ERROR: INCREASE NDIM5'
      ENDIF
C
C PACK
C
      IBNMX=NBINRM*NBINM
      IF(IBNMX.GT.NDIM0)THEN
        NBINRM=NDIM0/NBINM
        IF(NBINRM.LT.NBINM)THEN
          IBNMX=NBINM*NBINM
          WRITE(6,711)IBNMX
 711      FORMAT(' DIMENSION EXCEEDED IN SR.CROSSJ, INCREASE NDIM0 TO'
     X           ,' AT LEAST',I4)
          STOP 'ERROR: DIMENSION EXCEEDED IN SR.CROSSJ, INCREASE NDIM0'
        ELSE
          WRITE(6,721)NBINRM+1,IBNMX
 721      FORMAT(' REDUCING NTAR2 TO',I5,' BECAUSE OF DIMENSIONS.'
     X          /' INCREASE NDIM0 TO',I5,' IF YOU *REALLY* NEED NTAR2.')
        ENDIF
      ENDIF
C
      NBINP=NBINM
      IF(.NOT.BHYBRD)NBINP=MAX(NBINP,NBINRM)
C
      BPRNTP2=IPRINT.EQ.2
      BPRNTP1=IPRINT.EQ.1
      BPRNT0=IPRINT.EQ.0
      BPRNT1=IPRINT.GE.-1
      BPRNT2=IPRINT.GE.-2
C
C********
C ZEROISE
C********
C
      L0=0
      DO L=1,NBINM                          !INITIAL STATES
        DO M=1,NBINRM                       !FINAL PARENTS
          L0=L0+1
          DO NS=1,NSYM                      !2 SPIN SYSTEMS
            DO I=1,NDIM27
              DO J=1,JTMAX
                BN(J,I,NS,L0)=DZERO          !BUNDLE-N
              ENDDO
            ENDDO
          ENDDO
          IF(NLMAX.GT.0)THEN
            DO NS=1,NSYM
              DO I=1,NDIM27
                DO K=1,NDIM31
                  DO J=1,JTMAX
                    BNL(J,K,I,NS,L0)=DZERO   !BUNDLE-NL
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
C
        IF(BRSLF)THEN
          DO I=1,NDIM17
            DO J=1,JTMAX
              AN(J,I,L)=DZERO                !RESOLVED
            ENDDO
          ENDDO
        ENDIF
C
        DO I=1,NDIM27
          DO J=1,JTMAX
            ALFN(J,I,L)=DZERO                !TOTAL-N
          ENDDO
        ENDDO
        DO J=1,JTMAX
          ALF(J,L)=DZERO                     !TOTAL
        ENDDO
      ENDDO
C
C**************************************************************
C CASE OF BARE ION (OR PURE HYDROGENIC) NO on OR opn FILES READ
C**************************************************************
C
      IF(NZ.GT.0)THEN
        bare=.true.
        NE=1
        NZ0=NZ
        DZ=NZ*NZ
C        DZ0=DZ                             !TEST
        LV0=MAX(0,LMN)
        LV0=LV0-1
        LV00=LV0-1
        NCMX0=0                             !NOT USED
        WNP0=DZERO
        EIONMN=DZERO
        TOLB0=TOLB
        IB0=0
        IBN(0)=0
        IFLGL1=NBINM
        IFLGL2=NBINRM
        DO I=1,NBINM
          WNP(I)=DZERO
          LMP(I)=0
          DO J=1,5
            QSP(I,J)=MBLNK
          ENDDO
        ENDDO
        IF(IMESH.LE.0)THEN                  !DEFAULT
          MXETAB=MXTABC
          NENGP=MIN(MXETAB,NDIM33)
          IF(NP0.GT.0)NENGP=MIN(NENGP,NP0)
          DO N=1,NENGP
            ENERGP(N)=EDUMC(N)*DZ
          ENDDO
        ELSE                                !TEST FINER
          MXETAB=MXTABF
          NENGP=MIN(MXETAB,NDIM33)
          IF(NP0.GT.0)NENGP=MIN(NENGP,NP0)
          DO N=1,NENGP
            ENERGP(N)=EDUMF(N)*DZ
          ENDDO
        ENDIF
        IF(MXETAB.GT.NDIM33)THEN
          WRITE(6,*)'*** INCREASE NDIM33 TO AT LEAST',MXETAB
          STOP 'ERROR: *** INCREASE NDIM33'
        ENDIF
        IF(NENGP.LT.MXETAB)THEN
          WRITE(*,*)'*** WARNING: RR OF BARE ION, PI XSCTN EXTRAPOLATED
     XFROM TOO LOW AND ENERGY???'
          WRITE(6,*)'*** WARNING: RR OF BARE ION, PI XSCTN EXTRAPOLATED
     XFROM TOO LOW AND ENERGY???',ENERGP(NENGP),' INSTEAD OF' ,DZ*EDUMC(
     XMXETAB)
        ENDIF
C
        NENGT=NENGP
        NENGPT=NENGP
C
        WRITE(6,120)(ENERGP(M),M=1,NENGP)
        IF(BPRNT2)WRITE(6,90)
C
        IF(NR1.EQ.999)THEN
          IF(NCORE.GT.0)THEN
            NR1=NCORE+1
          ELSEIF(NMN.GT.0)THEN
            NR1=NMN
          ELSE
            NR1=1
          ENDIF
        ELSEIF(NR1.LE.0)THEN
          WRITE(6,*)'***ERROR, RR OF BARE ION, BUT NR1=',NR1
          STOP 'ERROR: NR1 INCORRECT FOR RR OF BARE ION'
        ENDIF
        INR1=NR1
C
        NRSOL=0
C
        GO TO 1001
      ENDIF
C
C******************************************************
C POSSIBLE UNIT NOS TO BE CHECKED FOR DATA: READ NV, LV
C******************************************************
C
      MR=70
      MRU=MR
      MRP=80
      MRPU=MRP
C
      IF(NBINRM.EQ.0)GO TO 200
C
C CHECK FOR TARGET CF DEF. CASE OF HYBRID RESOLUTION
C (ALSO USE CF'S ONLY WHEN INCONSISTENT IN CHANNEL DEF.)
C WE NOW REQUIRE COUPLING SCHEMES TO MATCH.
C
      NCFR=0
      BFORM=.FALSE.
      FILNAM='o_str'
      INQUIRE(FILE=FILNAM,EXIST=BEX)
      IF(BEX)THEN
        OPEN(MR,FILE=FILNAM)
        READ(MR,38,END=322)MDUM1,MDUM2
        BFORM=.TRUE.
        READ(MR,'(3X,A1)',END=322)LAB1
        IF(LAB1.EQ.'C')THEN
          F101='(I3,12X,I2,6X,I2,4X,50(I3,I2))'
        ELSE
          F101='(I5,10X,I2,6X,I2,4X,50(I3,I2))'
        ENDIF
        BACKSPACE(MR)
      ENDIF
 322  IF(.NOT.BFORM)THEN
        FILNAM='ou_str'
        INQUIRE(FILE=FILNAM,EXIST=BEX)
        IF(BEX)THEN
          OPEN(MRU,FILE=FILNAM,FORM='UNFORMATTED')
          READ(MRU,END=1000)MDUM1,MDUM2
        ELSE                                           !EXIT STAGE RIGHT
          IF(BHYBRD)THEN                               !NEEDED
cif bare target then
            bare=.true.                         !flag potential bare ion
            ncfr=1
            do j=1,10
              qsh(1,j)=mblnk
              qlh(1,j)=0
            enddo
            jkh(1)=1
            jjh(1)=1
            lmh(1)=0
            wnh(1)=0
            iws(1)=1
            iwl(1)=0
            iwj(1)=0
            wnp(1)=0
            lmp(1)=0
            do j=1,5
              qsp(1,j)=mblnk
            enddo
            go to 200
celse
c           WRITE(6,*)'NO TARGET DATA ON ',FILNAM,'-REQUIRED FOR HYBRID'
c           STOP 'ERROR: NO TARGET DATA ON FILE - REQUIRED FOR HYBRID'
cendif
          ELSE                     !ASSUME NOT NEEDED - NO NECOR+BUNDLED
            GO TO 200              !OR FOR CASE OF BAD LABEL (HOPEFULLY)
          ENDIF
        ENDIF
      ENDIF
      IF(BHYBRD.AND.MDUM2.GE.0)THEN                      !SUSPICIOUS...
        WRITE(6,*)'TARGET DATA HAS LV.GE.0???'
        STOP 'ERROR: TARGET DATA HAS LV.GE.0 - INVALID FOR HYBRID'
      ENDIF
      IF(BFORM)READ(MR,F101)NCFR,NZOLD,NEOLD              !CF HEADER
      IF(.NOT.BFORM)READ(MRU)NCFR,NZOLD,NEOLD
      NEOLD=NEOLD+1
      IF(NCFR.GT.NDIM14)THEN
        WRITE(6,*)'*** INCREASE NDIM14 CONFIGS FOR HYBRID TO',NCFR
        STOP 'ERROR: INCREASE NDIM14'
      ENDIF
      IF(BHYBRD)THEN
        IF(NCFR.GT.NBINRM)THEN
          WRITE(6,*)'*** WARNING: YOUR TARGET ',FILNAM,' FILE HAS',NCFR,
     X    ' CONFIGS BUT YOU HAVE ONLY REQUESTED',NBINRM,' FINAL PARENTS'
        ELSEIF(NCFR.LT.NBINRM)THEN!ASSUME TARGET CORRECT,QUIETLY REDUCE
          NBINRM=NCFR
          NBINR=NBINRM+1
        ENDIF
        IF(LCP(1).EQ.0)THEN
          WRITE(6,*)'*** ERROR: NEED TARGET CF NOS IN adasin FOR HYBRID'
     X             ,' OPERATION'
          STOP'*** ERROR: NEED TARGET CF NOS IN adasin FOR HYBRID MODE'
        ENDIF
      ENDIF
C
      DO I=1,NCFR                                        !CFGS
        IF(BFORM)READ(MR,179,END=1002)IN,NGR,MA0,MB0
     X                              ,(QS0(J),QL0(J),J=1,10)
        IF(.NOT.BFORM)READ(MRU,END=1002)IN,NGR,MA0,MB0
     X                              ,(QS0(J),QL0(J),J=1,10)
        IF(IN.LT.0)THEN
          IF(BHYBRD)THEN
            STOP 'INVALID TARGET FILE - REQUIRED FOR HYBRID'
          ELSE
            NCFR=0                                !HOPE WE DON'T NEED IT
            GO TO 200
          ENDIF
        ELSE
          JJH(I)=NGR
          WNH(I)=DZERO                             !INITIALIZE
          jkh(i)=0                                !checksum
        ENDIF
        DO 166 J=1,10
          QSH(I,J)=MBLNK
          IF(QL0(J).EQ.MBLNK)GO TO 166
          LMH(I)=J
          M=MOD(QS0(J),50)
          IF(M.GT.0)QSH(I,J)=LIT(M)
          DO K=1,NLIT
            IF(QL0(J).EQ.LIT(K))GO TO 199
          ENDDO
          QLH(I,J)=0
          GO TO 166
 199      QLH(I,J)=K
 166    ENDDO
      ENDDO
C
      IF(BHYBRD)THEN
        IF(BFORM)THEN                          !SKIP AA HEADERS
          DO I=1,3
            READ(MR,103,END=1002)
          ENDDO
        ELSEIF(.NOT.BFORM)THEN                 !THERE SHOULD BE NO RATES
          READ(MRU,END=1002)
          READ(MRU,END=1002)
        ENDIF
        IF(BFORM)READ(MR,121,END=1002)NENG,ECORE
        IF(.NOT.BFORM)READ(MRU,END=1002)NENG,ECORE
C
        IF(BFORM)READ(MR,106,END=1002)MTEST4
        IF(.NOT.BFORM)READ(MRU,END=1002)MTEST4
 106    FORMAT(21X,A4)
        IF(BFORM)BACKSPACE(MR)
        IF(.NOT.BFORM)BACKSPACE(MRU)
C
        BCAH=MTEST4.EQ.MBLNK4
        if(bcah.and.neng.ne.ncfr)stop 'confusion, is o_str CA or not?'
        IF(BCAH.NEQV.BCA)THEN
          WRITE(6,375)BCAH,BCA
 375    FORMAT('*** HYBRID ERROR: COUPLING SCHEME INCONSISTENCY BETWEEN'
     X        ,' o_str AND o1 ETC FILES; ONE IS CA THE OTHER NOT:',2L3)
          STOP '*** HYBRID ERROR: COUPLING SCHEME INCONSISTENCY'
        ENDIF
C
        IF(BFORM)READ(MR,105,END=1002)MTEST4
        IF(.NOT.BFORM)READ(MRU,END=1002)MTEST4
C
        BTEST=MTEST4.NE.MBLNK4
        IF(BTEST.NEQV.BIC)THEN
          WRITE(6,376)BTEST,BIC
 376    FORMAT('*** HYBRID ERROR: COUPLING SCHEME INCONSISTENCY BETWEEN'
     X        ,' o_str AND o1 ETC FILES; ONE IS IC THE OTHER NOT:',2L3)
          STOP '*** HYBRID ERROR: COUPLING SCHEME INCONSISTENCY'
        ENDIF
        IF(.NOT.BTEST.NEQV.BLS)THEN
          BTEST=.NOT.BTEST
          WRITE(6,377)BTEST,BLS
 377    FORMAT('*** HYBRID ERROR: COUPLING SCHEME INCONSISTENCY BETWEEN'
     X        ,' o_str AND o1 ETC FILES; ONE IS LS THE OTHER NOT:',2L3)
          STOP '*** HYBRID ERROR: COUPLING SCHEME INCONSISTENCY'
        ENDIF
C
        DO I=1,NENG                            !GET ENERGIES
          IF(BFORM)READ(MR,123,END=1002)MDUM,MT,IS,IL,IJ,IC,E
          IF(.NOT.BFORM)READ(MRU,END=1002)MDUM,MT,IS,IL,IJ,IC,E
          IS=IABS(IS)                       !CASE CORR.
          IF(BCAH)THEN
            IF(BFORM)THEN
              IS=100000*MT+IS               !AS 2I5
              IS=IS*(2*IL+1)                !UNLIKELY CASE 1 TERM PER CF
              if(is.ne.jjh(ic))stop 'is.ne.jjh test'
              JJH(IC)=IS                    !AS ONLY I5 IN CF READ
            ENDIF
            JKH(I)=IC                       !EO->SO
            WNH(IC)=-E-ECORE
          ELSE
            IF(BTEST)THEN
              JW=IJ+1
            ELSE
              JW=IS*(2*IL+1)
            ENDIF
            WNH(IC)=WNH(IC)+JW*E
            jkh(ic)=jkh(ic)+jw              !checksum
          ENDIF
        ENDDO
C
        IF(.NOT.BCAH)THEN                   !NEED EO INDEX
          DO I=1,NCFR
            WNH(I)=-WNH(I)/jkh(i)-ECORE
c            if(jkh(i).ne.JJH(I))then        !only if no corr.
c              write(6,*)'stat weight checksum failure:',i,jkh(i),jjh(i)
c              stop 'stat weight checksum failure'
c            endif
            JJH(I)=jkh(i)                   !CASE CORR.
          ENDDO
          CALL HPSRTI(NCFR,WNH,JKH)
        ENDIF
      ENDIF
C
 200  CONTINUE
C
C*****************************************
C INITIALIZE REGULAR UNITS FOR RR (NL, LV)
C*****************************************
C (DO NOT MIX FORMATTED AND FORMATTED BETWEEN O-AND OP-FILES)
C
      BFORM=.FALSE.
      BTWO=.FALSE.
      JFILE=1
      bnmesh=.false.
C
      IF(BLS)THEN
        O='ols'
        OP='opls'
      ENDIF
      IF(BCA)THEN
        O='oca'
        OP='opca'
      ENDIF
      IF(BIC)THEN
        O='oic'
        OP='opic'
      ENDIF
C
      IFILE=1
cpar      ifile=ifile+iam                                           !par
      IC1=IFILE/10
      IC2=IFILE-10*IC1
      IC0=ICHAR('0')
      IC1=IC1+IC0
      IC2=IC2+IC0
      INAM=CHAR(IC1)//CHAR(IC2)
      IF(IFILE.LT.10)THEN                            !LOOK FOR EX-SERIAL
        O1='o'//CHAR(IC2)
        O1U='o'//CHAR(IC2)//'u'
        OP1='op'//CHAR(IC2)
        OP1U='op'//CHAR(IC2)//'u'
      ELSE                                          !HAS NO LEADING ZERO
        O1='o'//INAM
        O1U='o'//INAM//'u'
        OP1='op'//INAM
        OP1U='op'//INAM//'u'
      ENDIF
C
 330  FILNAM=O1
      FILNAMP=OP1
      INQUIRE(FILE=FILNAM,EXIST=BEX)
      INQUIRE(FILE=FILNAMP,EXIST=BEXP)
      IF(BEX.NEQV.BEXP)THEN
        WRITE(6,1007)
 1007   FORMAT(/'*** ERROR: opn/u FILE EXISTS WITH NO CORRESPONDING'
     X  ,' on/u FILE')
        STOP '*** ERROR: opn FILE EXISTS WITH NO CORRESPONDING on FILE'
      ENDIF
      IF(BEX)THEN
        BFORM=.TRUE.
        OPEN(MR,FILE=FILNAM)
        READ(MR,38,END=332)MDUM1,MDUM2
        BFORM=.TRUE.
        READ(MR,'(3X,A1)',END=332)LAB1
        IF(LAB1.EQ.'C')THEN
          F101='(I3,12X,I2,6X,I2,4X,50(I3,I2))'
        ELSE
          F101='(I5,10X,I2,6X,I2,4X,50(I3,I2))'
        ENDIF
        BACKSPACE(MR)
        BACKSPACE(MR)
        OPEN(MRP,FILE=FILNAMP)
      ENDIF
 332  BFORMP=BFORM                                   !NO MIXED FILES NOW
      IF(.NOT.BFORM)THEN
        FILNAM=O1U
        FILNAMP=OP1U
        INQUIRE(FILE=FILNAM,EXIST=BEX)
        INQUIRE(FILE=FILNAMP,EXIST=BEXP)
        IF(BEX.NEQV.BEXP)THEN
          WRITE(6,1007)
          STOP '*** ERROR:opn FILE EXISTS WITH NO CORRESPONDING on FILE'
        ENDIF
        IF(BEX)THEN
          OPEN(MRU,FILE=FILNAM,FORM='UNFORMATTED')
          OPEN(MRPU,FILE=FILNAMP,FORM='UNFORMATTED')
        ELSE                                   !NOW LOOK FOR EX-PARALLEL
          IF(IFILE.GT.0)THEN                 !CHECK FOR SINGLE EXT FIRST
            IFILE=0
cpar            ifile=ifile+iam                                     !par
            IC1=IFILE/10
            IC2=IFILE-10*IC1
            IC0=ICHAR('0')
            IC1=IC1+IC0
            IC2=IC2+IC0
            INAM=CHAR(IC1)//CHAR(IC2)
            O1=O//INAM
            O1U=O//'u'//INAM
            OP1=OP//INAM
            OP1U=OP//'u'//INAM
cpar            ifile=-ifile                                        !par
            GO TO 330
          ELSEIF(JFILE.GT.0)THEN               !NOW CHECK FOR DOUBLE EXT
            write(6,*)
cpar     x     'padasrr not yet coded to handle double filename ext'!par
cpar            stop                                                !par
cpar     x       'padasrr not yet coded to handle double file ext'  !par
            BTWO=.TRUE.
            JFILE=0
cpar            jfile=jfile+jam                                     !par
            JC1=JFILE/10
            JC2=JFILE-10*JC1
            JC0=ICHAR('0')
            JC1=JC1+JC0
            JC2=JC2+JC0
            JNAM=CHAR(JC1)//CHAR(JC2)
            O1=O//JNAM//'.'//INAM
            O1U=O//'u'//JNAM//'.'//INAM
            OP1=OP//JNAM//'.'//INAM
            OP1U=OP//'u'//JNAM//'.'//INAM
cpar            jfile=-jfile                                        !par
            inquire(file='n-mesh.dat',exist=bex)      !require AS n-mesh
            if(bex)then
              bnmesh=.true.
            else
              write(6,*)
     x             ' AS nprocperl.gt.1 run requires an n-mesh.dat file'
              stop ' AS nprocperl.gt.1 run requires an n-mesh.dat file'
            endif
            GO TO 330
          ELSE                                 !EXIT STAGE RIGHT
            WRITE(6,*)'NO RATE INPUT DATA ON FILE o1 OR o1u!!!'
            WRITE(6,*)'NO XSCTN INPUT DATA ON FILE op1 OR op1u!!!'
            STOP 'ERROR: NO RATE/XSCTN INPUT DATA FOUND!!!'
          ENDIF
        ENDIF
      ENDIF
C
cpar      ifile=-1               !flag one file per proc            !par
C
C***********
C INITIALIZE
C***********
C
      BBIN=EI(1).GE.DZERO
      BPASS1=.TRUE.
      IF(JCFJ.LE.0)JCFJ=999
      BFAST=JCFJ.EQ.999.AND.EBDMIN.LT.-1.D29.AND.EBDMAX.GT.1.D29
     X      .AND.TOLB.LT.DZERO
      NZOLD=0
      NEOLD=0
      NEPOLD=0
      INIT=0
C
      IBX=MAX(NRSLMX,NLMAX)
      IF(IBX.GT.NDIM27)THEN
        WRITE(6,399)
        STOP 'ERROR: SR.CROSSJ TOO MANY N-STATES, INCREASE NDIM27'
      ENDIF
      IB0=IBX
      DO I=0,IB0
        IBN(I)=I
      ENDDO
c
      if(bnmesh)then
        open(80,file='n-mesh.dat')
 311    read(80,*,end=312)nv
        DO I=1,IB0
          IF(NV.EQ.IBN(I))GO TO 311
        ENDDO
        IF(NV.LT.IBN(IBX))THEN
          WRITE(6,390)NV,IBN(IBX)
          STOP 'ERROR: NV MUST HAVE NO GAPS UP TO NRSLMX'
        ENDIF
        IF(NV.LT.IBN(IB0))THEN
          WRITE(6,391)NV
          STOP 'ERROR: NV NOT FOUND'
        ENDIF
        IB0=IB0+1
        I=IB0
        IF(I.GT.NDIM27)THEN
          WRITE(6,399)
          STOP 'ERROR: SR.CROSSJ TOO MANY N-STATES, INCREASE NDIM27'
        ENDIF
        IBN(I)=NV
        IB00=MAX(IB00,I)                           !MAX (INDEX) READ
        go to 311
 312    close(80)
      endif
C
      NRSOL=0
      LV00=-1
C
C ENTRY POINT FOR NEW UNIT
C
 331  NV00=0
      NV0=100000
      LV0=-1
      MXORB=min(30,MXORB0)
C
 310  NV=0
      IF(BFORM)READ(MR,38,END=1000)NV,LV
      IF(.NOT.BFORM)READ(MRU,END=1000)NV,LV
      IF(BFORMP)READ(MRP,38,END=1002)NVP,LVP
      IF(.NOT.BFORMP)READ(MRPU,END=1002)NVP,LVP
  38  FORMAT(5X,I5,5X,I5)
C
      IF(NV.NE.NVP.OR.LV.NE.LVP)THEN
        WRITE(6,*)'MIS-MATCH BETWEEN on and opn FILES FOR NV,LV: '
     X            ,NV,LV,NVP,LVP
        STOP 'MIS-MATCH BETWEEN on and opn FILES FOR NV,LV: '
      ENDIF
C
      IF(NV.EQ.0.AND.NV0.EQ.100000)THEN
        NV=INR1
        LV=-1
      ENDIF
      IF(LV.GT.LMAXZ-1.AND.LV.NE.999)THEN
        WRITE(6,*)'INCREASE PARAMETER LMAXZ TO: ',LV+1
        STOP 'INCREASE PARAMETER LMAXZ'
      ENDIF
      IF(LV.LT.0.AND.NR1.LT.0)INIT=INR1
      IF(NV.GT.0.AND.INR1.NE.999.AND.LV.GE.0)THEN
        IF(NV.LT.INR1.AND.NV.GE.NMN)THEN
          WRITE(6,39)NV,NR1
  39      FORMAT(' ERROR IN CROSSJ: NV MUST BE .GE. ABS(NR1):',2I6)
          STOP ' ERROR IN CROSSJ: NV MUST BE .GE. ABS(NR1):'
        ENDIF
        IF(NV0.EQ.100000.AND.NR1.GT.0.AND.NV.GT.NR1)THEN
          IF(BLS)WRITE(6,388)NR1,NV-1,NR1
 388      FORMAT(' NOTE: TO OBTAIN NON-HYDROGENIC ENERGY LEVELS FOR N='
     X    ,I2,' TO',I2,' RE-RUN AUTOSTRUCTURE FROM NMIN=',I2)
          IF(BIC)WRITE(6,389)NR1,NV-1,NR1
 389      FORMAT(' NOTE: TO OBTAIN TERM LABELS AND NON-HYDROGENIC'
     X    ,' ENERGY LEVELS FOR N='
     X    ,I2,' TO',I2,' RE-RUN AUTOSTRUCTURE FROM NMIN=',I2)
        ENDIF
      ENDIF
      IF(LV.LT.0.AND.LV00.GE.0)THEN
        WRITE(6,*)'***ERROR: RE-ORDER INPUT FILES on(u) ETC SO THAT'
     X           ,'   EQUIVALENT ELECTRON FILES COME FIRST***'
        STOP '***ERROR: RE-ORDER INPUT FILES on(u)'
      ENDIF
C
      IF(LV0.LT.0.OR.LV.LT.0)THEN
        NRSOLZ=0
        WNP0=-1.0D0
        KFPM=0
        IF(INIT.EQ.0)EIONMN=DZERO
        IF(BBIN)EI(1)=1.0D0
      ENDIF
      IF(NV00.EQ.0)NV00=NV
      IF(NV0.EQ.100000)GO TO 70
      BEX=.FALSE.
      BPASS1=.FALSE.
      IF(NV.GT.0.AND.LV.EQ.LV0)GO TO 37
C
  91  IF(NV.EQ.0)GO TO 1000
      IF(LV.GT.LCUT.AND.LV.NE.999)GO TO 1000
C
C**************
C START A NEW L
C**************
C
  70  LV0=LV
      NV0=NV-1
      NV00=NV
      LV00=LV0-1
      IF(LV.GE.0)THEN
        MXLL=MAX(MXLL,LV)                        !+0 AS LV IS LOWER
      ELSE
        MXLL=MAX(MXLL,NV-1)
      ENDIF
C
      DO M=1,NDIM10         !INITIALIZE FOR T OMITTED BY ANY AS BUNDLE
        JK(M)=0
      ENDDO
C
C**************
C START A NEW N
C**************
C
  37  BINT=.FALSE.
      IF(NV.LE.NCUT.AND.NV.GE.NMN.AND.LV.GE.LMN)GO TO 85
      BINT=.TRUE.
      IF(LV.LT.LCUT.OR.NV.LT.NMN)GO TO 75
      LV=LV+1
      GO TO 91
  85  CONTINUE
  75  CONTINUE
      DO I=1,IB0
        IF(NV.EQ.IBN(I))GO TO 60
      ENDDO
      IF(NV.LT.IBN(IBX))THEN
        WRITE(6,390)NV,IBN(IBX)
 390    FORMAT('NV=',I5,' MUST HAVE NO GAPS UP TO NRSLMX=',I3)
        STOP 'ERROR: NV MUST HAVE NO GAPS UP TO NRSLMX'
      ENDIF
      IF(NV.LT.IBN(IB0))THEN
        WRITE(6,391)NV
 391    FORMAT('NV=',I5,' NOT FOUND - LIKELY CAUSE LV .GE. LAST ',
     X  'SEQUENTIAL N, REDUCE LMAX OR RERUN AUTOSTRUCTURE')
        STOP 'ERROR: NV NOT FOUND'
      ENDIF
      IB0=IB0+1
      I=IB0
      IF(I.GT.NDIM27)THEN
        WRITE(6,399)
 399    FORMAT(' SR.CROSSJ TOO MANY N-STATES, INCREASE NDIM27')
        STOP 'ERROR: SR.CROSSJ TOO MANY N-STATES, INCREASE NDIM27'
      ENDIF
      IBN(I)=NV
  60  IB=I                                       !CURRENT INDEX
      IB00=MAX(IB00,I)                           !MAX (INDEX) READ
      NV0=NV
C
C************************************
C READ HEADER, AND MAYBE ORBITAL CODE
C************************************
C
      DO I=1,MXORB
        QND(I)=0
      ENDDO
      IF(BFORM)THEN
 299    READ(MR,F101,END=1002) NCFD,NZ0D,NED,(QND(I),QLD(I),I=1,MXORB)
        if(kfpm.eq.0.and.qnd(mxorb).ne.0.and.mxorb.lt.mxorb0)then
          mxorb=mxorb+20
          mxorb=min(mxorb,mxorb0)
          backspace(mr)
          go to 299
        endif
      ELSE
       READ(MRU,END=1002,ERR=300)NCFD,NZ0D,NED,(QND(I),QLD(I),I=1,MXORB)
       GO TO 302
 300   IF(BEX)THEN                               !START OF A FILE
         REWIND(MRU)                             !SO REWIND
         READ(MRU)
         READ(MRU)
       ELSE
         STOP 'UNABLE TO READ ORBITAL HEADER...' !SHOULD NOT GET HERE
       ENDIF
      ENDIF
C
 302  NCF=NCFD                     !NOT EOF SO SAFE TO RELABEL
      IF(NCF.EQ.0)then
        write(6,*)'ncf=0, why??'
        stop 'ncf=0, why??'
cRETURN
      endif
c
      NZ0=NZ0D
      NE=NED
      NZ=NZ0-NE+1
      DZ=NZ*NZ
      DZ0=NZ0*NZ0
C
      IF(IFLAGB.EQ.0)THEN
        IFLAGB=1
C        TOLB0=TOLB
        IF(TOLB0.LT.DZERO)TOLB=MAX(1.5D-7,1.0D-9*DZ*NZ)
        TOLBE=TOLB
        IF(BFORM)TOLBE=MAX(TOLBE,2.D-6)
        DELTAF=ALFP*DELTAF/DZ
      ENDIF
C
      DO I=1,MXORB
        QN(I)=QND(I)
        QL(I)=QLD(I)
      ENDDO
C
      IF(NZOLD.NE.0.AND.NZ0.NE.NZOLD)THEN
        WRITE(6,*)'*** ERROR: DIFFERENT ELEMENTS ON TWO FILES, NZ='
     X            ,NZOLD,NZ0
        STOP '*** ERROR: DIFFERENT ELEMENTS ON TWO FILES'
      ENDIF
      NZOLD=NZ0
C
      IF(NEOLD.NE.0.AND.NE.NE.NEOLD)THEN
        if(lv.ge.0)then
          WRITE(6,*)'*** ERROR: DIFFERENT IONS ON TWO FILES, NE='
     X              ,NEOLD,NE
          STOP '*** ERROR: DIFFERENT IONS ON TWO FILES'
        else
          if(ne.ne.neold-1.or.ne.ne.neold+1)then
            WRITE(6,*)'*** WARNING: DIFFERENT CHARGES ON TWO FILES, NE='
     X                ,NEOLD,NE
          endif
        endif
      ENDIF
      NEOLD=NE
C
      IF(bare.AND.NEOLD.NE.1.and.neold.ne.3.and.neold.ne.11
     x                      .and.neold.ne.29)THEN
        WRITE(6,*)'NO TARGET DATA/FILE ON o_str - REQUIRED FOR HYBRID'
        STOP 'ERROR: NO TARGET DATA/FILE ON - REQUIRED FOR HYBRID'
      ENDIF
C
      DO I=1,MXORB                 !SHORT ORBITAL LIST
        IF(QN(I).LE.0)GO TO 301
      ENDDO
      I=MXORB+1
 301  MXORB=I-1
C
      IF(NCF.GT.NDIM14)THEN
        WRITE(6,136)NCF
 136  FORMAT(' DIMENSION EXCEEDED IN SR.CROSSJ, INCREASE NDIM14 TO',I5)
        STOP 'ERROR: DIMENSION EXCEEDED IN SR.CROSSJ, INCREASE NDIM14'
      ENDIF
C
C****************************
C READ PHOTOIONIZATION HEADER
C****************************
C
      IF(BFORMP)READ(MRP,201) NENGP0,NZ0P,NEP,EIONMNP
      IF(.NOT.BFORMP)READ(MRPU) NENGP0,NZ0P,NEP,EIONMNP
 201  FORMAT(I3,12X,I2,6X,I2,35X,F15.6)
C
      IF(NENGP0.EQ.0)GO TO 1000
      IF(NENGP0.GT.NDIM33)THEN
        WRITE(6,*)'INCREASE NDIM33 TO: ',NENGP0
        STOP 'ERROR: INCREASE NDIM33'
      ENDIF
C
      IF(NEPOLD.NE.0.AND.NEPOLD.NE.NENGP0)THEN
        WRITE(6,*)'*** ERROR: DIFFERENT PI ENERGIES ON TWO FILES, NEP='
     X            ,NEPOLD,NENGP0
        STOP '*** ERROR: DIFFERENT PI ENERGIES ON TWO FILES'
      ENDIF
      NEPOLD=NENGP0
C
      IF(NZ0.NE.NZ0P)THEN
        WRITE(6,*)'*** ERROR: DIFFERENT ELEMENTS ON O/OP FILES, NZ='
     X            ,NZ0P,NZ0
        STOP '*** ERROR: DIFFERENT ELEMENTS ON O/OP FILES'
      ENDIF
C
      IF(NE.NE.NEP)THEN
        WRITE(6,*)'*** ERROR: DIFFERENT IONS ON O/OP FILES, NE='
     X            ,NEOLD,NE
        STOP '*** ERROR: DIFFERENT IONS ON O/OP FILES'
      ENDIF
C
C******************************
C READ PHOTOIONIZATION ENERGIES
C******************************
C
      IF(BFORMP)THEN
        READ(MRP,202)(ENERGP0(I),I=1,NENGP0)
 202    FORMAT(5E15.5)
        READ(MRP,202)
        READ(MRP,202)
        READ(MRP,204)ICP,ITP
        BACKSPACE(MRP)
      ELSE
        READ(MRPU)(ENERGP0(I),I=1,NENGP0)
        READ(MRPU)
        READ(MRPU)
        READ(MRPU)ICP,ITP
        BACKSPACE(MRPU)
      ENDIF
C
      BDOWN=ITP.EQ.0                          !BUNDLED ALREADY DOWNWARDS
      IF(BDOWN)THEN
        IF(ICP.EQ.0)THEN
          IF(NBINRM.NE.0)IFLAGR=2
        ELSEIF(.NOT.BHYBRD)THEN
          IFLAGR=1
        ENDIF
        IF(IFLAGR.NE.0)GO TO 1000
      ENDIF
C
C******************************
C LOOK FOR A SUBSET OF ENERGIES
C******************************
C
      IF(BPASS1)THEN
        PMIN=ALF4*PMN0/DZ                             !SKIP
        NP00=MIN(IABS(NP0),NENGP0)
        IF(NENGP0.GT.2.AND.PFACT.GT.DZERO)THEN
          IF(BFORMP)THEN
            READ(MRP,204)ICP,ITP,IWP(1),JCP,JTP,L,PCS0(1),EI0,EC
            READ(MRP,202)(PCS0(I),I=1,NENGP0)
            NREC=(NENGP0-1)/5+2
            DO N=1,NREC
              BACKSPACE(MRP)
            ENDDO
          ELSE
            READ(MRPU)ICP,ITP,IWP(1),JCP,JTP,L,PCS0(1),EI0,EC
            READ(MRPU)(PCS0(I),I=1,NENGP0)
            BACKSPACE(MRPU)
            BACKSPACE(MRPU)
          ENDIF
          NENGP=1
          MENGP(NENGP)=1
          M=1
          DO N=2,NP00-1
            IF(1.2D0*ABS(PCS0(N)).LT.ABS(PCS0(M)))THEN
              NENGP=NENGP+1
              MENGP(NENGP)=N
              M=N
            ENDIF
          ENDDO
          NENGP=NENGP+1
          MENGP(NENGP)=NP00
          DO N=1,NENGP
            M=MENGP(N)
            ENERGP(N)=ENERGP0(M)
          ENDDO
        ELSE
          NENGP=NP00
          DO N=1,NENGP
            MENGP(N)=N
            ENERGP(N)=ENERGP0(N)
          ENDDO
        ENDIF
        NENGPT=NENGP
        NENGT=NENGP
        IF(NP0.LT.0)GO TO 304               !NO HIGH-E TOP-UP
C
C ADD-IN HIGHER ENERGIES HYDROGENICALLY
C
        IF(IMESH.LE.0)THEN                  !DEFAULT
          MXETAB=MXTABC
          DO N=1,MXETAB
            EDUM(N)=EDUMC(N)*DZ
          ENDDO
        ELSE                                !TEST FINER
          MXETAB=MXTABF
          DO N=1,MXETAB
            EDUM(N)=EDUMF(N)*DZ
          ENDDO
        ENDIF
        IF(MXETAB.GT.NDIM33)THEN
          WRITE(6,*)'*** INCREASE NDIM33 TO AT LEAST',MXETAB
          STOP 'ERROR: *** INCREASE NDIM33'
        ENDIF
        DO N=1,MXETAB
          IF(ENERGP(NENGP).LT.EDUM(N)*0.99D0)THEN
            NT=N
            NENGT=NENGP+MXETAB-N+1
            IF(NENGT+1.GT.NDIM33)THEN
              WRITE(6,*)'***INCREASE NDIM33 TO ',NENGT+1
              STOP 'ERROR: ***INCREASE NDIM33'
            ENDIF
            GO TO 303
          ENDIF
        ENDDO
  303   DO N=NENGP+1,NENGT
          ENERGP(N)=EDUM(NT)
          NT=NT+1
        ENDDO
        NENGPT=NENGT
        NENGT=NENGT+1            !ALLOWS REPLACEMENT BY COULOMB VIA PMIN
        ENERGP(NENGT)=ENERGP(NENGP)
C
  304   IF(ENERGP(NENGPT).LT.5*DZ/NR1**2)THEN
          WRITE(6,*)'*** MAX PI ENERGY TOO SMALL, INCREASE TO AT LEAST'
     X              ,5*DZ/NR1**2,' RYD'
          STOP 'ERROR: *** MAX PI ENERGY TOO SMALL ***'
        ENDIF
        IF(ENERGP(NENGPT).LT.10*DZ/NR1**2)THEN
          WRITE(6,*)'*** WARNING, CONSIDER INCREASING MAX PI ENERGY TO'
     X              ,10*DZ/NR1**2,' RYD'
        ENDIF
      ENDIF
C
C************************
C READ CONFIGURATION DATA
C************************
C
      DO 102 I=1,NCF
C
      IF(BFORM)READ(MR,179,END=1002)II,NGR,MA0,MB0,(QS0(J),QL0(J)
     X,J=1,10)
      IF(.NOT.BFORM)READ(MRU,END=1002)II,NGR,MA0,MB0,(QS0(J),QL0(J)
     X,J=1,10)
 179  FORMAT(2I5,2X,I3,I2,1X,10(I2,A1))
C
      IN=IABS(II)
C
C DECODE CONFIGURATIONS:
C   LMX(I) IS THE NO OF DISTINCT OPEN-SHELL ORBITALS IN CONFIG I.
C   QSB(I,J) IS THE OCCUPATION NO OF ORBITAL J IN CONFIG I.
C   QLB(I,J) IS THE ORBITAL NO OF ORBITAL J IN CONFIG I, J=1,LMX(I).
C   QS0,QL0 CONTAIN EISSNER SPECIFICATION OF CONFIG TO BE DECODED.
C
        DO 16 J=1,10
          QSB(I,J)=MBLNK
          IF(QL0(J).EQ.MBLNK)GO TO 16
          LMX(I)=J
          M=MOD(QS0(J),50)
          IF(M.GT.0)QSB(I,J)=LIT(M)
          DO K=1,NLIT
            IF(QL0(J).EQ.LIT(K))GO TO 19
          ENDDO
          QLB(I,J)=0
          GO TO 16
  19      QLB(I,J)=K
  16    CONTINUE
C
        J=LMX(I)
        IF(II.LT.0)THEN
          QSB(I,J)=MBLNK
          LMX(I)=J-1
          if(j.eq.1)then
            if(mb0.gt.0)ncmx0=max(ncmx0,qn(mb0))
            if(bhybrd)itarh(i)=1
            go to 102
          endif
        ENDIF
C
        IF(LV00.NE.LV0)THEN
          MM=QLB(I,LMX(I))
          IF(MM.GT.MXORB.OR.MM.EQ.0)THEN
            WRITE(6,*)'***ERROR, CF=',II,' USES ORBITAL NO=',MM
     X       ,' WHICH IS NOT DEFINED IN ORBITAL HEADER!!'
            STOP'***ERROR, NEED ORBITAL NOT DEFINED IN HEADER!!'
          ENDIF
C
          IF(BHYBRD)THEN                               !JUST MATCH BY CF
            J1=J
            IF(QSB(I,J1).NE.LIT(1).and.ii.gt.0)THEN
              J1M=J1
            ELSE
              J1M=J1-1
            ENDIF
            DO N0=1,NBINRM
              N=JKH(N0)                                !N=SYMMETRY ORDER
              IF(J1M.NE.LMH(N))GO TO 227
              DO J=1,J1-1
                IF(QSB(I,J).NE.QSH(N,J))GO TO 227
                IF(QLB(I,J).NE.QLH(N,J))GO TO 227
              ENDDO
              IF(QSB(I,J1).NE.LIT(1).and.ii.gt.0)THEN
                IF(QLB(I,J1).NE.QLH(N,J1M))GO TO 227
              ENDIF
              ITARH(I)=N0                               !N0=ENERGY ORDER
c          write(*,*)i,itarh(i)
              GO TO 226
 227          CONTINUE
            ENDDO
            ITARH(I)=9999
          ENDIF
C
C ATTEMPT TO DETERMINE MAX CORE N
 226      IF(QN(MM).NE.NV.OR.LV.LT.0)THEN   !SUITABLE
            NCMX0=MAX(NCMX0,QN(MM))
          ELSEif(lmx(i).gt.1)then
            M=MAX(1,LMX(I)-1)
            M=QLB(I,M)
            NCMX0=MAX(NCMX0,QN(M))   !ONLY APPROX
          ENDIF
        ENDIF
C
 102  CONTINUE
C
      IF(NCMX0.GT.IBN(IBX).AND.INIT.EQ.0)THEN
        WRITE(6,*)' INCREASE NRSLMX,NLMAX TO AT LEAST',NCMX0
        STOP 'ERROR: INCREASE NRSLMX,NLMAX'
      ENDIF
C
      IF(INR1.EQ.999)THEN       !DEFAULT
        NR1=MAX(NCMX0+1,NMN)
        INR1=NR1
        WRITE(6,*)' '
        WRITE(6,*)'*** NR1 RESET TO:',NR1
      ENDIF
      IF(NR1.EQ.0)THEN
c    BRAD=.FALSE. !COULD TEST IF LCUT SET TO SWITCH IT ON, STILL PARENTS
        WRITE(6,*)"SET NR1 .GT. DZERO TO GET PP'D HYDROGENIC DATA"
        STOP "ERROR: SET NR1 .GT. DZERO TO GET PP'D HYDROGENIC DATA"
      ENDIF
      IF(NR1.LT.0.AND.INIT.EQ.0)THEN
        WRITE(6,*)"*** WARNING!! ***"
        WRITE(6,*)"*** YOU HAVE SWITCHED-OFF PP'D HYDROGENIC DATA! ***"
        WRITE(*,*)"*** WARNING!! ***"
        WRITE(*,*)"*** YOU HAVE SWITCHED-OFF PP'D HYDROGENIC DATA! ***"
      ENDIF
C
C**************************
C SKIP AUTOIONIZATION RATES
C**************************
C
      IF(BFORM)READ(MR,103,END=1002)
      IF(.NOT.BFORM)READ(MRU,END=1002)NZTEST,NDUME
      IF(BFORM)READ(MR,103,END=1002)
 103  FORMAT(A1)
C
 111  IF(BFORM)READ(MR,112,END=1002)ICA,ITA
      IF(.NOT.BFORM)READ(MRU,END=1002)ICA,ITA
 112  FORMAT(5I5,5X,1PE15.5,2(0PF15.6))
C
      IF(ITA.GT.0) GO TO 111
C
      IF(BFORM)THEN
        BACKSPACE(MR)
        READ(MR,112,END=1002)ICA,ITA,IWA,JCA,JTA,AA,EC,EION
      ENDIF
      IF(.NOT.BFORM)THEN
        BACKSPACE(MRU)
        READ(MRU,END=1002)ICA,ITA,IWA,JCA,JTA,AA,EC,EION
      ENDIF
C
      EIONMN=EION          !NEED MIN(EIONMN,EION) IF EION NOT LOWEST...
C
C**************
C READ ENERGIES
C**************
C
      IF(BFORM)READ(MR,121,END=1002) NENG,ECORE
      IF(.NOT.BFORM)READ(MRU,END=1002) NENG,ECORE
 121  FORMAT(10X,I5,45X,F15.6)
C
      IF(NENG.EQ.0)THEN
        WRITE(6,*)'***ERROR: NO PI CROSS SECTIONS PRESENT FOR N .GE.',NV
        WRITE(6,*)'***APPLY EXPLICIT NCUT, OR RE-RUN AUTOSTRUCTURE'
        STOP '***ERROR: PI CROSS SECTIONS MISSING'
      ENDIF
C
      IF(BFORM)READ(MR,105,END=1002)MTEST4
      IF(.NOT.BFORM)READ(MRU,END=1002)MTEST4
 105  FORMAT(26X,A4)
C
      BTEST=MTEST4.NE.MBLNK4                               !IC=TRUE
      IF(BLS.AND.BTEST)THEN
        IF(BCA)THEN
          WRITE(6,371)
 371      FORMAT(' RUN INITIALIZED FOR CA BUT LS/IC DATA FOUND')
          STOP 'ERROR: RUN INITIALIZED FOR CA BUT LS/IC DATA FOUND'
        ELSE
          WRITE(6,370)
 370      FORMAT(' RUN INITIALIZED FOR LS BUT IC DATA FOUND')
          STOP 'ERROR: RUN INITIALIZED FOR LS BUT IC DATA FOUND'
        ENDIF
      ENDIF
      IF(BIC.AND..NOT.BTEST)THEN
        WRITE(6,374)
 374    FORMAT(' RUN INITIALIZED FOR IC BUT LS DATA FOUND')
        STOP 'ERROR: RUN INITIALIZED FOR IC BUT LS DATA FOUND'
      ENDIF
C
      IF(NENG.GT.NDIM13)THEN
        WRITE(6,369)NENG
 369  FORMAT(' NUMBER OF LEVELS EXCEEDS STORAGE,INCREASE NDIM13 TO',I6)
        STOP 'ERROR: NUMBER OF LEVELS EXCEEDS STORAGE,INCREASE NDIM13'
      ENDIF
C
      MFLAG=0
      KFLAG=NENG
C
      DO 122 I=1,NENG
C
      IF(BFORM)READ(MR,123,END=1002)IK(I),IT,SS(I),LL(I),JJ(I),LCF(I)
     X,ENERG(I)
      IF(.NOT.BFORM)READ(MRU,END=1002)IK(I),IT,SS(I),LL(I),JJ(I),LCF(
     XI),ENERG(I)
 123  FORMAT(5X,6I5,F15.6)
C
      M=IK(I)
      M=IABS(M)
      IF(M.LE.NDIM10)JK(M)=I
      MFLAG=MAX(MFLAG,M)
      K=IABS(LCF(I))
      IF(BCA.AND.BFORM)THEN
        SS(I)=100000*IT+SS(I)           !as write i10 read as 2i5
      ENDIF
C
      TE=ECORE+ENERG(I)
      IF(LCF(I).LT.0.AND.EIONMN.GE.DZERO)EIONMN=TE
C
C********************************
C SET-UP TARGET BINS AND INDEXING (ONLY DONE FOR A NEW UNIT).
C********************************
C
      IF(KFPM.GT.NBINP)GO TO 122
      IF(LCF(I).GT.0)GO TO 122
C
      IF(ENERG(I).GT.(WNP0+TOLB))THEN
C
        WNP0=ENERG(I)
        KFPM=KFPM+1
C
        IF(BBIN)THEN
          EII(KFPM)=TE
C
          IF(KFPM.LE.NBINM.AND.TOLR.GT.DZERO)THEN  !SET FINAL META RANGE
            T=EII(KFPM)-EII(1)
            IF(TOLR.LT.T)TOLR=T
          ENDIF
C
          IF(LV.GE.0)THEN           !CHECK TARGET ENERGIES
            T=TE-EII(1)
            T0=E1C(KFPM)
            IF(KFPM.EQ.1)THEN
              IF(IPRINT.GE.0)WRITE(6,372)TOLBE
 372          FORMAT(3X,'IE',10X,'E(N)',14X,'E(N+1)',2X,'TOLB='
     X              ,1PE10.3)
              T0=T0-E1C(1)
            ENDIF
            IF(ABS(T-T0).GT.TOLBE.AND.E1C(KFPM).NE.DZERO)THEN
              MMM=MSTAR
              IFLAGE=IFLAGE+1
            ELSE
              MMM=MBLNK
            ENDIF
            IF(IPRINT.GE.0)                            !.OR.MMM.NE.MBLNK
     X        WRITE(6,373)KFPM,T0,T,MMM
 373          FORMAT(I5,2F18.8,3X,A4)
          ENDIF
        ENDIF
C
        IF(KFPM.LE.NBINP)THEN
          WNP(KFPM)=-TE
          IF(LV.LT.0.AND.TEAPOT.GT.DZERO)WNP(KFPM)=WNP(KFPM)
     X                                -(WNP(1)+ENERG(1)+ECORE+TEAPOT)
          DO J=1,10
            QSP(KFPM,J)=QSB(K,J)
            QLP(KFPM,J)=QLB(K,J)
          ENDDO
          LMP(KFPM)=LMX(K)
        ELSEIF(BBIN.AND..NOT.BHYBRD)THEN
          IF(NR1.NE.0)THEN
            T=EII(KFPM)+QDT(QD0,NZ0,NE,INR1,0,0)
            IF(T.LT.EIONMN.AND.IK(I).GT.0)THEN   !ONLY SPEC AS RR 1-BODY
              WRITE(6,744) T,EIONMN,(EII(J),J=1,KFPM)
 744  FORMAT(' STRONG WARNING IN SR.CROSSJ, THERE MAYBE A FINAL-STATE'
     X,' WITH PARENT NOT SPECIFIED BY NTAR2'//' INCREASE NTAR2 AND/OR'
     X,' NDIM0'//'EBOUND=',F12.3,5X,'EION=',F12.3,' EBIN='/
     X (10F12.6))
c      STOP 'ERROR: INCREASE NTAR2 AND/OR NDIM0'
            ENDIF
          ENDIF
        ENDIF
C
      ELSEIF(KFPM.LE.NBINP)THEN                    !CHECK LABELLING
C
        IF(LMP(KFPM).NE.LMX(K))GO TO 381
        DO J=1,LMX(K)
          IF(QSP(KFPM,J).NE.QSB(K,J))GO TO 381
          IF(QLP(KFPM,J).NE.QLB(K,J))GO TO 381
        ENDDO
        GO TO 382                                  !WE ARE GOOD
c
C FLAG
 381    IF(NCFR.EQ.0)THEN
          WRITE(6,*)'LEVEL ',I,
     X      ' HAS AN INCONSISTENT TARGET CONFIG LABEL='
     X      ,K,', REDUCE NTAR2 TO',KFPM-1
          NFLAG2=MIN(KFPM-1,NFLAG2)
c        write(6,*)i,kfpm,k
c        write(6,383)((qsp(kfpm,j),qlp(kfpm,j)),j=1,lmp(kfpm)-1)
c        write(6,383)((qsb(k,j),qlb(k,j)),j=1,lmx(k)-1)
c 383    format(10(a2,i2,3x))
c        STOP ' INCONSISTENT TARGET CONFIG LABEL, REDUCE NTAR2'
        ELSE
          KFLAG=MIN(KFLAG,KFPM-1)
        ENDIF
C
      ENDIF
C
 382  CONTINUE
C
      WNP0=ENERG(I)           !ALLOW FOR ANY DRIFT OF CONTINUUM ENERGIES
C
 122  CONTINUE                  !END ENERGY LOOP READ
C
      IF(MFLAG.GT.NDIM10)THEN
        WRITE(6,368)MFLAG
 368  FORMAT(' NUMBER OF LEVELS EXCEEDS STORAGE,INCREASE NDIM10 TO',I6)
        STOP 'ERROR: NUMBER OF LEVELS EXCEEDS STORAGE,INCREASE NDIM10'
c      ELSE
c        MENG=MFLAG              !FOR CORRELATION LABELS
      ENDIF
C
C*****************************************
C END ENERGY READ AND INDEXING BY SYMMETRY
C*****************************************
C
C AS BUNDLED CASE MAY NOT HAVE ALL TARGETS, SO:
C
      KHOLD=KFPM
      IF(KFLAG.LT.NBINP)THEN
        IF(LCP(1).GT.0)THEN
          WRITE(6,*)
          WRITE(6,*)'COMMENT: USING TARGET CONFIG LABELS FROM FILE '
     X,'o_str'                    !or ou_str...
C          KFPM=KFLAG             !ASSUME TARGET CFS CORRECT TO HERE
          KFPM=0                  !MIGHT AS WELL REPLACE ALL...
        ELSE
          WRITE(6,*)'*** ERROR: NEED TARGET CF NOS IN adasin TO USE '
     X              ,'o_str'
          STOP '*** ERROR: NEED TARGET CF NOS IN adasin - see adasout'
        ENDIF
      ENDIF
C
      IF(KFPM.LT.NBINP)THEN
        DO KK=KFPM+1,NBINP
          K=LCP(KK)                                   !TARGET CF NO
          IF(K.GT.NCFR)GO TO 275
          IF(BHYBRD)THEN                              !NEED ENERGY
            TE=E1C(KK)+EIONMN
            IF(BBIN)EII(KK)=TE
            WNP(KK)=-TE
          ENDIF
          DO J=1,10
            QSP(KK,J)=QSH(K,J)
            QLP(KK,J)=QLH(K,J)
          ENDDO
          LMP(KK)=LMH(K)
        ENDDO
        IF(BHYBRD)THEN
          IF(KHOLD.LT.NBINM)GO TO 274                 !WE NEED INITIALS
          KFPM=NBINP
          IF(BBIN.AND.E1C(NBINP+1).GT.DZERO)THEN
            KFPM=KFPM+1
            EII(KFPM)=E1C(KFPM)+EIONMN
          ENDIF
        ELSE
          KFPM=KHOLD
        ENDIF
        GO TO 279
C
 274    WRITE(6,1378)NBINM,KFPM,TOLB
1378  FORMAT(//' *** ERROR: YOU HAVE REQUESTED NTAR1=',I2,' PARENTS BUT'
     X,' ONLY ',I2,' CAN BE DETERMINED FROM YOUR AUTOSTRUCTURE DATA.'/
     X/' *** POSSIBLE CAUSES:'//
     X'1.   THERE ARE NOT NTAR1 TARGET CONTINUA - '
     X,'CHECK YOUR AUTOSTRUCTURE TARGET.'/
     X'2.   CHECK AUTOSTRUCTURE NMETAR/J, TERMS/LEVELS *AND* EMXLS/IC.'/
     X'3.   LEVEL SPLITTING IS .LT. 1.D-6 BUT O1 IS BEING USED - '
     X,'SWITCHED TO UNFORMATTED DATA ON O1U.'/
     X'4.   TOLB IS TOO LARGE - TRY SETTING IT TO LESS THAN HALF OF '
     X,'THE SMALLEST LEVEL SPLITTING: TOLB=',1PE9.2)
        IF(LV.EQ.-1)WRITE(6,1379)
1379  FORMAT(//' *** MOST LIKELY CAUSE:'/
     X/'     SEPARATION BETWEEN TARGET NTAR1'
     X,' AND NTAR1+1 IS COMPARABLE WITH ENERGY DIFFERENCE'/
     X,'     BETWEEN USING N- AND (N+1)-ELECTRON ORBITAL POTENTIALS'/
     x/' *** WORKAROUND (IN AUTOSTRUCTURE):'//
     X'     EITHER SET EIXMLS/IC MANUALLY/EXPLICITLY IN SMINIM'/
     X'     OR INCREASE NMETAR/J TO REACH STATE WELL ABOVE LOWER ONE(S)'
     X)
      STOP'ERROR: NUMBER OF PARENTS REQUESTED EXCEEDS NO. FOUND ON FILE'
C
 275    WRITE(6,378)NBINP,KK-1,TOLB
 378    FORMAT(' ERROR: YOU HAVE REQUESTED NTAR2=',I2,' PARENTS BUT'
     X,' ONLY ',I2,' CAN BE DETERMINED FROM YOUR AUTOSTRUCTURE DATA.'/
     X/' POSSIBLE CAUSES:'/
     X'1.    THERE ARE NOT NTAR2 TARGET CONTINUA - '
     X,'CHECK YOUR AUTOSTRUCTURE TARGET.'/
     X,'2.   LEVEL SPLITTING IS .LT. 1.D-6 BUT O1 IS BEING USED - '
     X,'SWITCHED TO UNFORMATTED DATA ON O1U.'/
     X'3.   TOLB IS TOO LARGE - TRY SETTING IT TO LESS THAN HALF OF '
     X,'THE SMALLEST LEVEL SPLITTING: TOLB=',1PE9.2)
      STOP'ERROR: NUMBER OF PARENTS REQUESTED EXCEEDS NO. FOUND ON FILE'
      ENDIF
C
 279  IF(BBIN.AND.EI(1).GE.DZERO)THEN
        IF(KFPM.LE.NBINP)EII(KFPM+1)=0.8D0*EII(NBINP)+D1M18
        EI(1)=1.2D0*EII(1)
        DO M=1,NBINP
          EI(M+1)=EII(M+1)-0.5D0*(EII(M+1)-EII(M))
        ENDDO
      ENDIF
      KFPM=NENG
C
      IF(BPASS1)THEN
        WRITE(6,120)(ENERGP(M),M=1,NENGPT)
 120    FORMAT(//'  INCIDENT ELECTRON ENERGIES (RYD)'/(16X,1P8E15.5))
        IF(BPRNT2)WRITE(6,90)
  90    FORMAT(/4X,'N  L',5X,'EO ',8(4X,'E*CROSS(MB)'))
        IF(BPRNT0)WRITE(6,33)
  33    FORMAT(4X,'J',4X,'I',1X,8(4X,'E*CROSS(MB)'))
      ENDIF
C
      IF(BPRNT0.AND.NV.GE.0)WRITE(6,34)NV,LV
  34  FORMAT(I5,I3)
C
C********************************
C SET-UP RESOLVED DUMMY PARENTAGE
C********************************
C
      IF(BRSLF.AND.NR1.NE.0.AND.NRSOLZ.EQ.0)THEN
       if(nbinp.ne.nbinrm)then
         write(*,*)'nbinp,nbinrm mismatch?',nbinp,nbinrm
         stop 'nbinp,nbinrm mismatch?'
       endif
C
       DO KK=1,NBINP
C
C TERMS
C
        IF(BLS)THEN
        ISP=1
        IF(IWS(KK).GT.1)ISP=3
          L00=0
          LUPZ=LMAXZ
          IF(LCUT.GE.0)LUPZ=MIN(LCUT+1,LMAXZ)
          DO L=1,LUPZ                       !DUMMY
            IF(INIT.NE.0)L00=L
            L1=IABS(L-1-IWL(KK))
            L2=IABS(L-1+IWL(KK))
            DO LT=L1,L2
              DO IS=1,ISP,2
                NRSOLZ=NRSOLZ+1
                IF(NRSOLZ.LE.NDIM18)THEN
                  ITARZ(NRSOLZ)=KK
                  SSZ(NRSOLZ)=IWS(KK)+2-IS
                  LLZ(NRSOLZ)=LT
                  JJZ(NRSOLZ)=MAX(INIT,L00)
                  LMZ(NRSOLZ)=LMP(KK)+1
c        write(6,*)nrsolZ,kk,999,l-1,ssZ(nrsolZ),lt,lmZ(nrsolZ)
                  DO J=1,10
                    QSZ(NRSOLZ,J)=QSP(KK,J)
                    QLZ(NRSOLZ,J)=QLP(KK,J)
                  ENDDO
                  J=LMZ(NRSOLZ)
                  QSZ(NRSOLZ,J)=LIT(1)
                  QLZ(NRSOLZ,J)=L-1               !JUST L
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
C LEVELS
C
        IF(BIC)THEN
          L00=0
          LUPZ=LMAXZ
          IF(LCUT.GE.0)LUPZ=MIN(LCUT+1,LMAXZ)
          DO L=1,LUPZ                       !DUMMY
            IF(INIT.NE.0)L00=L
            JV1=IABS(2*L-3)
            JV2=IABS(2*L-1)
            DO JVT=JV1,JV2,2
              JMIN=IABS(JVT-IWJ(KK))
              JMAX=IABS(JVT+IWJ(KK))
              DO JT=JMIN,JMAX,2
                NRSOLZ=NRSOLZ+1
                IF(NRSOLZ.LE.NDIM18)THEN
                  ITARZ(NRSOLZ)=KK
                  SSZ(NRSOLZ)=MAX(INIT,L00)
                  LLZ(NRSOLZ)=0
                  JJZ(NRSOLZ)=JT
                  LMZ(NRSOLZ)=LMP(KK)+1
c        write(6,*)nrsolZ,kk,999,l-1,jvt,jjz(nrsolZ)
                  DO J=1,10
                    QSZ(NRSOLZ,J)=QSP(KK,J)
                    QLZ(NRSOLZ,J)=QLP(KK,J)
                  ENDDO
                  J=LMZ(NRSOLZ)
                  QSZ(NRSOLZ,J)=LIT(1)
                  QLZ(NRSOLZ,J)=L-1           !JUST L
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
C
        IF(NRSOLZ.GT.NDIM18)THEN
          WRITE(6,380)NRSOLZ
 380      FORMAT(' SR.CROSSJ: INCREASE NDIM18 TO:',I6)
          STOP 'ERROR: DIMENSION ERROR: INCREASE NDIM18'
        ENDIF
C
       ENDDO
C
      ENDIF    !END DUMMY
C
C******************
C DETERMINE PARENTS
C******************
C
      IFIRST=1
      DO 36 I=1,NENG
C
        IRSOL(I)=0
        KK=IABS(IK(I))
        ITAR(KK)=-1
        TE=ENERG(I)+ECORE
C
        IF(LCF(I).LT.0)THEN              !CONTINUUM
C
          DO M=IFIRST,NBINP
            IF(TE.GE.EI(M).AND.TE.LT.EI(M+1))THEN
              ITAR(KK)=M
              IF(M.GT.IFIRST)IFIRST=M
              GO TO 36
            ENDIF
          ENDDO
C
        ELSE                             !DISCRETE
C
         K=LCF(I)
         J1=LMX(K)
         M=QLB(K,J1)

c
         l=0
         if(lv.lt.0.or.qn(m).ne.nv)l=1                 !master
C
         IF(NBINRM.GT.0)THEN
C HYBRID
          IF(BHYBRD)THEN
C
            IF(ITARH(K).EQ.9999)THEN
              IF(TE.LT.EIONMN+TOLR.and.l.eq.0)THEN
                WRITE(6,775)I
                STOP 'ERROR: PARENT NOT FOUND: SEE adasout FOR CAUSES'
              ELSE
                if(l.ne.0)then                         !relax for master
                  bflagp=.true.
                  if(bprnt0)
     x               WRITE(6,*)' ***PARENT NOT DETERMINABLE FOR CF=',K
                else
                  WRITE(6,*)' ***PARENT NOT DETERMINABLE FOR CF=',K
                endif
              ENDIF
            ENDIF
C
C CA: SET-UP ALL PARENTS AS POTENTIAL FINAL RECOMBINED
C
          ELSEIF(BCA)THEN
            IF(QSB(K,J1).NE.LIT(1))THEN
              J1M=J1
            ELSE
              J1M=J1-1
            ENDIF
            DO NP=1,NBINRM
              IF(J1M.NE.LMP(NP))GO TO 228
              DO J=1,J1-1
                IF(QSB(K,J).NE.QSP(NP,J))GO TO 228
                IF(QLB(K,J).NE.QLP(NP,J))GO TO 228
              ENDDO
              IF(QSB(K,J1).NE.LIT(1))THEN
                IF(QLB(K,J1).NE.QLP(NP,J1M))GO TO 228
              ENDIF
              GO TO 230                                !WE HAVE A WINNER
  228         CONTINUE
            ENDDO
            IF(TE.LT.EIONMN+TOLR.and.l.eq.0)THEN
              WRITE(6,775)I
              STOP 'ERROR: PARENT NOT FOUND: SEE adasout FOR CAUSES'
            ELSE
              IF(NCFR.GT.0)                      !SHOULD BE ABLE TO FIND
     X           WRITE(6,*)' ***PARENT NOT DETERMINABLE FOR I=',I
              if(l.ne.0)then                           !relax for master
                bflagp=.true.
                if(bprnt0)
     x             WRITE(6,*)' ***PARENT NOT DETERMINABLE FOR CF=',K
              else
                WRITE(6,*)' ***PARENT NOT DETERMINABLE FOR CF=',K
              endif
            ENDIF
            NP=9999
  230       ITAR(KK)=NP
c      write(6,*)kk,np
          ENDIF
C
          IF(.NOT.BRSLF)GO TO 360
C
C USE DUMMY SET-UP TO DETERMINE PARENT
C
          IF(LV.GE.0.AND.NR1.NE.0.AND.QN(M).EQ.NV.OR.   !RYDBERG
     X       LV.LT.0.AND.NR1.LT.0.AND.QN(M).GE.INR1)THEN
            J1M=LMX(K)-1
            DO 128 N=1,NRSOLZ
C
              IF(BLS)THEN
                IF(INIT.EQ.0.AND.JJZ(N).NE.0)GO TO 128  !ALREADY MATCHED
                IF(INIT.NE.0.AND.JJZ(N).NE.QN(M))GO TO 128  !WRONG N
                IF(IABS(SS(I)).NE.SSZ(N))GO TO 128
                IF(LL(I).NE.LLZ(N))GO TO 128
              ENDIF
              IF(BIC)THEN
                IF(INIT.EQ.0.AND.SSZ(N).NE.0)GO TO 128  !ALREADY MATCHED
                IF(INIT.NE.0.AND.SSZ(N).NE.QN(M))GO TO 128  !WRONG N
                IF(JJ(I).NE.JJZ(N))GO TO 128
              ENDIF
C
              IF(LMX(K).NE.LMZ(N))GO TO 128
              IF(QL(QLB(K,J1)).NE.QLZ(N,J1))GO TO 128
              DO J=1,J1M
                IF(QSB(K,J).NE.QSZ(N,J))GO TO 128
                I1=QLB(K,J)
                I2=QLZ(N,J)
                IF(QL(I1).NE.QL(I2))GO TO 128
                IF(QN(I1).NE.QN(I2))GO TO 128
              ENDDO
C
              ITAR(KK)=ITARZ(N)
              IF(BLS)JJZ(N)=QN(M)+1        !TAG
              IF(BIC)SSZ(N)=QN(M)+1        !TAG
c              write(6,*)i,n,itarz(n),qn(m),ql(m)
              GO TO 36
 128        CONTINUE
C
            IF(TE.LT.EIONMN+TOLR)THEN
              IF(INIT.EQ.0)THEN
                WRITE(6,767)I
 767            FORMAT(' ERROR IN CROSSJ: THERE IS A RESOLVED FINAL',
     X        '  STATE I=',I5,', WITH A PARENT NOT SPECIFIED BY NTAR2.'
     X          //' POSSIBLE CAUSES:'
     X          //'   1/ DIMENSIONS, INCREASE NTAR2 AND/OR NDIM0.'
     X          //'   2/ CHECK ORDER OF TARGET SYMMETRIES IN adasin.')
                STOP 'ERROR: PARENT NOT FOUND: SEE adasout FOR CAUSES'
              ELSE
                WRITE(6,*)' ***USING CONTINUUM PARENT FOR I=',I
              ENDIF
            ENDIF
            ITAR(KK)=9999
C
          ELSE
C
C CASE LV.EQ.-1 (AS NR1.EQ.0 NOT ALLOWED NOW): BY ENERGY ALONE
C SO DO NOT REALLY NEED PARENT, BUT FLAG TO USE CONTINUUM PARENT FOR RR.
C
            ITAR(KK)=9999
c      write(6,*)-nv,lv,i,ITAR(KK)
          ENDIF
C
         ENDIF
C
C CASE LV.GE.0: NON-RYDBERG BOUND N+1, FLAG AS CORRELATION TO AVOID
C MULTIPLE COUNTING (USER HAS LEFT MXCCF.GT.0 - O.K. FOR DR, NOT RR.)
C
 360     IF(LV.GE.0.AND.QN(M).NE.NV)THEN
           IF(BHYBRD)ITARH(K)=-IABS(ITARH(K))
           IK(I)=-IABS(IK(I))                         !AS op MAY NOT BE
         ENDIF
C
        ENDIF
C
  36  ENDDO
C
      IF(.NOT.BRSLF)GO TO 100
C
C RESET (NOTE, DOES NOT WORK FOR CASE INIT.NE.0 - NEEDS RYD L)
C
      IF(BLS)THEN
        DO N=1,NRSOLZ
          JJZ(N)=0
        ENDDO
      ENDIF
      IF(BIC)THEN
        DO N=1,NRSOLZ
          SSZ(N)=0
        ENDDO
      ENDIF
C
C****************************************
C SET-UP INDEXING FOR NEW RESOLVED STATES
C****************************************
C
      DO 127 I=1,NENG
C
      IF(LCF(I).LT.0)GO TO 127
      IF(IK(I).LE.0)GO TO 127
      TE=ENERG(I)+ECORE
      IF(TE.GE.EIONMN+TOLR)GO TO 705
C
      KK=IABS(IK(I))
      K=LCF(I)
      J1=LMX(K)
      M=QLB(K,J1)
      IF(QN(M).LE.NRSLMX)THEN
        ITEST=0
        IF(QN(M).EQ.NV.AND.LV.GE.0)THEN
          ITEST=ITAR(KK)
        ELSE
          IF(ITAR(KK).GT.NBINM.AND.TE.GE.EIONMN)GO TO 127
          IF(LV.GE.0.AND.LCF(I).LE.JCFJ)THEN
            WRITE(6,*)'***DANGER: CONFIG ',K,' WILL BE MULTIPLY'
     X,' COUNTED. SET JCFJ TO EXCLUDE IT, OR REMOVE FROM DATASET'
            STOP 'ERROR: ***DANGER: CONFIG WILL BE MULTIPLY COUNTED'
          ENDIF
        ENDIF
C
        NRSOL=NRSOL+1
        IF(NRSOL.LE.NDIM17)THEN
c          write(6,*)nv,lv,i,itar(kk),nrsol
          QNV(NRSOL)=QN(M)
          IF(LV.LT.0)QNV(NRSOL)=-QNV(NRSOL)
          QLV(NRSOL)=QL(M)
          IRSOL(I)=NRSOL
          ITARR(NRSOL)=ITAR(KK)
          SSR(NRSOL)=-IABS(SS(I))    !TAG
          LLR(NRSOL)=LL(I)
          IF(BLS)JJR(NRSOL)=-1       !TAG
          IF(BIC)JJR(NRSOL)=JJ(I)
          WNR(NRSOL)=ENERG(I)+ECORE
          LMR(NRSOL)=LMX(K)
          DO J=1,10
            QSR(NRSOL,J)=QSB(K,J)
            QLR(NRSOL,J)=QLB(K,J)
          ENDDO
          IF(ITEST.GT.0)THEN
            QLR(NRSOL,LMR(NRSOL))=(NV*(NV-1))/2+QL(M)+1+MXORB0
            if(lv.ne.ql(m).and.lv.ne.999)stop 'lv .ne. ql'
          ENDIF
        ENDIF
      ENDIF
C
 127  CONTINUE
C
 705  IF(NRSOL.GT.NDIM17)THEN
        WRITE(6,379)NRSOL
        STOP'ERROR: TOO MANY RESOLVED RESONANCES,INCREASE NDIM17'
      ENDIF
C
C*********************
C SKIP RADIATIVE RATES
C*********************
C
 100  IF(BFORM)READ(MR,104,END=1002)NZTEST
      IF(.NOT.BFORM)READ(MRU,END=1002)NZTEST,NDUME
 104  FORMAT(66X,I2)
C
      IF(NZTEST.GE.1)THEN
C
        IF(BFORM)READ(MR,103,END=1002)
C
 131    IF(BFORM)READ(MR,132,END=1002) ICR,ITR
        IF(.NOT.BFORM)READ(MRU,END=1002)ICR,ITR
 132    FORMAT(6I5,1PE15.5,2(0PF15.6))
C
        IF(ITR.GT.0)GO TO 131
      ENDIF
C
C*************************************
C INITIALIZE FOR PHOTOIONIZATION XSCNS
C*************************************
C
      ITP0=-1
      IPP0=0
      L0=999999
      NRSOL0=0
      IF(BDOWN)THEN
        DE=DONE
        C4=DONE/D1M18                 !JUST CONVERT TO MB
      ELSE
        C4=ALF4
      ENDIF
C
      DO M=1,NENG
        IWP(M)=0
        DO N=1,NENGP
          PCS(N,M)=DZERO
        ENDDO
      ENDDO
C
      DO L=1,NBINM
        DO N=1,NENGP
          TC(N,L)=DZERO
          TCP(N,L)=DZERO
        ENDDO
      ENDDO
C
C***************************
C READ PHOTOIONIZATION XSCNS
C***************************
C
 133  IF(BFORMP)READ(MRP,204)ICP,ITP,IWP0,JCP,JTP,L,PCS0(1),EI0,EC
      IF(.NOT.BFORMP)READ(MRPU)ICP,ITP,IWP0,JCP,JTP,L,PCS0(1),EI0,EC
 204  FORMAT(6I5,E15.5,2F15.6)
C
      IF(JCP.EQ.0)THEN
        IPP=0
        GO TO 313
      ENDIF
C
      IF(NENGP0.GT.1)THEN
        IF(BFORMP)READ(MRP,202)(PCS0(I),I=1,NENGP0)
        IF(.NOT.BFORMP)READ(MRPU)(PCS0(I),I=1,NENGP0)
      ENDIF
      IFLGL2=MAX(IFLGL2,L)
C
      IF(BINT)GO TO 133                    !NV,LV NOT WANTED, SO SKIP
      IF(L.LE.0.OR.L.GT.NBINM)GO TO 133    !INITIAL STATE NOT METASTABLE
      IF(EI0.GT.EIONMN+TOLR)GO TO 133      !FINAL STATE NOT REQUIRED
C
      IF(.NOT.BDOWN)THEN
        J=JK(ITP)
        IF(IK(J).LE.0)GO TO 133            !CORRELATION (INTERNAL FLAG)
        IPP=ITP
      ELSE
        IF(ICP.GT.0)THEN
          IF(ITARH(ICP).LT.0)GO TO 133     !CORRELATION (INTERNAL FLAG)
          IPP=ICP
        ELSE
          IPP=1
        ENDIF
      ENDIF
C
      IF(.NOT.BFAST)THEN
        IF(ICP.GT.JCFJ)GO TO 133
        ITOLB0=ITAR(JTP)
        IF(L.NE.ITOLB0)THEN
           WRITE(88,*)JCP,JTP,ITOLB0,L,EC
           IF(ITOLB0.GT.0)L=ITOLB0
        ENDIF
        IF(EIONMNP-EI0.LT.EBDMIN.OR.EIONMNP-EI0.GT.EBDMAX)GO TO 133
      ENDIF
C
C**************************
C SUM PHOTOIONIZATION XSCNS
C**************************
C
      IF(ITP0.EQ.0)GO TO 313              !SINCE MUST STORE FIRST BUNDLE
      IF(L.GT.L0)GO TO 313
C
 203  IF(IWP(IPP).EQ.0)THEN       !NEW
        NRSOL0=NRSOL0+1
        IRSOL0(NRSOL0)=IPP
        IWP(IPP)=NRSOL0
      ENDIF
C
      N0=IWP(IPP)
      DO N=1,NENGP
        M=MENGP(N)
        PCS(N,N0)=PCS(N,N0)+ABS(PCS0(M))
      ENDDO
C
      L0=L
      IPP0=IPP
      IF(BDOWN)THEN
        JTP0=JTP
        ITP0=ITP
        EI00=EI0
      ENDIF
C
      GO TO 133
C
C**********************************
C STORE AS PARTIAL RR XSCNS * ECONT
C**********************************
C
 313  IF(IPP0.EQ.0)GO TO 414
      IPP0=IPP
      IFLGL1=MAX(IFLGL1,L0)
C
      DO N0=1,NRSOL0
C
        IPP=IRSOL0(N0)
C
        IF(.NOT.BDOWN)THEN                               !ITP=IPP HERE
          IJT=IPP
          KT=JK(IJT)
          ICP=LCF(KT)
          EI00=ENERG(KT)+ECORE                           !ABSOLUTE
          EC=-WNP(L0)
          DE0=EC-EI00
          BDROP=PCS(1,N0).LT.DELTAF
        ELSE                                             !NO CONT SUM
          IJT=JTP0
          KT=JK(IJT)
          ICP=IPP
        ENDIF
C
        MSS=IABS(SS(KT))                                 !FOR LS NSYS
        IF(BLS)TW=IABS(SS(KT))*(2*LL(KT)+1)              !ADJUST WEIGHT
        IF(BIC)TW=JJ(KT)+1                               !N+1 TO N
        TWW=C4*TW/IWT(L0)
        IF(BPRNTP2)WRITE(6,107)ICP,IJT,L0,(PCS(N,N0)/D1M18,N=1,NENGP)
C
        DO N=1,NENGP
          TCP(N,L0)=TCP(N,L0)+PCS(N,N0)/D1M18
          IF(.NOT.BDOWN)THEN
            DE=DE0+ENERGP(N)
            DE=DE*DE
          ENDIF
          PCS(N,N0)=TWW*DE*PCS(N,N0)                     !E*RRCS
          TC(N,L0)=TC(N,L0)+PCS(N,N0)
        ENDDO
        IF(BPRNT0)WRITE(6,107)ICP,IJT,L0,(PCS(N,N0),N=1,NENGP)
 107              FORMAT(3I5,8(2X,1PE13.4)/(15X,8(2X,1PE13.4)))
C
C RESOLVED
C
        MS=0
        I0=IRSOL(KT)
        IF(I0.GT.0)THEN                         !TERM/LEVEL RESOLVED
          MS=ITARR(I0)                          !ASSUME WANTED HERE
          IF(MS.LE.0)THEN
            WRITE(6,*)'***NO PARENT FOR RESOLVED CF, I:',ICP,IPP
            STOP 'ERROR: ***RESOLVED PARENT NOT FOUND'
          ENDIF
          IF(.NOT.BDROP)THEN
            SSR(I0)=IABS(SSR(I0))               !TAG RADIATION EXISTS
            DO N=1,NENGP                        !Only abs(ssr) in adasrr
              AN(N,I0,L0)=AN(N,I0,L0)+PCS(N,N0)
            ENDDO
          ENDIF
        ENDIF
C
C BUNDLED
C
c        if(ipp.ne.iabs(ik(kt)))stop 'ipp.ne.iabs(ik(kt))'   !test
c        if(ms.gt.0.and.ms.ne.itar(ipp))stop 'ms.ne.itar(ipp)'
        IF(BHYBRD)MS=ITARH(ICP)
        IF(MS.EQ.0)MS=ITAR(IPP)
        IF(MS.EQ.9999)MS=L0                      !NON-DETERMINABLE
C
        MN=QLB(ICP,LMX(ICP))
        ML=QL(MN)+1                              !=LV+1 NORMALLY
        MN=QN(MN)                                !=NV NORMALLY
C
        IF(MS*NBINRM.GT.NBINRM*NBINRM)THEN!SHOULD'VE ALREADY BEEN CAUGHT
          IF(MN.LT.INR1.OR.LV.LT.0.OR.NR1.EQ.0)THEN
            MS=L0                     !RAD INTO CORE, NO PARENTAGE
          ELSE
            WRITE(6,713)IPP
 713        FORMAT(' ERROR IN SR.CROSSJ, FINAL PARENT NOT FOUND, T=',I4
     X           ,'. INCREASE NTAR2/NDIM0')
            STOP 'ERROR IN SR.CROSSJ, FINAL PARENT NOT FOUND'
          ENDIF
        ENDIF
        IF(MS.LE.0)THEN
          WRITE(6,*)'***NO PARENT FOR BUNDLED CF, I:',ICP,IPP
          STOP 'ERROR: ***BUNDLED PARENT NOT FOUND'
        ENDIF
C
        NRX=MAX(NRX,MS)
        L00=(L0-1)*NBINRM+MS
C
        NSYS=1
        IF(NSYM.EQ.2.AND.MSS.GT.IWS(L0).AND.IWS(L0).GT.1   !LS NOT CA
     X      .AND.(.NOT.BRSLF.OR.IWS(MS).EQ.IWS(L0)))NSYS=2
        IF(MN.LE.NLMAX.AND.ML.LE.LLMAX+1)THEN              !BUNDLE-NL
          DO N=1,NENGP
            BNL(N,ML,IB,NSYS,L00)=BNL(N,ML,IB,NSYS,L00)+PCS(N,N0)
          ENDDO
        ENDIF
        IF(NBINRM.GT.0)THEN
          DO N=1,NENGP
            BN(N,IB,NSYS,L00)=BN(N,IB,NSYS,L00)+PCS(N,N0)  !BUNDLE-N
          ENDDO
        ENDIF
        IF(EI00.LT.EIONMN)THEN  !TRUE BOUND
          DO N=1,NENGP
            ALFN(N,IB,L0)=ALFN(N,IB,L0)+PCS(N,N0)          !FOR TOTAL
          ENDDO
        ENDIF
C
      ENDDO
C
      IF(IPP0.EQ.0)GO TO 414
C
      DO N=1,NENG               !RE-INITIALIZE
        IWP(N)=0
      ENDDO
      DO N0=1,NRSOL0
        DO N=1,NENGP
          PCS(N,N0)=DZERO
        ENDDO
      ENDDO
      NRSOL0=0
      IPP=IPP0
      GO TO 203
C
 414  CONTINUE
C
C**********************************************
C END PROCESSING PHOTOIONIZATION CROSS SECTIONS
C**********************************************
C
      IF(NV.GE.0)THEN
        LP=1
        IF(IPRINT.GE.0)LP=NBINM
c        LP=MAX(LP,IABS(IPRINT))
c        LP=MIN(LP,NBINM)
        IF(BPRNTP1)THEN
          DO L=1,LP
            IF(TCP(1,L).GT.DZERO)WRITE(6,35)NV,LV,L,(TCP(N,L),N=1,NENGP)
          ENDDO
        ENDIF
        IF(BPRNT1)THEN
          DO L=1,LP
            IF(TC(1,L).GT.DZERO)WRITE(6,35)NV,LV,L,(TC(N,L),N=1,NENGP)
  35        FORMAT(I5,I3,2X,I5,8(2X,1PE13.4)/(15X,8(2X,1PE13.4)))
          ENDDO
        ENDIF
      ENDIF
c      IF(TC(1).EQ.DZERO)NV00=0 !L=1
C
C  GO AND READ NEW NL BLOCK
C
      GO TO 310
C
C  ABORT
 1002 NV=0
      WRITE(6,1107)
 1107 FORMAT(/' ******WARNING, UNEXPECTED END OF DATA IN SR.CROSSJ !!!!'
     X,/)
C      GO TO 1001
C
 1000 CONTINUE
C
      IF(IFLAGR.EQ.2)THEN
        WRITE(6,*)' '
        WRITE(6,*)
     X          'AS DATA CONTAINS UNRESOLVED PHOTORECOMBINATION RATES'
        WRITE(6,*)'NO ADF48 POSSIBLE, TOTALS (NTAR2=0) ONLY'

      ELSEIF(IFLAGR.EQ.1)THEN
        WRITE(6,*)' '
        WRITE(6,*)
     X'AS DATA CONTAINS CONFIGURATION RESOLVED PHOTORECOMBINATION RATES'
        WRITE(6,*)
     X     'ADF48 FILE REQUIRES NTAR2.LT.0 AND THE TARGET o_str FILE'

      ENDIF
      IF(IFLAGR.NE.0)THEN
        WRITE(6,*)'***ERROR: RESOLUTION CONFLICT:'
        STOP '***ERROR: RESOLUTION CONFLICT - SEE adasout'
      ENDIF
C
C
C**********************
C  GO AND READ NEW FILE
C**********************
C
      CLOSE(MR)                !MRU=MR
c
cpar      if(ifile.lt.0)go to 1001                    !no more      !par
C
      IFILE=IFILE+1
      IC1=IFILE/10
      IC2=IFILE-10*IC1
      IC0=ICHAR('0')
      IC1=IC1+IC0
      IC2=IC2+IC0
      INAM=CHAR(IC1)//CHAR(IC2)
C
 500  IF(BFORM)THEN
        IF(O1.EQ.'o1')THEN                                    !EX-SERIAL
          IF(IFILE.LT.10)THEN
            FILNAM='o'//CHAR(IC2)
            FILNAMP='op'//CHAR(IC2)
          ELSE
            FILNAM='o'//INAM
            FILNAMP='op'//INAM
          ENDIF
        ELSE                                                !EX-PARALLEL
          IF(BTWO)THEN
            FILNAM=O//JNAM//'.'//INAM
            FILNAMP=OP//JNAM//'.'//INAM
          ELSE
            FILNAM=O//INAM
            FILNAMP=OP//INAM
          ENDIF
        ENDIF
        INQUIRE(FILE=FILNAM,EXIST=BEX)
        INQUIRE(FILE=FILNAMP,EXIST=BEXP)
        IF(BEX.NEQV.BEXP)THEN
          WRITE(6,1007)
          STOP '*** ERROR:opn FILE EXISTS WITH NO CORRESPONDING on FILE'
        ENDIF
        IF(BEX)THEN
          OPEN(MR,FILE=FILNAM)
          OPEN(MRP,FILE=FILNAMP)
        ENDIF
      ELSE
        IF(O1U.EQ.'o1u')THEN                                  !EX-SERIAL
          IF(IFILE.LT.10)THEN
            FILNAM='o'//CHAR(IC2)//'u'
            FILNAMP='op'//CHAR(IC2)//'u'
          ELSE
            FILNAM='o'//INAM//'u'
            FILNAMP='op'//INAM//'u'
          ENDIF
        ELSE                                                !EX-PARALLEL
          IF(BTWO)THEN
            FILNAM=O//'u'//JNAM//'.'//INAM
            FILNAMP=OP//'u'//JNAM//'.'//INAM
          ELSE
            FILNAM=O//'u'//INAM
            FILNAMP=OP//'u'//INAM
          ENDIF
        ENDIF
        INQUIRE(FILE=FILNAM,EXIST=BEX)
        INQUIRE(FILE=FILNAMP,EXIST=BEXP)
        IF(BEX.NEQV.BEXP)THEN
          WRITE(6,1007)
          STOP '*** ERROR:opn FILE EXISTS WITH NO CORRESPONDING on FILE'
        ENDIF
        IF(BEX)THEN
          OPEN(MRU,FILE=FILNAM,FORM='UNFORMATTED')
          OPEN(MRPU,FILE=FILNAMP,FORM='UNFORMATTED')
        ENDIF
      ENDIF
      IF(BEX)GO TO 331
C
      IF(BTWO.AND.IFILE.GT.0)THEN         !RESET IFILE & INCREMENT JFILE
        IFILE=0
ccpar        ifile=ifile+iam                                        !par
        IC1=IFILE/10
        IC2=IFILE-10*IC1
        IC0=ICHAR('0')
        IC1=IC1+IC0
        IC2=IC2+IC0
        INAM=CHAR(IC1)//CHAR(IC2)                     !SINCE EX-PARALLEL
        JFILE=JFILE+1
        JC1=JFILE/10
        JC2=JFILE-10*JC1
        JC0=ICHAR('0')
        JC1=JC1+JC0
        JC2=JC2+JC0
        JNAM=CHAR(JC1)//CHAR(JC2)
        GO TO 500
      ENDIF
C
 1001 CONTINUE                !HYDROGENIC RE-ENTRY POINT <<<<<<<<<<<<<<<
C
C**************************************************
C EVALUATE HYDROGENIC RADIATIVE RECOMBINATION RATES
C**************************************************
C (CAN'T USE WITH LV=999 UNLESS SET MAX L OF EXPLICIT PI)
C
      IBT=IB0      !=MAX N-MESH NEEDED=READ(=IB00) UNLESS IB00.LT.NRSLMX
      LV0=LV00+1
      IF(.NOT.BRAD)GO TO 160
      IF(LV0.EQ.999)GO TO 160
      IF(LV0.GE.LCUT.AND.NENGP.GE.NENGPT)GO TO 160
      IF(LCUT.GE.NDIM15)THEN
        WRITE(6,*)'***INCREASE NDIM15 TO AT LEAST',LCUT+1
        STOP 'ERROR: ***INCREASE NDIM15'
      ENDIF
C
      NENGP0=NENGP
      NENGP=NENGPT
C
      DO N=1,NTOP0
        IF(NTOP(N).GT.IBN(IB0).AND.NTOP(N).LE.NCUT)THEN
          IBT=IBT+1
          IBN(IBT)=NTOP(N)
        ENDIF
      ENDDO
      IF(IBT.GT.NDIM27)THEN
        WRITE(6,*)'*** DIMENSION ERROR, INCREASE NDIM27 TO',IBT
        STOP 'ERROR: ***INCREASE DIMENSION NDIM27'
      ENDIF
C
      IF(IREL.GT.0.AND.LLMAX0.EQ.-9999)LLMAX=MXLL   !NO H-LIKE BUNDLE-NL
C
      NRMXN=MAX(NRSLMX,NLMAX)                 !CASE NLMAX=0
      NRMXL=MAX(NRSLMX-1,LLMAX)               !    DITTO
C
      DO I=1,IBT
        DO L=1,NRMXL+1
          NLSWCH(L,I)=NENGP0
        ENDDO
        DO M=1,NENGT
          RNE(M,I)=DZERO
        ENDDO
        DO L=1,NRMXL+1
          DO M=1,NENGT
            RNLE(M,L,I)=DZERO
          ENDDO
        ENDDO
      ENDDO
      DO I=IBT,1,-1
        IF(INR1.GE.IBN(I))THEN
          IBT0=I
c         write(0,*)i,ibn(i),inr1
          GO TO 119
        ENDIF
      ENDDO
C
 119  DCZ=DZ                                  !DCZ=DZ0 TEST
      LVMIN=LV0+1
      LVMN=LVMIN
      LL0=LVMN
      IF(NENGT.GT.NENGP)LL0=MAX(0,LMN)
      NMIN=MAX(NR1,LL0+1)
      N15=MAX(15,NRMXN)
      LV0=LCUT
      IBL0=0
C
      DO I=IBT0,IBT                           !NR1 CASE NMIN.NE.I

        N=IBN(I)
        IF(N.GT.NCUT)GO TO 159
        IF(N.GE.NMIN)THEN
C
          IF(N.LE.NLMAX)IBL0=I
          LVMAX=MIN(LCUT,N-1)
          LPMAX=LVMAX+1
C
          IF(I.GT.IB00)THEN
            LVMN=MAX(0,LMN)
            LL0=MIN(LL0,LVMN)
          ENDIF
C
          TN=N*N
          DO M=1,NENGP
            DE=DZ/TN+ENERGP(M)
            E2=ENERGP(M)/DCZ
C
            CALL DIPOL(1,N,0,E2,LPMAX,CP,CM,JDUM)
C
            IF(NENGT.GT.NENGP)THEN                     !TOP-UP IN ENERGY
              IF(N.LE.N15)THEN
                LSWCH=MIN(LVMAX,NRMXL)
                IF(N.GT.NRMXN)LSWCH=0              !BUNDLE-N SWITCH ONLY
                DO L=0,LSWCH
                  LP=L+1
                  IF(M.LE.NLSWCH(LP,I))THEN
                    TL=L+L
                    TLP=LP+LP
                    T=TLP*CP(LP)*1.0D10**JDUM(LP)
                    IF(L.GT.0)T=T+TL*CM(L)*1.0D10**JDUM(L)
                    T=1.13953D-5*T*DE/(DCZ*DCZ)                    !SKIP
                    IF(T.LT.PMIN)THEN                   !/DE**(2*IPIG-2)
                      NLSWCH(LP,I)=M-1
                    ELSE
CE                                          !T=1.13953D-5*T*DE/(DCZ*DCZ)
                      IF(N.LE.NRMXN)RNLE(NENGT,LP,I)=T*DE**2
                    ENDIF
                  ENDIF
                ENDDO
              ELSE
                NLSWCH(1,I)=NLSWCH(1,I-1)
              ENDIF
            ENDIF
C
            IF(LL0.LE.LVMAX)THEN                    !TOP-UP IN L
              DO L=LL0,LVMAX
                LP=L+1
                TL=L+L
                TLP=LP+LP
                T=TLP*CP(LP)*1.0D10**JDUM(LP)
                IF(L.GT.0)T=T+TL*CM(L)*1.0D10**JDUM(L)
                T=1.13953D-5*T*DE**3/(DCZ*DCZ)
                IF(N.LE.NRMXN.AND.L.LE.NRMXL)THEN
                  MENGP0=NLSWCH(LP,I)
                  IF(M.GT.MENGP0.OR.L.GE.LVMN)RNLE(M,LP,I)=T
                ENDIF
                IF(M.GT.NLSWCH(1,I).OR.L.GE.LVMN)RNE(M,I)=RNE(M,I)+T
                IF(M.LE.NLSWCH(1,I).AND.NENGT.GT.NENGP)THEN
                  IF(L.EQ.LL0)RNE(NENGT,I)=DZERO
                  RNE(NENGT,I)=RNE(NENGT,I)+T
                ENDIF
              ENDDO
            ENDIF
C
            TC(M,1)=RNE(M,I)
          ENDDO                                     !END ENERGY LOOP
c
c         do l=0,lswch
c           write(6,*)'switch:  ',n,l,nlswch(l+1,i)+1
c         enddo
C
          IF(BPRNT1)WRITE(6,129)N,(TC(M,1),M=1,NENGP)
 129      FORMAT(I5,10X,8(2X,1PE13.4)/(15X,8(2X,1PE13.4)))
        ENDIF
C
      ENDDO                                         !END N-LOOP
C
 159  CONTINUE
C
C**************************************************
C TOP-UP IN A.M. AND ENERGY USING HYDROGENIC RATES
C**************************************************
C
      TS=DONE
C
C BUNDLE-N
C
      IF(IBN(IB0).LT.999)THEN
        IB0=IBT
        IF(.not.bare)THEN
          WRITE(6,*)
     X'***STRONG WARNING: HYDROGENIC N TOP-UP NOT RECOMMENDED FOR L=O-2'
          WRITE(0,*)
     X'***STRONG WARNING: HYDROGENIC N TOP-UP NOT RECOMMENDED FOR L=O-2'
        ENDIF
      ENDIF
C
CTOP-UP IN L
      DO K=1,NBINM                                  !METASTABLES ONLY
        NSYSM=1
        IF(NSYM.EQ.2.AND.IWS(K).GT.1)NSYSM=2        !LS NOT CA
        TOTS=2*IWS(K)
        K0=K
        IF(BHYBRD)K0=ITARH(LCP(K))
        L0=(K-1)*NBINRM+K0                          !NO CHANGE IN PARENT
        DO I=IBT0,IBT                               !NR1 CASE NMIN.NE.I
          N=IBN(I)
          TN=N*N
          EI0=-WNP(K)-DZ/TN
          IF(EI0.GT.EIONMN+TOLR)GO TO 151
          MENGP0=NLSWCH(1,I)
          IF(I.LE.IB0.AND.NBINRM.GT.0)THEN
            DO NS=1,NSYSM
              IF(NSYM.EQ.2)THEN                 !LS NOT CA
                IS=(IWS(K)+2*NS-3)
                IF(IS.EQ.0)IS=2                 !TEDIOUS
                TS=IS/TOTS                      !ADJUST DETAILED BALANCE
              ENDIF
              DO M=1,MENGP0
                BN(M,I,NS,L0)=BN(M,I,NS,L0)+RNE(M,I)*TS
              ENDDO
            ENDDO
          ENDIF
          IF(EI0.LT.EIONMN)THEN                   !TRUE BOUND
            DO M=1,MENGP0
              ALFN(M,I,K)=ALFN(M,I,K)+RNE(M,I)    !FOR TOTAL
            ENDDO
          ELSE
            IF(I.GT.IB0)GO TO 151
          ENDIF
        ENDDO
 151  ENDDO
C
CTOP-UP IN E
      IF(NENGT.GT.NENGP)THEN
        L0=0
        DO K=1,NBINM                            !METASTABLES ONLY
          NSYSM=1
          IF(NSYM.EQ.2.AND.IWS(K).GT.1)NSYSM=2  !LS NOT CA
          DO J=1,NBINRM                         !ALLOW CHANGE IN PARENT
            L0=L0+1
            DO I=IBT0,IB0                       !NR1 CASE NMIN.NE.I
              N=IBN(I)
              TN=N*N
              EI0=-WNP(J)-DZ/TN
              IF(EI0.GT.EIONMN+TOLR)GO TO 152
              MENGP0=NLSWCH(1,I)
              M1=MENGP0+1
              DO NS=1,NSYSM
                IF(BN(MENGP0,I,NS,L0).GT.DZERO)THEN
                  TT=BN(MENGP0,I,NS,L0)/RNE(NENGT,I)
                  DO M=M1,NENGP
                    BN(M,I,NS,L0)=RNE(M,I)*TT
                  ENDDO
c                  if(k.eq.j)write(78,*)k,ns,ibn(i),tt
                ENDIF
              ENDDO
            ENDDO
 152      ENDDO
          DO I=IBT0,IBT                         !FOR TOTALS
            N=IBN(I)
            TN=N*N
            EI0=-WNP(K)-DZ/TN
            IF(EI0.LT.EIONMN)THEN               !TRUE BOUND
              MENGP0=NLSWCH(1,I)
              M1=MENGP0+1
              IF(ALFN(MENGP0,I,K).GT.DZERO)THEN
                TT=ALFN(MENGP0,I,K)/RNE(NENGT,I)
                DO M=M1,NENGP
                  ALFN(M,I,K)=RNE(M,I)*TT
                ENDDO
c                write(79,*)k,n,tt
              ENDIF
            ELSE
              GO TO 153
            ENDIF
          ENDDO
 153    ENDDO
      ENDIF
C
C BUNDLE-NL
C
      LVMN=LVMIN
      IF(IBL0.GT.IB00)LVMN=MAX(0,LMN)
C
CTOP-UP IN L
      IF(NLMAX.GT.0.AND.LVMN.LE.LLMAX)THEN
        DO K=1,NBINM                            !METASTABLES ONLY
          NSYSM=1
          IF(NSYM.EQ.2.AND.IWS(K).GT.1)NSYSM=2        !LS NOT CA
          TOTS=2*IWS(K)
          K0=K
          IF(BHYBRD)K0=ITARH(LCP(K))
          L0=(K-1)*NBINRM+K0                    !NO CHANGE IN PARENT
          LVMN=LVMIN
          DO I=IBT0,IBL0                        !NR1 CASE NMIN.NE.I
            N=IBN(I)
            TN=N*N
            EI0=-WNP(K)-DZ/TN
            IF(EI0.GT.EIONMN+TOLR)GO TO 154
            LVMAX=MIN(LLMAX,N-1)
            IF(I.GT.IB00)LVMN=MAX(0,LMN)
            DO NS=1,NSYSM
              IF(NSYM.EQ.2)THEN                 !LS NOT CA
                IS=(IWS(K)+2*NS-3)
                IF(IS.EQ.0)IS=2                 !TEDIOUS
                TS=IS/TOTS                      !ADJUST DETAILED BALANCE
              ENDIF
              DO L=LVMN+1,LVMAX+1
                MENGP0=NLSWCH(L,I)
                DO M=1,MENGP0
                  BNL(M,L,I,NS,L0)=BNL(M,L,I,NS,L0)+RNLE(M,L,I)*TS
                ENDDO
              ENDDO
            ENDDO
          ENDDO
 154    ENDDO
C
CTOP-UP IN E
        IF(NENGT.GT.NENGP)THEN
          L0=0
          NSYSM=1
          IF(NSYM.EQ.2.AND.IWS(K).GT.1)NSYSM=2   !LS NOT CA
          DO K=1,NBINM                           !METASTABLES ONLY
            DO J=1,NBINRM                        !ALLOW CHANGE IN PARENT
              L0=L0+1
              DO I=IBT0,IBL0                     !NR1 CASE NMIN.NE.I
                N=IBN(I)
                TN=N*N
                EI0=-WNP(J)-DZ/TN
                IF(EI0.GT.EIONMN+TOLR)GO TO 155
                DO NS=1,NSYSM
                  LVMAX=MIN(LLMAX,N-1)
                  DO L=1,LVMAX+1
                    MENGP0=NLSWCH(L,I)
                    IF(BNL(MENGP0,L,I,NS,L0).GT.DZERO)THEN
                      TT=BNL(MENGP0,L,I,NS,L0)/RNLE(NENGT,L,I)
                      M1=MENGP0+1
                      DO M=M1,NENGP
                        BNL(M,L,I,NS,L0)=RNLE(M,L,I)*TT
                      ENDDO
c                      if(j.eq.k)write(77,*)k,n,l-1,tt
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
 155        ENDDO
          ENDDO
        ENDIF
      ENDIF
C
      IF(.NOT.BRSLF)GO TO 160
C
C RESOLVED (TOP-UP IN ENERGY ONLY)
C
      IF(NENGT.GT.NENGP)THEN
        DO N0=1,NRSOL
          L0=ITARR(N0)
          J=LMR(N0)
          II=QLR(N0,J)
          N=QN(II)
          IF(IBN(N).NE.N)STOP 'ERROR: MIS-MATCH ON N-INDEXING'
          L=QL(II)
          LP=L+1
          MENGP0=NLSWCH(LP,N)
          M1=MENGP0+1
          IF(L0.LE.NBINM)THEN
            IF(BLS)THEN
              TW=SSR(N0)*(2*LLR(N0)+1)
              TOTW=IWS(L0)*(2*IWL(L0)+1)*2*(2*L+1)
            ENDIF
            IF(BIC)THEN
              TW=JJR(N0)+1
              TOTW=(IWJ(L0)+1)*2*(2*L+1)
            ENDIF
            TW=TW/TOTW                           !TW FOR INFO ONLY
          ELSE
            TW=DONE
          ENDIF
          IF(RNLE(NENGT,LP,N).GT.DZERO)THEN
            DO K=1,NBINM                         !ALLOW CHANGE IN PARENT
              IF(AN(MENGP0,N0,K).GT.DZERO)THEN
                TT=AN(MENGP0,N0,K)/(RNLE(NENGT,LP,N)*TW)
                DO M=M1,NENGP
                  AN(M,N0,K)=RNLE(M,LP,N)*TW*TT
                ENDDO
c                if(k.eq.l0)write(76,*)n0,n,lp-1,tt
              ENDIF
            ENDDO
          ELSE
            DO K=1,NBINM                         !ALLOW CHANGE IN PARENT
              IF(AN(MENGP0,N0,K).GT.DZERO)THEN
                DO M=M1,NENGP
                  TT=ENERGP(MENGP0)/ENERGP(M)
                  TT=SQRT(TT)*TT**LP
                  AN(M,N0,K)=AN(MENGP0,N0,K)*TT
                ENDDO
c                if(k.eq.l0)write(76,*)n0,n,lp-1,tt
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF
C
C********************************************************
C SET-UP RESOLVED TERMS/LEVELS FOR HYDROGENIC PP'D STATES
C********************************************************
C
      NRSOL0=NRSOL+1
      LVMN=LVMIN
      IF(NRSLMX.GT.IB00)LVMN=MAX(0,LMN)
      LOW=LVMN+1
      IF(LOW.GT.NRSLMX)GO TO 160
      NLOW=MAX(INR1,LVMIN+1)
      NUP=NRSLMX
      IF(NCUT.GT.0)NUP=MIN(NCUT,NRSLMX)
C
C TERMS
C
      IF(BLS)THEN
       DO K=1,NBINM                      !INTIAL METASTABLES
        ISP=1
        IF(IWS(K).GT.1)ISP=3
        LVMN=LVMIN
        DO N=NLOW,NUP
          IN=(N*(N-1))/2 + MXORB0
          LUP=N
          IF(LCUT.GE.0)LUP=MIN(LCUT+1,N)
          IF(N.GT.IB00)LVMN=MAX(0,LMN)
          LOW=LVMN+1
          DO L=LOW,LUP
            E=-WNP(K)+QDT(QD0,NZ0,NE,N,L-1,0)
            IF(E.GT.EIONMN+TOLR)GO TO 471     !AUTOIONIZING
            IO=IN+L
            L1=IABS(L-1-IWL(K))
            L2=IABS(L-1+IWL(K))
            DO LT=L1,L2
              DO IS=1,ISP,2
                NRSOL=NRSOL+1
                IF(NRSOL.LE.NDIM17)THEN
                  ITARR(NRSOL)=K
                  QNV(NRSOL)=0
                  SSR(NRSOL)=IWS(K)+2-IS
                  LLR(NRSOL)=LT
                  JJR(NRSOL)=-1          !TAG
                  WNR(NRSOL)=E
                  LMR(NRSOL)=LMP(K)+1
                  DO J=1,10
                    QSR(NRSOL,J)=QSP(K,J)
                    QLR(NRSOL,J)=QLP(K,J)
                  ENDDO
                  J=LMR(NRSOL)
                  QSR(NRSOL,J)=LIT(1)
                  QLR(NRSOL,J)=IO
                ENDIF
              ENDDO
            ENDDO
          ENDDO
 471      IF(L.EQ.LOW)GO TO 470
        ENDDO
 470    CONTINUE
       ENDDO
      ENDIF
C
C LEVELS
C
      IF(BIC)THEN
       DO K=1,NBINM                      !INTIAL METASTABLES
        LVMN=LVMIN
        DO N=NLOW,NUP
          IN=(N*(N-1))/2 + MXORB0
          LUP=N
          IF(LCUT.GE.0)LUP=MIN(LCUT+1,N)
          IF(N.GT.IB00)LVMN=MAX(0,LMN)
          LOW=LVMN+1
          DO L=LOW,LUP
            E=-WNP(K)+QDT(QD0,NZ0,NE,N,L-1,0)
            IF(E.GT.EIONMN+TOLR)GO TO 473     !AUTOIONIZING
            IO=IN+L
            JV1=IABS(2*L-3)
            JV2=IABS(2*L-1)
            DO JVT=JV1,JV2,2
              JMIN=IABS(JVT-IWJ(K))
              JMAX=IABS(JVT+IWJ(K))
              DO JT=JMIN,JMAX,2
                NRSOL=NRSOL+1
                IF(NRSOL.LE.NDIM17)THEN
                  ITARR(NRSOL)=K
                  QNV(NRSOL)=0
                  SSR(NRSOL)=-1          !TAG
                  LLR(NRSOL)=-1          !FLAG
                  JJR(NRSOL)=JT
                  WNR(NRSOL)=E
                  LMR(NRSOL)=LMP(K)+1
                  DO J=1,10
                    QSR(NRSOL,J)=QSP(K,J)
                    QLR(NRSOL,J)=QLP(K,J)
                  ENDDO
                  J=LMR(NRSOL)
                  QSR(NRSOL,J)=LIT(1)
                  QLR(NRSOL,J)=IO
                ENDIF
              ENDDO
            ENDDO
          ENDDO
 473      IF(L.EQ.LOW)GO TO 472
        ENDDO
 472    CONTINUE
       ENDDO
      ENDIF
C
      IF(NRSOL.GT.NDIM17)THEN
        WRITE(6,379)NRSOL
 379    FORMAT(' SR.CROSSJ: INCREASE NDIM17 TO:',I6)
        STOP 'ERROR: DIMENSION ERROR: INCREASE NDIM17'
      ENDIF
C
C RESOLVED
C
      DO N0=NRSOL0,NRSOL
        L0=ITARR(N0)
        J=LMR(N0)
        II=QLR(N0,J)
        N=QN(II)
        L=QL(II)
        LP=L+1
        IF(BLS)THEN
          TW=SSR(N0)*(2*LLR(N0)+1)
          TOTW=IWS(L0)*(2*IWL(L0)+1)*2*(2*L+1)
        ENDIF
        IF(BIC)THEN
          TW=JJR(N0)+1
          TOTW=(IWJ(L0)+1)*2*(2*L+1)
        ENDIF
        TW=TW/TOTW                         !ADJUST DETAILED BALANCE
        DO M=1,NENGP
          AN(M,N0,L0)=RNLE(M,LP,N)*TW       !TT=1. HERE OF COURSE
        ENDDO
      ENDDO
C
C******
C FLAGS
C******
C
 160  IF(TOLB.GT.TOLB0)WRITE(6,137)TOLB
 137  FORMAT(/' *** ATTN: TOLB HAS BEEN RESET TO =',1PE10.2,' RYD'/)
C
      IF(IFLAGE.NE.0.and.iprint.ge.0)WRITE(6,1006)IFLAGE
 1006 FORMAT(//'NOTE: ',I4,' UNIT5 TARGET ENERGIES DID NOT MATCH WITH'
     X,' THOSE PRESENT IN THE RATE FILE, SEE ABOVE "***" !'//)
c
      if(bflagp)then
        if(bprnt0)WRITE(0,*)'*** MASTER PARENT(S) NOT DETERMINABLE...'
        WRITE(6,*)'*** MASTER PARENT(S) NOT DETERMINABLE...'
      endif
C
      IF(NFLAG2.LT.NBINM)THEN                     !PARENT PROBLEMS
        WRITE(6,1008)NFLAG2
 1008   FORMAT('*** PRIOR PROBLEMS DECODING TARGET, TRY AND REDUCE'
     X,' NTAR2 TO LAST SAFE PARENT',I5)
       STOP'*** PARENT PROBLEMS - SEE adasout FOR DETAILS, REDUCE NTAR2'
      ENDIF
C
      IF(IFLGL1.NE.NBINM)THEN
        WRITE(0,140)NBINM,IFLGL1
        WRITE(6,140)NBINM,IFLGL1
 140    FORMAT(/' *** NOTE, NTAR1=',I5,' BUT LAST PI IS TO LEVEL',I5//)
      ENDIF
      IF(.NOT.BHYBRD.AND.IFLGL2.NE.NBINRM)THEN
        WRITE(0,141)NBINRM,IFLGL2
        WRITE(6,141)NBINRM,IFLGL2
 141    FORMAT(/' *** NOTE, NTAR2=',I5,' BUT LAST PI IS TO LEVEL',I5//)
      ENDIF
C
C
C***********************************
C WRITE RR DATA TO ADAS FORMAT ADF48
C *****            ****        *****
C***********************************
C
      IF(JTHETA.GE.0)THEN
        OPEN(10,FILE='adf48'//nam0)           !A SINGLE TYPE 3 OR 5 FILE
      ELSE
        JTEMP=0                                    !SET FOR TYPE-5 FIRST
        OPEN(10,FILE='adf48_5'//nam0)             !FIRSTLY A TYPE 5 FILE
      ENDIF
C
      IREL=ABS(IREL)                     !AS ONLY .GT.0 APPPLIES JUTTNER
C BUT DON'T MULTIPLY CROSS SECTION BY JUTTNER
      IF(JTHETA.LE.0.AND.IREL.GT.0)IREL=-1    !.LT. WILL TURN IT BACK ON
C
C******************
C INITIAL TEMPS ETC
C******************
C
 2000 CONTINUE                   !RE-ENTRY POINT IF WRITING TYPE 5 AND 3
C
      IF(JTEMP.EQ.0)JTMAX=NENGPT
      IF(JTEMP.GT.0)JTMAX=JTEMP
C
      IF(JTEMP.GT.0)THEN
        IF(TEMP(1).LE.DZERO)THEN
          DO J=1,JTEMP
            TEMP(J)=DZ*THTIC(J)
            TEMP(J)=TEMP(J)/CONRYK
          ENDDO
        ENDIF
        IF(ENERGP(1).GT.TEMP(1)/100)THEN
          WRITE(6,*)'*** RE-SET MIN PI ENERGY TO DZERO ***'
          STOP 'ERROR: *** RE-SET MIN PI ENERGY TO DZERO ***'
        ENDIF
        IF(LMAX.EQ.-77)THEN                !DEFAULT
          T=TEMP(1)*CONRYK
          IF(T.LT.9*DZ)LCUT=500
          IF(T.LT.0.9*DZ)LCUT=900
          IF(T.LT.0.09*DZ)THEN
            WRITE(6,*)'*** MIN TEMP TOO LOW TO CONVERGE RR TOTAL'
            STOP 'ERROR: *** MIN TEMP TOO LOW TO CONVERGE RR TOTAL'
          ENDIF
        ENDIF
        IF(NMAX.EQ.-66)THEN                !DEFAULT
          T=TEMP(1)*CONRYK
          IF(T.LT.9*DZ)NCUT=9999
          IF(T.LT.0.9*DZ)NCUT=45000
          IF(T.LT.0.09*DZ)THEN
            WRITE(6,*)'*** MIN TEMP TOO LOW TO CONVERGE RR TOTAL'
            STOP 'ERROR: *** MIN TEMP TOO LOW TO CONVERGE RR TOTAL'
          ENDIF
        ENDIF
      ENDIF
C
C APPLY RELATIVISTIC (JUTTNER) CORRECTION TO DISTRIBUTION
C
      IF(IREL.GT.0)THEN
C
        NU=2
        MU=4*NU*NU
C
        DO J=1,JTEMP
C
          THETA=TEMP(J)*DALF/2                            !KT(a.u.)/C**2
          T8=THETA/8
C START-OFF WITH ABRAMOWITZ & STEGUN 9.7.2, AVOIDS SMALL THETA OVERFLOW
          KSUM=-5/LOG10(THETA)
          KSUM=KSUM+1
          KSUM=MIN(KSUM,10)
          T=DONE
          SUM=DONE
          DO K=1,KSUM
            T=T*(MU-(2*K-1)**2)*T8/K
            SUM=SUM+T
c            write(*,*)k,t,sum,done/sum
          ENDDO
C IF NOT CONVERGED, SAFE TO EVALUATE EXPLICITLY NOW AS THETA LARGE
          IF(ABS(T/SUM).GT.D1M4)THEN
c            write(*,*)temp(j),theta,sum,t,t/sum
            TT=DONE/THETA
            FX=SQRT(PI*THETA/2)/(BESSK(2,TT)*EXP(TT))
            FREL(J)=FX
          ELSE
            FREL(J)=DONE/SUM
          ENDIF
c          write(*,*)j,ksum,theta,t/sum,done/sum,fx
C
        ENDDO
C
      ELSE
C
        DO J=1,JTMAX
          FREL(J)=DONE
        ENDDO
C
      ENDIF
C
      IF(NBINRM.EQ.0)GO TO 2001               !SKIP ADF48 PARTIAL WRITES
C
C TARGET/PARENT INFO
C
      NTEST=1+NENGP/2
      IF(BLS)LAB4='/LS/'
      IF(BIC)LAB4='/IC/'
      IF(BCA)LAB4='/CA/'
      IF(NE.LE.NDIM20)THEN
        LAB2=LSQ(NE)
      ELSE
        LAB2='**'
      ENDIF
      IF(JTEMP.GT.0)THEN
        ITYPE=3
      ELSE
        ITYPE=5
      ENDIF
      WRITE(10,21)LAB2,NZ0,LAB4,ITYPE
      WNP0=WNP(1)*DKCM
      IF(.NOT.BHYBRD)THEN
        IF(BCA)THEN
          WRITE(10,253)WNP0,NBINM,NRX
        ELSE
          IF(NRX.LT.100)THEN
            IF(BLS)WRITE(10,212)WNP0,NBINM,NRX
            IF(BIC)WRITE(10,2332)WNP0,NBINM,NRX
          ELSE
            IF(BLS)WRITE(10,213)WNP0,NBINM,NRX
            IF(BIC)WRITE(10,2333)WNP0,NBINM,NRX
          ENDIF
          IF(BLS)WRITE(10,23)
          IF(BIC)WRITE(10,233)
        ENDIF
      ELSE
        IF(BCA)THEN
          WRITE(10,254)WNP0,NBINM
        ELSE
          IF(BLS)WRITE(10,223)WNP0,NBINM
          IF(BIC)WRITE(10,224)WNP0,NBINM
        ENDIF
        NRX0=NRX
        NRX=NBINM
      ENDIF
C
C PARENT INDEXING
C
      IWF=26
      DO M=1,NRX
        WNPM=-WNP(M)*DKCM+WNP0
        IF(BLS)THEN
          TW=IWS(M)*(2*IWL(M)+1)
          TW=(TW-DONE)/DTWO
        ENDIF
        IF(BIC)THEN
          TW=IWJ(M)
          TW=TW/DTWO
        ENDIF
        DO J=1,10
          QS0(J)=MBLNK
          QL0(J)=MBLNK
          IF(J.LE.LMP(M))THEN
            K=QLP(M,J)
            QS0(J)=LIT(QN(K))
            L=MIN((QL(K)+1),NLAB)
            QL0(J)=LABL(L)
          ENDIF
        ENDDO
        J1=MAX(5,LMP(M))
        J0=J1-4
        IF(M.GT.NBINM)THEN
          MP=LABL(NLAB)
        ELSE
          MP=MBLNK
        ENDIF
        IF(BCA)THEN
          WRITE(10,27)M,(QS0(J),QL0(J),QSP(M,J),J=J0,J1),TW,WNPM,MP
        ELSE
          WRITE(10,26)M,(QS0(J),QL0(J),QSP(M,J),J=J0,J1),IWS(M),IWL(M)
     X               ,TW,WNPM,MP
        ENDIF
      ENDDO
C
      IF(BHYBRD)THEN
        M=JKH(1)
        if(bare)then
          wnh0=wnp0
          wnh(1)=wnh0/dkcm
        else
          WNH0=WNH(M)*DKCM
        endif
        NRX=NRX0
        WRITE(10,254)WNH0,NRX
        DO M0=1,NRX
          M=JKH(M0)
          WNHM=-WNH(M)*DKCM+WNH0
          TW=JJH(M)-1
          TW=TW/DTWO
          DO J=1,10
            QS0(J)=MBLNK
            QL0(J)=MBLNK
            IF(J.LE.LMH(M))THEN
              K=QLH(M,J)
              QS0(J)=LIT(QN(K))
              L=MIN((QL(K)+1),NLAB)
              QL0(J)=LABL(L)
            ENDIF
          ENDDO
          J1=MAX(5,LMH(M))
          J0=J1-4
C          IF(M0.GT.NBINM)THEN
C            MP=LABL(NLAB)
C          ELSE
            MP=MBLNK
C          ENDIF
          WRITE(10,27)M0,(QS0(J),QL0(J),QSH(M,J),J=J0,J1),TW,WNHM,MP
        ENDDO
      ENDIF
C
C RESOLVED INDEXING
C
      IF(NRSOL.GT.0)THEN
C
      if(.not.bskp)then                                 !for double pass
      NSKP=0
      DO M=1,NRSOL
        JVR(M)=M
        IF(BLS.AND.JJR(M).EQ.0)JVR(M)=0 !UNMATCHED SO ASSUME UNWANTED
        IF(BIC.AND.SSR(M).EQ.0)JVR(M)=0 !UNMATCHED SO ASSUME UNWANTED
        IF(BHYBRD)THEN
          IF(BLS)THEN
            WNR(M)=WNR(M)/SSR(M)
            JJR(M)=-JJR(M)
          ENDIF
          IF(BIC)THEN
            WNR(M)=WNR(M)/JJR(M)
            JJR(M)=JJR(M)-1                                  !BACK TO 2J
          ENDIF
        ENDIF
        IF(ITARR(M).LE.0)THEN
          JVR(M)=0
        ELSE
          DO L=1,NBINM
            IF(AN(NTEST,M,L).GT.DZERO)GO TO 44
          ENDDO
          JVR(M)=0
        ENDIF
  44    CONTINUE
c        IF(WNR(M).GE.EIONMN)THEN                           !FINAL CATCH
c          IF(ITARR(M).GT.NBINM0.and.ITARR(M).LT.9999)JVR(M)=0
c        ENDIF
        IF(JVR(M).EQ.0)THEN
          NSKP=NSKP+1
          WNR(M)=DZERO
        ENDIF
      ENDDO
c      write(0,*)'nskp=',nskp
C
      CALL HPSRTI(NRSOL,WNR,JVR(1))
C
      NRSOL=NRSOL-NSKP
c
c      do m0=1,nrsol+nskp
c        m=jvr(m0)
c        write(6,*)m,itarr(m),(qsr(m,i),qlr(m,i),i=1,lmr(m)),wnr(m)
c      enddo
      endif
C
      IF(NRSOL.GT.0)THEN
        M=JVR(1)
        WNR0=WNR(M)
      ELSE
        WNR0=DZERO
      ENDIF
      E00=WNR0
      WNR0=-WNR0*DKCM
      IF(BLS)LAB2='LS'
      IF(BCA)LAB2='CA'
      IF(BIC)LAB2='IC'
      IF(BHYBRD)THEN
        IF(BLS)WRITE(10,251)LAB2,WNR0,NRSOL
        IF(BIC)WRITE(10,241)LAB2,WNR0,NRSOL
      ELSE
        IF(BCA)THEN
          WRITE(10,250)LAB2,WNR0,NRSOL
        ELSE
          IF(BLS)WRITE(10,24)LAB2,WNR0,NRSOL
          IF(BIC)THEN
            IF(NRSOL.LT.10000)THEN
              WRITE(10,2401)LAB2,WNR0,NRSOL
            ELSE
              WRITE(10,2402)LAB2,WNR0,NRSOL
            ENDIF
            WRITE(10,240)
          ENDIF
        ENDIF
      ENDIF
C
      DO M0=1,NRSOL
        M=JVR(M0)
        IF(WNR(M).GE.EIONMN)THEN
          MP=LABL(NLAB)
        ELSE
          MP=MBLNK
        ENDIF
        WNRM=(WNR(M)-E00)*DKCM
        ISSR=IABS(SSR(M))
        IF(BLS)THEN
          TW=ISSR*(2*LLR(M)+1)
          TW=(TW-DONE)/DTWO
        ENDIF
        IF(BIC)THEN
          TW=IABS(JJR(M))
          TW=TW/DTWO
        ENDIF
        DO 48 J=1,10
          QS0(J)=MBLNK
          QL0(J)=MBLNK
C          IF(J-LMR(M))46,47,48              !REPLACEMENT UGLY, BUT SAFE
          IF(J.GT.LMR(M))GO TO 48
          IF(J.LT.LMR(M))GO TO 46
  47      CONTINUE
          IF(QNV(M).EQ.0)GO TO 46
          K=IABS(QNV(M))
          QS0(J)=LIT(K)
          L=MIN((QLV(M)+1),NLAB)
          QL0(J)=LABL(L)
          GO TO 48
  46      K=QLR(M,J)
          IF(QN(K).GT.NRSLMX)THEN
            WRITE(6,*)'***ERROR: ORBITAL CONFUSION BETWEEN onu FILES:'
     X                ,' M0,M,N,L='
            WRITE(6,*)M0,M,QN(K),QL(K)
            STOP'ERROR: ***ERROR: ORBITAL CONFUSION BETWEEN onu FILES'
          ENDIF
          QS0(J)=LIT(QN(K))
          L=MIN((QL(K)+1),NLAB)
          QL0(J)=LABL(L)
  48    CONTINUE
        J1=MAX(5,LMR(M))
        J0=J1-4
        IF(BRSLF)THEN                  !IF(.NOT.BHYBRD.AND..NOT.BCA)
          WRITE(10,29)M0,ITARR(M),(QS0(J),QL0(J),QSR(M,J),J=J0,J1)
     X         ,ISSR,CLIT(LLR(M)),TW,WNRM,MP
        ELSEIF(BCA.OR.BIC)THEN
          WRITE(10,30)M0,ITARR(M),(QS0(J),QL0(J),QSR(M,J),J=J0,J1)
     X                           ,TW,WNRM,MP
        ELSE                          !LS
          WRITE(10,31)M0,ITARR(M),(QS0(J),QL0(J),QSR(M,J),J=J0,J1)
     X              ,CLIT(JJR(M)),TW,WNRM,MP
        ENDIF
      ENDDO
C
      ENDIF
C
C BUNDLED-NL INDEXING
C
      IF(NLMAX.GT.0)THEN
        MXLL=LLMAX             !NOT MIN(MXLL,LLMAX) SINCE WE HAVE TOP-UP
c        mxll=9
        NLREP=0
        DO I=1,IB0
          N=IBN(I)
          IF(N.GT.NLMAX)GO TO 43
          NLREP=NLREP+MIN(N,MXLL+1)
        ENDDO
  43    IBL0=I-1
        WRITE(10,14)NLREP
        IL=0
        M1=1
        DO I=1,IBL0
          N=IBN(I)
          L2=MIN(N,MXLL+1)
          DO L1=1,L2
            IL=IL+1
            WRITE(10,17)IL,N,L1-1
          ENDDO
        ENDDO
      ENDIF
C
C BUNDLED-N INDEXING
C
      IF(NLMAX.GT.0)THEN
        LAB5='INREP'
      ELSE
        LAB5='IREP '
      ENDIF
      WRITE(10,13)IB0,LAB5
C
      DO I=1,IB0
        WRITE(10,15)I,IBN(I)
      ENDDO
C
      WRITE(10,1)
      IF(JTEMP.NE.0)THEN
        INDX='INDX TE='
      ELSE
        INDX='INDX EE='
      ENDIF
C
      JTW10=MIN(JTMAX,10)
      JTW20=MIN(JTMAX,20)
C
C PARTIAL RR DATA
C
      L0=0
      DO L=1,NBINM             !LOOP OVER INITIAL PARENTS
C
        WRITE(10,1)
        LP=IWL(L)+1
        IF(BCA)THEN
          TW=(IWS(L)-DONE)/DTWO
          WRITE(10,53)L,TW
        ELSEIF(BLS)THEN
          WRITE(10,9)L,IWS(L),LABL(LP),IWS(L)
        ELSEIF(BIC)THEN
          TW=IWJ(L)
          TW=TW/DTWO
          WRITE(10,52)L,IWS(L),LABL(LP),TW
        ENDIF
C
        IF(JTEMP.GT.0)THEN
          WRITE(10,12)INDX,(TEMP(J)*CONRYK,J=1,JTW20)
          IF(JTMAX.GT.JTW20)WRITE(10,6)(TEMP(J)*CONRYK,J=JTW20+1,JTMAX)
        ELSE
          WRITE(10,612)INDX,(ENERGP(J),J=1,JTW20)
          IF(JTMAX.GT.JTW20)WRITE(10,5)(ENERGP(J),J=JTW20+1,JTMAX)
        ENDIF
C
C RESOLVED DATA
C
        IF(NRSOL.GT.0)THEN
          DO M0=1,NRSOL
            M=JVR(M0)
            IF(AN(NTEST,M,L).GT.DZERO)THEN
              II=QLR(M,LMR(M))
              LX=QL(II)
              NX=QN(II)
              DO J=1,NENGP
                PCS0(J)=AN(J,M,L)                      !AS AN() REAL*4
              ENDDO
              CALL CONVOL(JTEMP,TEMP,ICON,TKAPPA,XDRY,NENG37,E37,F37
     X                   ,ENERGP,PCS0,NENGP,LX,NX,NZ)
              WRITE(10,3)M0,(FREL(J)*PCS0(J),J=1,JTW10)
              IF(JTMAX.GT.JTW10)
     X        WRITE(10,6)(FREL(J)*PCS0(J),J=JTW10+1,JTMAX)
            ENDIF
          ENDDO
          WRITE(10,1)
        ENDIF
C
C BUNDLED DATA
C
        DO 708 M=1,NBINRM      !LOOP OVER FINAL PARENTS
C
          L0=L0+1
          IF(M.GT.NRX)GO TO 708
C
          DO K=1,NSYM
            DO I=1,IB0
              IF(BN(NTEST,I,K,L0).GT.DZERO)GO TO 704
            ENDDO
          ENDDO
          GO TO 708
C
  704     NSYSM=1
          IF(NSYM.EQ.2)THEN                                   !LS NOT CA
            IF(BRSLF)THEN
              IF(ABS(IWS(L)-IWS(M)).GT.2)GO TO 708
              IF(IWS(L).GT.1.AND.IWS(L).EQ.IWS(M))NSYSM=2
              LM=IWL(M)+1
              IF(NBINRM.LT.100)THEN
                WRITE(10,716)M,IWS(M),LABL(LM),IWS(M),NSYSM
              ELSE
                WRITE(10,726)M,IWS(M),LABL(LM),IWS(M),NSYSM
              ENDIF
            ELSE
              M0=JKH(M)
              TW=JJH(M0)-1
              TW=TW/DTWO
              IF(IWS(L).GT.1)NSYSM=2
              WRITE(10,727)M,TW,NSYSM
            ENDIF
          ENDIF
          IF(BCA)THEN
            IF(BHYBRD)THEN
              M0=JKH(M)
              TW=JJH(M0)-1
            ELSE
              TW=IWS(M)-DONE
            ENDIF
            TW=TW/DTWO
            WRITE(10,18)M,TW
          ELSEIF(BIC)THEN
            IF(BRSLF)THEN
              LM=IWL(M)+1
              TW=IWJ(M)
              TW=TW/DTWO
              IF(NBINRM.LT.100)THEN
                WRITE(10,2)M,IWS(M),LABL(LM),TW
              ELSE
                WRITE(10,20)M,IWS(M),LABL(LM),TW
              ENDIF
            ELSE
              M0=JKH(M)
              TW=JJH(M0)-1
              TW=TW/DTWO
              WRITE(10,18)M,TW
            ENDIF
          ENDIF
C
C LOOP-OVER 2 SPIN-SYSTEMS IF LS, ELSE 1.
C
          DO K=1,NSYSM
            IF(NSYM.EQ.2)THEN                                 !LS NOT CA
              ISPTL=IWS(L)+2*K-3
              ISPTM=IWS(M)+2*K-3
              ISPT=MAX(ISPTL,ISPTM)
              IF(ISPT.EQ.0)ISPT=2
              WRITE(10,10)K,ISPT
            ENDIF
C
C BUNDLED-NL
C
            IF(NLMAX.GT.0)THEN
              WRITE(10,7)
              IL=0
              DO I=1,IBL0
                N=IBN(I)
                L2=MIN(N,MXLL+1)
                DO L1=1,L2
                  IL=IL+1
                  IF(BNL(NTEST,L1,I,K,L0).GT.DZERO)THEN
                    DO J=1,NENGP
                      PCS0(J)=BNL(J,L1,I,K,L0)         !AS BNL() REAL*4
                    ENDDO
                    CALL CONVOL(JTEMP,TEMP,ICON,TKAPPA,XDRY,NENG37
     X                         ,E37,F37,ENERGP,PCS0,NENGP,L1-1,N,NZ)
                    WRITE(10,3)IL,(FREL(J)*PCS0(J),J=1,JTW10)
                    IF(JTMAX.GT.JTW10)
     X              WRITE(10,6)(FREL(J)*PCS0(J),J=JTW10+1,JTMAX)
                  ENDIF
                ENDDO
              ENDDO
              WRITE(10,1)
            ENDIF
C
C BUNDLED-N
C
            WRITE(10,11)LAB5
            DO I=1,IB0
              IF(BN(NTEST,I,K,L0).GT.DZERO)THEN
                DO J=1,NENGP
                  PCS0(J)=BN(J,I,K,L0)                 !AS BN() REAL*4
                ENDDO
                CALL CONVOL(JTEMP,TEMP,ICON,TKAPPA,XDRY,NENG37,E37,F37
     X                     ,ENERGP,PCS0,NENGP,0,IBN(I),NZ)
                WRITE(10,3)I,(FREL(J)*PCS0(J),J=1,JTW10)
                IF(JTMAX.GT.JTW10)
     X          WRITE(10,6)(FREL(J)*PCS0(J),J=JTW10+1,JTMAX)
              ENDIF
            ENDDO
            WRITE(10,1)
          ENDDO
C
 708    CONTINUE         !END LOOP OVER FINAL PARENTS
C
      ENDDO              !END LOOP OVER INITIAL PARENTS
C
C END OF (PARTIAL) WRITES FOR ADF48
C *********************************
C
C NOW SUM TOTAL FROM GROUND+METASTABLES AND WRITE TO UNIT6
C AND END OF ADF48 (NOT REQUIRED BY ADAS BUT FOR CONVENIENCE OF OTHERS).
C
 2001 CONTINUE                           !RE-ENRTY POINT FOR NO PARTIALS
C
      if(bskp)go to 2002                               !no need to re-do
c
c      write(6,*)'bsp0,inag,ncp,nlag,imesh=',bsp0,inag,ncp,nlag,imesh
C
C FIRST SUM SEQUENTIAL N
C
      N1=1
      DO I=1,IBT
        IF(IBN(I).GT.N1)GO TO 710
        N1=N1+1
        DO L=1,NBINM
          DO J=1,NENGP
            ALF(J,L)=ALF(J,L)+ALFN(J,I,L)
          ENDDO
        ENDDO
      ENDDO
C
  710 CONTINUE
c
c      write(6,*)'sum seq to',n1
c      write(6,*)(alf(j,1),j=1,nengp)
C
C INTERPOLATE TO NVINT (=100, DEFAULT) OR NCUT, IF SPECIFIED.
C
      I=I+1
      i11=mod(ibt-i,2)
      IS=I-i11
      N1=IBN(IS-2)+1+i11
      DO I=IS,IBT,2
        T1=IBN(I-2)
        T2=IBN(I-1)
        T3=IBN(I)
C        N1=IBN(I-2)+1
        N3=IBN(I)
c      write(0,*)ibn(i-2),ibn(i-1),ibn(i)
        DO N=N1,N3
          IF(N.GT.NCUT)GO TO 737
          IF(N.GE.NMN)THEN
c      write(0,*)'n=',n
            TN=N
            S1=(T2-TN)*(T3-TN)/((T2-T1)*(T3-T1))
            S2=(T1-TN)*(T3-TN)/((T1-T2)*(T3-T2))
            S3=(T1-TN)*(T2-TN)/((T1-T3)*(T2-T3))
            DO L=1,NBINM
              DO J=1,NENGP
                TS=S1*ALFN(J,I-2,L)+S2*ALFN(J,I-1,L)+S3*ALFN(J,I,L)
                IF(TS.GT.DZERO)ALF(J,L)=ALF(J,L)+TS         !JUST INCASE
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        IF(N3.GE.NVINT)GO TO 712
        N1=N3+1
      ENDDO
C
      IF(N3*NCUT.GT.N3*N3)THEN
        WRITE(6,215)NCUT,N3
        WRITE(0,215)NCUT,N3
      ENDIF
C
      GO TO 737
C
  712 CONTINUE
c
c      write(6,*)'plus interp to',nvint
c      write(6,*)(alf(j,1),j=1,nengp)
C
C SIMPSON'S RULE TO N=IBN(IBT)
C
      IS=I+2
      DO I=IS,IBT,2
        T1=IBN(I-2)*IBN(I-2)
        T3=IBN(I)*IBN(I)
        H=(T3-T1)/(T1*T3)
        H=H/DTWELV
        T=IBN(I-2)
        T1=T1*T*H
        T2=IBN(I-1)
        T2=T2*T2*T2
        T2=T2*DFOUR*H
        T=IBN(I)
        T3=T3*T*H
        DO L=1,NBINM
          DO J=1,NENGP
            ALF(J,L)=ALF(J,L)+T1*ALFN(J,I-2,L)+T2*ALFN(J,I-1,L)
     X                       +T3*ALFN(J,I,L)
          ENDDO
        ENDDO
      ENDDO
C
      DO L=1,NBINM
        DO J=1,NENGP
          TI=-ALFN(J,IS-2,L)/DTWO
          TF=DZERO
          IF(I-2.EQ.IBT)TF=ALFN(J,IBT,L)/DTWO
          ALF(J,L)=ALF(J,L)+TI+TF
        ENDDO
      ENDDO
c
c      write(6,*)'plus simp to',nint(t)
c      write(6,*)(alf(j,1),j=1,nengp)
C
 737  CONTINUE
C
C WRITE CROSS SECTION TOTALS
C
      DO L=1,NBINM
        WRITE(14,220)-(WNP(1)-WNP(L))*DKCM
        DO J=1,NENGP
          WRITE(14,219)ENERGP(J),ALF(J,L)              !E,E*SIGMA_RR(MB)
        ENDDO
      ENDDO
c
 2002 continue
C
C CONVERT TO RATE COEFF.
C
      DO L=1,NBINM
        DO J=1,NENGP
          PCS0(J)=ALF(J,L)                                 !AS ALF IS *4
        ENDDO
        CALL CONVOL(JTEMP,TEMP,ICON,TKAPPA,XDRY,NENG37,E37,F37
     X             ,ENERGP,PCS0,NENGP,IZERO,NCMX0+1,NZ)
        DO J=1,JTMAX
          ALF(J,L)=FREL(J)*PCS0(J)
        ENDDO
      ENDDO
C
C WRITE TOTALS TO UNIT6 AND ADF48
C
      IF(BLS)WRITE(6,729)NZ0,NE-1,(IWS(M)*(2*IWL(M)+1),M=1,NBINM)
      IF(BIC)WRITE(6,729)NZ0,NE-1,(IWJ(M)+1,M=1,NBINM)
C
      NMXW=MIN(NBINM,9)
      IF(JTEMP.GT.0)THEN
        WRITE(6,730)
        F732='(1PE10.2,1X,(10E10.2))'
        DO J=1,JTMAX
          WRITE(6,F732)TEMP(J)*CONRYK,(ALF(J,L),L=1,NBINM)
        ENDDO
        IF(NBINRM.NE.0)THEN
          WRITE(10,730)
          DO J=1,JTMAX
            WRITE(10,F732)TEMP(J)*CONRYK,(ALF(J,L),L=1,NMXW)
          ENDDO
        ENDIF
      ELSE
        WRITE(6,731)
        F732='(1PE10.2,1X,(10E10.2))'                   !N.B. DR 1PE10.3
        DO J=1,JTMAX
          WRITE(6,F732)ENERGP(J),(ALF(J,L),L=1,NBINM)
        ENDDO
        IF(NBINRM.NE.0)THEN
          WRITE(10,731)
          DO J=1,JTMAX
            WRITE(10,F732)ENERGP(J),(ALF(J,L),L=1,NMXW)
          ENDDO
        ENDIF
      ENDIF
C
      IF(NBINRM.NE.0)THEN
        WRITE(10,770)
        IF(JTEMP.GT.0)THEN
          IF(INAG.NE.0)THEN
            WRITE(10,*)
     X           'C*** WARNING *** WARNING *** WARNING *** WARNING ***'
            WRITE(10,*)
     X           'C INAG.NE.0: HIGH- N/T PARTIALS MAY NOT BE CONVERGED'
            WRITE(10,*)
     X           'C*** WARNING *** WARNING *** WARNING *** WARNING ***'
          ENDIF
          IF(ICON.GT.0)THEN
            IF(ICON.EQ.1)WRITE(10,772)TKAPPA
            IF(ICON.EQ.2)WRITE(10,773)XDRY
            IF(ICON.EQ.3)WRITE(10,774)                        !NUMERICAL
          ENDIF
          IF(IREL.GT.0)THEN
            WRITE(10,775)
            IF(ICON.NE.0)WRITE(6,776)
          ENDIF
        ENDIF
        WRITE(10,1005)(COD(I),I=2,20)
        WRITE(10,790)NAME,DATE
      ENDIF
C
      IF(JTHETA.LT.0)THEN                    !LOOP BACK UP AND CONVOLUTE
        JTHETA=-JTHETA
        JTEMP=JTHETA
        IREL=ABS(IREL)
        bskp=.true.
        CLOSE(10)
        OPEN(10,FILE='adf48_3'//nam0)            !SECONDLY A TYPE 3 FILE
        GO TO 2000
      ENDIF
C
      CLOSE(10)
C
      RETURN
C
C
   1  FORMAT(1X)
   2  FORMAT(85X,'--------------------------'/85X,'PRTF=',I2,2X,
     X'LVLPRT=',' (',I1,A1,1X,F4.1,')')
  18  FORMAT(85X,'--------------------------'/85X,'PRTF=',I2,2X,
     X'CFGPRT=',' (',F7.1,')')
  20  FORMAT(85X,'--------------------------'/85X,'PRTF=',I3,1X,
     X'LVLPRT=',' (',I1,A1,1X,F4.1,')')
   3  FORMAT(I6,5X,1P,10E10.2)
   5  FORMAT(11X,1P,10E10.2)                            !N.B. DR is 10.3
   6  FORMAT(11X,1P,10E10.2)
   7  FORMAT(3X,'ILREP'/3X,'-----')
   9  FORMAT(/1X,'--------------------------------'/
     X1X,'PRTI=',I2,2X,'TRMPRT=',' (',I1,A1,')',2X,'SPNPRT=',I2/)
  10  FORMAT(94X,'-----------------'/94X,'SYS=',I2,2X,'SPNSYS=',I2)
  11  FORMAT(3X,A5/3X,'-----')
  12  FORMAT(3X,A8,10(1PE10.2)/3X,'---- ---',10E10.2/)
 612  FORMAT(3X,A8,1P,10E10.2/3X,'---- ---',10E10.2)    !N.B. DR is 10.3
  13  FORMAT(/3X,'N-SHELL INDEXING',37X,'NREP=',I3
     X/3X,'----------------'/3X,A5,3X,'N'
     X/3X,'-----',3X,'-')
  14  FORMAT(/3X,'NL-SHELL INDEXING',21X,'NLREP=',I3
     X/3X,'-----------------'/3X,'ILREP',3X,'N',3X,'L'
     X/3X,'-----',3X,'-',3X,'-')
  15  FORMAT(I6,I6,9X,8(1PE10.2),(/21X,8E10.2))
  17  FORMAT(I6,I6,I4,5X,8(1PE10.2),(/21X,8E10.2))
  21  FORMAT("SEQ='",A2,"'",5X,"NUCCHG=",I2,50X,A4,5X,"ITYPE=",I1)
 212  FORMAT(/3X,'PARENT TERM INDEXING',17X,'BWNP=',F12.1,3X,'NPRNTI=',
     XI2,3X,'NPRNTF=',I2)
 213  FORMAT(/3X,'PARENT TERM INDEXING',17X,'BWNP=',F12.1,3X,'NPRNTI=',
     XI2,3X,'NPRNTF=',I3)
  23  FORMAT(
     X3X,'--------------------'/3X,'INDP',9X,'CODE',17X,'S L   WI',8X,
     X'WNP'/3X,'----',9X,'----',17X,'- -   --',2X,'----------')
  24  FORMAT(/3X,A2,' RESOLVED TERM INDEXING',12X,'BWNR=',F12.1,3X,'NTRM
     X=',I4/3X,'-------------------------'/3X,'INDX',2X,'INDP',3X,'CODE'
     X,17X,'S L   WJ',8X,'WNR'
     X/3X,'----',2X,'----',3X,'----',17X
     X,'- -   --',2X,'----------')
 251  FORMAT(/3X,A2,' AVERAGED CONFG INDEXING',11X,'BWNR=',F12.1,3X
     X,'NCFG=',I4/3X,'--------------------------'/3X,'INDX',2X,'INDP',3X
     X,'CODE',17X,'S     WJ',8X,'WNR'/3X,'----',2X,'----',3X,'----',17X
     X,'- ------',2X,'----------')
  26  FORMAT(I6,10X,5(A1,A1,A1,1X),'(',I1,')',I1,'(',F4.1,')',F11.1,A1)
  27  FORMAT(I6,10X,5(A1,A1,A1,1X),'(',F8.1,')',F11.1,A1)
  29  FORMAT(I6,I6,4X,5(A1,A1,A1,1X),'(',I1,')',A1,'(',F4.1,')',F11.1,A1
     X)
  30  FORMAT(I6,I6,4X,5(A1,A1,A1,1X),'(',F8.1,')',F11.1,A1)
  31  FORMAT(I6,I6,4X,5(A1,A1,A1,1X),'(',A1,F7.1,')',F11.1,A1)
  52  FORMAT(/1X,'--------------------------'/
     X1X,'PRTI=',I2,2X,'LVLPRT=',' (',I1,A1,1X,F4.1,')'/)
  53  FORMAT(/1X,'--------------------------'/
     X1X,'PRTI=',I2,2X,'CFGPRT=',' (',F7.1,')'/)
 215  FORMAT(/' *** WARNING:YOUR NMAX=',I4,' DOES NOT MATCH A'
     X,' REPRESENTATIVE N-VALUE.'/5X,'CUT-OFF APPLIED AT N=',I4/)
 219  FORMAT(1PE12.4,3E14.4)
 220  FORMAT('#',1PE16.6)
 223  FORMAT(/3X,'PARENT TERM INDEXING',17X,'BWNP=',F12.1,4X,'NPRNT=',I2
     X/3X,'--------------------'/3X,'INDP',9X,'CODE',17X,'S L   WI',8X,
     X'WNP'/3X,'----',9X,'----',17X,'- -   --',2X,'----------')
 224  FORMAT(/3X,'PARENT LEVEL INDEXING',16X,'BWNP=',F12.1,4X,'NPRNT=',
     XI2/3X,'---------------------'/3X,'INDP',9X,'CODE',17X,'S L   WI',
     X8X,'WNP'/3X,'----',9X,'----',17X,'- -   --',2X,'----------')
2332  FORMAT(/3X,'PARENT LEVEL INDEXING',16X,'BWNP=',F12.1,3X,'NPRNTI=',
     XI2,3X,'NPRNTF=',I2)
2333  FORMAT(/3X,'PARENT LEVEL INDEXING',16X,'BWNP=',F12.1,3X,'NPRNTI=',
     XI2,3X,'NPRNTF=',I3)
 233  FORMAT(
     X3X,'---------------------'/3X,'INDP',9X,'CODE',17X,'S L   WI',8X,'
     XWNP'/3X,'----',9X,'----',17X,'- -   --',2X,'----------')
 241  FORMAT(/3X,A2,' AVERAGED CONFG INDEXING',11X,'BWNR=',F12.1,3X
     X,'NCFG=',I4/3X,'--------------------------'/3X,'INDX',2X,'INDP',3X
     X,'CODE',17X,'      WJ',8X,'WNR'/3X,'----',2X,'----',3X,'----',17X
     X,'--------',2X,'----------')
2401  FORMAT(/3X,A2,' RESOLVED LEVEL INDEXING',11X,'BWNR=',F12.1,3X
     X,'NLVL=',I4)
2402  FORMAT(/3X,A2,' RESOLVED LEVEL INDEXING',11X,'BWNR=',F12.1,3X
     X,'NLVL=',I5)
 240  FORMAT(3X,'--------------------------'/3X,'INDX',2X,'INDP',3X
     X,'CODE',17X,'S L   WJ',8X,'WNR'/3X,'----',2X,'----',3X,'----',17X
     X,'- -   --',2X,'----------')
 250  FORMAT(/3X,A2,' RESOLVED CONFG INDEXING',11X,'BWNR=',F12.1,3X
     X,'NCFG=',I4/3X,'--------------------------'/3X,'INDX',2X,'INDP',3X
     X,'CODE',17X,'      WJ',8X,'WNR'/3X,'----',2X,'----',3X,'----',17X
     X,'--------',2X,'----------')
 253  FORMAT(/3X,'PARENT CONFG INDEXING',16X,'BWNP=',F12.1,3X,'NPRNTI=',
     XI2,3X,'NPRNTF=',I2/
     X3X,'---------------------'/3X,'INDP',9X,'CODE',17X,'      WI',8X,'
     XWNP'/3X,'----',9X,'----',17X,'--------',2X,'----------')
 254  FORMAT(/3X,'PARENT CONFG INDEXING',16X,'BWNP=',F12.1,4X,'NPRNT=',
     XI2/3X,'---------------------'/3X,'INDP',9X,'CODE',17X,'      WI',
     X8X,'WNP'/3X,'----',9X,'----',17X,'--------',2X,'----------')
 716  FORMAT(/1X,'-----------------------------------------'/
     X1X,'PRTF=',I2,2X,'TRMPRT=',' (',I1,A1,')',2X,'SPNPRT=',I2,2X
     X,'NSYS=',I2)
 726  FORMAT(/1X,'-----------------------------------------'/
     X1X,'PRTF=',I3,1X,'TRMPRT=',' (',I1,A1,')',2X,'SPNPRT=',I2,2X
     X,'NSYS=',I2)
 727  FORMAT(/1X,'------------------------------------'/
     X1X,'PRTF=',I2,2X,'CFGPRT=',' (',F7.1,')',2X,'NSYS=',I2)
 729  FORMAT('Z=',I2,1X,'N=',I2,1X,'W:',10(I3))
 730  FORMAT(//'    T(K) ',4X,'ALFT( 1)',2X,'ALFT( 2)',2X,'ALFT( 3)'
     X,2X,'ALFT( 4)',2X,'ALFT( 5)',2X,'ALFT( 6)'
     X,2X,'ALFT( 7)',2X,'ALFT( 8)',2X,'ALFT( 9)'
     X/4X,'----',3X,9(2X,'--------'))
 731  FORMAT(//'   E(RYD)',4X,'E*SIG(1)',2X,'E*SIG(2)',2X,'E*SIG(3)'
     X,2X,'E*SIG(4)',2X,'E*SIG(5)',2X,'E*SIG(6)'
     X,2X,'E*SIG(7)',2X,'E*SIG(8)',2X,'E*SIG(9)'
     X/3X,'------',2X,9(2X,'--------'))
C 732  FORMAT(1PE10.2,1X,(10E10.2))
 770  FORMAT('C',110('-')/'C')
 772  FORMAT('C *** RATE COFFICIENTS FOR A KAPPA DISTRIBUTION, WITH K='
     X,F6.1/'C')
 773  FORMAT('C *** RATE COFFICIENTS FOR A DRUYVESTEYN DISTRIBUTION,'
     X,' WITH X=',F6.1/'C')
 774  FORMAT('C *** RATE COFFICIENTS FOR A NUMERICAL DISTRIBUTION...')
 775  FORMAT('C     JUTTNER RELATIVISTIC CORRECTION APPLIED TO THE',
     X ' DISTRIBUTION'/'C')
 776  FORMAT(/'C *** ATTENTION: YOU ARE APPLYING A JUTTNER RELATIVISTIC'
     X      ,' CORRECTION TO A NON-MAXWELLIAN DISTRIBUTION...'/'C')
 790  FORMAT('C'/'C',1X,A30/'C',1X,A30/'C',110('-'))
 1005 FORMAT('C',19A4)
C
      END
C
C***********************************************************************
C
      SUBROUTINE CONVOL(JTEMP,TEMP,ICON,TKAPPA,XDRY,NENG37,E37,F37
     X                 ,ENERGP,PCS0,NENGP,LX,NX,NZ)
C
C-----------------------------------------------------------------------
C
CNRB
C  CONVOLUTE RR CROSS SECTION WITH ELECTRON ENERGY DISTRIBUTION ICON.
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (NDIM33=46) !PI XSCTNS
      PARAMETER (IMAX=9)    !MAX NO SETS OF GAUSS-LAGUERRE PTS INAG.GT.0
      PARAMETER (MMAX=100) !MAX NO OF ITERTNS OF MIDPOINT RULE INAG.LT.0
C      PARAMETER (NMAX=2)  !SET=2 TO TEST CONVRGNCE OF NEWTON-COTES STEP
C     DIMENSION OF D (BELOW) SHOULD BE MAX(IMAX,MMAX,NMAX)
C
      PARAMETER (DZERO=0.0D0)
      PARAMETER (DONE=1.0D0)
      PARAMETER (DEIGHT=8.0D0)
      PARAMETER (DNINE=9.0D0)
      PARAMETER (BIG=1.D30)
C
      PARAMETER (T2ORTPI=1.1283792D0)     !2/sqrt(pi)
      PARAMETER (CONQ=21.876870D-11)       !a_0/tau_0
      PARAMETER (CONRYK=1.5789D5) !1.578885D5)      !RYDBERGS TO KELVIN
C
      LOGICAL BFIRST
C
      DIMENSION TEMP(*),ENERGP(*),PCS0(*)
      DIMENSION D(0:MMAX)
      DIMENSION N(IMAX),WLAG(180),XLAG(180) !see data
      DIMENSION DUM(NDIM33)
C
      DIMENSION E37(NENG37,*),F37(NENG37,*)

c     dimension tpow(ndim33)
C
      COMMON /CSPLYN/SP1(NDIM33),SP2(NDIM33),SP3(NDIM33),SP4(NDIM33),JSP
      COMMON /CINT/ENERG(NDIM33),PCS(NDIM33),BSP,NENG,LXTRP
C
      COMMON /NRBTST/BSP0,INAG,NCP,NLAG,IMESH
C
      DATA BFIRST/.TRUE./
      DATA (N(I),I=1,IMAX)/16,20,24,32,48,64,96,128,180/
C
      IF(JTEMP.EQ.0)RETURN
C
      IF(INAG.GT.IMAX)THEN
        WRITE(6,*)'***CONVOL, TOO MANY G-L POINTS REQUESTED:',INAG,IMAX
        STOP 'ERROR: ***CONVOL, TOO MANY G-L POINTS REQUESTED:'
      ENDIF
      IF(INAG.LT.-MMAX)THEN
        WRITE(6,*)'***CONVOL, TOO MANY MID-POINTS REQUESTED:',INAG,MMAX
        STOP 'ERROR: ***CONVOL, TOO MANY MID-POINTS REQUESTED:'
      ENDIF
C
      LXTRP=LX
      LP=MAX(1,LXTRP+1)
      BSP=BSP0
      NENG=NENGP
C
      DO M=1,NENG
        ENERG(M)=ENERGP(M)
        PCS(M)=PCS0(M)
      ENDDO
C
      IF(BSP0.NE.DZERO)THEN                 !SPLINE SCALED CROSS SECTION
        IF(BSP0.LT.DZERO)THEN
          DE=DONE
          IF(LXTRP.GE.0)DE=SQRT(ENERG(NENG))
          ECRIT=NZ
          ECRIT=ECRIT/NX
          ECRIT=ECRIT*ECRIT
C
          BSP=-BSP0
          IF(PCS(1)*PCS(NENG).NE.DZERO)THEN
            BSP=BSP*ECRIT*NX*PCS(NENG)*DE*ENERG(NENG)**LP
     X          /(PCS(1)-PCS(NENG))
c           BSP=BSP*ECRIT*PCS(1)/PCS(NENG)/NX/NX
          ENDIF
        ENDIF
        DE=DONE
c        pold=1
c        eold=1
        DO M=1,NENG
          IF(LXTRP.GE.0)DE=SQRT(ENERG(M))
c          tpow(m)=log(pold/pcs(m))/(log(energ(m)/eold))
c          pold=pcs(m)
c          eold=energ(m)
          PCS(M)=PCS(M)*(BSP+DE*ENERG(M)**LP)
        ENDDO
c        write(73,111)nx,lx,(tpow(m),m=1,neng,3)
c        write(73,111)nx,lx,(pcs(m),m=1,neng,3)
c111     format(2i5,30(1pe10.2))
        CALL SPLYN(NENG,ENERG,PCS,3,DZERO,3,DZERO,SP1,SP2,SP3,SP4,DUM)
      ENDIF
C
      DO J=1,JTEMP
C
        TEMPE=TEMP(J)
        B=DONE/TEMPE
        TE=T2ORTPI/(TEMPE*SQRT(TEMPE))         !AS MIDEXP/GAUSS-LAG OMIT
        D(0)=BIG
C
        IF(INAG.GT.0)GO TO 1
        IF(INAG.LT.0)GO TO 2
C
C NEWTON-COTES (BODE)
C
        TE=DONE     !as absorb 2/e*sqrt(pi*e) in general fbar distrib.
C
        NMAX=1
        IF(NCP.LT.0.AND.NCP.EQ.MOD(NCP,10))NMAX=2
        JT=J
        IF(IMESH.GT.0)JT=-J
        DO I=1,NMAX
          CALL BODE(JT,TEMPE,ICON,TKAPPA,XDRY,NENG37,E37(1,J),F37(1,J)
     X             ,I,NCP,NX,LX,NZ,D(I))
          IF(ABS(D(I)/D(I-1)-DONE).LT.0.01)GO TO 3
          IF(D(I).EQ.DZERO)GO TO 3
        ENDDO
        I=NMAX
        IF(NMAX.GT.1)THEN
          IF(BFIRST)THEN
            WRITE(6,105)TEMP(J)*CONRYK
            BFIRST=.FALSE.
            WRITE(6,104)(JJ,D(JJ),JJ=1,I)
          ENDIF
          WRITE(66,104)(JJ,D(JJ),JJ=1,I)
        ENDIF
        GO TO 3
C
C GAUSS-LAGUERRE
C
    1   DO I=1,INAG
          JSP=0
          D(I)=GAUSSQ(N(I),WLAG,XLAG,B)
          IF(ABS(D(I)/D(I-1)-DONE).LT.0.01)GO TO 3
          IF(D(I).EQ.DZERO)GO TO 3
        ENDDO
        I=INAG
        IF(BFIRST)THEN
          WRITE(6,102)TEMP(J)*CONRYK
          BFIRST=.FALSE.
          WRITE(6,100)(N(JJ),D(JJ),JJ=1,I)
        ENDIF
        WRITE(66,100)(N(JJ),D(JJ),JJ=1,I)
        GO TO 3
C
C MID-EXP RULE
C
   2    OST=DZERO
        DO I=1,-INAG
          JSP=0
          CALL MIDEXP(I,B,ST)
          S=(DNINE*ST-OST)/DEIGHT
          D(I)=S
          IF(ABS(D(I)/D(I-1)-DONE).LT.0.005)GO TO 3
          IF(D(I).EQ.DZERO)GO TO 3
          OST=ST
        ENDDO
        I=-INAG
        IF(BFIRST)THEN
          WRITE(6,103)TEMP(J)*CONRYK
          BFIRST=.FALSE.
          WRITE(6,101)(JJ,D(JJ),JJ=1,I)
        ENDIF
        WRITE(66,101)(JJ,D(JJ),JJ=1,I)
C
C
   3    IF(D(I).LT.DZERO)D(I)=DZERO
        PCS0(J)=CONQ*TE*D(I)
C
      ENDDO
C
      RETURN
C
 100  FORMAT(' SR.CONVOL: GAUSS-LAGUERRE FAILURE; N(I),D(I)='
     X,10(I4,1PE10.3))
 101  FORMAT(' SR.CONVOL: MIDEXP FAILURE; I,D(I)=',10(I4,1PE10.3))
 102  FORMAT(' SR.CONVOL: GAUSS-LAUGUERRE FIRST FAILS AT TEMP (K)='
     X,1PE10.2/)
 103  FORMAT(' SR.CONVOL: MID-POINT RULE FIRST FAILS AT TEMP (K)='
     X,1PE10.2/)
 104  FORMAT(' SR.CONVOL: NEWTON-COTES FAILURE; I,D(I)='
     X,10(I4,1PE10.3))
 105  FORMAT(' SR.CONVOL: NEWTON-COTES FIRST FAILS AT TEMP (K)='
     X,1PE10.2/)
C
      END
C
C***********************************************************************
C
      SUBROUTINE BODE(JT0,TEMPE,ICON,TKAPPA,XDRY,NENG37,E37,F37
     X               ,J2,ICP,NX,LX,NZ,D)
C
C-----------------------------------------------------------------------
C
C N R BADNELL
C INPUT:
C   JT:  TEMPERATURE INDEX
CTEMPE: = KT
C ICON: CONVOLUTION DISTRIBUTION, DEFINED BY KAPPA OR XDRY
C   J2: = 1 NORMALLY (SINGLE PASS) =2 A RETURN TO CHECK CONVERGENCE
C  NCP: NO. POINTS OF NEWTON-COTES FORMULA, 3 OR 5.
C   NX: PRINCIPAL QUANTUM NO. OF LOWER STATE
C   LX: ORBITAL ANG. MOM. OF LOWER STATE (OR LOWEST OF WHEN BUNDLED)
C   NZ: RESIDUAL CHARGE OF ELECTRON TARGET
C OUTPUT:
C   D: QUADRATURE RESULT OF DISTRIBUTION INTEGRAL AT INPUT TEMP.
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (MAXK=1000)
C
      PARAMETER (IZERO=0)
C
      PARAMETER (DZERO=0.0D0)
      PARAMETER (DONE=1.0D0)
      PARAMETER (DTWO=2.0D0)
      PARAMETER (DTHREE=3.0D0)
      PARAMETER (DTEN=1.0D1)
C
      PARAMETER (T2ORTPI=1.1283792D0)     !2/sqrt(pi)
      PARAMETER (CONRYK=1.5789D5) !1.578885D5)      !RYDBERGS TO KELVIN
C
      LOGICAL BFIRST
C
      DIMENSION E37(*),F37(*)
      DIMENSION FST(0:MAXK),XK(0:MAXK),IWT(0:4)
C
      SAVE
      DATA BFIRST/.TRUE./,M00/-999/
      DATA THT0/1.0D1/ !LOWEST ADAS SCALED TEMP AS SET IN CROSSJ - SYNC!
C
      TEMP10=10*TEMPE
C
      IF(J2.GT.1)BFIRST=.TRUE.
C
C INITIALIZE NEWTON-COTES
C
      IF(BFIRST)THEN
        BFIRST=.FALSE.
        MOLD=-999
        NCP=MOD(ABS(ICP),10)
        NCPM=NCP-1
        IF(NCP.EQ.3)THEN
          IH0=10                  !11 FOR 3-PT DOUBLED STEP LENGTH
          IH0=IH0+J2              !12 for 3-pt standard length (steps*2)
          IWT(0)=1
          IWT(1)=4
          IWT(2)=1
          FFFF=DONE/3
        ELSEIF(NCP.EQ.5)THEN
          IH0=11                  !12 FOR 5-PT STANDARD STEP LENGTH
          IH0=IH0+J2              !13 for 5-pt halved step (steps*2)
          IWT(0)=7
          IWT(1)=32
          IWT(2)=12
          IWT(3)=32
          IWT(4)=7
          FFFF=DTWO/45
        ELSE
          WRITE(6,*)'SR.BODE: ONLY 3- OR 5-PT NEWTON-COTES ALLOWED:',NCP
          STOP 'ERROR: SR.BODE: ONLY 3- OR 5-PT NEWTON-COTES ALLOWED:'
        ENDIF
        TMULT=LOG(DTEN)/DTHREE    !ADAS DEFAULT TEMP-SPACING (APPROX)
      ENDIF
C
C SET ENERGY MESH (BASED-ON BURGESS MEM.R.ASTR.SOC. V69,1 (1964))
C
      JT=IABS(JT0)
      IF(JT0.LT.0)JT=JT/2+1
      IF(JT.EQ.1)THEN
        MOLD=-999
        TEMP0=TEMPE
        T0=NZ*NZ*THT0/(TEMPE*CONRYK)
        T=LOG(T0)/LOG(DTWO)
        M00=NINT(T)               !=0 IF DEFAULT ADAS TEMP(1) USED
        M00=MAX(M00,0)
        JM=0
      ELSEIF(M00.EQ.-999)THEN
        WRITE(6,*)'SR.BODE: FIRST CALL MUST BE WITH TEMP(JT=1), JT=',JT
        STOP 'ERROR: SR.BODE: FIRST CALL MUST BE WITH TEMP(1)'
      ELSE
        T=LOG(TEMPE/TEMP0)/TMULT
        JM=NINT(T)                !=JT IF DEFAULT ADAS TEMP(JT) USED
      ENDIF
C
      JM=JM+1
      M0=M00
      IF(NX.EQ.1)THEN
        MX=26
        M0=M0+MAX(4-JM,0)
      ELSEIF(NX.EQ.2)THEN
        MX=25
        M0=M0+MAX(3-JM,0)
      ELSEIF(NX.LT.5)THEN
        MX=24
        M0=M0+MAX(2-JM,0)
      ELSEIF(NX.LT.10)THEN
        MX=23
      ELSEIF(NX.LT.17)THEN
        MX=22
      ELSEIF(NX.LT.81)THEN
        MX=21
      ELSEIF(NX.LT.181)THEN
        MX=22
      ELSEIF(NX.LT.400)THEN
        MX=23
      ELSEIF(NX.LT.900)THEN
        MX=24
      ELSE
        MX=25
      ENDIF
      mx=mx+icp/10                !test further in/out -/+
c      m0=m0+1                     !test finer near threshold
      H0=NZ*NZ
      H0=H0/NX
      H0=H0/2**(IH0+M0+1)
      MM=MX+M0
C
C STORE INTEGRAND (E*SIGMA_RR)
C
      IF(M0.NE.MOLD)THEN
c        if(nx.gt.20)write(64,*)nx
        B=DONE/TEMPE
        H=H0
        XK(0)=DZERO
        K=0
        FST(K)=GLNAG(B,XK(K))
        DO M=1,MM
          H=H+H
          DO J=1,J2
            DO N=1,NCPM
              K=K+1
              XK(K)=XK(K-1)+H
              FST(K)=GLNAG(B,XK(K))
c        if(nx.gt.20)
c     x write(64,50)k,xk(k),fst(k),fst(k)*xk(k),fst(k)*xk(k)*sqrt(xk(k))
c  50    format(i5,4d14.5)
            ENDDO
          ENDDO
        ENDDO
        KMAX=K
        if(mod(kmax,ncpm).ne.0)stop 'array index error'
        IF(KMAX.GT.MAXK)THEN
          WRITE(6,*)'SR.BODE: INCREASE MAXK TO ', KMAX
          STOP 'ERROR: SR.BODE: INCREASE MAXK'
        ENDIF
      ENDIF
C
C APPLY QUADRATURE
C
      H=H0
      D=DZERO
      DT=DZERO
      DO K=NCPM,KMAX,NCPM*J2
        H=H+H
        D0=DZERO
        DO J=1,J2
          N1=(1-J)*NCPM
          N2=N1+NCPM
          DO N=N1,N2
            D0=D0+IWT(N-N1)*FST(K-N)
     X           *FBAR(XK(K-N),TEMPE,ICON,TKAPPA,XDRY,NENG37,E37,F37)
          ENDDO
        ENDDO
        D=D+D0*H
        IF(XK(K).GT.TEMP10)GO TO 100
      ENDDO
C
C ADD-IN HIGH-E TOP-UP
C
      IF(NX.GT.20)THEN
c        write(0,*)'nx=',nx
        LP=MAX(1,LX+1)
        X=XK(KMAX)/TEMPE
        ET=ENINT(LP,.TRUE.,X)             !LXTRP.GE.0
        IF(ET.LE.DZERO)THEN
          IF(ET.LT.DZERO)WRITE(6,*)'EN_INT FAILED TO CONVERGE:',LP,X
          DT=DZERO
          GO TO 100
        ENDIF
        DT=T2ORTPI*ET*FST(KMAX)*X/SQRT(TEMPE)
        IF(ICON.GT.0)THEN                        !COMPENSATE FOR NON-MAX
          X0=TEMPE
          FF=FBAR(X0,TEMPE,ICON,TKAPPA,XDRY,NENG37,E37,F37)     !NON-MAX
          FM=FBAR(X0,TEMPE,IZERO,TKAPPA,XDRY,NENG37,E37,F37)        !MAX
          R0=FF/FM
          X2=3*X0/2
          FF=FBAR(X2,TEMPE,ICON,TKAPPA,XDRY,NENG37,E37,F37)     !NON-MAX
          FM=FBAR(X2,TEMPE,IZERO,TKAPPA,XDRY,NENG37,E37,F37)        !MAX
          R2=FF/FM
C
C LINEAR C1*X+C0 (NOT GOOD)
C          C0=(X2*R0-X0*R2)/(X2-X0)
C          C1=(R2-R0)/(X2-X0)
c          et1=enint(lp-1,.true.,x)
c          write(0,*)et1
C          ET1=(EXP(-X)-(2*LP-1)*ET/2)/X
C          DT0=C0*DT
C          DT1=C1*DT*XK(KMAX)*ET1/ET
c          write(10,'(1p,6d12.3)')x0,r0,c0,x2,r2,c1
c          write(10,'(1p,3d12.3," sum=",d12.3)')dt,dt0,dt1,dt0+dt1
C          DT=DT0+DT1
C
C QUADRATIC C2*X**2+C1*X+C0 (GOOD ENOUGH)
          X1=LOG(X0)
          T=LOG(X2)-X1
          X1=X1+T/2
          X1=EXP(X1)
          FF=FBAR(X1,TEMPE,ICON,TKAPPA,XDRY,NENG37,E37,F37)     !NON-MAX
          FM=FBAR(X1,TEMPE,IZERO,TKAPPA,XDRY,NENG37,E37,F37)        !MAX
          R1=FF/FM
          C0=X1*X2*(X1-X2)*R0+X0*X2*(X2-X0)*R1+X0*X1*(X0-X1)*R2
          C1=X2*X2*(R0-R1)+X1*X1*(R2-R0)+X0*X0*(R1-R2)
          C2=(X1-X2)*R0+(X2-X0)*R1+(X0-X1)*R2
          T=(X0-X1)*(X0-X2)*(X1-X2)
          C0=C0/T
          C1=C1/T
          C2=C2/T
          DT0=C0*DT
          T=DT/ET
          ET1=(EXP(-X)-(2*LP-1)*ET/2)/X
          ET2=(EXP(-X)-(2*LP-3)*ET1/2)/X
          T=T*XK(KMAX)
          DT1=C1*ET1*T
          DT2=C2*ET2*T*XK(KMAX)
C
CTEST C2/X**2+C1/X+C0
C          ET1=2*(EXP(-X)-X*ET)/(2*LP+1)
C          ET2=2*(EXP(-X)-X*ET1)/(2*LP+3)
C          T=T/XK(KMAX)
C          DT1=C1*ET1*T
C          DT2=C2*ET2*T/XK(KMAX)
c
c          write(10,'(1p,9d12.3)')x0,r0,c0,x1,r1,c1,x2,r2,c2
c          write(10,'(1p,5d12.3," sum=",d12.3)')
c     x                            d*ffff,dt,dt0,dt1,dt2,dt0+dt1+dt2
          DT=DT0+DT1+DT2
        ENDIF
      ENDIF
C
 100  D=D*FFFF
      D=D+DT
C
      MOLD=M0
      IF(J2.GT.1)BFIRST=.TRUE.
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION GLNAG(B,X)
C
C-----------------------------------------------------------------------
C
CNRB
C FUNCTION TO EVALUATE E*SIGMA E.G. FOR ENERGY CONVOLUTION INTEGRALS
C USED BY NEWTON-COTES (SR.BODE) AS WELL AS ORIGINAL GAUSS-LAGUERRE.
C (B=1/KT IS ONLY USED NOW TO FLAG EXCESSIVE EXTRAPOLATION IN EENRGY.)
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (NDIM33=46)   !PI XSCTNS
C
      PARAMETER (DZERO=0.0D0)
      PARAMETER (DONE=1.0D0)
C
      COMMON /CSPLYN/SP1(NDIM33),SP2(NDIM33),SP3(NDIM33),SP4(NDIM33),JSP
      COMMON /CINT/ENERG(NDIM33),PCS(NDIM33),BSP,NENG,LXTRP
C
      COMMON /NRBTST/BSP0,INAG,NCP,NLAG,IMESH
C
      DATA IFLAG/1/
      DATA NOLD/1/,XOLD/1.D99/
C
      NP=NLAG/2
C
      LP=MAX(1,LXTRP+1)
      DX=DONE
      IF(LXTRP.GE.0)DX=SQRT(X)
      TT=DZERO
C
      IF(X.LT.ENERG(NENG))THEN                              !INTERPOLATE
C
        IF(X.GT.ENERG(1))THEN
          IF(BSP0.NE.DZERO)THEN                                  !SPLINE
            TT=SPVAL(NENG,ENERG,PCS,SP1,SP2,SP3,SP4,X,JSP)
            TT=TT/(BSP+DX*X**LP)
            GO TO 2
          ENDIF
C LAGRANGE
          IF(X.LT.XOLD)NOLD=1
          XOLD=X
          DO N=NOLD,NENG                    !FIND POSITION OF X IN ENERG
            IF(ENERG(N).GT.X)THEN
              N2=N+NP-1
              N1=N-NP
              IF(N1.LE.0)THEN
                N2=NLAG
                N1=1
              ELSEIF(N2.GT.NENG)THEN
                N2=NENG
                N1=N2-NLAG+1
              ENDIF
              NOLD=N
              GO TO 1
            ENDIF
          ENDDO
          NOLD=NENG
          N2=NENG
          N1=N2-NLAG+1
C WEIGHTS
    1     EP=DZERO
          IF(N1.GT.1.AND.PCS(N1)*PCS(N2).NE.DZERO)THEN
C                   .AND.PCS(N2).LT.PCS(N1)/1.5
            EP=-LOG(PCS(N1)/PCS(N2))/LOG(ENERG(N1)/ENERG(N2))
          ENDIF
          DO N=N1,N2
            DD0=DONE
            DO M=N1,N2
              IF(N.NE.M)THEN
                DD0=DD0*(X-ENERG(M))
                DD0=DD0/(ENERG(N)-ENERG(M))
              ENDIF
            ENDDO
            IF(ENERG(N).GT.DZERO)DD0=DD0*ENERG(N)**EP
            TT=TT+DD0*PCS(N)
          ENDDO
          TT=TT/X**EP
C
        ELSE      !ASSUME CONSTANT NEAR THRESH (ALREADY TESTED ENERG(1))
C
          TT=PCS(1)
          IF(BSP0.NE.DZERO)THEN                                  !SPLINE
            DE=DONE
            IF(LXTRP.GE.0)DE=SQRT(ENERG(1))
            TT=TT/(BSP+DE*ENERG(1)**LP)
          ENDIF
C
        ENDIF
C
      ELSE                                                  !EXTRAPOLATE
C
        IF(B*X.LT.5.0.AND.IFLAG.EQ.1)THEN
          WRITE(6,100)ENERG(NENG)
          IFLAG=0
        ENDIF
C
C E*RR=TT/E**(L+1.5)
C
        DE=DONE
        IF(LXTRP.GE.0)DE=SQRT(ENERG(NENG))
        TT=PCS(NENG)
        IF(BSP0.NE.DZERO)TT=TT/(BSP+DE*ENERG(NENG)**LP)          !SPLINE
        TT=TT*DE*ENERG(NENG)**LP
C
C E*RR=(A1+A2/X)/X**(L+1.5)
C
c        NENGM=NENG-1
c        DE1=DONE
c        DE2=DONE
c        IF(LXTRP.GE.0)DE1=SQRT(ENERG(NENGM))
c        IF(LXTRP.GE.0)DE2=SQRT(ENERG(NENG))
c        E1=DE1*ENERG(NENGM)**LP
c        E2=DE2*ENERG(NENG)**LP
c        P1=PCS(NENGM)*E1
c        IF(BSP0.NE.DZERO)P1=P1/(BSP+DE1*ENERG(NENG-1)**LP)      !SPLINE
c        P2=TT                                      !PCS(NENG)*E2
c        A2=ENERG(NENGM)*ENERG(NENG)*(P1-P2)/(ENERG(NENG)-ENERG(NENGM))
c        A1=P1-A2/ENERG(NENGM)
c        TT=A1+A2/X
C
        TT=TT/(DX*X**LP)
      ENDIF
C
   2  GLNAG=TT
C
      RETURN
 100  FORMAT(' SR.GLNAG: NOTE, PCS EXTRAPOLATED USING ENERG(NENG)='
     X       ,1PE9.1/)
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION FBAR(EEE,TEMPE,ICON,KAPPA,XDRY,NENG37,E37,F37)
C
C-----------------------------------------------------------------------
C
C EVALUATES THE ELECTRON ENERGY DISTRIBUTION FUNCTION FBAR AT ENERGY EEE
C AT AN EFFECTIVE TEMPERATURE TEMPE FOR MAXWELLIAN, KAPPA, DRUYVESTEYN &
C NUMERICAL.
C N.B. AS "MUST", OMITS THE FACTOR SQRT(EEE) - HENCE FBAR.
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 KAPPA,KAPPAP
C
      PARAMETER (DZERO=0.0D0)
      PARAMETER (DONE=1.0D0)
      PARAMETER (D0PT5=0.5D0)
      PARAMETER (D1PT5=1.5D0)
      PARAMETER (D2PT5=2.5D0)
      PARAMETER (D1M5=1.0D-5)
C
      PARAMETER (T2ORTPI=1.1283792D0)                        !2/sqrt(pi)
C
      PARAMETER (ARGFCT=57.0D0)
      PARAMETER (EXPFCT=650.0D0)
      PARAMETER (YPARL=80.0D0)
      PARAMETER (YPAR=20)                                 !MAX (E/KT)**X
C
      PARAMETER (NLAG=4)
      PARAMETER (NLAGH=NLAG/2)
C
      DIMENSION E37(*),F37(*)
C
      DATA FRACK/-1/,FRACX/-1/            !SINCE DISTRIB DOES NOT CHANGE
      DATA NOLD/1/,EOLD/1.D99/,TOLD/1.D99/
C
      SAVE T53
C
      IF(ICON.EQ.0)THEN                                      !MAXWELLIAN
        TSQRTE=SQRT(TEMPE)
        TE=EEE/TEMPE
        IF(TE.LT.75.0D0)THEN
          QUP=EXP(-TE)/(TSQRTE*TEMPE)
        ELSE
          QUP=DZERO
        ENDIF
      ELSEIF(ICON.EQ.1)THEN                                       !KAPPA
        IF(FRACK.LT.DZERO)THEN
          KAPPAP=KAPPA+DONE
          IF(KAPPAP.LT.ARGFCT)THEN
            FRACK=GAMA(KAPPAP)/(GAMA(KAPPA-D0PT5))
          ELSE
            FRACK=KAPPA*SQRT(KAPPA)
          ENDIF
        ENDIF
        TEMPK=TEMPE*(KAPPA-D1PT5)          !=KAPPA*CHARACTERISTIC ENERGY
        TSQRTK=SQRT(TEMPK)
        TE=EEE/TEMPK
        T=(DONE+TE)
        IF(KAPPAP*LOG(T).LT.EXPFCT)THEN
          QUP=DONE/T**KAPPAP
          QUP=FRACK*QUP/(TSQRTK*TEMPK)
        ELSE
          QUP=DZERO
        ENDIF
      ELSEIF(ICON.EQ.2)THEN                                 !DRUYVESTEYN
        IF(FRACX.LT.DZERO)THEN
          C=DONE/T2ORTPI
          T5=D2PT5/XDRY
          T3=D1PT5/XDRY
          T5=GAMA(T5)
          T3=GAMA(T3)
          T53=T5/T3
          FRACX=C*XDRY/T3
        ENDIF
        TEMPX=3*TEMPE/(2*T53)                !=CHARACTERISTIC ENERGY/T53
        TSQRTX=SQRT(TEMPX)
        TE=EEE/TEMPX
        IF(XDRY*LOG(TE+D1M5).LT.YPARL)THEN       !AS EEE CAN BE ZERO NO
          T=TE**XDRY
        ELSE
          T=2*YPAR
        ENDIF
        IF(T.LT.YPAR)THEN
          QUP=EXP(-T)
          QUP=FRACX*QUP/(TSQRTX*TEMPX)
        ELSE
          QUP=DZERO
        ENDIF
C
      ELSEIF(ICON.EQ.3)THEN                                   !NUMERICAL
        QUP=0
        N1=1
        N2=0
        IF(EEE.LT.E37(1).OR.EEE.GE.E37(NENG37))then
c        if(eee.lt.e37(1))go to 1 !test high-t sensitivity to low-e tail
c fill-in maxwellian test-case - no need to renorm.
          tsqrte=sqrt(tempe)
          te=eee/tempe
          if(te.lt.75.0d0)then
            qup=t2ortpi*exp(-te)/(tsqrte*tempe)
          else
            qup=dzero
          endif
c need non-maxwellian extrapolate option, esp. high-e for small kappa?
c Unlikely to need to renorm, except kappa-like->3/2?
          GO TO 1
        endif
        IF(EEE.LE.EOLD.OR.TEMPE.LT.TOLD)NOLD=1
        EOLD=EEE
        DO N=NOLD,NENG37
          IF(EEE.LT.E37(N))THEN
            N2=N+NLAGH-1
            N1=N-NLAGH
            IF(N1.LE.0)THEN
              N2=NLAG
              N1=1
            ELSEIF(N2.GT.NENG37)THEN
              N2=NENG37
              N1=N2-NLAG+1
            ENDIF
            NOLD=N
            GO TO 1
          ENDIF
        ENDDO
        NOLD=NENG37
        N2=NENG37
        N1=N2-NLAG+1
   1    CONTINUE
        DO N=N1,N2
          DD0=DONE
          DO M=N1,N2
            IF(N.NE.M)THEN
              DD0=DD0*(EEE-E37(M))
              DD0=DD0/(E37(N)-E37(M))
            ENDIF
          ENDDO
c          dd0=dd0*sqrt(e37(n))               !interp full f(e) e.g. SOL
          QUP=QUP+DD0*F37(N)
        ENDDO
c        if(eee.ne.dzero)qup=qup/sqrt(eee)        !if interp full f
        QUP=QUP/T2ORTPI
      ENDIF
C
      FBAR=T2ORTPI*QUP                            !*SQRT(EEE) FOR FULL F
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE ADF37(IPRINT,J,TEMPE37,ICON,NENG37,E37,F37)
C
C CHECK DISTRIBUTION NORMALIZATION AND TRY AND SET EFFECTIVE TEMPS
C TO HELP LATER RR QUADRATURE (SEE SR.BODE). MORE COSMETIC FOR DR.
C NORMALLY CALLED WITH ICON=3, AS ANALYTIC SHOULD BE KNOWN/EXACT.
C BUT MAY HAVE ALTERNATE NUMERICAL ICON FLAGGED IN THE FUTURE.
C
C KNOWN ISSUES FOR NUMERICAL DISTRIBUTIONS:
C
C CONTRIBUTION FROM THE TAILS I.E.OUTSIDE THE RANGE OF ENERGY TABULATION
C CURRENTLY, NO EXTRAPOLATION IS CARRIED-OUT - A MAXWELLIAN IS COMMENTED
C OUT FOR THE ADF37 TEST CASE TO ILLUSTRATE ITS EFFECT/IMPORTANCE.
C
C INTERPOLATION OF THE DISTRIBUTION F(E) OR F(E)/SQRT(E) ETC (WHERE F(E)
C IS NORMED TO UNITY). THE FORMER IS NECESSARY FOR THE ADF37 SOL EXAMPLE
C (NORM IS WILDLY OFF OTHERWISE) BUT THE LATTER IS MUCH MORE ACCURATE
C AT HIGH-E FOR THE MAXWELLIAN CASE. WOULD EXPECT THE LATTER TO BE
C PREFERABLE. THE SOL EXAMPLE HAS VERY FEW POINTS AT LOW E, SO SUSPECT.
C
C SEE FN.FBAR FOR CURRENT EXTRAPOLATION & INTERPOLATION AND DEVELOPMENT.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (IZERO=0)
      PARAMETER (DONE=1.0D0)
C
      PARAMETER (TNORM=1.0D-3)       !RENORM F(E) IF ERROR EXCEEDS TNORM
      PARAMETER (TMEAN=0.1D-0)     !RECALC E*F(E) IF ERROR EXCEEDS TMEAN
      PARAMETER (ITMAX=5)                 !FOR ITERATION ON TEMPE=E*F(E)
c
c      parameter (conryk=1.5789d5) !1.578885d5)      !rydbergs to kelvin
C
      DIMENSION E37(*),F37(*)
C
      DATA XDUM/0/                            !NOT USED
C
        FMAX=-1
        NFMAX=0
        DO N=1,NENG37
          F=F37(N)*SQRT(E37(N))          !AS F37 IS FBAR NOW
          IF(F.GT.FMAX)THEN
            FMAX=F
            NFMAX=N
          ENDIF
        ENDDO
        TEMPE37=2*E37(NFMAX)               !CASE MAXWELLIAN
C
      NSTEP=4001
C
      ITER=0
   2  TEMPE=TEMPE37
      ITER=ITER+1
C
      IF(IPRINT.GE.1)WRITE(0,"(/'T_EFF=',1P,D12.5)")TEMPE
      IF(IPRINT.GE.0)WRITE(6,"(/'T_EFF=',1P,D12.5)")TEMPE
C
      EMAX=1000*TEMPE
      EMIN=TEMPE/1000
C
      EE=EMIN
      EMAX=LOG(EMAX)
      EMIN=LOG(EMIN)
      DE=(EMAX-EMIN)/(NSTEP-1)
C
      SUMFF=0
      SUMFM=0
      SUMFFE=0
      SUMFME=0
      EE=0  !MAX(0,EE-EXP(DE))
      FF=0
      FM=0
      E=EMIN-DE
      DO N=1,NSTEP
        E0=EE
        FF0=FF
        FM0=FM
        E=E+DE
        EE=EXP(E)
        TSQRTE=SQRT(EE)
        FF=FBAR(EE,TEMPE,ICON,XDUM,XDUM,NENG37,E37,F37)*TSQRTE
        FM=FBAR(EE,TEMPE,IZERO,XDUM,XDUM,NENG37,E37,F37)*TSQRTE
        IF(IPRINT.GE.2)
     X     WRITE(77,'(I5,1P,3D12.3)')N,EE,FM/TSQRTE,FF/TSQRTE
        H=EE-E0
        SUMFF=SUMFF+H*(FF+FF0)
        SUMFM=SUMFM+H*(FM+FM0)
        SUMFFE=SUMFFE+H*(FF*EE+FF0*E0)
        SUMFME=SUMFME+H*(FM*EE+FM0*E0)
      ENDDO
C
      SUMFF=SUMFF/2
      SUMFM=SUMFM/2
C
      IF(IPRINT.GE.1)
     XWRITE(0,"('MAX NORM=',1P,D12.5,3X,'DIST NORM=',D12.5)")SUMFM,SUMFF
      IF(IPRINT.GE.0)
     XWRITE(6,"('MAX NORM=',1P,D12.5,3X,'DIST NORM=',D12.5)")SUMFM,SUMFF
C
      SUMFFE=2*SUMFFE/6
      SUMFME=2*SUMFME/6
C
      IF(IPRINT.GE.1)
     XWRITE(0,"('MAXWELL TE  =',1P,D12.5,3X,'DISTRIBN TE  =',D12.5)")
     XSUMFME,SUMFFE                                             !*conryk
      IF(IPRINT.GE.0)
     XWRITE(6,"('MAXWELL TE  =',1P,D12.5,3X,'DISTRIBN TE  =',D12.5)")
     XSUMFME,SUMFFE                                             !*conryk
C
      TEMPE37=SUMFFE
C
      RAT=SUMFFE/TEMPE
      IF((ABS(RAT-DONE).GT.TMEAN.OR.ITER.EQ.1).AND.ITER.LT.ITMAX)GO TO 2
C
      IF(ABS(SUMFF-DONE).GT.TNORM)THEN
        IF(IPRINT.GE.0)WRITE(0,1000)J,SUMFF
        IF(IPRINT.GE.-1)WRITE(6,1000)J,SUMFF
        DO N=1,NENG37
          F37(N)=F37(N)/SUMFF
        ENDDO
      ENDIF
C
      RETURN
C
 1000 FORMAT(/'*** ATTENTION: RENORMALIZING ADF37 DISTRIBUTION FOR '
     X      ,'TEMPERATURE INDEX',I4,' BY FACTOR',F6.3)
C
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION GAMA(X)
C
C-----------------------------------------------------------------------
C
C  FN.GAMA EVALUATES THE GAMMA FUNCTION OF ARGUMENT X,
C  USING LANCZOS' FORMULA.
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (NMAX=6)
      PARAMETER (ARGFCT=120.0D0)                    !=D.P., 57. FOR S.P.
C
      PARAMETER (DZERO=0.0D0)
      PARAMETER (DONE=1.0D0)
      PARAMETER (DHALF=0.5D0)
C
      DIMENSION COEF(6)
C
      DATA COEF,STP,SER0 /76.18009172947146D0,-86.50532032941677D0,
     &                    24.01409824083091D0, -1.231739572450155D0,
     &                    .1208650973866179D-2, -.5395239384953D-5,
     &                    2.5066282746310005D0, 1.000000000190015D0/
C
      IF(X.GT.ARGFCT)THEN
        WRITE(6,*)'***ERROR IN FN.GAMA: ARGUMENT TOO LARGE (OVERFLOW) ='
     X            ,X
        STOP '***ERROR IN FN.GAMA: ARGUMENT TOO LARGE (OVERFLOW)'
      ENDIF
C
      IF(X.LE.DZERO)THEN
        WRITE(6,*)'***ERROR IN FN.GAMA: ARGUMENT IS NOT POSITIVE =',X
        STOP '***ERROR IN FN.GAMA: ARGUMENT  IS NOT POSITIVE'
      ENDIF
C
      G=NMAX-1
      TE=X+DHALF
      TMP=G+TE
      TMP=EXP(-TMP)*TMP**TE
C
      SER=SER0
      Y=X
C
      DO J=1,NMAX
        Y=Y+DONE
        SER=SER+COEF(J)/Y
      ENDDO
C
      GAMA=TMP*STP*SER/X
C
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION ENINT(N,BHALF,X)
C
C-----------------------------------------------------------------------
C
C N R BADNELL, ADAPTED FROM NUMERICAL RECIPES AND ABRAMOWITZ & STEGUN
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      LOGICAL BHALF
C
      PARAMETER (DZERO=0.0D0)
      PARAMETER (DONE=1.0D0)
      PARAMETER (DTWO=2.0D0)
      PARAMETER (BIG=1.D30)
      PARAMETER (TOLE=1.D-5)
      PARAMETER (MAXIT=100)
      PARAMETER (EULER=.5772156649D0)
C
C GENERATE INTEGER (BHALF.EQ.FALSE.) EXPONENTIAL INTEGRAL ORDER N
C GENERATE HALF-INTEGER (BHALF.EQ.TRUE.) EXPONENTIAL INTEGRAL ORDER N+.5
C
      NM=N-1
      IF(N.LT.0.OR.X.LT.DZERO.OR.X.EQ.DZERO.AND.(N.EQ.0.OR.N.EQ.1))THEN
        WRITE(6,*)'ILLEGAL ARGUMENT IN ENINT: N,X=',N,X
        STOP 'ERROR: ENINT ILLEGAL ARGUMENT'
      ENDIF
C
      H=DZERO
      IF(BHALF)H=DONE/DTWO
C
      IF(N.EQ.0.AND..NOT.BHALF)THEN
        EXPINT=EXP(-X)/X
      ELSEIF(X.EQ.DZERO)THEN
        EXPINT=DONE/(NM+H)
      ELSEIF(X.GT.DONE)THEN              !CONTINUED FRACTION
        BX=X+N+H
        CX=BIG
        DX=DONE/BX
        HX=DX
        DO I=1,MAXIT
          AX=-I*(NM+H+I)
          BX=BX+DTWO
          DX=DONE/(AX*DX+BX)
          CX=BX+AX/CX
          DEL=CX*DX
          HX=HX*DEL
          IF(ABS(DEL-DONE).LT.TOLE)THEN
             EXPINT=HX*EXP(-X)
             GO TO 100
          ENDIF
        ENDDO
c        write(6,*)'ENINT continued fraction not converged'
        EXPINT=-DONE
      ELSE                               !SERIES
        FACT=DONE
        IF(BHALF)THEN
          PI=ACOS(-DONE)
          FACH=-SQRT(PI*X)/(N-H)
        ENDIF
        IF(NM.LE.0)THEN
          IF(BHALF)THEN
            IF(NM.EQ.0)THEN
              EXPINT=FACH+DTWO
            ELSE
              EXPINT=SQRT(PI/X)-DTWO
            ENDIF
          ELSE
            EXPINT=-LOG(X)-EULER
          ENDIF
        ELSE
          EXPINT=DONE/(NM+H)
        ENDIF
        DO I=1,MAXIT
          FACT=-FACT*X/I
          IF(BHALF)FACH=-FACH*X/(I-H)
          IF(I.NE.NM)THEN
            DEL=-FACT/(I-NM-H)
          ELSE
            IF(BHALF)THEN
              DEL=FACH+DTWO*FACT
            ELSE              !DIGAMMA
              PSI=-EULER
              DO II=1,NM
                PSI=PSI+DONE/II
              ENDDO
              DEL=FACT*(-LOG(X)+PSI)
            ENDIF
          ENDIF
          EXPINT=EXPINT+DEL
          IF(ABS(DEL).LT.ABS(EXPINT)*TOLE)GO TO 100
        ENDDO
c        write(6,*)'ENINT series not converged'
        EXPINT=-DONE
      ENDIF
C
  100 ENINT=EXPINT
C
      END
C
C***********************************************************************
C
      SUBROUTINE MIDEXP(N,B,S)
C
C-----------------------------------------------------------------------
C
CNRB
C ADAPTATION OF NUMERICAL RECIPES ROUTINE,
C FOR APPLICATION TO MAXWELLIAN INTEGRATION AT HIGH TEMP.
C THEN B=1/KT.
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      IF(N.EQ.1)THEN
        X=0.5D0
        T=-LOG(X)/B
        S=GLNAG(B,T)
        S=S/B
      ELSE
        IT=3**(N-2)
        TNM=IT
        DEL=1/(3*TNM)
        DDEL=DEL+DEL
        X=DEL/2
        T=-LOG(X)/B
        SUM=0
        DO J=1,IT
          SUM=SUM+GLNAG(B,T)
          X=X+DDEL
          T=-LOG(X)/B
          SUM=SUM+GLNAG(B,T)
          X=X+DEL
          T=-LOG(X)/B
        ENDDO
        SUM=SUM/B
        S=(S+SUM/TNM)/3
      ENDIF
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION GAUSSQ(NQ,W,X,B)
C
C-----------------------------------------------------------------------
C
C  (C) Copr. 1986-92 Numerical Recipes Software .37.
C
C GAUSSIAN QUADRATURE DRIVER
C
C NQ: NO OF QUADRATURE POINTS.
C W: WEIGHTS
C X: DZEROES
C LAGUERRE CASE: B DEFINED BY EXP(-E*B), I.E. B=1/KT FOR MAXWELLIAN.
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      PARAMETER (DZERO=0.0D0)
      DIMENSION X(NQ), W(NQ)
C
      CALL GAULAG(X,W,NQ,DZERO)
C
      SUM=DZERO
      DO N=1,NQ
      T=X(N)/B
      SUM=SUM+W(N)*GLNAG(B,T)
      ENDDO
C
      GAUSSQ=SUM/B
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE gaulag(x,w,n,alf)
C
C-----------------------------------------------------------------------
C
C  (C) Copr. 1986-92 Numerical Recipes Software .37.
C MODIFIED NRB 26/09/01 FOR REAL*8 AND OPT CASE ALF=0.
C The x(i) are in ascending order.
C
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
      INTEGER n,MAXIT
      real*8 alf,w(n),x(n)
      REAL*8 EPS
      PARAMETER (EPS=3.D-12,MAXIT=20)
CU    USES gammln
      INTEGER i,its,j
      real*8 ai,gammln,TN
      REAL*8 p1,p2,p3,pp,z,z1
C
      do 13 i=1,n
        if(i.eq.1)then
          z=(1.d0+alf)*(3.d0+.92d0*alf)/(1.d0+2.4d0*n+1.8d0*alf)
        else if(i.eq.2)then
          z=z+(15.d0+6.25d0*alf)/(1.d0+.9d0*alf+2.5d0*n)
        else
          ai=i-2
          z=z+((1.d0+2.55d0*ai)/(1.9d0*ai)+1.26d0*ai*alf/
     X(1.d0+3.5d0*ai))*(z-x(i-2))/(1.d0+.3d0*alf)
        endif
        do 12 its=1,MAXIT
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
11        continue
          pp=(n*p1-(n+alf)*p2)/z
          z1=z
          z=z1-p1/pp
          if(abs(z-z1).le.EPS)goto 1
12      continue
        WRITE(6,*) 'WARNING: too many iterations in gaulag',n,z,z1
1       x(i)=z
        TN=N
        w(i)=-1.0d0/(pp*n*p2)
        IF(ALF.NE.0.d0)w(i)=W(I)*exp(gammln(alf+n)-gammln(TN))
13    continue
      return
      END
C
C***********************************************************************
C
      FUNCTION gammln(xx)
C
C-----------------------------------------------------------------------
C
C  (C) Copr. 1986-92 Numerical Recipes Software .37.
C
C-----------------------------------------------------------------------
C
      real*8 gammln,xx
      INTEGER j
      REAL*8 ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     &  24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     &  -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
C
C***********************************************************************
C
       SUBROUTINE DIPOL(JSW,N1,N2,E2,LMAX,CP,CM,JC)
C
C-----------------------------------------------------------------------
C
C  ALAN BURGESS DAMTP CAMBRIDGE, MODS BY NRB.
C
C  SR.DIPOL CALCULATES SQUARES OF HYDROGENIC DIPOLE LENGTH RADIAL MATRIX
C  ELEMENTS FOR BOUND-BOUND OR BOUND-FREE TRANSITIONS.
C
C  BOUND STATES ARE NORMALISED TO UNITY.
C  FREE STATES ARE NORMALISED TO ASYMPTOTIC AMPLITUDE K**(-0.5).
C
C  N.B. DIPOLE ACCELERATION MATRIX ELEMENT = (E12**2/4Z) * DIPOLE LENGTH
C  WHERE E12 = - N1**(-2) + N2**(-2)  FOR BOUND-BOUND
C            = - N1**(-2) + E2        FOR BOUND-FREE
C          Z = REDUCED CHARGE
C  INPUT:
C   FOR BOUND-BOUND,SET JSW=NEGATIVE
C                     N1,N2=PRINCIPAL QUANTUM NUMBERS OF STATES
C                      LMAX=RANGE OF ANGULAR MOMENTUM QUANTUM NUMBERS
C   FOR BOUND-FREE, SET JSW=POSITIVE
C                       N1=BOUND STATE PRINCIPAL QUANTUM NUMBER
C                       E2=FREE STATE ENERGY IN RYDBERGS (=K**2)
C
C  OUTPUT:
C   VECTOR CP(L),L=1,LMAX,CONTAINS SQUARED MATRIX ELEMENTS FOR ANGULAR
C                         MOMENTUM TRANSITIONS FROM L-1 TO L,
C   VECTOR CM(L),L=1,LMAX,CONTAINS SQUARED MATRIX ELEMENTS FOR ANGULAR
C                         MOMENTUM TRANSITIONS FROM L TO L-1,
C               IN BOTH CASES THE TRANSITION IS FROM LOWER TO HIGHER
C               ENERGY, INDEPENDANT OF THE SIGN OF N1-N2 FOR BOUND-BOUND
C               CASES. IF N1=N2 THEN CP(L)=CM(L).
C   VECTOR JC(L),L=1,LMAX WILL USUALLY BE DZERO AND MAY THEN BE IGNORED,
C               BUT FOR EXTREME INPUT VALUES THERE IS POSSIBILITY OF
C               OVER OR UNDERFLOW OF CP(L) OR CM(L),IN WHICH CASE THE
C               OUTPUT VALUES OF CP(L) AND CM(L) SHOULD BE MULTIPLIED
C               BY (1.0D10)**JC(L) TO OBTAIN TRUE VALUES.
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (DZERO=0.0D0)
      PARAMETER (DONE=1.0D0)
C      PARAMETER (PI=3.14159265359D0)
      PARAMETER (S1=1.0D10)
      PARAMETER (S2=1.0D-10)
      PARAMETER (TEST1=1.0D-20)
      PARAMETER (TEST2=1.0D20)
      PARAMETER (TEST3=0.044D0)
      PARAMETER (TEST4=0.1D0)
      PARAMETER (TEST5=300.0D0)
      PARAMETER (TEST6=1.0D-30)
      PARAMETER (TEST7=1.0D30)
C
      DIMENSION CP(LMAX),CM(LMAX),JC(LMAX)
C
      PI=ACOS(-DONE)
C
      N=N1
      E=E2
      IF(JSW.LE.0)THEN
        EN2=N2
        N3=N2
        IF(N1.EQ.N2)GO TO 59
        IF(N2.LT.N1)THEN
          N=N2
          EN2=N1
          N3=N1
        ENDIF
        E=-DONE/(EN2*EN2)
      ENDIF
C
      EN=N
      ENN=EN*EN
      E1=-DONE/ENN
      JMAX=LMAX
      C1=DONE
      C2=DZERO
      JS=0
      L=N+1
      IF(N.LE.LMAX)THEN
        CP(N)=DONE
        CM(N)=DZERO
        JC(N)=0
        JMAX=N-1
        DO I=L,LMAX
          CP(I)=DZERO
          CM(I)=DZERO
          JC(I)=0
        ENDDO
      ENDIF
C
    9 L=L-1
      IF(L.GT.1)THEN
        EL=L
        ELL=EL*EL
        T1=DONE+ELL*E1
        T2=DONE+ELL*E
        T3=L+L-1
        T4=DONE/(T3+DONE)
        T5=(T3*T1*C2+T2*C1)*T4
        C1=(T1*C2+T3*T2*C1)*T4
        C2=T5
   11   IF(C1*C1.GT.TEST2)THEN
          C1=S2*C1
          C2=S2*C2
          JS=JS+1
          GO TO 11
        ENDIF
        IF(L.LE.LMAX+1)THEN
          CP(L-1)=C1
          CM(L-1)=C2
          JC(L-1)=JS
        ENDIF
        GO TO 9
      ENDIF
C
      JS=0
      T=4
      T=DONE/(T*EN*ENN)
      IF(JSW.LE.0)THEN                          !JSW.LT.0
        ENN2=EN2*EN2
        T1=4
        T1=T1*ENN*ENN2/(ENN2-ENN)
        T1=T1*T1
        T=T*T1*T1/(EN2*ENN2)
        IF(N3.LE.30)THEN
          T=T*((EN2-EN)/(EN2+EN))**(N3+N3)
          GO TO 34
        ENDIF
        E21=E/E1
        IF(E21.LE.TEST4)THEN
          T2=DZERO
          DO J=1,11
            T3=2*(11-J)+1
            T2=DONE/T3+T2*E21
          ENDDO
          T2=T2+T2
        ELSE
          T3=EN/EN2
          T2=LOG((DONE+T3)/(DONE-T3))/T3
        ENDIF
        T2=T2+T2
        T1=T1*EXP(-T2)
C
      ELSE                                      !JSW.GT.0
C
        T1=4
        T1=T1*ENN/(DONE+ENN*E)
        T1=T1*T1
        T=T*T1*T1
        IF(E.LT.TEST3)THEN
          T3=2
          T=T*(PI/T3)
        ELSE
          T4=SQRT(E)
          IF(T4.LE.TEST5)THEN
            T3=(PI+PI)/T4
            T3=DONE-EXP(-T3)
            T3=DONE/T3
          ELSE
            T4=PI/T4
            T3=3
            T3=(DONE+T4+T4*T4/T3)/(T4+T4)
          ENDIF
          T2=2
          T=T*(PI*T3/T2)
        ENDIF
C
        T4=ENN*E
        IF(T4.LE.TEST4)THEN
          T2=DZERO
          DO J=1,11
            T3=2*(11-J)+1
            T2=DONE/T3-T2*T4
          ENDDO
        ELSE
          T3=SQRT(T4)
          T2=ATAN(T3)/T3
        ENDIF
        T2=T2+T2
        T2=T2+T2
        T1=T1*EXP(-T2)
      ENDIF
C                                               !ALL JSW
   34 DO J=1,N
        TJ=J+J
        T2=TJ*(TJ-DONE)
        T2=T2*T2
        T=T*T1/T2
   35   IF(T.LE.TEST1)THEN
          T=T*S1
          JS=JS-1
          GO TO 35
        ENDIF
   37   IF(T.GE.TEST2)THEN
          T=T*S2
          JS=JS+1
          GO TO 37
        ENDIF
      ENDDO
      J=0
C
   40 J=J+1
      IF(J.LE.JMAX)THEN
        TJ=J
        TJ=TJ*TJ
        T1=DONE+TJ*E1
        T2=DONE+TJ*E
        T3=CP(J)
        T3=T2*T*T3*T3
        T4=CM(J)
        T4=T1*T*T4*T4
        L1=JC(J)+JC(J)+JS
C
   42   IF(L1.LT.0)THEN
          IF(T4.GT.TEST6)THEN
            L1=L1+1
            T3=T3*S2
            T4=T4*S2
            GO TO 42
          ENDIF
        ELSEIF(L1.GT.0)THEN
          IF(T3.LT.TEST7)THEN
            L1=L1-1
            T3=T3*S1
            T4=T4*S1
            GO TO 42
          ENDIF
        ENDIF
C
        CP(J)=T3
        CM(J)=T4
        JC(J)=L1
        T=T*T1*T2
   48   IF(T.GT.TEST2)THEN
          T=T*S2
          JS=JS+1
          GO TO 48
        ENDIF
        GO TO 40
      ENDIF
C
      IF(N.LE.LMAX)THEN
        T2=DONE+ENN*E
        T3=CP(N)
        T3=T2*T*T3*T3
        L1=JC(N)+JC(N)+JS
C
   52   IF(L1.LT.0)THEN
          IF(T3.GT.TEST6)THEN
            L1=L1+1
            T3=T3*S2
            GO TO 52
          ENDIF
        ELSEIF(L1.GT.0)THEN
          IF(T3.LT.TEST7)THEN
            L1=L1-1
            T3=T3*S1
            GO TO 52
          ENDIF
        ENDIF
C
        CP(N)=T3
        JC(N)=L1
      ENDIF
C
      RETURN
C
   59 JMAX=LMAX
      IF(N.LE.LMAX)THEN
        DO L=N,LMAX
          CP(L)=DZERO
          CM(L)=DZERO
          JC(L)=0
        ENDDO
        JMAX=N-1
      ENDIF
      T1=9
      T2=4
      T3=(T1/T2)
      T1=EN2*EN2
      T2=T1*T3
      DO J=1,JMAX
        TJ=J
        JC(J)=0
        T=T2*(T1-TJ*TJ)
        CP(J)=T
        CM(J)=T
      ENDDO
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE HPSRTI (N,A,IP)
C
C-----------------------------------------------------------------------
C
C IMPLICIT HEAPSORT  BY *MAGNITUDE* OF
C INPUT:  VECTOR A, LENGTH N.
C OUTPUT: DOWN-ORDERED POINTER IN IP, A IS UNCHANGED.
C        (UP-ORDERED CAN BE OBTAINED BY CHANGING .LT. TO .GT. AS BELOW).
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(*),IP(*)
C
      DO I=1,N
        IP(I)=I
      ENDDO
C
      IF(N.LT.2)RETURN
C
      L=N/2+1
      IT=N
C
  1   IF(L.GT.1)THEN
        L=L-1
        IPT=IP(L)
      ELSE
        IPT=IP(IT)
        IP(IT)=IP(1)
        IT=IT-1
        IF(IT.EQ.1)THEN
          IP(1)=IPT
          RETURN
        ENDIF
      ENDIF
      I=L
      J=L+L
C
  2   IF(J.LE.IT)THEN
        IF(J.LT.IT)THEN
          IF(abs(A(IP(J+1))).lT.abs(A(IP(J))))J=J+1  !.lt. down, .gt .up
        ENDIF
        IF(abs(A(IP(J))).lT.abs(A(IPT)))THEN         !.lt. down, .gt .up
          IP(I)=IP(J)
          I=J
          J=J+J
        ELSE
          J=IT+1
        ENDIF
        GO TO 2
      ENDIF
      IP(I)=IPT
      GO TO 1
C
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION QDT(QD,NZ0,NE,N,L,KAPPA)
C
C-----------------------------------------------------------------------
C
C NRB:
C EVALUATES ONE-ELECTRON ENERGY WITH NON-DZERO QUANTUM DEFECT
C
C   : QD0, UNIVERSAL QUANTUM DEFECT GIVEN BY
C         QD0*(NE**1.67-1)/(Z0**.67*Z**.33*(1+L**3))
C         CURRENT VALUE IN FUNCTION QDT IS QD0=0.182
C
C KAPPA= 0 NON-RELATIVISTIC
C      =-1 KAPPA-AVERAGE RELATIVISTIC
C      = L RELATIVISTIC FOR J=L-0.5
C      =-L-1 RELATIVISTIC FOR J=L+0.5
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (DZERO=0.0D0)
      PARAMETER (DONE=1.0D0)
      PARAMETER (DTWO=2.0D0)
      PARAMETER (QD0=0.182D0)
      PARAMETER (ALF2=5.325D-5)
C
      COMMON /QDTS/QDTS(0:30),NQDT
C
      IF(N.LE.0)THEN
        QD=DZERO
        QDT=DZERO
        RETURN
      ENDIF
C
      TZ0=NZ0
      NZ=NZ0-NE+1
      TZ=NZ
      IF(L.LT.0.OR.NE.LE.1)THEN
        QD=DZERO
      ELSE
        IF(NQDT.GT.L)THEN
          QD=QDTS(L)
        ELSE
          TL=L**3+1
          TE=NE
          QD=QD0*(TE**1.667D0-DONE)/(TZ0**0.667D0*TZ**0.333D0*TL)
        ENDIF
      ENDIF
      TN=N
      T3=TN*TN*TN
      TN=TN-QD
C
      QDT=-(TZ/TN)**2
      IF(KAPPA.EQ.0)RETURN
C
      IF(KAPPA.EQ.-1)THEN
        ESO=DZERO
      ELSE
        IF(L.EQ.0)THEN
          WRITE(6,*)'*** FN.QDT ERROR: L=0 FOR KAPPA.NE.-1'
          ESO=DZERO
        ELSE
          ESO=ALF2*NZ**4*KAPPA/(T3*L*(L+1)*(2*L+1))
        ENDIF
      ENDIF
C
      IF(L.EQ.0)THEN
        ED=ALF2*NZ**4/T3
      ELSE
        ED=DZERO
      ENDIF
C
      EM=-ALF2*NZ**4*(4*N/(L+DONE/DTWO)-3)/(4*N*T3)
C
      QDT=QDT+ESO+ED+EM
C
      RETURN
      END
C
C***********************************************************************
C
       SUBROUTINE SPLYN(N,X,Y,I1,E1,I2,E2,A,B,C,D,S)
C
C-----------------------------------------------------------------------
C
C ALAN BURGESS, D.A.M.T.P. CAMBRIDGE.
C CUBIC SPLINE FITTING TO THE DATA POINTS (X(J),Y(J)),J=1,2...N
C WITH KNOTS AT X(J),J=2,3...(N-1),
C IN THE FORM Y(X)=A(J)+Z*(B(J)+Z*(C(J)+Z*D(J))),
C FOR X IN THE RANGE (X(J),X(J+1)),
C WHERE Z=X-(X(J)+X(J+1))/2.
C ONE OF THE FOLLOWING END CONDITIONS MUST BE CHOSEN FOR EACH END:
C (1)SPECIFIED END FIRST DERIVATIVES; SET I1=1, E1=(DY/DX)(X=X(1))
C                                         I2=1, E2=(DY/DX)(X=X(N))
C (2)SPECIFIED END 2ND DERIVATIVES; SET I1=2, E1=((D/DX)**2)Y (X=X(1))
C                                       I2=2, E2=((D/DX)**2)Y (X=X(N))
C (3)END 2ND DERIVATIVE =NEXT-TO-END 2ND DERIVATIVE; SET I1=3, I2=3
C                                            (NO NEED TO SET E1,E2)
C (4)3RD DERIVATIVE CONTINUOUS AT FIRST AND LAST KNOTS; SET I1=4, I2=4
C                                              (NO NEED TO SET E1,E2).
C N.B. THE CHOSEN CONDITIONS NEED NOT BE THE SAME FOR THE TWO ENDS.
C  INPUT: N,X(J),Y(J) (J=1,2...N),I1,E1,I2,E2.
C  OUTPUT: A(J),B(J),C(J),D(J) (J=1,2...(N-1))
C          S(J) (J=1,2...N), THE SECOND DERIVATIVE OF Y.
C
C-----------------------------------------------------------------------
C
       IMPLICIT REAL*8 (A-H,O-Z)
C
       DIMENSION A(N),B(N),C(N),D(N),X(N),Y(N),S(N)
C
       H1=X(2)-X(1)
       T1=(Y(2)-Y(1))/H1
C
C       GO TO (1,2,3,4),I1
      IF(I1.EQ.1)THEN
       B(1)=H1+H1
       C(1)=H1
       D(1)=6.0 *(T1-E1)
      ELSEIF(I1.EQ.2)THEN
       B(1)=1.0
       C(1)=0.0
       D(1)=E1
      ELSEIF(I1.EQ.3)THEN
       B(1)=1.0
       C(1)=-1.0
       D(1)=0.0
      ELSEIF(I1.EQ.4)THEN
       B(1)=1.0
       C(1)=0.0
       D(1)=0.0
      ELSE
        WRITE(MW6,*)'SR.SPLYN: ILLEGAL VALUE FOR I1:',I1
        STOP 'SR.SPLYN: ILLEGAL VALUE FOR I1'
      ENDIF
C
C    5  GO TO (6,7,8,9),I2
      IF(I2.EQ.1)THEN
       H2=X(N)-X(N-1)
       A(N)=H2
       B(N)=H2+H2
       D(N)=6.0 *(E2-(Y(N)-Y(N-1))/H2)
      ELSEIF(I2.EQ.2)THEN
       A(N)=0.0
       B(N)=1.0
       D(N)=E2
      ELSEIF(I2.EQ.3)THEN
       A(N)=-1.0
       B(N)=1.0
       D(N)=0.0
      ELSEIF(I2.EQ.4)THEN
       A(N)=0.0
       B(N)=1.0
       D(N)=0.0
      ELSE
        WRITE(MW6,*)'SR.SPLYN: ILLEGAL VALUE FOR I2:',I2
        STOP 'SR.SPLYN: ILLEGAL VALUE FOR I2'
      ENDIF
C
   10  N1=N-1
       DO 11 J=2,N1
       H2=X(J+1)-X(J)
       T2=(Y(J+1)-Y(J))/H2
       A(J)=H1
       B(J)=2.0 *(H1+H2)
       C(J)=H2
       D(J)=6.0 *(T2-T1)
       H1=H2
       T1=T2
   11  CONTINUE
C
C       IF(I1-4)13,12,13
      IF(I1.EQ.4)THEN
       A(2)=0.0
       H1=X(2)-X(1)
       H2=X(3)-X(2)
       B(2)=B(2)+H1+H1*H1/H2
       C(2)=C(2)-H1*H1/H2
      ENDIF
C
C   13  IF(I2-4)15,14,15
      IF(I2.EQ.4)THEN
       C(N-1)=0.0
       T1=X(N-1)-X(N-2)
       T2=X(N)-X(N-1)
       A(N-1)=A(N-1)-T2*T2/T1
       B(N-1)=B(N-1)+T2+T2*T2/T1
      ENDIF
C
       CALL TRIMAT(N,A,B,C,D,S)
C
C       IF(I1-4)17,16,17
       IF(I1.EQ.4)S(1)=((H1+H2)*S(2)-H1*S(3))/H2
C   17  IF(I2-4)19,18,19
       IF(I2.EQ.4)S(N)=((T1+T2)*S(N-1)-T2*S(N-2))/T1
C
       DO 20 J=1,N1
       H1=X(J+1)-X(J)
       T1=0.25 *H1*H1
       D(J)=(S(J+1)-S(J))/(6.0 *H1)
       C(J)=0.25 *(S(J+1)+S(J))
       B(J)=(Y(J+1)-Y(J))/H1-T1*D(J)
       A(J)=0.5 *(Y(J+1)+Y(J))-T1*C(J)
   20  CONTINUE
C
       RETURN
      END
C
C***********************************************************************
C
       REAL*8 FUNCTION SPVAL(N,X,Y,A,B,C,D,X1,J1)
C
C-----------------------------------------------------------------------
C
C ALAN BURGESS, D.A.M.T.P. CAMBRIDGE.
C EVALUATES  CUBIC SPLINE FIT TO Y(X) AT X=X1, WHERE X1 LIES IN THE
C INTERVAL (X(J1),X(J1+1)), GIVEN THE SPLINE COEFFICIENTS
C A(J),B(J),C(J),D(J) PRODUCED BY SUBROUTINE SPLYN.
C J1 NEED NOT BE SET AS INPUT, BUT EXECUTION MAY BE QUICKER IF IT IS.
C  INPUT: X(J), Y(J), J=1,2...N  : NOTE, Y(J) NOT NEEDED. (NRB)
C         A(J),B(J),C(J),D(J), J=1,2...N-1
C         X1
C         J1 (OPTIONAL).
C  OUTPUT: SPVAL=Y(X1).
C
C-----------------------------------------------------------------------
C
       IMPLICIT REAL*8 (A-H,O-Z)
C
       DIMENSION X(N),Y(N),A(N),B(N),C(N),D(N)
C
      IF (X(1).LT.X(N))THEN
       I1=1
       I2=N
       I=1
      ELSE
       I1=N
       I2=1
       I=-1
      ENDIF
C
       IF((J1-I1)*(J1-I2+I).GT.0)J1=I1
C
C    2  IF(X1-X(J1))3,9,6
    2 IF(X1.LT.X(J1))THEN
       IF (J1.EQ.I1)GO TO 9
       J1=J1-I
       GO TO 2
      ENDIF
C
    6 IF(X1.GT.X(J1+I))THEN
       IF (J1.EQ.(I2-I))GO TO 9
       J1=J1+I
       GO TO 6
      ENDIF
C
    9  J2=J1+I
       Z=X1-0.5*(X(J2)+X(J1))
       IF(J2.GT.J1)J2=J1
C
       SPVAL=A(J2)+Z*(B(J2)+Z*(C(J2)+Z*D(J2)))
C
       RETURN
      END
C
C***********************************************************************
C
       SUBROUTINE TRIMAT(N,A,B,C,D,Y)
C
C-----------------------------------------------------------------------
C
C ALAN BURGESS, D.A.M.T.P. CAMBRIDGE.
C SOLUTION OF TRI-DIAGONAL MATRIX EQUATION BY FORWARD AND BACKWARD PASS
C (SEE 'MODERN COMPUTING METHODS', PAGES 97,98).
C  INPUT: N, (THE MATRIX IS N BY N)
C         A(J), J=2,3...N (SUB-DIAGONAL ELEMENTS)
C         B(J), J=1,2...N (DIAGONAL ELEMENTS)
C         C(J), J=1,2...N-1 (SUPER-DIAGONAL ELEMENTS)
C         D(J), J=1,2...N (RIGHT-HAND SIDE).
C  OUTPUT: N,A,B AND C ARE PRESERVED, D IS OVERWRITTEN
C          Y(J), J=1,2...N CONTAINS THE SOLUTION VECTOR.
C
C-----------------------------------------------------------------------
C
       IMPLICIT REAL*8 (A-H,O-Z)
C
       DIMENSION A(N),B(N),C(N),D(N),Y(N)
C
       Y(1)=B(1)
       DO 1 J=2,N
       T=A(J)/Y(J-1)
       Y(J)=B(J)-T*C(J-1)
       D(J)=D(J)-T*D(J-1)
    1  CONTINUE
       Y(N)=D(N)/Y(N)
       DO 2 I=2,N
       J=N-I+1
       Y(J)=(D(J)-C(J)*Y(J+1))/Y(J)
    2  CONTINUE
       RETURN
      END
C
C***********************************************************************
C
C-----------------------------------------------------------------------
C
C A COLLECTION OF INTEGER ORDER MODIFIED BESSEL FUNCTION ROUTINES.
C
C-----------------------------------------------------------------------
C
C***********************************************************************
C
      REAL*8 FUNCTION BESSK(N,X)
      IMPLICIT REAL*8(A-H,O-Z)
C
C-----------------------------------------------------------------------
C
C MODIFIED BESSEL FUNCTION K(N,X) GENERATED BY UP RECURRENCE RELATION:
C
C K(N+1,X)=(2N/X)*K(N,X)+K(N-1,X)                          (A&S: 9.6.26)
C
C-----------------------------------------------------------------------
C
      IF(N.LT.2)STOP 'BAD ARGUMENT N IN BESSK'
C
      TOX=2.0D0/X
      BKM=BESSK0(X)
      BK=BESSK1(X)
      DO J=1,N-1
        BKP=BKM+J*TOX*BK
        BKM=BK
        BK=BKP
      ENDDO
      BESSK=BK
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION BESSK0(X)
C
C-----------------------------------------------------------------------
C
C MODIFIED BESSEL FUNCTION K(0,X) GENERATED FROM ABRAMOWITZ & STEGUN
C POLYNOMIAL FITS 9.8.5 AND 9.8.6
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0,
     X                      0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1,
     X                -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
C
      IF (X.LE.2.0D0) THEN
        Y=X*X/4.0D0
        BESSK0=(-LOG(X/2.0D0)*BESSI0(X))+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*
     X         (P6+Y*P7))))))
      ELSE
        Y=(2.0D0/X)
        BESSK0=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*
     X          Q7))))))
      ENDIF
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION BESSK1(X)
      IMPLICIT REAL*8(A-H,O-Z)
C
C-----------------------------------------------------------------------
C
C MODIFIED BESSEL FUNCTION K(1,X) GENERATED FROM ABRAMOWITZ & STEGUN
C POLYNOMIAL FITS 9.8.7 AND 9.8.8
C
C-----------------------------------------------------------------------
C
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,0.15443144D0,-0.67278579D0,
     X       -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1,
     X                0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
C
      IF (X.LE.2.0D0) THEN
        Y=X*X/4.0D0
        BESSK1=(LOG(X/2.0D0)*BESSI1(X))+(1.0D0/X)*(P1+Y*(P2+Y*(P3+Y*
     X         (P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=2.0D0/X
        BESSK1=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*
     X          Q7))))))
      ENDIF
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION BESSI(N,X)
      IMPLICIT REAL*8(A-H,O-Z)
C
C-----------------------------------------------------------------------
C
C MODIFIED BESSEL FUNCTION I(N,X) GENERATED BY DOWN RECURRENCE RELATION:
C
C I(N+1,X)=-(2N/X)*I(N,X)+I(N-1,X)                         (A&S: 9.6.26)
C
C-----------------------------------------------------------------------
C
      PARAMETER (BIGNO=1.0D10)
      PARAMETER (BIGNI=1.0D-10)
C
      PARAMETER (IACC=40)
C
      IF (N.LT.2)STOP 'BAD ARGUMENT N IN BESSI'
C
      IF (X.EQ.0.0D0) THEN
        BESSI=0.0D0
      ELSE
        TOX=2.0D0/ABS(X)
        BIP=0.0D0
        BI=1.0D0
        BESSI=0.0D0
        M=2*((N+INT(SQRT(DFLOAT(IACC*N)))))
        DO J=M,1,-1
          BIM=BIP+DFLOAT(J)*TOX*BI
          BIP=BI
          BI=BIM
          IF (ABS(BI).GT.BIGNO) THEN
            BESSI=BESSI*BIGNI
            BI=BI*BIGNI
            BIP=BIP*BIGNI
          ENDIF
          IF (J.EQ.N) BESSI=BIP
        ENDDO
        BESSI=BESSI*BESSI0(X)/BI
        IF (X.LT.0.0D0.AND.MOD(N,2).EQ.1) BESSI=-BESSI
      ENDIF
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION BESSI0(X)
      IMPLICIT REAL*8(A-H,O-Z)
C
C-----------------------------------------------------------------------
C
C MODIFIED BESSEL FUNCTION I(0,X) GENERATED FROM ABRAMOWITZ & STEGUN
C POLYNOMIAL FITS 9.8.1 AND 9.8.2
C
C-----------------------------------------------------------------------
C
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,
     X        1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     X       0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     X      0.2635537D-1,-0.1647633D-1,0.392377D-2/
C
      IF (ABS(X).LT.3.75D0)THEN
        Y=(X/3.75D0)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=ABS(X)
        Y=3.75D0/AX
        BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*
     X         (Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
C
      RETURN
      END
C
C***********************************************************************
C
      REAL*8 FUNCTION BESSI1(X)
      IMPLICIT REAL*8(A-H,O-Z)
C
C-----------------------------------------------------------------------
C
C MODIFIED BESSEL FUNCTION I(1,X) GENERATED FROM ABRAMOWITZ & STEGUN
C POLYNOMIAL FITS 9.8.3 AND 9.8.4
C
C-----------------------------------------------------------------------
C
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     X         0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     X       -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     X                   -0.2895312D-1,0.1787654D-1,-0.420059D-2/
C
      IF (ABS(X).LT.3.75D0) THEN
        Y=(X/3.75D0)**2
        BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=ABS(X)
        Y=3.75D0/AX
        BESSI1=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*
     X         (Q7+Y*(Q8+Y*Q9))))))))
        IF(X.LT.0.0D0)BESSI1=-BESSI1
      ENDIF
C
      RETURN
      END
