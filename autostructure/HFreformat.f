!----------------------------------------------------------------------
!
!          w-file   --->  orbital.inp    (for RMATRIX1, AUTOSTRUCTURE)
!
!----------------------------------------------------------------------
!
!     using interpolation Lagrange formular of order 5
!
!     and dirrect determination of derivaties from finite-difference
!     formular
!---------------------------------------------------------------------- 


      Implicit real *8 (A-H,O-Z)

      Parameter (mr=220,mfun=99)

      Character AF*12,EL*4,AT*6,TR*6

      COMMON
     */RFF/ H,EH,RHO,Z,NO,ND,R(mr),RR(mr),R2(mr),X(mr)
     */FUN/ P(mr,mfun),YK(mr),N(mfun),L(mfun),KN(mfun),MX(mfun),
     *      AZ(mfun),EI(mfun),ZETA(mfun),EL(mfun),NWF

      PARAMETER (MZLR2=  20)
      PARAMETER (MZNPT= 800)
      PARAMETER (MXNIX=MZNPT/16)
      COMMON /MESHS/ RX(MZNPT),HINT,IHX(MXNIX),IRX(MXNIX),NIX
      Dimension PX(MZNPT),QX(MZNPT),PY(MZNPT),QY(MZNPT)

      Data ONE/1.0/,ZERO/0.0/


!-----------------------------------------------------------------------
!                                                           read w-file:

c     write(*,'('' Enter the name of w-file: '')')

c     Read(*,*) AF ! only works for f90
c     Read(*,'(a)') AF ! works also for f77
cSGZ  Read(*,*) AF
      DATA AF/'wfn.out'/
      Call Read_fun(1,AF,0,0,AT,TR)

      AF='orbital.inp'
      Open(1,file=AF,status='unknown')
cSGZ  Open(1,file=AF)

!-----------------------------------------------------------------------
!                                                         R-matrix mesh:
      write(*,'('' Enter the bound radius ....: '')')
      Read(*,*) RA

      write(*,'('' Enter NRANG2 (usually 31) .: '')')
      Read(*,*) NRANG2

      write(*,'('' Enter LRANG2 (usually 31) .:'')')
      Read(*,*) LRANG2

      Call MESH(RA,Z,NRANG2,LRANG2)

      NTOT=IRX(NIX)
      NC=NTOT+1

!----------------------------------------------------------------------
!                                           first line of orbital.inp:

      write(*,'('' Enter the number of electrons: '')')
      Read(*,*) Nelc
      NZ=Z+0.1
      mps = 0                  ! ???
      write(1,3020) -9,nwf,mps,RX(NTOT),NC,Nelc,NZ,AT,TR
      write(*,3020) -9,nwf,mps,RX(NTOT),NC,Nelc,NZ,AT,TR
 3020 FORMAT(3i5,F8.2,i5,i4,i4,33x,2A6)

!----------------------------------------------------------------------
!                                         radial mesh and potensial(?): 
!                                         see eq.80 in Berr.et.al(1995)
      zzz=Z**(ONE/3)
      Do i=1,NTOT
       PX(i)= (EXP(-ZZZ*RX(i))-ONE)*(Nelc-1)+Z
      End do

      write(1,3000) -8,1,0.0,Z,2,RX(1),PX(1),AT,TR
      Do i=2,NTOT-2,2
       j=i+1
       write(1,3000) -8,i+1,RX(i),PX(i),j+1,RX(j),PX(j),AT,TR
      End do
      i=NC
      write(1,3000) -8,i,RX(i-1),PX(i-1),1,0.0,Z,AT,TR
 3000 FORMAT(i5,2(I4,2E14.7),2A6)

!----------------------------------------------------------------------
!                                                      radial orbitals:
      Do k=1,NWF

       CALL DIFF(k)
       QZ=-YK(4)/R2(4)/R(4)**(L(k)+2)
       Do I=1,MX(k)
        PX(I)= P(i,k)*R2(I)
        QX(I)=-YK(I)/R2(I)/R(I)
       End do
       Call LAGRN(5,MX(k),NTOT,R,RX,PX,PY)
       Call LAGRN(5,MX(k),NTOT,R,RX,QX,QY)

       RM=R(MX(K)-2)
       Do i=1,NTOT
        if(RX(i).gt.RM) PY(i)=0.0
        if(RX(i).gt.RM) QY(i)=0.0
       End do

       NOC=NTOT/2+1
       X1=RX(ntot-1)
       X2=RX(ntot-2)
       P1=PY(ntot-1)
       P2=PY(ntot-2)
       Q1=QY(ntot-1)
       Q2=QY(ntot-2)
       EPS=ZERO
       if(N(k)*N(k)*ABS(Q1).lt.Z*Z*ABS(P1))
     : EPS = (X1*Q1/P1-X2*Q2/P2)/(X1-X2)
       write(1,3010) -7,k,N(k),L(k),Z,NELC,0.0,NOC,EPS,AT,TR,EL(k)
 3010 FORMAT(2I5,2X,2I3,f4.0,i3,F12.6,i6,F12.6,1x,2a6,4x,a4)


      write(1,3001) -6,1,AZ(k),QZ,2,PY(1),QY(1),EL(k)
      Do i=2,NTOT-2,2
       j=i+1
       write(1,3001) -6,i+1,PY(i),QY(i),j+1,PY(j),QY(j),EL(k)
      End do
      i=NC
      write(1,3001) -6,i,PY(i-1),QY(i-1),1,AZ(k),QZ,EL(k)
 3001 FORMAT(i5,2(I4,2E14.7),A7)

      write(*,'(a7,e15.5)') EL(k),PY(NC-1)
cSGZ  write(*,'(a7,f10.5)') EL(k),PY(NC-1)
      End do

      write(1,3000) 0
      End



C--------------------------------------------------------------------
C        R e a d _ f u n
C--------------------------------------------------------------------

      Subroutine Read_fun(nf,AF,kz,k,AT,TR)

c     read radial function from file AF (unit=nf)
c     k - serial number for given set

      Implicit real *8 (A-H,O-Z)
      Parameter (mr=220,mfun=99)
      Character EL*4, ELF4*4, AF*12, AT*6,TR*6, EL3*3, EL4*4
      COMMON
     */RFF/ H,EH,RHO,Z,NO,ND,R(mr),RR(mr),R2(mr),X(mr)
     */FUN/ P(mr,mfun),YK(mr),N(mfun),L(mfun),KN(mfun),MX(mfun),
     *      AZ(mfun),EI(mfun),ZETA(mfun),EL(mfun),NWF

      OPEN(nf,file=AF,status='OLD',form='UNFORMATTED')
      Rewind(nf)
      READ(nf,END=2) AT,TR,EL3,MM,Z
      Call ZRFF

      nwf = 0
      rewind(nf)
    1 nwf=nwf+1
      if(nwf.gt.mfun) Stop ' Read_fun: NWF > mfun'
      READ(nf,END=2) AT,TR,EL3,MX(nwf),zz,EI(nwf),ZETA(nwf),AZ(nwf),
     *         (P(j,nwf),j=1,MX(nwf))
      if(zz.ne.Z) Stop ' Read_fun: unconsistent Z'
      EL4 = EL3
      Call EL4_nlk(EL4,N(nwf),L(nwf),KN(nwf))
      if(k.ge.0)  KN(nwf)=KN(nwf)+k
      EL(nwf)=ELF4(N(nwf),L(nwf),KN(nwf))
      Do j=MX(nwf)+1,NO
        P(j,nwf)=0.0
      End do
      go to 1
    2 Close(nf)
      nwf=nwf-1

      Return
      End


C--------------------------------------------------------------------
C        Z R F F
C--------------------------------------------------------------------

      Subroutine ZRFF
c     prepare of logarithmic radial scale
      Implicit real *8 (A-H,O-Z)
      Parameter (mr=220)
      COMMON
     */RFF/  H,EH,RHO,Z,NO,ND,R(mr),RR(mr),R2(mr),X(mr)
 
      NO=mr
      ND=NO-2
      RHO=-4.0
      H=1./16.
      EH=EXP(-H)

      DO I=1,NO
       R(I)=EXP(RHO+(I-1)*H)/Z
       RR(I)=R(I)*R(I)
       R2(I)=SQRT(R(I))
      END DO

      Return
      End


C----------------------------------------------------------------------
C        E L 4 _ n l k
C----------------------------------------------------------------------

      Subroutine EL4_nlk(EL,n,l,k)
      Character EL*4

      if(EL(3:4).eq.'  ') then
       EL(3:4)=EL(1:2)
       EL(1:2)='  '
      elseif(EL(4:4).eq.' ') then
       EL(4:4)=EL(3:3)
       EL(3:3)=EL(2:2)
       EL(2:2)=EL(1:1)
       EL(1:1)=' '
      end if

      if(ichar(EL(3:3)).gt.57) then
       read(EL(1:2),'(i2)') n
       l=LA(EL(3:3))
       k=ICHAR(EL(4:4))-ICHAR('1')+1
      elseif(ichar(EL(4:4)).gt.57) then
       read(EL(2:3),'(i2)') n
       l=LA(EL(4:4))
       k=0
      else
       Stop ' EL4_nlk: unknown format for EL'
      end if

      Return
      End


C----------------------------------------------------------------------
C        E L F 4
C----------------------------------------------------------------------

      Character *4 Function ELF4(n,l,k)
      Character AL,EL*4

      write(EL(1:3),'(i2,a1)') n,AL(l,1)
      if(k.eq.0) EL(4:4)=' '
      if(k.gt.0) then
       i=k+ICHAR('1')-1
       EL(4:4)=CHAR(i)
      end if

      ELF4=EL
      Return
      End


C--------------------------------------------------------------------
C        A L
C--------------------------------------------------------------------

      CHARACTER FUNCTION AL(L,K)
C     provides some spectroscopic symbols
      CHARACTER AS*11,AB*11,AN*16
      DATA AS/'spdfghiklm*'/,
     *     AB/'SPDFGHIKLM*'/,
     *     AN/'0123456789ABCDEF'/
      I=L+1
      IF(K.EQ.5.OR.K.EQ.6) I=(L-1)/2+1
      AL='?'
      IF(I.GE.1.AND.I.LE.11) THEN
        IF(K.EQ.1) AL=AS(I:I)
        IF(K.EQ.2) AL=AB(I:I)
        IF(K.EQ.3) AL=AN(I:I)
        IF(K.EQ.5) AL=AS(I:I)
        IF(K.EQ.6) AL=AB(I:I)
      END IF
      IF(K.EQ.4) THEN
        if(L.eq.-1.OR.L.eq.0) AL='o'
        if(L.eq.+1          ) AL='e'
      END IF
      IF(K.EQ.7) THEN
       if(L.le.0) AL='-'
       if(L.gt.0) AL='+'
       END IF
      RETURN
      END


C--------------------------------------------------------------------
C       L A
C--------------------------------------------------------------------

      Function LA(a)
      Character a, SET*20
      Data SET/'spdfghiklmSPDFGHIKLM'/
      i = INDEX (SET,a)
      if(i.le.10) then
       la=i-1
      elseif(i.le.20) then
       la=i-11
      else
       la=10
      end if
      Return
      End



C---------------------------------------------------------------------
C        D I F F
C---------------------------------------------------------------------

      SUBROUTINE DIFF(I)
C     stores L{F(i)} in the array YK.
C     L{F(i)} = (DD + 2z*r - (l+1/2)^2) [F(i)>
C     L{F(i)} = r**3/2 * L{P(i)}
C     L{P(i)} = (DD + 2z/r - l(l+1)/rr) [P(i)>
      Implicit real *8 (A-H,O-Z)
      Parameter (mr=220,mfun=99)
      Character EL*4
      COMMON
     */RFF/ H,EH,RHO,Z,NO,ND,R(mr),RR(mr),R2(mr),X(mr)
     */FUN/ P(mr,mfun),YK(mr),N(mfun),L(mfun),KN(mfun),MX(mfun),
     *      AZ(mfun),EI(mfun),ZETA(mfun),EL(mfun),NWF
      Data D1/1.0/
C---------------------------------------------------------------------
      MM = MX(I) - 3
      FL = L(I)
      in = L(i)+1
      TZ = Z + Z
      C  = (FL+.5)**2
      HH = 180.*H*H
      DO  K =  4,MM
      YK(K) = (2.*(P(K+3,i)+P(K-3,i)) - 27.*(P(K+2,i)+P(K-2,i)) +
     *       270.*(P(K+1,i)+P(K-1,i)) - 490.*P(K,i) )/HH +
     *             P(K,i)*( TZ*R(K) - C )
      End do
C
C  *****  because of the possibility of extensive cancellation near the
C  *****  origin, search for the point where the asymptotic behaviour
C  *****  begins and smooth the origin.
C  *****  LP(i,r) = O(r^l+1) --> LF(i,r) = O(r^l+5/2)
C
      LEXP = L(I) + 2
      Y1 = YK(4)/R2(4)/R(4)**LEXP
      Y2 = YK(5)/R2(5)/R(5)**LEXP
      Do 10 k = 4,100
       KP = K+2
       Y3 = YK(KP)/R2(KP)/R(KP)**LEXP
       IF (Y2 .EQ. 0.0) go to 10
       IF (ABS(Y1/Y2 - D1) .LT. .05 .AND. ABS(Y3/Y2 - D1) .LT. .05)
     *        GO TO 2
       Y1 = Y2
       Y2 = Y3
   10 Continue
      WRITE (6,1)  I
    1 FORMAT(' ASYMPTOTIC REGION NOT FOUND FOR FUNCTION NUMBER',I3)
      STOP
C
C  *****  asymptotic region has been found
C
    2 KP = K
      KM = KP - 1
      DO K = 1,KM
       YK(K) = Y1*R2(K)*R(K)**LEXP
      End do
c
      MM = MM + 1
      YK(MM) = (-(P(MM+2,i)+P(MM-2,i)) + 16.*(P(MM+1,i)+P(MM-1,i))
     *          -30.*P(MM,i))/(12.*H*H) +
     *          P(MM,i)*(TZ*R(MM) - C )
      MM = MM + 1
      YK(MM) = (P(MM+1,i) + P(MM-1,i) - 2.*P(MM,i))/(H*H) +
     *          P(MM,i)*(TZ*R(MM) - C )
      MM = MM+1
      DO K =MM,NO
       YK(K) = 0.0
      End do
c
      Return
      END


C----------------------------------------------------------------------
      SUBROUTINE MESH(RA,Z,NRANG2,LRANG2)
C----------------------------------------------------------------------
C
C      AUTOMATICALLY GENERATES THE INTEGRATION MESH,
C      ON THE BASIS OF THE NUCLEAR BEHAVIOUR OF BOUND ORBITALS,
C      THE NUMBER OF CONTINUUM ORBITAL LOOPS AND THE CURRENT ARRAY SIZES
C
C-----------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (MZLR2=  20)
      PARAMETER (MZNPT= 800)

      PARAMETER (MXNIX=MZNPT/16)

      COMMON /MESHS/ RX(MZNPT),HINT,IHX(MXNIX),IRX(MXNIX),NIX

      PARAMETER (ZERO=0.0D0)

      DATA MINFAC,MAXFAC/1,2/,MSHDIM/16/

C-----------------------------------------------------------------------
C
C     NCORSE = NRANG2*MSHDIM
      NCORSE = (NRANG2+(LRANG2-2)/2)*MSHDIM
      IF (NCORSE.LT.96) NCORSE = 96
C
C      CALCULATE THE COARSEST MESH, HMAX, AND THE MESH REQUIRED
C      NEAR THE ORIGIN. CALCULATE HMIN = HMAX/2**M WHERE M+1 IS
C      THE NUMBER OF STEP SIZES
C
      HMAX = RA/NCORSE
      HINNER = 0.025D0/Z
      DELTA = HINNER/5.0D0
      M = 0
      HMIN = HMAX
   10 CONTINUE
      IF (HMIN.GT.HINNER+DELTA) THEN
       M = M + 1
       HMIN = HMIN/2
       GOTO 10
      ENDIF

      HINT = HMIN
      NIX = M + 1
      IF (NIX.GT.MXNIX) then
       write(*,*) ' MESH:  NIX > MXNIX',nix,mxnix
       Stop
      End if

C      SET UP THE IHX ARRAY

      IH = 1
      DO I = 1,NIX
        IHX(I) = IH
        IH = IH + IH
      END DO

C      CONSIDER SEPARATELY M .LT. 4  AND  M .GE. 4
C      NA IS THE NUMBER OF STEPS AT EACH STEP SIZE
C      IT IS A MULTIPLE OF 16 AND CAN TAKE VALUES
C      FROM 16*MAXFAC DOWN TO 16*MINFAC

      MPOW2 = 2**M
      NAFAC = MAXFAC
      IF (M.GE.4) THEN

   30   CONTINUE
        NA = 16*NAFAC
        NTOT = NCORSE + (M-1)*NA + NA/8
        IF (NTOT.GE.MZNPT) THEN
          NAFAC = NAFAC - 1
          IF (NAFAC.GE.MINFAC) GOTO 30
        ENDIF

C      SET UP IRX ARRAY

        IRX(2) = NA + NA
        IRX(3) = NA + IRX(2)
        IRX(4) = NA + NA/8 + IRX(3)
        IA = 5

      ELSE

   40   CONTINUE
        NA = 16*NAFAC
        NTOT = NCORSE + M*NA - NA* (MPOW2-1)/MPOW2
        IF (NTOT.GE.MZNPT) THEN
          NAFAC = NAFAC - 1
          IF (NAFAC.GE.MINFAC) GOTO 40
        ENDIF

        IA = 2

      ENDIF

C     FILL IRX ARRAY

      IRX(1) = NA
      DO 50 I = IA,NIX - 1
        IRX(I) = NA + IRX(I-1)
   50 CONTINUE

      IF (NAFAC.GE.MINFAC) GOTO 60
      NA = (((1-M)*MINFAC*16+MZNPT)*MPOW2-MINFAC*16)*NRANG2*HMIN/RA
      WRITE (*,3000) NA
 3000 FORMAT (/
     :' TO SATISFY INTEGRATION MESH CONDITIONS NRANG2 SHOULD BE REDUCED
     :BELOW',I3/' RECOMPILE IF THIS IS UNDESIRABLE:')
      NTOT = NTOT - 2

C      NUMBER OF STEPS AT EACH STEP SIZE MUST BE EVEN

   60 CONTINUE
      IF (MOD(NTOT,2).NE.0) THEN
        NTOT = NTOT + 1
        IRX(NIX-1) = IRX(NIX-1) + 2
      ENDIF

      IRX(NIX) = NTOT

C      PERFORM CHECK

      RVAL = ZERO
      IR = 0
      Do I = 1,NIX
       j1=1
       if(i.gt.1) j1=IRX(i-1)
       j2=IRX(i)
       Do J = J1,J2
        RX(j) = RVAL + HINT*(J-IR)*IHX(I)
       End do
       RVAL = RVAL + HINT* (IRX(I)-IR)*IHX(I)
       IR = IRX(I)
      End do

      Open(2,file='rrr',status='unknown')
cSGZ  Open(2,file='rrr')
      WRITE(2,'(/'' RVAL ='',E14.7)') RVAL
      write(2,'(/'' NIX='',i5)') NIX
      write(2,'( '' IHX='',12I5)') (IHX(i),i=1,NIX)
      write(2,'( '' IRX='',12I5)') (IRX(i),i=1,NIX)
      write(2,'(5F10.5)') (RX(i),i=1,NTOT)

      END


C----------------------------------------------------------------------
C        L A G R N
C----------------------------------------------------------------------

      SUBROUTINE LAGRN(K,N,N1,R,R1,F,F1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(N),F(N),R1(N1),F1(N1),X(10),Y(10)
      K1=K/2+1
      JX=K1
      N2=N-K1+1
      IF(K.EQ.K/2*2) N2=N2+1
      DO 5 I1=1,N1
      XX=R1(I1)
      DO 1 J=JX,N2
      A=ABS(R(J)-XX)
      B=ABS(R(J+1)-XX)
      IF(B.GT.A) GO TO 2
    1 CONTINUE
      J=N2
    2 JX=J
      DO 3 I=1,K
      J=JX-K1+I
      X(I)=R(J)
    3 Y(I)=F(J)
      DO 4 I=1,K
      S=XX-X(I)
      DO 4 J=1,K
    4 IF(I.NE.J) Y(J)=Y(J)*S/(X(J)-X(I))
      F1(I1)=0.0
      DO 5 I=1,K
    5 F1(I1)=F1(I1)+Y(I)
      RETURN
      END
