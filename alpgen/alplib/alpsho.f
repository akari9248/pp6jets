*-- AUTHOR :    MICHELANGELO MANGANO
C----------------------------------------------------------------------
      SUBROUTINE UPVETO(IPVETO)
C----------------------------------------------------------------------
C     SUBROUTINE TO IMPLEMENT THE MLM JET MATCHING CRITERION
C     USING FASTJET TO FORM CLUSTERS WITH ANTI-KT
C----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ITMP
C...GUP EVENT COMMON BLOCK
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &              IDUP(MAXNUP),ISTUP(MAXNUP),MOTHUP(2,MAXNUP),
     &              ICOLUP(2,MAXNUP),PUP(5,MAXNUP),VTIMUP(MAXNUP),
     &              SPINUP(MAXNUP)
C...HEPEVT COMMONBLOCK.
      INTEGER NMXHEP,NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
      PARAMETER (NMXHEP=4000)
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      DOUBLE PRECISION PHEP,VHEP
      SAVE /HEPEVT/
      INCLUDE 'alpsho.inc'
      INTEGER IPVETO

C     JET VARIABLES                             
      INTEGER NJMAX,NCLUS,NCJET
      DOUBLE PRECISION PI,PCLUS,PCJET,ETJET,ETAJ,PHIJ
     $     ,DRJ 
      PARAMETER (PI=3.141593D0)
      PARAMETER (NJMAX=1000)
      COMMON/CLUS/PCLUS(4,NJMAX),PCJET(4,NJMAX),ETJET(NJMAX)
     $     ,ETAJ(NJMAX),PHIJ(NJMAX),DRJ(NJMAX,NJMAX),NCLUS
     $     ,NCJET
C     PARTONIC VARIABLES FOR FASTJET
      INTEGER NFJPART,NMAX
      PARAMETER (NMAX=NJMAX)
      DOUBLE PRECISION PPART(4,NMAX)
      INTEGER IDFJ(NMAX),IJMAP(NJMAX)
C      COMMON/PART/PPART(4,NMAX),NFJPART
C     
      DOUBLE PRECISION PSERAP
      INTEGER K(NJMAX),KP(NJMAX),KPJ(NJMAX)
C LOCAL VARIABLES
      INTEGER I,J,IHEP,NMATCH,JRMIN
      DOUBLE PRECISION ETAJET,PHIJET,DELR,DPHI,DELRMIN
      DOUBLE PRECISION P(4,20),PT(20),ETA(20),PHI(20)
C HEAVY QUARK MATCHING
      INTEGER IHVQ(20),NHVQ,NMJET,ID
      DOUBLE PRECISION ETAHVQ(20),PHIHVQ(20)
      INTEGER IEND,INORAD
      COMMON/SHVETO/IEND,INORAD(MAXNUP)
C HEAVY QUARK VETOES
      INTEGER NC,NB
C DEBUGGING OPTIONS
      INTEGER IDBG
      PARAMETER (IDBG=2)
      DOUBLE PRECISION PTPART,PTJETS,ETAPART,ETAJETS
      INTEGER NMMAX
      COMMON/MTCHDBG/PTPART(20),PTJETS(20),ETAPART(20),ETAJETS(20)
     $     ,NMMAX
C-LIZA
      INTEGER IVETO,ISVETO,NJVETO
      COMMON/MTCHSTAT/IVETO,ISVETO,NJVETO
C JETS FROM RADIATION OF LIGHT PARTONS:
C ISVETO=1 EVENT FAILS SINCE THERE ARE FEWER RECONSTRUCTED JETS THAN PARTONS
C     NJVETO=NUMBER OF RECONSTRUCTED JETS
C ISVETO=2 EVENT FAILS SINCE NOT ALL PARTONS MATCH A JET
C     NJVETO=NUMBER OF MATCHED PARTONS/JETS
C ISVETO=3 EVENT FAILS SINCE THERE ARE ADDITIONAL JETS
C     NJVETO= NUMBER OF JETS
C ISVETO=4 EVENT FAILS SINCE MATCHED JETS ARE NOT THE HARDEST (INCLUSIVE CASE)
C     NJVETO = HARDEST JET NOT MATCHED (1=HIGHEST ET, 2=2ND...
C 
C JETS FROM RADIATION OF HEAVY QUARKS, AT LARGE ANGLE W.R.T. HQ:
C ISVETO=5 EVENT FAILS SINCE THERE ARE ADDITIONAL JETS
C     NJVETO= NUMBER OF UNMATCHED JETS
C ISVETO=6 EVENT FAILS SINCE MATCHED JETS ARE NOT THE HARDEST (INCLUSIVE CASE)
C     NJVETO=NUMBER OF ADDITIONAL JETS
C LIZA-
C
      DOUBLE PRECISION ETMIN  !,ETMAX
      DOUBLE PRECISION TINY
      PARAMETER (TINY=1D-3)
      INTEGER ICOUNT,NTMP
      DATA ICOUNT/0/
C
      IPVETO=0
C-LIZA
      IVETO=IPVETO
      ISVETO=0
      NJVETO=0
C LIZA-
C     VETO ON EXTRA HEAVY QUARKS FROM SHOWER OR MATRIX ELEMENT
      IF(IHVYV.EQ.1) THEN
        CALL HVQVETO(NPFST,NPLST,NJLAST,NC,NB)
        IF(NC.LE.NCMAX.AND.NB.LE.NBMAX) GOTO 5
        IF(NC.GT.NCMAX) THEN
          NVETOC(NC)=NVETOC(NC)+1
        ELSEIF(NB.GT.NBMAX) THEN 
          NVETOB(NB)=NVETOB(NB)+1
        ENDIF
        GOTO 999
      ENDIF
 5    IF(ICKKW.EQ.0) RETURN
      IF(IHRD.EQ.8.OR.IHRD.EQ.13.OR.IHRD.EQ.20) THEN
        WRITE(*,*) 'JET MATCHING FOR HARD PROCESS ',IHRD
     $       ,' NOT IMPLEMENTED, STOP'
        STOP
      ENDIF
      IF(NLJETS.EQ.0.AND.IEXC.EQ.0) RETURN
C CHECK FOR EVENT ERROR OR ZERO WGT
      I=0
C     HERWIG/PYTHIA SPECIFIC
      CALL ALSHER(I)
      IF(I.EQ.1) RETURN
C
C     INITIALIZE DEBUG IF NEEDED
      IF(IDBG.EQ.1) THEN
        WRITE(1,*) ' '
        WRITE(1,*) 'NEW EVENT '
        WRITE(1,*) 'PARTONS'
      ENDIF
      IF(IDBG.EQ.2) THEN
        DO I=1,20
          PTPART(I)=0D0
          ETAPART(I)=0D0
          PTJETS(I)=0D0
          ETAJETS(I)=0D0
        ENDDO
        NMMAX=0
      ENDIF
C
C     RECONSTRUCT PARTON-LEVEL EVENT
C     START FROM THE PARTONIC SYSTEM
      DO I=1,NLJETS
        IHEP=I+NJSTART
        DO J=1,4
          P(J,I)=PUP(J,IHEP)
        ENDDO
        PT(I)=SQRT(P(1,I)**2+P(2,I)**2)
        ETA(I)=-LOG(TAN(0.5D0*ATAN2(PT(I)+TINY,P(3,I))))
        PHI(I)=ATAN2(P(2,I),P(1,I))
        IF(IDBG.EQ.1) THEN
          WRITE(1,*) PT(I),ETA(I),PHI(I)
        ENDIF
      ENDDO
      IF(NLJETS.GT.0) CALL ALPSOR(PT,NLJETS,KP,2)  
      IF(IDBG.EQ.2) THEN
        DO I=1,NLJETS
          PTPART(I)=PT(KP(NLJETS-I+1))
          ETAPART(I)=ETA(KP(NLJETS-I+1))
        ENDDO
      ENDIF
C     
C     DISPLAY EVENT SEEN BY UPVETO:
C      IF(IDBG.EQ.1) THEN
C        DO I=1,NHEP
C          WRITE(1,111) I,ISTHEP(I),IDHEP(I),JMOHEP(1,I),JMOHEP(2,I)
C     $         ,PHEP(1,I),PHEP(2,I),PHEP(3,I)
C        ENDDO
C 111  FORMAT(5(I4,1X),3(F12.5,1X))
C      ENDIF
C     DISPLAY PYTHIA EVENT:
C      CALL PYLIST(7)    ! PYTHIA USER PROCESS EVENT DUMP
C      CALL PYLIST(2)    ! PYTHIA FULL EVENT DUMP
C      CALL PYLIST(5)    ! PYTHIA HEPEVT DUMP
C
C     RECONSTRUCT SHOWERED JETS:
C CONVERT MOMENTA TO FASTJET FORMAT
      CALL FJHEVT_M(NPFST,NPLST,NJLAST,PPART,IDFJ,NFJPART)
C FASTJET CLUSTERING, RETURNS FASTJET JETS
      CALL FJCLUS(RCLUS,-1d0,PPART,NFJPART,PCLUS,NCLUS)
C APPLY PT/ETA CUTS
      CALL FJSORT(ETCLUS,ETACLMAX,PCLUS,NCLUS,PCJET,ETJET
     $     ,ETAJ,PHIJ,DRJ ,NCJET, IJMAP )
c
      IF(NCJET.GT.0) CALL ALPSOR(ETJET,NCJET,K,2) 
      IF(IDBG.EQ.1) THEN
        WRITE(1,*) 'JETS'
        DO I=1,NCJET
          J=K(NCJET+1-I)
          WRITE(1,*) ETJET(J),ETAJ(J),PHIJ(J)
        ENDDO
      ENDIF
      IF(IDBG.EQ.2) THEN
        NMMAX=NCJET
        DO I=1,NCJET
          J=K(NCJET+1-I)
          PTJETS(I)=ETJET(J)
          ETAJETS(I)=ETAJ(J)
        ENDDO
      ENDIF
C     ANALYSE ONLY EVENTS WITH AT LEAST NLJETS-RECONSTRUCTED JETS
      IF(NCJET.LT.NLJETS) THEN
        NVETONOJ=NVETONOJ+1
C EVENT FAILS SINCE THERE ARE FEWER RECONSTRUCTED JETS THAN PARTONS
        ISVETO=1
C NJVETO HERE REPRESENT THE NUMBER OF RECONSTRUCTED JETS
        NJVETO=NCJET
        GOTO 999
      ENDIF
C     ASSOCIATE PARTONS AND JETS, USING MIN(DELR) AS CRITERION
      NMATCH=0
      DO I=1,NCJET
        KPJ(I)=0
      ENDDO
      ETMIN=1D10
      DO I=1,NLJETS
        DELRMIN=1D5
        DO 110 J=1,NCJET
          IF(KPJ(J).NE.0) GO TO 110
          ETAJET=ETAJ(J)
          PHIJET=PHIJ(J)
          DPHI=ABS(PHI(KP(NLJETS-I+1))-PHIJET)
          IF(DPHI.GT.PI) DPHI=2.*PI-DPHI
          DELR=SQRT((ETA(KP(NLJETS-I+1))-ETAJET)**2+(DPHI)**2)
          IF(DELR.LT.DELRMIN) THEN
            DELRMIN=DELR
            JRMIN=J
          ENDIF
 110    CONTINUE
        IF(DELRMIN.LT.1.5*RCLUS) THEN
          NMATCH=NMATCH+1
          KPJ(JRMIN)=I
          ETMIN=MIN(ETMIN,ETJET(JRMIN))
C     ASSOCIATE PARTONS AND MATCHED JETS:
          IF(IDBG.EQ.2) THEN
c            PTJETS(I)=ETJET(JRMIN)
c            ETAJETS(I)=ETAJ(JRMIN)
          ENDIF
C          WRITE(*,*) 'PARTON-JET',I,' BEST MATCH:',K(NCJET+1-JRMIN)
C     $           ,DELRMIN
        ENDIF
      ENDDO
      IF(NMATCH.LT.NLJETS) THEN
        NVETOMCH(NMATCH)=NVETOMCH(NMATCH)+1
C EVENT FAILS SINCE NOT ALL PARTONS MATCH A JET
        ISVETO=2
C NJVETO HERE REPRESENT THE NUMBER OF MATCHED PARTONS/JETS
        NJVETO=NMATCH
        GOTO 999
      ENDIF
C     REJECT EVENTS WITH LARGER JET MULTIPLICITY FROM EXCLUSIVE SAMPLE
      IF(NCJET.GT.NLJETS.AND.IEXC.EQ.1) THEN
        NVETOEXC(NCJET-NLJETS)=NVETOEXC(NCJET-NLJETS)+1
C EVENT FAILS SINCE THERE ARE ADDITIONAL JETS
        ISVETO=3
C NJVETO HERE REPRESENT THE NUMBER OF JETS
        NJVETO=NCJET
        GOTO 999
      ENDIF
C     VETO EVENTS WHERE MATCHED JETS ARE SOFTER THAN NON-MATCHED ONES
      IF(IEXC.NE.1) THEN
        NTMP=0
        J=NCJET
        DO I=1,NLJETS
C     KPJ(K(J)) IS THE PARTON THAT MATCHES THE J-TH SOFTEST JET. HERE
C     K(NCJET) IS THE HARDEST JET, K(1) IS THE SOFTEST
          IF(KPJ(K(J)).EQ.0) THEN
C EVENT FAILS SINCE MATCHED JETS ARE NOT THE HARDEST
            ISVETO=4
C NJVETO HERE REPRESENT THE HARDEST JET NOT MATCHED (1=HIGHEST ET, 2=2ND...
            IF(NTMP.EQ.0) NJVETO=NCJET-J+1
            NTMP=NTMP+1
          ENDIF
          J=J-1
        ENDDO
        IF(NTMP.GT.0) THEN
          NVETOEXC(NTMP)=NVETOEXC(NTMP)+1
          GOTO 999
        ENDIF
      ENDIF
C
C     ADDITIONAL TREATMENT FOR HVQ EVENTS: VETO/ACCEPT JETS EMITTED BY HVQ'S
      IF(IHRD.LE.2.OR.IHRD.EQ.6.OR.IHRD.EQ.7.OR.IHRD.EQ.10.OR.IHRD.EQ.15
     $   .OR.IHRD.EQ.16  ) THEN
C     RECOSTRUCT POSSIBLE JETS FROM RADIATION OFF THE HEAVY QUARKS
C     CONVERT MOMENTA TO FASTJET FORMAT
        CALL FJHEVT_HVQ(NPFST,NPLST,NJLAST,PPART,IDFJ,NFJPART)
C     FASTJET CLUSTERING, RETURNS FASTJET JETS
        CALL FJCLUS(RCLUS,-1d0,PPART,NFJPART,PCLUS,NCLUS)
C     APPLY PT/ETA CUTS
        CALL FJSORT(ETCLUS,ETACLMAX,PCLUS,NCLUS,PCJET,ETJET
     $     ,ETAJ,PHIJ,DRJ ,NCJET )
c
C      IF(NCJET.GT.0) CALL ALPSOR(ETJET,NCJET,K,2) 
C        CALL CALINI_M
C     RECOSTRUCT POSSIBLE JETS FROM RADIATION OFF THE TOP QUARKS
C        CALL CALDEL_HVQ(NPFST,NPLST,NJLAST)
C        CALL GETJET_M(RCLUS,ETCLUS,ETACLMAX)


C     IF NO EXTRA JET: ACCEPT EVENT
        IF(NCJET.EQ.0) RETURN
C     IF EXTRA JETS, REMOVE THOSE LYING WITHIN DRJMIN OF A B/C QUARK, TO
C     ALLOW THE SHOWER TO GOVERN THE DEVELOPMENT OF A B/C JET.
C     START BY FLAGGING VETOED Q AND QBAR OBJECTS:
        NHVQ=0
        DO I=1,NHEP
          ID=IDHEP(I)
          IF(INORAD(I).EQ.1.AND.ABS(ID).LE.5.AND.ABS(ID)
     $         .GE.4) THEN
            NHVQ=NHVQ+1
            IHVQ(NHVQ)=I
            ETAHVQ(NHVQ)=PSERAP(PHEP(1,I))
            PHIHVQ(NHVQ)=ATAN2(PHEP(2,I),PHEP(1,I))
          ENDIF
        ENDDO
        NMJET=NCJET
        DO I=1,NCJET
          ETAJET=ETAJ(I)
          PHIJET=PHIJ(I)
          DO J=1,NHVQ
            DPHI=ABS(PHIHVQ(J)-PHIJET)
            IF(DPHI.GT.PI) DPHI=ABS(DPHI-2*PI)
            DELR=SQRT(DPHI**2+(ETAJET-ETAHVQ(J))**2)
            IF(DELR.LT.DRJMIN) THEN
              NMJET=NMJET-1
              ETJET(I)=0D0
              GOTO 200
            ENDIF
          ENDDO
 200      CONTINUE
        ENDDO
C     IF NO JETS UNMATCHED TO HEAVY QUARKS: ACCEPT EVENT
        IF(NMJET.EQ.0) RETURN
C     IF UNMATCHED JETS AND IEXC=1: REJECT EVENT
        IF(IEXC.EQ.1) THEN
          NVETOHVY(NMJET)=NVETOHVY(NMJET)+1
C EVENT FAILS SINCE THERE ARE ADDITIONAL JETS
          ISVETO=5
C NJVETO HERE REPRESENT THE NUMBER OF UNMATCHED JETS
          NJVETO=NMJET
          GOTO 999
        ENDIF
C     IF JETS AND IEXC=0: CHECK THAT JETS ARE SOFTER THAN MATCHED ONES
        NTMP=0
        DO I=1,NCJET
          IF(ETJET(I).GT.ETMIN) NTMP=NTMP+1
        ENDDO
        IF(NTMP.GT.0) THEN
          NVETOHVY(NTMP)=NVETOHVY(NTMP)+1
C EVENT FAILS SINCE THERE ARE JETS HARDER THAN THE MATCHED ONES
          ISVETO=6
C NJVETO HERE REPRESENTS THE NUMBER OF ADDITIONAL JETS
          NJVETO=NCJET
          GOTO 999
        ENDIF
      ENDIF
      RETURN
C     HERWIG/PYTHIA TERMINATION:
 999  CALL ALSHEN
      IPVETO=1
c flag the vent as rejected
      IVETO=IPVETO
c uncomment line below allow the event to arrive to the analysis stage
c      IPVETO=0
      END


C----------------------------------------------------------------------
      SUBROUTINE UPINIT
C----------------------------------------------------------------------
C     HERWIG/PYTHIA UNIVERSAL EVENT INITIALITION ROUTINE
C----------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'alpsho.inc'
      data nvetoc/100*0/,nvetob/100*0/,nvetomch/11*0/,nvetoexc/10*0/
     $     ,nvetohvy/10*0/,nvetonoj/0/
      CHARACTER *3 CSHO
C--   GUP Event common block
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &     IDUP(MAXNUP),ISTUP(MAXNUP),MOTHUP(2,MAXNUP),
     &     ICOLUP(2,MAXNUP),PUP(5,MAXNUP),VTIMUP(MAXNUP),
     &     SPINUP(MAXNUP)
C--   GUP common block
      INTEGER MAXPUP
      PARAMETER(MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON /HEPRUP/ IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &     IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),
     &     XMAXUP(MAXPUP),LPRUP(MAXPUP)
C     CALSIM AND JET VARIABLES
c      INTEGER NCY,NCPHI,NJMAX,JETNO,NCJET
c      DOUBLE PRECISION YCMIN,YCMAX,PI,ET,DELPHI,CPHCAL,SPHCAL,DELY,
c     &     CTHCAL,STHCAL,PCJET,ETJET
c      PARAMETER (NCY=100)
c      PARAMETER (NCPHI=60,PI=3.141593D0)
c      COMMON/CALOR_M/DELY,DELPHI,ET(NCY,NCPHI),
c     $     CTHCAL(NCY),STHCAL(NCY),CPHCAL(NCPHI),SPHCAL(NCPHI),YCMIN
c     $     ,YCMAX
C     
      INTEGER IEND,INORAD
      COMMON/SHVETO/IEND,INORAD(MAXNUP)
c     
C     LOCAL VARIABLES
      CHARACTER*70 STDUMMY
      INTEGER IPR
      INTEGER MAXLEN
      PARAMETER (MAXLEN=200)
      CHARACTER*(MAXLEN) STRING
      INTEGER I,NTMP,IHEPMIN,IERR,IBEG,ISTRLH
      DOUBLE PRECISION TMP
      DOUBLE PRECISION PBEAM1,PBEAM2
C     USER ACCESS TO MATCHING PARAMETERS
      INTEGER IUSRMAT
      PARAMETER (IUSRMAT=1)
C     
      WRITE(*,*) 'INPUT NAME OF FILE CONTAINING EVENTS'
      WRITE(*,*) '(FOR "file.unw" or "file.lhe" ENTER "file")'
      READ(*,*) FILENAME
C
      CALL ALSHCD(CSHO)
C     OPEN A LOG FILE
      IF(CSHO.EQ.'HER') THEN
        CALL STRCATH(FILENAME,'.her-log',TMPSTR)
      ELSE 
        CALL STRCATH(FILENAME,'.pyt-log',TMPSTR)
      ENDIF
      CALL GETUNIT(NUNITOUT)
      OPEN(UNIT=NUNITOUT,FILE=TMPSTR,STATUS='UNKNOWN')
C CHECK IF OLD FORMAT
      CALL STRCATH(FILENAME,'_unw.par',TMPSTR)
      CALL GETUNIT(NUNITINI)
      OPEN(UNIT=NUNITINI,FILE=TMPSTR,STATUS='OLD',ERR=101)
      WRITE(NUNITOUT,*) 'READING PARAMS FROM ',TMPSTR(1:istrlH(tmpstr)),
     $     ' UNW FILE'
      CALL STRCATH(FILENAME,'.unw',TMPSTR)
      CALL GETUNIT(NUNIT)
      OPEN(UNIT=NUNIT,FILE=TMPSTR,STATUS='OLD',ERR=101)
      WRITE(NUNITOUT,*) 'READING EVENTS FROM ',TMPSTR(1:istrlH(tmpstr)),
     $     ' UNW FILE'
      ILHE=0
      GOTO 200
C IF WE ARE HERE, IT'S NEW FORMAT
 101  CALL STRCATH(FILENAME,'.lhe',TMPSTR)
      CALL GETUNIT(NUNITINI)
      NUNIT=NUNITINI                  ! read events from same file as params
      OPEN(UNIT=NUNITINI,FILE=TMPSTR,STATUS='OLD',ERR=102)
      WRITE(NUNITOUT,*) 'READING PARAMS/EVENTS FROM '
     $     ,TMPSTR(1:istrlH(tmpstr)),' LHE FILE'
C SKIP FIRST TWO LINES OF LHE FILE HEADER
      READ(NUNITINI,*) 
      READ(NUNITINI,*)
      ILHE=1
      GOTO 200
 102  WRITE(*,*) 'NO DATA CORRESPONDING TO THIS FILE, STOP'
      WRITE(NUNITOUT,*) 'NO DATA CORRESPONDING TO THIS FILE, STOP'
      STOP
 200  CONTINUE
C read in parameters
C     START READING FILE
      DO I=1,10000
        READ(NUNITINI,'(A)') STDUMMY
        IF(STDUMMY(1:4).EQ.'****') GOTO 300
        WRITE(*,*) STDUMMY
        WRITE(NUNITOUT,*) STDUMMY
      ENDDO
C     
C     READ IN INPUT PARAMETERS
 300  READ(NUNITINI,*) IHRD
C     
      READ(NUNITINI,*) MC,MB,MT,MW,MZ,MH
      DO I=1,1000
        READ(NUNITINI,*,ERR=310) NTMP,TMP
        PARVAL(NTMP)=TMP
      ENDDO
 310  CONTINUE
C
C     WRITE PARAMETER VALUES INHERITED FROM PARTON LEVEL GENERATION
      CALL AHSPAR
      PBEAM1=DBLE(EBEAM)
      PBEAM2=DBLE(EBEAM)
      IH1=1
      READ(NUNITINI,*) AVGWGT,ERRWGT
      READ(NUNITINI,*) UNWEV,TOTLUM
      IF(ILHE.EQ.0) CLOSE(NUNITINI)
      WRITE(NUNITOUT,*) " "
      WRITE(NUNITOUT,*) "INPUT CROSS SECTION (PB):",AVGWGT," +/-"
     $     ,ERRWGT
      WRITE(NUNITOUT,*) "NUMBER OF INPUT EVENTS:",UNWEV
      WRITE(NUNITOUT,*) "INTEGRATED LUMINOSITY:",TOTLUM
C     FROM NOW ON, PROCESS THE INFORMATION READ IN, TO COMPLETE SETTING
C     UP THE GUP COMMON
C     
      IF(ILHE.EQ.0) THEN
C--   SET UP THE BEAMS
C--   ID'S OF BEAM PARTICLES
        IDBMUP(1) = 2212
        IF(IH2.EQ.1) THEN
          IDBMUP(2) = 2212
        ELSEIF(IH2.EQ.-1) THEN
          IDBMUP(2) =-2212
        ELSE
          WRITE(*,*) 'BEAM 2 NOT PROPERLY INITIALISED, STOP'
          STOP
        ENDIF
        EBMUP(1) = ABS(PBEAM1)
        EBMUP(2) = ABS(PBEAM2)
C--   PDF'S FOR THE BEAMS; WILL BE EVALUATED USING THE NDNS VARIABLE
C     READ IN EARLIER
        PDFGUP(1) = -1
        PDFGUP(2) = -1
        PDFSUP(1) = -1
        PDFSUP(2) = -1
C--   WHAT DO DO WITH THE WEIGHTS(WE ARE GENERATING UNWEIGHTED EVENTS)
        IDWTUP = 3
C--   ONLY ONE PROCESS
        NPRUP  = 1
C--   CROSS SECTION
        XSECUP(1) = avgwgt
C--   ERROR ON THE CROSS SECTION
        XERRUP(1) = errwgt
C--   MAXIMUM WEIGHT
        XMAXUP(1) = avgwgt
      ELSE
C...Loop until finds line beginning with "<init>" or "<init ". 
        IERR=1
        
 400    READ(NUNITINI,'(A200)',END=510,ERR=510) STRING
        IBEG=0
 410    IBEG=IBEG+1
C...  Allow indentation.
        IF(STRING(IBEG:IBEG).EQ.' '.AND.IBEG.LT.MAXLEN-5) GOTO 410 
        IF(STRING(IBEG:IBEG+5).NE.'<init>'.AND.
     &       STRING(IBEG:IBEG+5).NE.'<init ') GOTO 400
        READ(NUNITINI,*,END=510,ERR=510) IDBMUP(1),IDBMUP(2),EBMUP(1),
     &       EBMUP(2),PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),IDWTUP
     $       ,NPRUP
        WRITE(*,*) IDBMUP(1),IDBMUP(2),EBMUP(1),
     &       EBMUP(2),PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),IDWTUP
     $       ,NPRUP
        DO 500 IPR=1,NPRUP
          READ(NUNITINI,*,END=510,ERR=510) XSECUP(IPR),XERRUP(IPR),
     &         XMAXUP(IPR),LPRUP(IPR)
          WRITE(*,*) XSECUP(IPR),XERRUP(IPR),
     &         XMAXUP(IPR),LPRUP(IPR)
 500    CONTINUE
        IERR=0
 510    IF(IERR.EQ.1) THEN
          WRITE(*,*) 'CANNOT READ INITIALIZATION INFORMATION, STOP'
          STOP
        ENDIF
C SKIP </init> line of header
        READ(NUNITINI,*) 
      ENDIF
C--   HERWIG/PYTHIA SPECIFIC SETTIGS
      CALL ALSHIN(I)
      LPRUP(1) = I
C
C     ASSIGN DEFAULT VALUES FOR SHOWER GENERATION, AND ALLOW USER TO MODIFY
      CALL AHSHPA
C     CALORIMETER ETA RANGE
C      YCMAX=ETACLMAX+RCLUS
C      YCMIN=-YCMAX
C
C     CONVERT PDF TYPES (NOT FOR LHAPDF SETS)
      IF(NDNS.LT.500) CALL PDFCONVH(NDNS,NTMP,PDFTYP)
C     SETUP MATCHING PARAMETERS
C     DEFINE RANGE FOR PARTONS TO BE USED IN MATCHING
      DO I=1,MAXNUP
        INORAD(I)=0
      ENDDO
      IF(CSHO.EQ.'HER') THEN
        NPFST=149
        NPLST=149
C     HERWIG: ALL SHOWERS ORIGINATE FROM IHEP=6
        IEND=6
C     HERWIG: HEPEVT EVENT RECORD FOR FINAL STATE STARTS AT 7=6+1
        IHEPMIN=6
      ELSE
        NPFST=1
        NPLST=1
C     PYTHIA: ALL SHOWERS ORIGINATE FROM IHEP=O
        IEND=0
C     PYTHIA: HEPEVT EVENT RECORD FOR FINAL STATE STARTS AT 1=0+1
        IHEPMIN=0
        IDPRUP=661
      ENDIF
      IF(IHRD.LE.2) THEN
C     NLJETS=NPART-6
        NJSTART=4
        NJLAST=155
C     DO NOT INCLUDE IN MATCHING RADIATION ORIGINATING FROM THE W
        INORAD(IHEPMIN+2+NLJETS+1)=1
C     DO NOT INCLUDE IN MATCHING RADIATION ORIGINATING FROM HEAVY QUARK
C     PAIR
        INORAD(IHEPMIN+1)=1
        INORAD(IHEPMIN+2)=1
      ELSEIF(IHRD.LE.4) THEN
C     NLJETS=NPART-4
        NJSTART=2
        NJLAST=155
C     DO NOT INCLUDE IN MATCHING RADIATION ORIGINATING FROM THE W
        IF(CSHO.EQ.'HER'.AND.NLJETS.EQ.0) THEN
C     ALLOW RADIATION FROM IHEPMIN+1 IN HERWIG WHEN NJET=0
          INORAD(IHEPMIN+NLJETS+1)=0
        ELSE 
          INORAD(IHEPMIN+NLJETS+1)=1
        ENDIF
      ELSEIF(IHRD.EQ.5) THEN
C     NLJETS=NPART-3*(NW+NZ)-NH-2
        NJSTART=2
        NJLAST=155
C     DO NOT INCLUDE IN MATCHING RADIATION ORIGINATING FROM THE W
        IF(CSHO.EQ.'HER'.AND.NLJETS.EQ.0.AND.(NW+NZ+NH+NPH).EQ.1) THEN
C     ALLOW RADIATION FROM IHEPMIN+1 IN HERWIG WHEN NJET=0 AND NW+NZ+NH
C     +NPH=1
          INORAD(IHEPMIN+NLJETS+1)=0
        ELSE 
          DO I=1,NW+NZ+NH+NPH
            INORAD(IHEPMIN+NLJETS+I)=1
          ENDDO
        ENDIF
      ELSEIF(IHRD.EQ.6) THEN
C     DO NOT INCLUDE IN MATCHING RADIATION ORIGINATING FROM HEAVY QUARK
C     PAIR
C     NLJETS=NPART-8  (IHVY.EQ.6)     NLJETS=NPART-4  (IHVY.LT.6)
        INORAD(IHEPMIN+1)=1
        INORAD(IHEPMIN+2)=1
        NJSTART=4
        NJLAST=155
      ELSEIF(IHRD.EQ.7) THEN
C     DO NOT INCLUDE IN MATCHING RADIATION ORIGINATING FROM HEAVY QUARK
C     PAIR
C     NLJETS=NPART-8  (IHVY.EQ.6)     NLJETS=NPART-4  (IHVY.LT.6)
        INORAD(IHEPMIN+1)=1
        INORAD(IHEPMIN+2)=1
        INORAD(IHEPMIN+3)=1
        INORAD(IHEPMIN+4)=1
        NJSTART=6
        NJLAST=155
      ELSEIF(IHRD.EQ.9) THEN
C     NLJETS=NPART-2
        NJSTART=2
        NJLAST=155
      ELSEIF(IHRD.EQ.10) THEN
C     NLJETS=NPART-4
        NJSTART=3
        NJLAST=155
C     DO NOT INCLUDE IN MATCHING RADIATION ORIGINATING FROM THE W
        INORAD(IHEPMIN+NLJETS+1+1)=1
C     DO NOT INCLUDE IN MATCHING RADIATION ORIGINATING FROM THE CHARM
        INORAD(IHEPMIN+1)=1
      ELSEIF(IHRD.EQ.11) THEN
C     NLJETS=NPART-2-NPH
        NJSTART=2
        NJLAST=155
C     DO NOT INCLUDE IN MATCHING THE HARD PHOTONS
        DO I=1,NPH
          INORAD(IHEPMIN+NLJETS+I)=1
        ENDDO
      ELSEIF(IHRD.EQ.12) THEN
C     NLJETS=NPART-2-NH
        NJSTART=2
        NJLAST=155
C     DO NOT INCLUDE IN MATCHING THE HIGGS DECAY PRODUCTS
        NH=1
        IF(NLJETS+NH.GT.1) THEN
          DO I=1,NH
            INORAD(IHEPMIN+NLJETS+I)=1
          ENDDO
        ELSE
          IF(CSHO.EQ.'HER'.AND.NLJETS.EQ.0) THEN
C     ALLOW RADIATION FROM IHEPMIN+1 IN HERWIG WHEN NJET=0
            INORAD(IHEPMIN+1)=0
          ELSE 
            INORAD(IHEPMIN+NLJETS+1)=1
          ENDIF
        ENDIF
      ELSEIF(IHRD.EQ.14) THEN
C     NLJETS=NPART-4-NPH
        NJSTART=2
        NJLAST=155
C     DO NOT INCLUDE IN MATCHING RADIATION ORIGINATING FROM THE W
        IF(CSHO.EQ.'HER'.AND.NLJETS.EQ.0.AND.NPH.EQ.0) THEN
C     ALLOW RADIATION FROM IHEPMIN+1 IN HERWIG WHEN NJET=0
          INORAD(IHEPMIN+NLJETS+1)=0
        ELSE 
          INORAD(IHEPMIN+NLJETS+1)=1
        ENDIF
      ELSEIF(IHRD.EQ.15) THEN
C     NLJETS=NPART-6-NPH
        NJSTART=2
        NJLAST=155
C     DO NOT INCLUDE IN MATCHING RADIATION ORIGINATING FROM THE W
        INORAD(IHEPMIN+2+NLJETS+NPH+1)=1
C     DO NOT INCLUDE IN MATCHING RADIATION ORIGINATING FROM HEAVY QUARK
C     PAIR
C     NOTICE THAT HERE THE LIGHT JETS PRECEDE THE QQ PAIR
        INORAD(IHEPMIN+NLJETS+1)=1
        INORAD(IHEPMIN+NLJETS+2)=1
      ELSEIF(IHRD.EQ.16) THEN
C     DO NOT INCLUDE IN MATCHING RADIATION ORIGINATING FROM HEAVY QUARK
C     PAIR
C     NLJETS=NPART-8  (IHVY.EQ.6)     NLJETS=NPART-4  (IHVY.LT.6)
        INORAD(IHEPMIN+1)=1
        INORAD(IHEPMIN+2)=1
        NJSTART=4
        NJLAST=155
C     DO NOT INCLUDE IN MATCHING THE HARD PHOTONS
        DO I=1,NPH
          INORAD(IHEPMIN+2+NLJETS+I)=1
        ENDDO
      ELSEIF(IHRD.EQ.17) THEN
C     NLJETS=NPART-4-NPH
        NJSTART=2
        NJLAST=155
C     DO NOT INCLUDE IN MATCHING RADIATION ORIGINATING FROM THE Z
        IF(CSHO.EQ.'HER'.AND.NLJETS.EQ.0.AND.NPH.EQ.0) THEN
C     ALLOW RADIATION FROM IHEPMIN+1 IN HERWIG WHEN NJET=0
          INORAD(IHEPMIN+NLJETS+1)=0
        ELSE 
          INORAD(IHEPMIN+NLJETS+1)=1
        ENDIF
      ENDIF
C
      END

*--   Author :    Michelangelo Mangano
C----------------------------------------------------------------------
      SUBROUTINE UPEVNT
C     DECK  ID>, UPEVNT.
*     CMZ :-        -13/02/02  07.20.46  by  Peter Richardson
*--   Author :    Michelangelo Mangano
      implicit none
      INCLUDE 'alpsho.inc'
      IF(ILHE.EQ.0) THEN
        CALL UPAUNW
      ELSEIF(ILHE.EQ.1) THEN
        CALL UPALHE
      ELSE
        WRITE(*,*) 'UNRECOGNIZED DATA FORMAT, STOP'
        STOP
      ENDIF
      END
C----------------------------------------------------------------------
      SUBROUTINE UPAUNW
C----------------------------------------------------------------------
c     Puts Alpgen event from .unw format into GUPI common block HEPEU
c----------------------------------------------------------------------
      implicit none
      INCLUDE 'alpsho.inc'
C--   GUP Event common block
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &     IDUP(MAXNUP),ISTUP(MAXNUP),MOTHUP(2,MAXNUP),
     &     ICOLUP(2,MAXNUP),PUP(5,MAXNUP),VTIMUP(MAXNUP),
     &     SPINUP(MAXNUP)
C--   GUP Run common block
      INTEGER MAXPUP
      PARAMETER(MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON /HEPRUP/ IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &     IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),
     &     XMAXUP(MAXPUP),LPRUP(MAXPUP)
c     
c     local variables
      INTEGER INIT
      DATA INIT/0/
      CHARACTER *3 CSHO
      INTEGER MAXPAR
      PARAMETER (MAXPAR=100)
      INTEGER NEV,IPROC,IFL(MAXPAR),NEVPROC
      DATA NEVPROC/0/
      REAL SQ,SP(3,MAXPAR),SM(MAXPAR),SWGTRES
      INTEGER I,IUP,IWCH,IST
      REAL *8 TMP,WGTRES
C     LOCAL VARIABLES FOR TOP DECAYS
      INTEGER IT(2),ITB(2),IW,IWDEC,IBUP,IWUP,IJETS
C     LOCAL VARIABLES TOP HIGGS DECAYS
      INTEGER IH
C     LOCAL VARIABLES FOR GAUGE BOSON  DECAYS
      INTEGER IVSTART,IVEND,NVB
C UPDATE MAXIMUM NMBER OF ALLOWED ERRORS
      CALL ALSHER(I)
C
      IST=0
      NEVPROC=NEVPROC+1
      IF(NEVMAX.GT.0.AND.NEVPROC.GT.NEVMAX) GOTO 500
C     INPUT EVENT NUMBER, PROCESS TYPE, N PARTONS, SAMPLE'S AVERAGE
C     WEIGHT AND QSCALE
      READ(NUNIT,2,END=500,ERR=501) NEV,IPROC,NPART,SWGTRES,SQ
 2    FORMAT(I8,1X,I4,1X,I2,2(1X,E12.6))
C     FLAVOUR, COLOUR AND Z-MOMENTUM OF INCOMING PARTONS
      READ(NUNIT,8) IFL(1),ICOLUP(1,1),ICOLUP(2,1),SP(3,1)
      READ(NUNIT,8) IFL(2),ICOLUP(1,2),ICOLUP(2,2),SP(3,2)
C     FLAVOUR, COLOUR, 3-MOMENTUM AND MASS OF OUTGOING PARTONS
      DO I=3,NPART
        READ(NUNIT,9) IFL(I),ICOLUP(1,I),ICOLUP(2,I),SP(1,I),SP(2,I)
     $       ,SP(3,I),SM(I)
      ENDDO
 8    FORMAT(I8,1X,2(I4,1X),F10.3)
 9    FORMAT(I8,1X,2(I4,1X),4(1X,F10.3))
c
C     START PROCESSING INPUT DATA
C     
C     SCALES AND WEIGHTS
      SCALUP=DBLE(SQ)
      IF(IDWTUP.EQ.3) THEN
	XWGTUP=DBLE(SWGTRES) !AVGWGT
      ELSE
        WRITE(*,*) 'ONLY UNWEIGHTED EVENTS ACCEPTED AS INPUT, STOP'
        STOP
      ENDIF
c     
c---  incoming lines
      do 100 i=1,2
        iup=i
        idup(iup)=ifl(i)
        istup(iup)=-1
        mothup(1,iup)=0
        mothup(2,iup)=0
        pup(1,iup)=0.
        pup(2,iup)=0.
        pup(3,iup)=dble(Sp(3,iup))
        pup(4,iup)=abs(pup(3,iup))
        pup(5,iup)=0d0
 100  continue
c---  outgoing lines
      do 110 i=3,npart
        iup=i
        idup(iup)=ifl(i)
        istup(iup)=1
        mothup(1,iup)=1
        mothup(2,iup)=2
        pup(1,iup)=dble(Sp(1,i))
        pup(2,iup)=dble(Sp(2,i))
        pup(3,iup)=dble(Sp(3,i))
        pup(5,iup)=dble(Sm(i))
        tmp=(pup(5,iup)**2+pup(1,iup)**2+pup(2,iup)**2+pup(3,iup)**2)
        pup(4,iup)=sqrt(tmp)
 110  continue
c     
      nup=npart
c---  set up colour structure labels
      Do iup=1,nup
        if(icolup(1,iup).ne.0) icolup(1,iup)=icolup(1,iup)+500
        if(icolup(2,iup).ne.0) icolup(2,iup)=icolup(2,iup)+500
      Enddo
c
c     and now consider assignements specific to individual hard
C     processes
c     
c---  W/Z/gamma b bbar + jets, or W/Z + jets
      if (ihrd.le.4.or.ihrd.eq.10.OR.IHRD.EQ.14.OR.IHRD.EQ.15.OR.IHRD.EQ
     $     .17) then
        iwch=0
        if((ihrd.ne.1.and.ihrd.ne.2.and.ihrd.ne.15).or.
     $     (ihrd.eq.1.and.itdec.eq.0).or.
     $     (ihrd.eq.2.and.itdec.eq.0).or.
     $     (ihrd.eq.15.and.itdec.eq.0)) then
c     here only stable quarks are involved
           do iup=nup-1,nup
              mothup(1,iup)=nup+1
              mothup(2,iup)=0
              if(ihrd.ne.2) iwch=iwch-mod(idup(iup),2)
c     electron+nubar -> 11 + (-12) => -(1)+0 = -1  => W-
c     positron+nu    -> -11+ 12    => -(-1)+0 = -1 => W+
c     u dbar -> 2 -1  => 0 -(-1) = 1 => W+
c     c dbar -> 4 -1  => W+
c etc.
           enddo
           iup=nup+1
           If (iwch.gt.0) then
              idup(iup)=24
           Elseif (iwch.lt.0) then
              idup(iup)=-24
           Else
              idup(iup)=23
           Endif
           istup(iup)=2
           mothup(1,iup)=1
           mothup(2,iup)=2
           tmp=pup(4,iup-2)+pup(4,iup-1)
           pup(4,iup)=tmp
           tmp=tmp**2
           do i=1,3
              pup(i,iup)=pup(i,iup-2)+pup(i,iup-1)
              tmp=tmp-pup(i,iup)**2
           enddo
           pup(5,iup)=sqrt(tmp)
           nup=nup+1
           icolup(1,nup)=0
           icolup(2,nup)=0
        else
c
c first: the W from lepton-neutrino is reconstructed
c
           ijets= 0
           do iup=1,nup
              if(abs(idup(iup)).ne.6) then
                 ijets= ijets+1
              else
                 ijets= ijets-2
                 goto 150
              endif
           enddo
 150       continue
           do iup=nup-4-1,nup-4
              mothup(1,iup)=nup+1
              mothup(2,iup)=0
              if(ihrd.ne.2) iwch=iwch-mod(idup(iup),2)
c     electron+nubar -> 11 + (-12) => -(1)+0 = -1  => W-
c     positron+nu    -> -11+ 12    => -(-1)+0 = -1 => W+
c     u dbar -> 2 -1  => 0 -(-1) = 1 => W+
c     c dbar -> 4 -1  => W+
c etc.
           enddo
           iup=nup+1
           If (iwch.gt.0) then
              idup(iup)=24
           Elseif (iwch.lt.0) then
              idup(iup)=-24
           Else
              idup(iup)=23
           Endif
           istup(iup)=2
           mothup(1,iup)=1
           mothup(2,iup)=2
           tmp=pup(4,iup-4-2)+pup(4,iup-4-1)
           pup(4,iup)=tmp
           tmp=tmp**2
           do i=1,3
              pup(i,iup)=pup(i,iup-4-2)+pup(i,iup-4-1)
              tmp=tmp-pup(i,iup)**2
           enddo
           pup(5,iup)=sqrt(tmp)
           nup=nup+1
           icolup(1,nup)=0
           icolup(2,nup)=0
c
c now reconstruct W's from t and tbar
c
           if(ihrd.ne.15) then
              istup(3)=2
              istup(4)=2
              if(ifl(3).eq.6) then
                 it(1)=3
                 itb(1)=4
              else
                 it(1)=4
                 itb(1)=3
              endif
           else
              istup(3+ijets)=2
              istup(4+ijets)=2
              if(ifl(3+ijets).eq.6) then
                 it(1)=3 + ijets
                 itb(1)=4 + ijets
              else
                 it(1)=4 + ijets
                 itb(1)=3 + ijets
              endif
           endif
c     reconstruct W's from decay products
           do iw=1,2
              iwdec=nup-1-5+2*iw
              iwup=nup+iw
              ibup=iwup+2
              iwch=0
              do iup=iwdec,iwdec+1
                 mothup(1,iup)=iwup
                 mothup(2,iup)=0
                 iwch=iwch-mod(idup(iup),2)
c     electron+nubar -> 11 + (-12) = -1 => W-
c     d + ubar -> 1 + (-2) = -1 => W-
c     positron+nu    -> -11+ 12    =  1 => W+
c     u + dbar -> 2 + (-1) = 1 => W+
              enddo
              If (iwch.gt.0) then
                 idup(iwup)=24
                 idup(ibup)=5
                 mothup(1,iwup)=it(1)
                 mothup(2,iwup)=0
                 mothup(1,ibup)=it(1)
                 mothup(2,ibup)=0
              Elseif (iwch.lt.0) then
                 idup(iwup)=-24
                 idup(ibup)=-5
                 mothup(1,iwup)=itb(1)
                 mothup(2,iwup)=0
                 mothup(1,ibup)=itb(1)
                 mothup(2,ibup)=0
              Endif
              istup(iwup)=2
              istup(ibup)=1
c     reconstruct W momentum
              tmp=pup(4,iwdec)+pup(4,iwdec+1)
              pup(4,iwup)=tmp
              tmp=tmp**2
              do i=1,3
                 pup(i,iwup)=pup(i,iwdec)+pup(i,iwdec+1)
                 tmp=tmp-pup(i,iwup)**2
              enddo
              pup(5,iwup)=sqrt(tmp)
c     reconstruct b momentum
              tmp=pup(4,mothup(1,iwup))-pup(4,iwup)
              pup(4,ibup)=tmp 
              tmp=tmp**2
              do i=1,3
                 pup(i,ibup)=pup(i,mothup(1,iwup))-pup(i,iwup)
                 tmp=tmp-pup(i,ibup)**2
              enddo
              pup(5,ibup)=sqrt(tmp)
              icolup(1,iwup)=0
              icolup(2,iwup)=0
              icolup(1,ibup)=icolup(1,mothup(1,iwup))
              icolup(2,ibup)=icolup(2,mothup(1,iwup))
           enddo
c     stop
           nup=nup+4

c           print*,'  '
c           print*,'  '
c           do iup=1,nup
c              print*,iup,'idup(iup)= ',idup(iup)
c              print*,iup,'istup(iup)= ',istup(iup)
c              print*,iup,'mothup(1,iup)= ',mothup(1,iup)
c              print*,iup,'mothup(2,iup)= ',mothup(2,iup)
c           enddo
c           print*,'  '
c           print*,'  '



        endif

c---  nW + mZ + kH + jets
      elseif (ihrd.eq.5) then
c     find first gauge bosons
        ivstart=0
        ivend=0
        do i=1,npart
          if(abs(idup(i)).eq.24.or.idup(i).eq.23) then
            istup(i)=2
            if(ivstart.eq.0) ivstart=i
            ivend=i+1
          endif
        enddo
        nvb=ivend-ivstart
c     decay products pointers, starting from the end
        do i=1,nvb
          mothup(1,npart-2*i+2)=ivend-i
          mothup(1,npart-2*i+1)=ivend-i
          mothup(2,npart-2*i+2)=0
          mothup(2,npart-2*i+1)=0
        enddo
c---  t tbar + jets and t tbar + gamma + jets
c     t tb jets f fbar f fbar W+ b W- bbar
      elseif( (ihrd.eq.6.or.ihrd.eq.16) .and.abs(ifl(3)).eq.6) then
         if(itdec.ne.0) then
c         reset top status codes
            istup(3)=2
            istup(4)=2
            if(ifl(3).eq.6) then
               it(1)=3
               itb(1)=4
            else
               it(1)=4
               itb(1)=3
            endif
c         reconstruct W's from decay products
            do iw=1,2
               iwdec=nup-5+2*iw
               iwup=nup+iw
               ibup=iwup+2
               iwch=0
               do iup=iwdec,iwdec+1
                  mothup(1,iup)=iwup
                  mothup(2,iup)=0
                  iwch=iwch-mod(idup(iup),2)
c         electron+nubar -> 11 + (-12) = -1 => W-
c         d + ubar -> 1 + (-2) = -1 => W-
c         positron+nu    -> -11+ 12    =  1 => W+
c         u + dbar -> 2 + (-1) = 1 => W+
               enddo
               If (iwch.gt.0) then
                  idup(iwup)=24
                  idup(ibup)=5
                  mothup(1,iwup)=it(1)
                  mothup(2,iwup)=0
                  mothup(1,ibup)=it(1)
                  mothup(2,ibup)=0
               Elseif (iwch.lt.0) then
                  idup(iwup)=-24
                  idup(ibup)=-5
                  mothup(1,iwup)=itb(1)
                  mothup(2,iwup)=0
                  mothup(1,ibup)=itb(1)
                  mothup(2,ibup)=0
               Endif
               istup(iwup)=2
               istup(ibup)=1
c         reconstruct W momentum
               tmp=pup(4,iwdec)+pup(4,iwdec+1)
               pup(4,iwup)=tmp
               tmp=tmp**2
               do i=1,3
                  pup(i,iwup)=pup(i,iwdec)+pup(i,iwdec+1)
                  tmp=tmp-pup(i,iwup)**2
               enddo
               pup(5,iwup)=sqrt(tmp)
c         reconstruct b momentum
               tmp=pup(4,mothup(1,iwup))-pup(4,iwup)
               pup(4,ibup)=tmp 
               tmp=tmp**2
               do i=1,3
                  pup(i,ibup)=pup(i,mothup(1,iwup))-pup(i,iwup)
                  tmp=tmp-pup(i,ibup)**2
               enddo
               pup(5,ibup)=sqrt(tmp)
               icolup(1,iwup)=0
               icolup(2,iwup)=0
               icolup(1,ibup)=icolup(1,mothup(1,iwup))
               icolup(2,ibup)=icolup(2,mothup(1,iwup))
            enddo
c         stop
            nup=nup+4
         endif
c---  b bbar t tbar + jets and t tbar t tbar + jets
      elseif (ihrd.eq.7) then !and.abs(ifl(5)).eq.6) then
         if(abs(ifl(3)).eq.6.and.abs(ifl(5)).eq.5) then
            write(*,*) 'abs(ifl(3))= 6 and abs(ifl(5))= 5 not possible'
            stop
         endif
         if(abs(ifl(3)).eq.5) then
c       b bbar t tbar + jets
            if(itdec.ne.0) then
c         reset top status codes
               istup(5)=2
               istup(6)=2
               if(ifl(5).eq.6) then
                  it(1)=5
                  itb(1)=6
               else
                  it(1)=6
                  itb(1)=5
               endif
c         reconstruct W's from decay products
               do iw=1,2
                  iwdec=nup-5+2*iw
                  iwup=nup+iw
                  ibup=iwup+2
                  iwch=0
                  do iup=iwdec,iwdec+1
                     mothup(1,iup)=iwup
                     mothup(2,iup)=0
                     iwch=iwch-mod(idup(iup),2)
c         electron+nubar -> 11 + (-12) = -1 => W-
c         d + ubar -> 1 + (-2) = -1 => W-
c         positron+nu    -> -11+ 12    =  1 => W+
c         u + dbar -> 2 + (-1) = 1 => W+
                  enddo
                  If (iwch.gt.0) then
                     idup(iwup)=24
                     idup(ibup)=5
                     mothup(1,iwup)=it(1)
                     mothup(2,iwup)=0
                     mothup(1,ibup)=it(1)
                     mothup(2,ibup)=0
                  Elseif (iwch.lt.0) then
                     idup(iwup)=-24
                     idup(ibup)=-5
                     mothup(1,iwup)=itb(1)
                     mothup(2,iwup)=0
                     mothup(1,ibup)=itb(1)
                     mothup(2,ibup)=0
                  Endif
                  istup(iwup)=2
                  istup(ibup)=1
c     reconstruct W momentum
                  tmp=pup(4,iwdec)+pup(4,iwdec+1)
                  pup(4,iwup)=tmp
                  tmp=tmp**2
                  do i=1,3
                     pup(i,iwup)=pup(i,iwdec)+pup(i,iwdec+1)
                     tmp=tmp-pup(i,iwup)**2
                  enddo
                  pup(5,iwup)=sqrt(tmp)
c     reconstruct b momentum
                  tmp=pup(4,mothup(1,iwup))-pup(4,iwup)
                  pup(4,ibup)=tmp 
                  tmp=tmp**2
                  do i=1,3
                     pup(i,ibup)=pup(i,mothup(1,iwup))-pup(i,iwup)
                     tmp=tmp-pup(i,ibup)**2
                  enddo
                  pup(5,ibup)=sqrt(tmp)
                  icolup(1,iwup)=0
                  icolup(2,iwup)=0
                  icolup(1,ibup)=icolup(1,mothup(1,iwup))
                  icolup(2,ibup)=icolup(2,mothup(1,iwup))
               enddo
c         stop
               nup=nup+4
            endif
         elseif(abs(ifl(3)).eq.6) then
c       t tbar t tbar + jets
            if(itdec.ne.0) then
c         reset top status codes
               istup(3)=2
               istup(4)=2
               istup(5)=2
               istup(6)=2
               if(ifl(3).eq.6) then
                  it(1)=3
                  itb(1)=4
               else
                  it(1)=4
                  itb(1)=3
               endif
               if(ifl(5).eq.6) then
                  it(2)=5
                  itb(2)=6
               else
                  it(2)=6
                  itb(2)=5
               endif
c         reconstruct W's from decay products
               do iw=1,4
                  iwdec=nup-9+2*iw
                  iwup=nup+iw
                  ibup=iwup+4
                  iwch=0
                  do iup=iwdec,iwdec+1
                     mothup(1,iup)=iwup
                     mothup(2,iup)=0
                     iwch=iwch-mod(idup(iup),2)
c         electron+nubar -> 11 + (-12) = -1 => W-
c         d + ubar -> 1 + (-2) = -1 => W-
c         positron+nu    -> -11+ 12    =  1 => W+
c         u + dbar -> 2 + (-1) = 1 => W+
                  enddo
                  If (iwch.gt.0) then
                     idup(iwup)=24
                     idup(ibup)=5
                     if(iw.lt.3) then
                        mothup(1,iwup)=it(1)
                        mothup(2,iwup)=0
                        mothup(1,ibup)=it(1)
                        mothup(2,ibup)=0
                     else
                        mothup(1,iwup)=it(2)
                        mothup(2,iwup)=0
                        mothup(1,ibup)=it(2)
                        mothup(2,ibup)=0
                     endif
                  Elseif (iwch.lt.0) then
                     idup(iwup)=-24
                     idup(ibup)=-5
                     if(iw.lt.3) then
                        mothup(1,iwup)=itb(1)
                        mothup(2,iwup)=0
                        mothup(1,ibup)=itb(1)
                        mothup(2,ibup)=0
                     else
                        mothup(1,iwup)=itb(2)
                        mothup(2,iwup)=0
                        mothup(1,ibup)=itb(2)
                        mothup(2,ibup)=0
                     endif
                  endif
                  istup(iwup)=2
                  istup(ibup)=1
c     reconstruct W momentum
                  tmp=pup(4,iwdec)+pup(4,iwdec+1)
                  pup(4,iwup)=tmp
                  tmp=tmp**2
                  do i=1,3
                     pup(i,iwup)=pup(i,iwdec)+pup(i,iwdec+1)
                     tmp=tmp-pup(i,iwup)**2
                  enddo
                  pup(5,iwup)=sqrt(tmp)
c     reconstruct b momentum
                  tmp=pup(4,mothup(1,iwup))-pup(4,iwup)
                  pup(4,ibup)=tmp 
                  tmp=tmp**2
                  do i=1,3
                     pup(i,ibup)=pup(i,mothup(1,iwup))-pup(i,iwup)
                     tmp=tmp-pup(i,ibup)**2
                  enddo
                  pup(5,ibup)=sqrt(tmp)
                  icolup(1,iwup)=0
                  icolup(2,iwup)=0
                  icolup(1,ibup)=icolup(1,mothup(1,iwup))
                  icolup(2,ibup)=icolup(2,mothup(1,iwup))
               enddo
c         stop
               nup=nup+8
            endif
         endif
c---  H t tbar + jets
c     H t tb jets f fbar f fbar W+ b W- bbar
      elseif (ihrd.eq.8.and.abs(ifl(4)).eq.6) then
         if(itdec.ne.0) then
c         reset top status codes
            istup(4)=2
            istup(5)=2
            if(ifl(4).eq.6) then
               it(1)=4
               itb(1)=5
            else
               it(1)=5
               itb(1)=4
            endif
c         reconstruct W's from decay products
            do iw=1,2
               iwdec=nup-5+2*iw
               iwup=nup+iw
               ibup=iwup+2
               iwch=0
               do iup=iwdec,iwdec+1
                  mothup(1,iup)=iwup
                  mothup(2,iup)=0
                  iwch=iwch-mod(idup(iup),2)
c         electron+nubar -> 11 + (-12) = -1 => W-
c         d + ubar -> 1 + (-2) = -1 => W-
c         positron+nu    -> -11+ 12    =  1 => W+
c         u + dbar -> 2 + (-1) = 1 => W+
               enddo
               If (iwch.gt.0) then
                  idup(iwup)=24
                  idup(ibup)=5
                  mothup(1,iwup)=it(1)
                  mothup(2,iwup)=0
                  mothup(1,ibup)=it(1)
                  mothup(2,ibup)=0
               elseif (iwch.lt.0) then
                  idup(iwup)=-24
                  idup(ibup)=-5
                  mothup(1,iwup)=itb(1)
                  mothup(2,iwup)=0
                  mothup(1,ibup)=itb(1)
                  mothup(2,ibup)=0
               endif
               istup(iwup)=2
               istup(ibup)=1
c         reconstruct W momentum
               tmp=pup(4,iwdec)+pup(4,iwdec+1)
               pup(4,iwup)=tmp
               tmp=tmp**2
               do i=1,3
                  pup(i,iwup)=pup(i,iwdec)+pup(i,iwdec+1)
                  tmp=tmp-pup(i,iwup)**2
               enddo
               pup(5,iwup)=sqrt(tmp)
c         reconstruct b momentum
               tmp=pup(4,mothup(1,iwup))-pup(4,iwup)
               pup(4,ibup)=tmp 
               tmp=tmp**2
               do i=1,3
                  pup(i,ibup)=pup(i,mothup(1,iwup))-pup(i,iwup)
                  tmp=tmp-pup(i,ibup)**2
               enddo
               pup(5,ibup)=sqrt(tmp)
               icolup(1,iwup)=0
               icolup(2,iwup)=0
               icolup(1,ibup)=icolup(1,mothup(1,iwup))
               icolup(2,ibup)=icolup(2,mothup(1,iwup))
            enddo
c         stop
            nup=nup+4
         endif
c---  SINGLE TOP
c     Input: T  
c     output: jets t b w f fbar t b w f fbar 
      elseif (ihrd.eq.13) then
         if(itdec.ne.0) then
            nw=1
            if(itopprc.ge.3) nw=2
c     assign mass to the incoming bottom quark, if required
            DO I=1,2
               IF(ABS(IFL(I)).EQ.5) THEN
                  IUP=I
                  PUP(5,IUP)=mb
                  PUP(4,IUP)=SQRT(PUP(3,IUP)**2+PUP(5,IUP)**2)
               ENDIF
            ENDDO
            istup(3)=2
            it(1)=0
            itb(1)=0
            if(ifl(3).eq.6) then
               it(1)=3
            elseif(ifl(3).eq.-6) then
               itb(1)=3
            else
               write(*,*) 'wrong assumption about top position, stop'
               stop
            endif
c     
c     TOP DECAY
c     reconstruct W's from decay products
c     
c     iwdec: 1st W decay product.
            if(nw.eq.1) then
               iwdec=nup-1
            elseif(nw.eq.2) then
               iwdec=nup-3
            endif
c     put W and b  at the end
            iwup=nup+1
            ibup=iwup+1
c     
            iwch=0
            do iup=iwdec,iwdec+1
               mothup(1,iup)=iwup
               mothup(2,iup)=0
               iwch=iwch-mod(idup(iup),2)
c     electron+nubar -> 11 + (-12) = -1 => W-
c     d + ubar -> 1 + (-2) = -1 => W-
c     positron+nu    -> -11+ 12    =  1 => W+
c     u + dbar -> 2 + (-1) = 1 => W+
            enddo
            If (iwch.gt.0) then
               idup(iwup)=24
               idup(ibup)=5
               mothup(1,iwup)=it(1)
               mothup(2,iwup)=0
               mothup(1,ibup)=it(1)
               mothup(2,ibup)=0
            Elseif (iwch.lt.0) then
               idup(iwup)=-24
               idup(ibup)=-5
               mothup(1,iwup)=itb(1)
               mothup(2,iwup)=0
               mothup(1,ibup)=itb(1)
               mothup(2,ibup)=0
            Endif
            istup(iwup)=2
            istup(ibup)=1
c     reconstruct W momentum
            tmp=pup(4,iwdec)+pup(4,iwdec+1)
            pup(4,iwup)=tmp
            tmp=tmp**2
            do i=1,3
               pup(i,iwup)=pup(i,iwdec)+pup(i,iwdec+1)
               tmp=tmp-pup(i,iwup)**2
            enddo
            pup(5,iwup)=sqrt(tmp)
c     reconstruct b momentum
            tmp=pup(4,mothup(1,iwup))-pup(4,iwup)
            pup(4,ibup)=tmp 
            tmp=tmp**2
            do i=1,3
               pup(i,ibup)=pup(i,mothup(1,iwup))-pup(i,iwup)
               tmp=tmp-pup(i,ibup)**2
            enddo
c     write(*,*) (pup(i,ibup),i=1,4),sqrt((tmp))
            pup(5,ibup)=sqrt(tmp)
            icolup(1,iwup)=0
            icolup(2,iwup)=0
            icolup(1,ibup)=icolup(1,mothup(1,iwup))
            icolup(2,ibup)=icolup(2,mothup(1,iwup))
c     
            nup=nup+2
            if(nw.eq.2) then
c     
c     W DECAY
c     
c     iwdec: 1st W decay product. 
               iwdec=nup-3
c     iwup: location of the W in the event record
               iwup=nup-6
               iwch=0
               do iup=iwdec,iwdec+1
                  mothup(1,iup)=iwup
                  mothup(2,iup)=0
                  iwch=iwch-mod(idup(iup),2)
c     electron+nubar -> 11 + (-12) = -1 => W-
c     d + ubar -> 1 + (-2) = -1 => W-
c     positron+nu    -> -11+ 12    =  1 => W+
c     u + dbar -> 2 + (-1) = 1 => W+
               enddo
               istup(iwup)=2
               icolup(1,iwup)=0
               icolup(2,iwup)=0
            endif
         endif
c---  hvy N production 
      elseif (ihrd.eq.19) then
         if (idup(3)*idup(4).lt.0) then
          nup= nup+3
C     the s-channel W or Z
          iwch=0
          iwup=7
          do iup=1,2
             iwch=iwch-mod(idup(iup),2)
c     electron+nubar -> 11 + (-12) = -1 => W-
c     d + ubar -> 1 + (-2) = -1 => W-
c     positron+nu    -> -11+ 12    =  1 => W+
c     u + dbar -> 2 + (-1) = 1 => W+
          enddo
          If (iwch.gt.0) then
             idup(iwup)=24
          Elseif (iwch.lt.0) then
             idup(iwup)=-24
          elseif (iwch.eq.0) then
             idup(iwup)= 23
          else
             print*,'HVYN: Wrong EW boson assignement in alpsho.f, stop'
             stop
          endif
          istup(iwup)= +2
          mothup(1,7)= 1
          mothup(2,7)= 2
          icolup(1,7)= 0
          icolup(2,7)= 0
          do i=1,4
             pup(i,7)= pup(i,3)+pup(i,4)+pup(i,5)+pup(i,6)
          enddo 
          pup(5,7)= sqrt(dabs(pup(4,7)**2-pup(1,7)**2
     &         -pup(2,7)**2-pup(3,7)**2))
C the Heavy N
          istup(8)= +2
          idup(8) = 0  
          mothup(1,8)= 7
          mothup(2,8)= 7
          icolup(1,8)= 0
          icolup(2,8)= 0
          do i=1,4
            pup(i,8)= pup(i,4)+pup(i,5)+pup(i,6)
          enddo 
          pup(5,8)= sqrt(dabs(pup(4,8)**2-pup(1,8)**2
     &                       -pup(2,8)**2-pup(3,8)**2))
C the second W,or Z or H
          iwch=0
          iwup=9
          do iup=5,6
            iwch=iwch-mod(idup(iup),2)
c     electron+nubar -> 11 + (-12) = -1 => W-
c     d + ubar -> 1 + (-2) = -1 => W-
c     positron+nu    -> -11+ 12    =  1 => W+
c     u + dbar -> 2 + (-1) = 1 => W+
          enddo
          If (iwch.gt.0) then
            idup(iwup)=24
          Elseif (iwch.lt.0) then
            idup(iwup)=-24
          elseif(nz.eq.1) then
            idup(iwup)= +23
          elseif (nh.eq.1) then
            idup(iwup)= +25
          else
            print*,'error2 in alpsho.f!'
            stop
          endif
          istup(iwup)= +2
          mothup(1,9)= 8
          mothup(2,9)= 8
          icolup(1,9)= 0
          icolup(2,9)= 0
          do i=1,4
            pup(i,9)= pup(i,5)+pup(i,6)
          enddo 
          pup(5,9)= sqrt(dabs(pup(4,9)**2-pup(1,9)**2
     &                       -pup(2,9)**2-pup(3,9)**2))
          mothup(1,3)= 7 
          mothup(2,3)= 7
          mothup(1,4)= 8
          mothup(2,4)= 8
          mothup(1,5)= 9
          mothup(2,5)= 9
          mothup(1,6)= 9
          mothup(2,6)= 9
        else
          nup= nup+2
C the s-channel W
          iwch=0
          iwup=7
          do iup=1,2
            iwch=iwch-mod(idup(iup),2)
c     electron+nubar -> 11 + (-12) = -1 => W-
c     d + ubar -> 1 + (-2) = -1 => W-
c     positron+nu    -> -11+ 12    =  1 => W+
c     u + dbar -> 2 + (-1) = 1 => W+
          enddo
          If (iwch.gt.0) then
            idup(iwup)=24
          Elseif (iwch.lt.0) then
            idup(iwup)=-24
          else
            print*,'error3 in alpsho.f!'
            stop
          endif
          istup(iwup)= +2
          mothup(1,7)= 1
          mothup(2,7)= 2
          icolup(1,7)= 0
          icolup(2,7)= 0
          do i=1,4
            pup(i,7)= pup(i,3)+pup(i,4)+pup(i,5)+pup(i,6)
          enddo 
          pup(5,7)= sqrt(dabs(pup(4,7)**2-pup(1,7)**2
     &                       -pup(2,7)**2-pup(3,7)**2))
C the second W
          iwch=0
          iwup=8
          do iup=5,6
            iwch=iwch-mod(idup(iup),2)
c     electron+nubar -> 11 + (-12) = -1 => W-
c     d + ubar -> 1 + (-2) = -1 => W-
c     positron+nu    -> -11+ 12    =  1 => W+
c     u + dbar -> 2 + (-1) = 1 => W+
          enddo
          If (iwch.gt.0) then
            idup(iwup)=24
          Elseif (iwch.lt.0) then
            idup(iwup)=-24
          else
            print*,'error4 in alpsho.f!'
            stop
          endif
          istup(iwup)= +2
          mothup(1,8)= 7
          mothup(2,8)= 7
          icolup(1,8)= 0
          icolup(2,8)= 0
          do i=1,4
            pup(i,8)= pup(i,5)+pup(i,6)
          enddo 
          pup(5,8)= sqrt(dabs(pup(4,8)**2-pup(1,8)**2
     &                       -pup(2,8)**2-pup(3,8)**2))
          mothup(1,3)= 7 
          mothup(2,3)= 7
          mothup(1,4)= 7
          mothup(2,4)= 7
          mothup(1,5)= 8
          mothup(2,5)= 8
          mothup(1,6)= 8
          mothup(2,6)= 8
        endif  
      endif
c     herwig debugging:
c      call HWUPUP
      return
c     
c     
c     end of file
 500  ist=1
c     error reading file
 501  if(ist.eq.0) ist=2
C     RESET CROSS-SECTION INFORMATION FOR END OF RUN AND FINALIZE
      IF(IST.GT.0) THEN
        CALL AHVSTA(NEVPROC-1)
        CALL ALSFIN
      ENDIF
      close(Nunit)
      close(NunitOut)
      END

C----------------------------------------------------------------------
      SUBROUTINE UPALHE
C----------------------------------------------------------------------
c     Puts Alpgen event from .lhe format into GUPI common block HEPEU
c----------------------------------------------------------------------
      implicit none
      INCLUDE 'alpsho.inc'
C--   GUP Event common block
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &     IDUP(MAXNUP),ISTUP(MAXNUP),MOTHUP(2,MAXNUP),
     &     ICOLUP(2,MAXNUP),PUP(5,MAXNUP),VTIMUP(MAXNUP),
     &     SPINUP(MAXNUP)
C--   GUP Run common block
      INTEGER MAXPUP
      PARAMETER(MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON /HEPRUP/ IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &     IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),
     &     XMAXUP(MAXPUP),LPRUP(MAXPUP)
c     
c     local variables
      INTEGER INIT
      DATA INIT/0/
      INTEGER IST,I,J
      INTEGER NEV,IPROC,NEVPROC
      DATA NEVPROC/0/
C...Lines to read in assumed never longer than 200 characters. 
      INTEGER IBEG,MAXLEN
      PARAMETER (MAXLEN=200)
      CHARACTER*(MAXLEN) STRING
C...Format for reading lines.
      CHARACTER*6 STRFMT
      DATA STRFMT/'(A200)'/
C
      IST=0
      NEVPROC=NEVPROC+1
      IF(NEVMAX.GT.0.AND.NEVPROC.GT.NEVMAX) GOTO 500
      CALL ALSHER(I)
C...Loop until finds line beginning with "<event>" or "<event ". 
  100 READ(NUNIT,STRFMT,END=500,ERR=501) STRING
      IBEG=0
  110 IBEG=IBEG+1
C...Allow indentation.
      IF(STRING(IBEG:IBEG).EQ.' '.AND.IBEG.LT.MAXLEN-6) GOTO 110 
      IF(STRING(IBEG:IBEG+6).NE.'<event>'.AND.
     &STRING(IBEG:IBEG+6).NE.'<event ') GOTO 100
C
C...Read first line of event info.
      READ(NUNIT,*,END=500,ERR=501) NUP,IDPRUP,XWGTUP,SCALUP,
     &AQEDUP,AQCDUP
C...Read NUP subsequent lines with information on each particle.
      DO 200 I=1,NUP
        READ(NUNIT,*,END=500,ERR=501) IDUP(I),ISTUP(I),
     &  MOTHUP(1,I),MOTHUP(2,I),ICOLUP(1,I),ICOLUP(2,I),
     &  (PUP(J,I),J=1,5),SPINUP(I)
c     &  (PUP(J,I),J=1,5),VTIMUP(I),SPINUP(I)
  200 CONTINUE

C     SCALES AND WEIGHTS
C	XWGTUP=DBLE(SWGTRES) !AVGWGT
c     herwig debugging:
c      call HWUPUP
      RETURN
c     
c     
c     end of file
 500  ist=1
c     error reading file
 501  if(ist.eq.0) ist=2
C     RESET CROSS-SECTION INFORMATION FOR END OF RUN AND FINALIZE
      IF(IST.GT.0) THEN
        CALL AHVSTA(NEVPROC-1)
        CALL ALSFIN
      ENDIF
      close(Nunit)
      close(NunitOut)
      END

      SUBROUTINE AHVSTA(NEVPROC)
C     VETO STATISTICS
      IMPLICIT NONE
      include 'alpsho.inc'
      INTEGER NEVPROC,I,N,NTOT,NQTOT
      NTOT=0
      IF(IHVYV.EQ.1) then
        WRITE(NUNITOUT,*) '**** CHARM VETOES:'
C
        DO I=1,100
          N=NVETOC(I)
          NTOT=NTOT+N
          IF(N.GT.0) WRITE(NUNITOUT,*) 'NCH= ',I,' VETOED EVNTS:',N
     $         ,' VETOED FRAC:',REAL(N)/REAL(NEVPROC)
        ENDDO
C     
        WRITE(NUNITOUT,*) '**** BOTTOM VETOES:'
        DO I=1,100
          N=NVETOB(I)
          NTOT=NTOT+N
          IF(N.GT.0) WRITE(NUNITOUT,*) 'NBOT=',I,' VETOED EVNTS:',N
     $         ,' VETOED FRAC:',REAL(N)/REAL(NEVPROC)
        ENDDO
        NQTOT=NTOT
      ENDIF
C
      WRITE(NUNITOUT,*) '**** NJET<NPART VETOES:'
      N=NVETONOJ
      NTOT=NTOT+N
      IF(N.GT.0) WRITE(NUNITOUT,*) ' VETOED EVNTS:',N
     $     ,' VETOED FRAC:',REAL(N)/REAL(NEVPROC)
C
      WRITE(NUNITOUT,*) '**** MATCHING VETOES:'
      DO I=0,10
        N=NVETOMCH(I)
        NTOT=NTOT+N
        IF(N.GT.0) WRITE(NUNITOUT,*) 'NMATCH=',I,'< NPARTON'
     $       ,' VETOED EVNTS:',N,' VETOED FRAC:',REAL(N)/REAL(NEVPROC)
      ENDDO
C
      WRITE(NUNITOUT,*) '**** EXTRA-JET VETOES:'
      DO I=1,10
        N=NVETOEXC(I)
        NTOT=NTOT+N
        IF(N.GT.0) WRITE(NUNITOUT,*) 'EXTRA JETS=',I
     $       ,' VETOED EVNTS:',N,' VETOED FRAC:',REAL(N)/REAL(NEVPROC)
      ENDDO
C
      WRITE(NUNITOUT,*) '**** EXTRA-JETS FROM HEAVY QUARKS VETOES:'
      DO I=1,10
        N=NVETOHVY(I)
        NTOT=NTOT+N
        IF(N.GT.0) WRITE(NUNITOUT,*) 'EXTRA JETS=',I
     $       ,' VETOED EVNTS:',N,' VETOED FRAC:',REAL(N)/REAL(NEVPROC)
      ENDDO
      IF(NTOT.GT.0) WRITE(NUNITOUT,*) 'TOTAL VETOED EVENTS=',NTOT
     $       ,' VETOED FRAC:',REAL(NTOT)/REAL(NEVPROC)
      IF(NQTOT.GT.0) WRITE(NUNITOUT,*) 'FRACTION SURVIVING HF CUTS='
     $     ,REAL(NQTOT)/REAL(NEVPROC-NTOT+NQTOT)
      END

CDECK  ID>, STRNUM.
*CMZ :-        -13/02/02  07.20.46  by  Peter Richardson
*-- Author :    Michelangelo Mangano
C----------------------------------------------------------------------
      subroutine strnumH(string,num)
C----------------------------------------------------------------------
c- writes the number num on the string string starting at the blank
c- following the last non-blank character
C----------------------------------------------------------------------
      character * (*) string
      character * 20 tmp
      l = len(string)
      write(tmp,'(i15)')num
      j=1
      dowhile(tmp(j:j).eq.' ')
        j=j+1
      enddo
      ipos = istrlH(string)
      ito = ipos+1+(15-j)
      if(ito.gt.l) then
         write(*,*)'error, string too short'
         write(*,*) string
         stop
      endif
      string(ipos+1:ito)=tmp(j:)
      end

      function istrlH(string)
c returns the position of the last non-blank character in string
      character * (*) string
      i = len(string)
      dowhile(i.gt.0.and.string(i:i).eq.' ')
         i=i-1
      enddo
      istrlH= i
      end
CDECK  ID>, STRCATH.
*CMZ :-        -13/02/02  07.20.46  by  Peter Richardson
*-- Author :    Michelangelo Mangano
C----------------------------------------------------------------------
      subroutine strcatH(str1,str2,str)
C----------------------------------------------------------------------
c concatenates str1 and str2 into str. Ignores trailing blanks of str1,str2
C----------------------------------------------------------------------
      character *(*) str1,str2,str
      l1=istrlH(str1)
      l2=istrlH(str2)
      l =len(str)
      if(l.lt.l1+l2) then
        write(*,*) str1,str2
          write(*,*) 'error: l1+l2>l in strcatH'
          write(*,*) 'l1=',l1,' str1=',str1
          write(*,*) 'l2=',l2,' str2=',str2
          write(*,*) 'l=',l
          stop
      endif
      if(l1.ne.0) str(1:l1)=str1(1:l1)
      if(l2.ne.0) str(l1+1:l1+l2)=str2(1:l2)
      if(l1+l2+1.le.l) str(l1+l2+1:l)= ' '
      end
c-------------------------------------------------------------------
      subroutine pdfconvH(nin,nout,type)
c-------------------------------------------------------------------
c converts ALPHA convention for PDF namings to hvqpdf conventions
      implicit none
      integer nin,nout
      character*25 type
      character*25 pdftyp(20,2)
      data pdftyp/
c cteq sets
     $     'CTEQ4M ','CTEQ4L ','CTEQ4HJ',
     $     'CTEQ5M ','CTEQ5L ','CTEQ5HJ',
     $     'CTEQ6M ','CTEQ6L ',12*' ',
C MRST SETS
     $     'MRST99 ',
     $     'MRST01; as=0.119','MRST01; as=0.117','MRST01; as=0.121'
     $     ,'MRST01J; as=0.121','MRST02LO',14*' '/
      integer pdfmap(20,2)
      data pdfmap/
     $   81,83,88,   101,103, 104,   131,133, 12*0,
     $  111,  185,186,187,188,189,   14*0/
c
      nout=pdfmap(mod(nin ,100),1+nin /100)
      type=pdftyp(mod(nin ,100),1+nin /100)
      
      end

      SUBROUTINE HVQVETO(NPFST,NPLST,NJLAST,NC,NB)
C...  REJECT EVENTS WITH UNWANTED HVQ'S FROM SHOWER/PL GENERATION
      INTEGER NMXHEP,NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
      PARAMETER (NMXHEP=4000)
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      DOUBLE PRECISION PHEP,VHEP
      SAVE /HEPEVT/
      LOGICAL FOUND
      INTEGER NPFST,NPLST,NJLAST,NC,NB,IST
C
      FOUND=.FALSE.
      NC=0
      NB=0
C      write(*,*) 'enter, nhep=',nhep,' istlo..=',istlo,isthi,istop
      DO 200 IHEP=1,NHEP
        IST=ISTHEP(IHEP)
        IF (IST.EQ.NJLAST) FOUND=.TRUE.
        IF (IST.GE.NPFST.AND.IST.LE.NPLST.AND..NOT.FOUND) THEN
          ID=ABS(IDHEP(IHEP))
C          write(*,*) ID
          IF (ID.EQ.4) NC=NC+1
          IF (ID.EQ.5) NB=NB+1
        ENDIF
 200  CONTINUE
      IF(NC.GT.NCMAX.OR.NB.GT.NBMAX) IFLAG=1
      END

C-----------------------------------------------------------------------
      SUBROUTINE FJHEVT_M(ISTLO,ISTHI,ISTOP,P,ID,NPART)
C     LABEL ALL PARTICLES WITH STATUS BETWEEN ISTLO AND ISTHI (UNTIL A
C     PARTICLE WITH STATUS ISTOP IS FOUND) AS FINAL-STATE, CALL CALSIM_M
C     AND THEN PUT LABELS BACK TO NORMAL
C
C     THIS VERSION LEAVES OUT PARTICLES THAT POINT BACK TO MOTHERS TO BE
C     LEFT OUT OF MATCHING
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER MAXNUP
      PARAMETER(MAXNUP=500)
      INTEGER IEND,INORAD
      COMMON/SHVETO/IEND,INORAD(MAXNUP)
C...HEPEVT commonblock.
      INTEGER NMXHEP,NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
      PARAMETER (NMXHEP=4000)
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      DOUBLE PRECISION PHEP,VHEP
      SAVE /HEPEVT/
      INTEGER ISTOLD(NMXHEP),IHEP,IST,ISTLO,ISTHI,ISTOP,IMO
      INTEGER NPART
      INTEGER NMAX
      PARAMETER (NMAX=1000)
      DOUBLE PRECISION P(4,NMAX)
      INTEGER ID(NMAX)
      LOGICAL FOUND
c      write(3,*) 'new event',nevhep
      FOUND=.FALSE.
      DO 10 IHEP=1,NHEP
        IST=ISTHEP(IHEP)
        ISTOLD(IHEP)=IST
        IF (IST.EQ.ISTOP) FOUND=.TRUE.
        IF (IST.GE.ISTLO.AND.IST.LE.ISTHI.AND..NOT.FOUND) THEN
C     FOUND A RADIATED PARTON, CHECK MOTHER
          IMO=IHEP
 1        IMO=JMOHEP(1,IMO)
          IF(IMO.EQ.IEND) THEN
C     PARENTHOOD OK
            IST=1
c            write(3,*) ihep,ist
            GOTO 9
          ENDIF
          IF(INORAD(IMO).EQ.1) THEN
C     PARTON COMES FROM A VETOED MOTHER
            IST=0
            GOTO 9
          ELSE
C     CHECK GRANDMOTHER
            GOTO 1
          ENDIF
        ELSE
          IST=0
        ENDIF
 9      ISTHEP(IHEP)=IST
 10   CONTINUE
C      CALL FJHEVT(ISTLO,ISTHI,ISTOP,P,NPART)
      CALL FJHEVT(1,1,200,P,ID,NPART)
      DO 20 IHEP=1,NHEP
        ISTHEP(IHEP)=ISTOLD(IHEP)
 20   CONTINUE
      END
C-----------------------------------------------------------------------
      SUBROUTINE FJHEVT_HVQ(ISTLO,ISTHI,ISTOP,P,ID,NPART)
C     LABEL ALL PARTICLES WITH STATUS BETWEEN ISTLO AND ISTHI (UNTIL A
C     PARTICLE WITH STATUS ISTOP IS FOUND) AS FINAL-STATE, CALL CALSIM_M
C     AND THEN PUT LABELS BACK TO NORMAL
C
C     THIS VERSION KEEPS ONLY ALL IST=1 PARTICLES REJECTED BY CALDEL AS
C     DAUGHTERS OF VETOED HEAVY-QUARK MOTHERS: JETS COMPLEMENTARY TO
C     THOSE RECONSTRUCTED BY CALDEL
C-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER MAXNUP
      PARAMETER(MAXNUP=500)
      INTEGER IEND,INORAD
      COMMON/SHVETO/IEND,INORAD(MAXNUP)
C...HEPEVT commonblock.
      INTEGER NMXHEP,NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
      PARAMETER (NMXHEP=4000)
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      DOUBLE PRECISION PHEP,VHEP
      SAVE /HEPEVT/
      INTEGER ISTOLD(NMXHEP),IHEP,IST,ISTLO,ISTHI,ISTOP,IMO
      INTEGER IDMOTH,IDDAUG
C
      INTEGER NPART
      INTEGER NMAX
      PARAMETER (NMAX=1000)
      DOUBLE PRECISION P(4,NMAX)
      INTEGER ID(NMAX)
      LOGICAL FOUND
      FOUND=.FALSE.
C TAG W-B FROM TOP DECAYS
      DO IHEP=1,NHEP
        IMO=JMOHEP(1,IHEP)
        IDMOTH=ABS(IDHEP(IMO))
        IDDAUG=ABS(IDHEP(IHEP))
        IF(IDMOTH.EQ.6.AND.(IDDAUG.EQ.5.OR.IDDAUG.EQ.24)) INORAD(IHEP)=2
      ENDDO
C
      DO 10 IHEP=1,NHEP
        IST=ISTHEP(IHEP)
        ISTOLD(IHEP)=IST
        IF (IST.EQ.ISTOP) FOUND=.TRUE.
        IF (IST.GE.ISTLO.AND.IST.LE.ISTHI.AND..NOT.FOUND) THEN
C     FOUND A RADIATED PARTON, CHECK MOTHER
          IMO=IHEP
 1        IMO=JMOHEP(1,IMO)
          IF(IMO.EQ.IEND) THEN
C     PARENTHOOD OK, VETO
            IST=0
            GOTO 9
          ENDIF
          IF(INORAD(IMO).EQ.1) THEN
            IDMOTH=ABS(IDHEP(IMO))
            IDDAUG=ABS(IDHEP(IHEP))
C     VERIFY IT'S A HEAVY QUARK -- LEAVE OUT GAUGE BOSON DECAYS
            IF(IDMOTH.GE.4.AND.IDMOTH.LE.6) THEN
C     PARTON COMES FROM A VETOED MOTHER, KEEP UNLESS IT IS THE EVOLVED
C     MOTHER
              IF(IDMOTH.NE.IDDAUG) THEN
                IST=1
                GOTO 9
              ELSE
                IST=0
              ENDIF
            ELSE
C     NOT A HEAVY QUARK MOTHER, LEAVE OUT OF JET RECONSUTRCTION
              IST=0
            ENDIF
          ELSEIF(INORAD(IMO).EQ.2) THEN
C     IT'S A TOP DECAY PRODUCT
            IST=0
          ELSE
C     GO CHECK GRANDMOTHER
            GOTO 1
          ENDIF
        ELSE
          IST=0
        ENDIF
 9      ISTHEP(IHEP)=IST
 10   CONTINUE
C      CALL FJHEVT(ISTLO,ISTHI,ISTOP,P,NPART)
      CALL FJHEVT(1,1,200,P,ID,NPART)
      DO 20 IHEP=1,NHEP
        ISTHEP(IHEP)=ISTOLD(IHEP)
 20   CONTINUE
      END

C-----------------------------------------------------------------------
      FUNCTION PSERAP(P)
C     PSEUDO-RAPIDITY (-LOG TAN THETA/2)
C-----------------------------------------------------------------------
      DOUBLE PRECISION PSERAP,P(3),PT,PL,TINY,THETA
      PARAMETER (TINY=1D-3)
      PT=SQRT(P(1)**2+P(2)**2)+TINY
      PL=P(3)
      THETA=ATAN2(PT,PL)
      PSERAP=-LOG(TAN(0.5*THETA))
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE ALPSOR(A,N,K,IOPT)
C-----------------------------------------------------------------------
C     Sort A(N) into ascending order
C     IOPT = 1 : return sorted A and index array K
C     IOPT = 2 : return index array K only
C-----------------------------------------------------------------------
      DOUBLE PRECISION A(N),B(5000)
      INTEGER N,I,J,IOPT,K(N),IL(5000),IR(5000)
      IF (N.GT.5000) then
        write(*,*) 'Too many entries to sort in alpsrt, stop'
        stop
      endif
      if(n.le.0) return
      IL(1)=0
      IR(1)=0
      DO 10 I=2,N
      IL(I)=0
      IR(I)=0
      J=1
   2  IF(A(I).GT.A(J)) GOTO 5
   3  IF(IL(J).EQ.0) GOTO 4
      J=IL(J)
      GOTO 2
   4  IR(I)=-J
      IL(J)=I
      GOTO 10
   5  IF(IR(J).LE.0) GOTO 6
      J=IR(J)
      GOTO 2
   6  IR(I)=IR(J)
      IR(J)=I
  10  CONTINUE
      I=1
      J=1
      GOTO 8
  20  J=IL(J)
   8  IF(IL(J).GT.0) GOTO 20
   9  K(I)=J
      B(I)=A(J)
      I=I+1
      IF(IR(J)) 12,30,13
  13  J=IR(J)
      GOTO 8
  12  J=-IR(J)
      GOTO 9
  30  IF(IOPT.EQ.2) RETURN
      DO 31 I=1,N
  31  A(I)=B(I)
      END
C-----------------------------------------------------------------------
      subroutine getunit(n)
      implicit none
      integer n,i
      logical yes
      do i=10,100
        inquire(unit=i,opened=yes)
        if(.not.yes) goto 10
      enddo
      write(*,*) 'no free units to write to available, stop'
      stop
 10   n=i
      end

c-------------------------------------------------------------------
      SUBROUTINE AHSPAR
c     set list of parameters types and assign default values
c-------------------------------------------------------------------
      implicit none
      include 'alpsho.inc'
      double precision chvalue
      character chparam*8
      integer n,iunit,i,j,aluisl,itmp
c     beam parameters
c
      chpar(2)='ih2'
      chpdes(2)='Select pp (1) or ppbar (-1) collisions'
      partyp(2)=1
      ih2=int(parval(2))
c
      chpar(3)='ebeam'
      chpdes(3)='beam energy in CM frame (e.g. 7000 for LHC)'
      partyp(3)=0
      ebeam=parval(3)
c
      chpar(4)='ndns'
      chpdes(4)='parton density set'
c currently available:
c ndns= 1      2      3       4      5      6       7      8
c pdf = cteq4m cteq4l cteq4hj cteq5m cteq5l cteq5hj cteq6m cteq6l
c ndns= 101    102        103        104        105       
c pdf = mrst99 mrst2002-1 mrst2002-2 mrst2002-3 mrst2002-4
      partyp(4)=1
      ndns=int(parval(4))
c
      chpar(7)='ickkw'
      chpdes(7)
     $     ='CKKW scale option: set to 1 to enable jet-parton matching'
      partyp(7)=1
      ickkw=int(parval(7))
c
      chpar(10)='njets'
      chpdes(10)='number of light jets'
      partyp(10)=1
      NLJETS=int(parval(10))
c
      chpar(11)='ihvy'
      chpdes(11)='heavy flavour type for procs like WQQ, ZQQ, 2Q, etc'/
     $     /'(4=c, 5=b, 6=t)'
      partyp(11)=1
      ihvy=int(parval(11))
c
      chpar(12)='ihvy2'
      chpdes(12)='2nd heavy flavour type for procs like 4Q'
      partyp(12)=1
      ihvy2=int(parval(12))
c
      chpar(13)='nw'
      chpdes(13)='number of W bosons'
      partyp(13)=1
      nw=int(parval(13))
c
      chpar(14)='nz'
      chpdes(14)='number of Z bosons'
      partyp(14)=1
      nz=int(parval(14))
c
      chpar(15)='nh'
      chpdes(15)='number of H bosons'
      partyp(15)=1
      nh=int(parval(15))
c
      chpar(16)='nph'
      chpdes(16)='number of photons'
      partyp(16)=1
      nph=int(parval(16))
c
c     masses
      chpar(20)='mc'
      chpdes(20)='charm mass'
c     the charm quark is considered massless unless explicitly requested
C     (e.g. in processes like W c cbar or c cbar)
c     it is 0 in W c processes
      partyp(20)=0
      parval(20)=mc
c
      chpar(21)='mb'
      chpdes(21)='bottom mass'
      partyp(21)=0
      parval(21)=mb
c
      chpar(22)='mt'
      chpdes(22)='top mass'
      partyp(22)=0
      parval(22)=mt
c
      chpar(23)='mh'
      chpdes(23)='higgs mass'
      partyp(23)=0
      parval(23)=mh
c
c     pt cuts
      chpar(30)='ptjmin'
      chpdes(30)='minimum pt for light jets'
      partyp(30)=0
      ptjmin=parval(30)
c
      chpar(31)='ptbmin'
      chpdes(31)='ptmin for bottom quarks (in procs with explicit b)'
      partyp(31)=0
      ptbmin=parval(31)
c
      chpar(32)='ptcmin'
      chpdes(32)='ptmin for charm quarks (in procs with explicit c)'
      partyp(32)=0
      ptcmin=parval(32)
c
      chpar(33)='ptlmin'
      chpdes(33)='minimum pt for charged leptons'
      partyp(33)=0
      ptlmin=parval(33)
c
      chpar(34)='metmin'
      chpdes(34)='minimum missing et'
      partyp(34)=0
      metmin=parval(34)
c
      chpar(35)='ptphmin'
      chpdes(35)='minimum pt for photons'
      partyp(35)=0
      ptphmin=parval(35)
c
      chpar(36)='ptcen'
      chpdes(36)='min pt for central jet in VBF 3-jet final states'/
     $ /' (used if irapgap = 1 and njets= 3)'
      partyp(36)=0
c
c     eta cuts
      chpar(40)='etajmax'
      chpdes(40)='max|eta| for light jets'
      partyp(40)=0
      etajmax=parval(40)
c
      chpar(41)='etabmax'
      chpdes(41)='max|eta| for b quarks (in procs with explicit b)'
      partyp(41)=0
      etabmax=parval(41)
c
      chpar(42)='etacmax'
      chpdes(42)='max|eta| for c quarks (in procs with explicit c)'
      partyp(42)=0
      etacmax=parval(42)
c
      chpar(43)='etalmax'
      chpdes(43)='max abs(eta) for charged leptons'
      partyp(43)=0
      etalmax=parval(43)
c
      chpar(44)='etaphmax'
      chpdes(44)='max abs(eta) for photons'
      partyp(44)=0
      etaphmax=parval(44)
c
      chpar(45)='irapgap'
      chpdes(45)='enable central rap-gap in VBF >=2 jet events'
      partyp(45)=1
c
      chpar(46)='etagap'
      chpdes(46)='min rap for 2 "fwd" jets in VBF >=2 jet events'/
     $ /' (used if irapgap = 1)'
      partyp(46)=0
c
c     isolation cuts
      chpar(50)='drjmin'
      chpdes(50)='min deltaR(j-j), deltaR(Q-j) [j=light jet, Q=c/b]'
      partyp(50)=0
      drjmin=parval(50)
c
      chpar(51)='drbmin'
      chpdes(51)='min deltaR(b-b) (procs with explicit b)'
      partyp(51)=0
      drbmin=parval(51)
c
      chpar(52)='drcmin'
      chpdes(52)='min deltaR(c-c) (procs with explicit charm)'
      partyp(52)=0
      drcmin=parval(52)
c
      chpar(55)='drlmin'
      chpdes(55)='min deltaR between charged lepton and light jets'
      partyp(55)=0
      drlmin=parval(55)
c
      chpar(56)='drphjmin'
      chpdes(56)='min deltaR between photon and light jets'
      partyp(56)=0
      drphjmin=parval(56)
c
      chpar(57)='drphlmin'
      chpdes(57)='min deltaR between photon and charged lepton'
      partyp(57)=0
      drphlmin=parval(57)
c
      chpar(58)='drphmin'
      chpdes(58)='min deltaR between photons'
      partyp(58)=0
      drphmin=parval(58)
c
c     dilepton cuts
      chpar(61)='mllmin'
      chpdes(61)='min dilepton inv mass'
      partyp(61)=0
      mllmin=parval(61)
c
      chpar(62)='mllmax'
      chpdes(62)='max dilepton inv mass'
      partyp(62)=0
      mllmax=parval(62)
c
      chpar(101)='itdec'
      chpdes(101)='forces top decays, with spin-correlations'
      partyp(101)=1
      itdec=int(parval(101))
c
      chpar(102)='itopprc'
      chpdes(102)='Selection of single-top process'
      partyp(102)=1
      itopprc=int(parval(102))
c
      end


c-------------------------------------------------------------------
      SUBROUTINE AHSHPA
c     set list of parameters types and assign default values
c-------------------------------------------------------------------
      implicit none
      include 'alpsho.inc'
      double precision chvalue
      character chparam*8,nchar*10
      integer i,j,k,aluisl,itmp,iflag
      nchar='0123456789'
C
      do i=1,nSHparam
        chSHpar(i)='***'
        chSHpdes(i)='parameter not assigned'
        do itmp=0,1
          SHparuse(i,itmp)=0
        enddo
      enddo
c
      chSHpar(1)='nevmax'
      chSHpdes(1)='max number of events to process, -1 for full file'
      SHpartyp(1)=1
      SHparval(1)=-1
      do i=0,1
        SHparuse(1,i)=1
      enddo
c
      chSHpar(2)='iexc'
      chSHpdes(2)
     $     ='exclusive (iexc=1) / inclusive (iexc=0) extra-jet veto'
      SHpartyp(2)=1
      SHparval(2)=1
      SHparuse(2,0)=0
      SHparuse(2,1)=1
c
      chSHpar(3)='ncmax'
      chSHpdes(3)='maximum number of charm quarks allowed after shower'
      SHpartyp(3)=1
      SHparval(3)=100
      do i=0,1
        SHparuse(3,i)=1
      enddo
c
      chSHpar(4)='nbmax'
      chSHpdes(4)='maximum number of bottom quarks allowed after shower'
      SHpartyp(4)=1
      SHparval(4)=100
      do i=0,1
        SHparuse(4,i)=1
      enddo
      ihvyv=0
c
      chSHpar(10)='etclus'
      chSHpdes(10)='min et for jet clusters used in matching'
      SHpartyp(10)=0
      SHparval(10)=MAX(PTJMIN+5,1.2*PTJMIN)
      SHparuse(10,0)=0
      SHparuse(10,1)=1
c
      chSHpar(11)='rclus'
      chSHpdes(11)='cone size for jet clusters used in matching'
      SHpartyp(11)=0
      SHparval(11)=DRJMIN
      SHparuse(11,0)=0
      SHparuse(11,1)=1
c
      chSHpar(12)='etaclmax'
      chSHpdes(12)='eta max for jet clusters used in matching'
      SHpartyp(12)=0
      SHparval(12)=ETAJMAX
      SHparuse(12,0)=0
      SHparuse(12,1)=1
c
      do k=1,10
        j=90+k
        chSHpar(j)='par'//nchar(k:k)
        chpdes(j)='auxiliary parameter (real*8)'
        SHpartyp(j)=0
        SHparval(j)=0
        do i=0,1
          SHparuse(j,i)=1
        enddo
      enddo
C
      do i=1,nSHparam
        SHparlen(i)=aluisl(chSHpar(i))
      enddo
c
      IEXC=int(SHPARVAL(2))
      ETCLUS=SHPARVAL(10)
      RCLUS=SHPARVAL(11)
      ETACLMAX=SHPARVAL(12)      
      NCMAX=int(SHPARVAL(3))
      NBMAX=int(SHPARVAL(4))
c
      IF(ICKKW.EQ.1) THEN
        WRITE(6,*) ' '
        WRITE(6,*) 'JET PARAMETERS FOR MATCHING:'
        WRITE(6,111) 'ET>',ETCLUS,' ETACLUS<',ETACLMAX,' RCLUS=',RCLUS
 111    FORMAT(3(A,F9.3,1X))


C      WRITE(6,*) 'DR(PARTON-JET)<',1.5*RCLUS
        IF(NLJETS.EQ.0) THEN
          WRITE(6,*)
     $         'FOR NLJETS=0 MAKE SURE THAT THESE ARE CONSISTENT WITH'
          WRITE(*,*) 'THE PARAMETERS USED FOR NLJETS>0'
        ENDIF
      ENDIF

c
 1    write(6,*) ' '
      write(6,*) 'Input parameters to replace defaults:  type and value'
      write(6,*)
     $    '(input ''print 1'' to display the list of parameter types and
     $ their current values)'
      write(6,*)
     $ '(input ''eoi 1'' to terminate the input sequence)'
      write(6,*)
     $ '(from the terminal, ''ctrl-D'' to terminate the input sequence)'
 2    read(5,*,end=3,err=10) chparam,chvalue
 10   if(chparam(1:3).eq.'eoi') goto 3
      if(chparam(1:1).eq.'*') goto 2
      itmp=aluisl(chparam)
c      call AHfpar(chparam(1:itmp),chvalue,iflag)
      call AHfpar(chparam,chvalue,iflag)
      if(iflag.le.1) then
        goto 2
      elseif(iflag.eq.2) then
        goto 1
      elseif(iflag.eq.3) then
        write(6,*) 'param ',chparam(1:itmp),' not available/changeable'
        goto 1
      endif
 3    continue
C     assign final set of parameters
      NEVMAX=int(SHPARVAL(1))
      IEXC=int(SHPARVAL(2))
      IF(NCMAX.NE.SHPARVAL(3).OR.NBMAX.NE.SHPARVAL(4)) IHVYV=1
      NCMAX=int(SHPARVAL(3))
      NBMAX=int(SHPARVAL(4))
      ETCLUS=SHPARVAL(10)
      RCLUS=SHPARVAL(11)
      ETACLMAX=SHPARVAL(12)
C     WRITE PARAMETERS IN LOG FILE
      WRITE(NUNITOUT,*) '==========================='
      WRITE(NUNITOUT,*) 'SHOWER/MATCHING PARAMETERS:'
      do i=1,nSHparam
        IF(SHparuse(i,iCKKW).EQ.1) THEN
          WRITE(NUNITOUT,*) CHSHPAR(I),'=',SHPARVAL(I)
        ENDIF
      ENDDO
      do i=0,9
        xpar(i)=SHparval(91+i)
      enddo
      end

c-------------------------------------------------------------------
      subroutine AHfpar(chparam,chvalue,iflag)
c     deposit parameter values into common blocks and store in fname.par
c-------------------------------------------------------------------
      implicit none
      include 'alpsho.inc'
      double precision chvalue
      character chparam*8,chtmp*8
      integer i,iflag,ipar,itmp,aluisl
      iflag=0
      itmp=aluisl(chparam)
      if(chparam(1:5).eq.'print') then
CCCCCMODIFIED
C        call ahppar(int(chvalue))
        call ahppar
        iflag=2
        return
      elseif(chparam(1:3).eq.'***') then
        iflag=1
        return
      else
        do i=1,nSHparam
          chtmp=chSHpar(i)
          if(chparam(1:itmp).eq.chtmp(1:itmp)) then
            ipar=i
            if(SHparuse(ipar,ickkw).eq.0) then
              iflag=3
              return
            endif
            goto 100
          endif
        enddo
      endif
      iflag=3
      return
 100  SHparval(ipar)=chvalue
c      write(NunitOut,*) chSHpar(ipar),chvalue
      end


      function aluisl(string)
c returns the position of the last non-blank character in string
      integer aluisl
      character * (*) string
      i = len(string)
      dowhile(i.gt.0.and.string(i:i).eq.' ')
         i=i-1
      enddo
      aluisl = i
      end

      subroutine AHppar
      implicit none
      include 'alpsho.inc'
      integer iunit,i,itmp,aluisl
c print definitions and current values of parameters
      iunit=6
      write(iunit,*) '------'
      do i=1,nSHparam
        if(chSHpar(i).ne.'***'.and.SHparuse(i,ickkw).eq.1) then
          itmp=aluisl(chSHpdes(i))
          write(iunit,*) '------'
          write(iunit,*) chSHpdes(i)(1:itmp),':'
          itmp=aluisl(chSHpar(i))
          if(SHpartyp(i).eq.0) then
            write(iunit,*) chSHpar(i)(1:itmp),'=',SHparval(i)
          else
            write(iunit,*) chSHpar(i)(1:itmp),'=',int(SHparval(i))
          endif
        endif
      enddo
      end
