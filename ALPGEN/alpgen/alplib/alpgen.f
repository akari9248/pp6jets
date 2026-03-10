c     Version to allow reading in a standard event, and perform the kt
C     clustering and alphas rescaling
C-----------------------------------------------------
      program alpgen
      character*5 nversion
      common/version/nversion
      nversion='3.1'
C-------------------------------------------------------------------------
C 
C driver for multi-parton matrix element generator based on ALPHA
C
C--------------------------------------------------------------------------
C
      call alsprc
c
c     Routine location: XXX.f
c     Purpose: 
c        -  assign the hard process code for selected process XXX
c
C setup running defaults:
c
      call alsdef
c
c     Routine location: alpgen.f
c     Purpose: 
c        -  initialise the default generation parameters (e.g. beam
c           type, energy, PDF set) 
c        -  initialise default mass and couplings for particles
c
c
C setup  event generation options, bookkeeping, etc
      call alhset
c
c     Routine location: XXX.f
c     Purpose: 
c        - setup event generation options, bookkeeping, etc, for
c              the specific hard process XXX
c        - write run information on stat, unwpar files

c initialise alpha parameters
      call alinit
c     Routine location: alpgen.f
c     Purpose: 
c        - evaluate parameter-dependent quantities (e.g. Higgs width)
c	 - fill the apar array of ALPHA (which contains all parameters
c	   required by ALPHA
c
c
C setup internal bookkeeping histograms
      call alsbkk
c
c     Routine location: alpgen.f
c     Purpose: initialise histograms common to all processes, and needed
C     for internal purposes
c
c
C setup user histograms
      call alshis
c
c     Routine location: XXXusr.f
c     Purpose: initialise histograms
c
c
c
C setup integration grids, including optimization if required
c
      call alsgrd
c
c     Routine location: XXX.f
c     Purpose: setup integration grid variables
      call aligrd
c
c     Routine location: alpgen.f
c     Purpose: initialise grid with warm-up iterations, if required
c
c
C generate events
c
      Call alegen
c
c     Routine location: alpgen.f
c     Purpose: generates events, calling in a standad format the
c          the process-depepdent phase-space and flavour-selection
c          routines, contained in XXX.f

C finalise histograms
c
      call alfhis
c
c     Routine location: XXXusr.f
c     Purpose: finalize analysis and histograms
      call alfbkk
c
c     Routine location: alpgen.f
c     Purpose: finalize internal histograms
c
c
      end

c-------------------------------------------------------------------
      subroutine alsbkk
c     setup weight bookeeping histograms
c-------------------------------------------------------------------
      implicit none
      double precision logwgt,wgtdis,wmin,wmax,wbin
      common/wgtbkk/logwgt(1000),wgtdis(0:1001),wmin,wmax,wbin
      integer i
      wmin=log10(1d-20)
      wmax=log10(1d20)
      wbin=(wmax-wmin)/1000d0
      do i=1,1000
        logwgt(i)=wmin+wbin*(i-0.5)
        wgtdis(i)=0d0
      enddo
c under- and overflow bin
      wgtdis(0)=0d0
      wgtdis(1001)=0d0
      call mbook(190,'rewgt factors',0.02e0,0.e0,2.e0)
      end

c-------------------------------------------------------------------
      subroutine alhbkk(wgt)
c     online bookeeping of weight distributions
c-------------------------------------------------------------------
      implicit none
      double precision wgt
      double precision logwgt,wgtdis,wmin,wmax,wbin
      common/wgtbkk/logwgt(1000),wgtdis(0:1001),wmin,wmax,wbin
      integer i
      double precision lwgt
      if(wgt.le.0.d0) return
      lwgt=log10(wgt)
      i=int((lwgt-wmin)/wbin)+1
      if(i.lt.0) then
c underflow bin
        wgtdis(0)=wgtdis(0)+wgt
      elseif(i.le.1000) then
        wgtdis(i)=wgtdis(i)+wgt
      else
c overflow bin
        wgtdis(1001)=wgtdis(1001)+wgt
      endif
      end

c-------------------------------------------------------------------
      subroutine alfbkk
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      double precision logwgt,wgtdis,wmin,wmax,wbin
      common/wgtbkk/logwgt(1000),wgtdis(0:1001),wmin,wmax,wbin
c locals
      double precision sigwgt(0:1001),newwmax(5)
      integer nev,iproc
      real Sq,Savgwgt
      integer i,j
      real  xnorm
      character*80 tmpstr
      integer aluisl,nfile1,nfile2

c find maxwgt thresholds for 5-4-3-2-1% of total weight
      if(imode.eq.1) then
        sigwgt(1001)=wgtdis(1001)
        do i=1000,0,-1
          sigwgt(i)=sigwgt(i+1)+wgtdis(i)
        enddo
        do i=1,1000
          do j=1,5
            if(sigwgt(i)/sigwgt(0).gt.0.01*float(j)) then
              newwmax(j)=logwgt(i)
            endif
          enddo
        enddo
        write(niopar,'(5(e12.6,1x))')(10d0**(newwmax(i)+0.5*wbin),i=1,5)
        close(niopar)
      endif

c     
      if(imode.eq.2) then
        close(niounw)
      endif
c      open(unit=99,file=topfile,err=999,status='old')
c      call aluend(99)
c compatbility with latest gfortran
      open(unit=99,file=topfile,access='append',err=999,status='old')
      if(imode.le.1) then
         xnorm=sngl(avgwgt/totwgt)
      elseif(imode.eq.2) then
         xnorm=1e0/real(unwev)
      else
         write(6,*) 'imode type not allowed, stop'
         stop
      endif
c
      do i=160,200
         if(i.ne.192) call mopera(i,'F',i,i,xnorm,1.e0)
         call mfinal(i)
      enddo 
c
c
      call mtop(190,99,'rewgt ','dN/d rewgt ','LIN')
      close(99)
 999  return
      end

c-------------------------------------------------------------------
      subroutine alsdef
c     setup default run parameters
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
c common declarations
c local variables
      character*50 tmpstr
      integer aluisl,l
      double precision dummy,ranset,ranset2,ccab
      integer i
      data imaxiter/1,4,6,10,17,28,45,49,49,10*49/
c
c list processes
      do i=1,nprocs
        chprc(i)='unavailable'
      enddo
      chprc(1)='wqq'
      chprc(2)='zqq'
      chprc(3)='wjet'
      chprc(4)='zjet'
      chprc(5)='vbjet'
      chprc(6)='2Q'
      chprc(7)='4Q'
      chprc(8)='QQh'
      chprc(9)='Njet'
      chprc(10)='wcjet'
      chprc(11)='phjet'
      chprc(12)='hjet'
      chprc(13)='top'
      chprc(14)='wphjet'
      chprc(15)='wphqq'
      chprc(16)='2Qph'
      chprc(17)='zphjet'
      chprc(18)='zphqq'
      chprc(19)='hvyN'
      chprc(20)='4ljet'
c modified_b
c      nactprc=16
      chprc(nexternal)='external'
      nactprc=nprocs
c modified_e
c fixed parameters:
      resc=1
c     number of calls for the spin/colour event-by-event average
      navg=1
c     beam 1   (ih=1: proton;  ih=-1: pbar)
      ih1= 1
c     lepton masses 
      mlep(1)=0.511d-3
      mlep(2)=0.10566d0
      mlep(3)=1.77682d0
c     ckm mixings
      ccab=0.975
      ccab2=ccab**2
      scab2=1-ccab2
c     EW parameters:
      mw=80.385d0
      mz=91.1876d0
      stw=sqrt(0.231d0)
      aem=1.d0/128.89
      gfermi=1.1663787d-5
C--   To ensure the gauge invariance of the ALPHA calculations,
c     we shall use the LO relations between mW, mZ, weak couplings,
c     GFermi and alpha(em).  Which parameters can be input, and which
c     are calculated, is governed by the option IEWOPT, which can be
c     reassigned in the ALSUSR routine. The implementation of the option
c     is contained in the routine ALINIT
c
      iewopt=3
c eliminate propagating photons in (w)phjet proceses:
      if(ihrd.eq.2.or.ihrd.eq.4.or.ihrd.eq.11.or.
     +   ihrd.eq.14.or.ihrd.eq.15.or.ihrd.eq.17.or.ihrd.eq.18) resc=1d3
c scale setting parameters:
      iqopt=1
      qfac=1d0
      ktfac=1d0
      ickkw=0
c clustering parameters
c     1: pclu(1:4)=p1(1:4)+p2(1:4)  2: pclu(1:3)=p1(1:3)+p2(1:3), mclu=0
      mrgopt=1
      cluopt=1
c
c     assign Default PARameter values for usr-accessible parameters
      call aldpar(1)
c
C     MANDATORY INPUTS:
      write(6,*) 'Input RUN generation mode:'
      write(6,*) '0: generate weighted events, no evt dumps to file'
      write(6,*)
     $   '1: generate wgtd events, write to file for later unweighting'
      write(6,*)
     $   '2: read events from file for unweighting or processing'
      write(6,*)
     $ 'or documentation modes:'
      write(6,*)
     $ '3: print parameter options and defaults, then stop'
      write(6,*)
     $ '4: write to par.list parameter options and defaults, then stop'
      write(6,*)
     $ '5: write to prc.list complete list of processes, parameter ',
     $ '   options and defaults, scale choices, PDF, etc., then stop'
      read(5,*) imode
c
c documentation modes:
      if(imode.ge.3.and.imode.le.5) then
        call alppar(imode)
        stop
      endif

      write(6,*) 'Input string labeling output and input files'
      write(6,*) '(e.g. w2j to output files w2j.stat, etc.)'
      read(5,*) fname
c     define input files
      call alsinp
c     
      if(imode.eq.2) then
c read in parameters from fname.par and from input
        call alrpar
      elseif(imode.lt.2) then
c     read in parameters from input
        call alrusr
      else
        write(6,*) 'With imode=',imode,' we should not be here, stop'
        stop
      endif
c     save parameter values to relative variables in common blocks
      call alspar
c     define output files
      write(*,*) 'ILHE=',ilhe
      call alsout
c     deposit random number generator seeds
      dummy= ranset(iseed)
      if(imode.eq.1) then
        dummy= ranset2(iseed)
      elseif(imode.eq.2) then
        dummy= ranset2(iseed2)
      endif
      end
      

c-------------------------------------------------------------------
      subroutine alrusr
c     input datacards superseding default parameters
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
CC      double precision chvalue
      character chvalue*20
      double precision tmp
      character chparam*8
      integer iflag,i
      integer aluisl,itmp,ntmp
      character*8 tmpstrfix
      character*8 blnk/'        '/
      character*1 blnk1/' '/
c     
c     select options for grid selection/generation
      write(6,*) ' '
      write(6,*) 'To generate new grid input 0, '
      write(6,*) 'To use grid generated by the warmup iterations of',
     +     ' the previous run, input 1,'
      write(6,*) 'To use grid generated by the event generation of',
     +     ' the previous run, input 2:'
      read(5,*) igrid
      write(6,*) 'Input N(events)/iteration and N(iterations) for the 
     + warmup iterations:'
      read(5,*) nopt,niter
      write(6,*) 'Input number evts to generate:'
      read(5,*) maxev
c     input scale options
      write(6,*) ' '
      call alhsca(ihrd,6)
c     
c     
c     
 1    write(6,*) ' '
      write(6,*) 'Input parameters to replace defaults:  type and value'
      write(6,*)
     $    '(input ''print 1'' to display the list of parameter types and
     $ their current values)'
      write(6,*)
     $     '(input ''print 2'' to write the list to file par.list)'
      write(6,*)
     $     '(input ''ctrl-D'' to terminate the input sequence)'
 2    read(5,*,end=3,err=10) chparam,chvalue
 10   if(chparam(1:3).eq.'eoi') goto 3
      if(chparam(1:1).eq.'*') goto 2
      itmp=aluisl(chparam)
C      if(chparam(1:3).eq.'new') call alfpar_external
C     $(chparam(1:itmp),chvalue,iflag)
c      tmpstrfix=chparam(1:itmp)//blnk
      tmpstrfix=chparam(1:itmp)
      do i=itmp+1,8
         tmpstrfix=tmpstrfix//blnk1
      enddo
      call alfpar(tmpstrfix,chvalue,iflag)
      if(iflag.le.1) then
        goto 2
      elseif(iflag.eq.2) then
        goto 1
      elseif(iflag.eq.3) then
        write(6,*) 'param ',chparam(1:itmp),
     +   ' not available/changeable for this process or imode=',imode
        goto 1
      endif
 3    write(niopar,*) 'eoi',1
      return
      end

c-------------------------------------------------------------------
      subroutine alhsca(ihrd,iunit)
c     defines the hard process scales
c-------------------------------------------------------------------
      implicit none
      integer ihrd,iunit
c     qsq scale; 
      write(iunit,*)
     $     'Options for Factorization/renormalization scale Q:'
      write(iunit,*) 'iqopt=0 => Q=qfac'
c modified_b
      call alhsca_external(ihrd,iunit)
c modified_e
      if(ihrd.eq.1) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m_W^2+ sum_jets(m_tr^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*mW'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m_W^2+ pt_W^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum_jets(m_tr^2)}'
        write(iunit,*) 'iqopt=5 => Q=qfac*HT'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- m_tr^2=m^2+pt^2, summed over heavy quarks and light jets'
      elseif(ihrd.eq.2) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m0^2+ sum_jets(m_tr^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*m0'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m0^2 + pt_Z^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum_jets(m_tr^2)}'
        write(iunit,*) 'iqopt=5 => Q=qfac*HT'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- m_tr^2=m^2+pt^2, summed over heavy quarks and light jets'
        write(iunit,*) '- m0=mll' 
      elseif(ihrd.eq.3) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m_W^2+ sum(pt_jet^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*mW'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m_W^2+ pt_W^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum(pt_jet^2)}'
        write(iunit,*) 'iqopt=5 => Q=qfac*HT'
      elseif(ihrd.eq.4) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m0^2+ sum(pt_jet^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*m0'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m0^2 + pt_Z^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum(pt_jet^2)}'
        write(iunit,*) 'iqopt=5 => Q=qfac*HT'
        write(iunit,*) 'where m0=mll' 
      elseif(ihrd.eq.5) then
        write(iunit,*) 
     $ 'iqopt=1 => Q=qfac*sqrt{sum(mV)^2+sum(pt_photons^2+pt_jet^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sum(mV)'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{shat}'
      elseif(ihrd.eq.6) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{sum(m_tr^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt(x1*x2*S)'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- m_tr^2=m^2+pt^2, summed over heavy quarks and light jets'
      elseif(ihrd.eq.7) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{sum(m_tr^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt(x1*x2*S)'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- m_tr^2=m^2+pt^2, summed over heavy quarks and light jets'
      elseif(ihrd.eq.8) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{mh^2+sum(m_tr^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt(x1*x2*S)'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- m_tr^2=m^2+pt^2, summed over heavy quarks and light jets'
      elseif(ihrd.eq.9) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{sum(pt_jet^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt(x1*x2*s)'
      elseif(ihrd.eq.10) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m_W^2+ sum(pt_jet^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*mW'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m_W^2+ pt_W^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum(pt_jet^2)}'
      elseif(ihrd.eq.11) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{sum(pt_ph^2+pt_jets^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt{pt_jets^2}'
      elseif(ihrd.eq.12) then 
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{nh*mh^2+pt_jets^2}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt{shat}'
      elseif(ihrd.eq.13) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt(mw^2+sum(m_tr^2))'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt{shat}'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- m_tr^2=m^2+pt^2, summed over heavy quarks and light jets'
        write(iunit,*) 
     $ '- the mw^2 term is present only for processes with a W'
      elseif(ihrd.eq.14) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m_W^2+ sum(pt^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*mW'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m_W^2+ pt_W^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum(pt^2)}'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- pt^2 is summed over photons and jets'
      elseif(ihrd.eq.15) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m_W^2+ sum(pt^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*mW'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m_W^2+ pt_W^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum(pt^2)}'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- pt^2 is summed over photons, heavy quarks and light jets'
      elseif(ihrd.eq.16) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{sum(m_tr^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt(x1*x2*S)'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- m_tr^2=m^2+pt^2, summed over heavy quarks and light jets'
      elseif(ihrd.eq.17) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m_Z^2+ sum(pt^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*mZ'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m_Z^2+ pt_Z^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum(pt^2)}'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- pt^2 is summed over photons and jets'
      elseif(ihrd.eq.18) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m_Z^2+ sum(pt^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*mZ'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m_Z^2+ pt_Z^2}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum(pt^2)}'
        write(iunit,*) 'where:'
        write(iunit,*) 
     $ '- pt^2 is summed over photons heavy quarks and jets'
      elseif(ihrd.eq.20) then
        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m_4l^2 + sum(ptj^2)}'
        write(iunit,*) 'iqopt=2 => Q=qfac*sqrt{sum(ptj^2)}'
        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m_4l^2+ sum(pt_l^2)}'
        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum(ptj^2)}'
      endif
c
c      write(iunit,*) ' '
c      write(iunit,*) 'Default is iqopt=1, qfac=1'
c      write(iunit,*) ' '
      if(ihrd.le.6.or.ihrd.eq.9.or.ihrd.eq.10.or.ihrd.eq.11.or.ihrd.eq
     $     .12.or.ihrd.eq.14.or.ihrd.eq.15.or.ihrd.eq.16
     $     .or.ihrd.eq.17.or.ihrd.eq.18.or.ihrd.eq.20) then
        write(iunit,*) 'To select CKKW scale input ''ickkw 1'' '
        write(iunit,*) '(mandatory for later use of jet matching)'
        write(iunit,*) 'In imode=2 select clustering option ''cluopt'':'
        write(iunit,*) 'cluopt=1: kperp propto pt(cluster) (default)'
        write(iunit,*) 'cluopt=2: kperp propto mt(cluster)'
        write(iunit,*) 'kperp is then rescaled by ktfac'
      endif
      end

c-------------------------------------------------------------------
      subroutine alrpar
c     read in parameters different from defaults for imode=2
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
CC      double precision chvalue
      character chvalue*20
      character chparam*8
      integer aluisl,iflag,itmp
      character*50 tmpstr
      integer l
C first read in parameters from the imode=1 run
      itmp=0
 1    read(niopar,*,end=2,err=10) chparam,chvalue
 10   if(chparam(1:3).eq.'eoi') goto 2
      if(chparam(1:3).eq.'***') goto 1
      if(itmp.eq.0) then
        itmp=1
        write(6,*) ' '
        write(6,*) 'read in generation parameters:'
      endif
      call alfpar(chparam,chvalue,iflag)
      if(iflag.le.1) then
        goto 1
      elseif(iflag.eq.3) then
        write(6,*) 'param ',chparam,' not recognised, stop'
        stop
      else
        write(6,*) 'unrecognized status after parameter input, stop'
        stop
      endif
 2    continue
c Restore defaults of imode=2 params that depend on imode=1 inputs:
c No option for Z decays if Z->nunu
      if((ihrd.eq.2.or.ihrd.eq.4.or.ihrd.eq.17.or.ihrd.eq.18)
     +   .and.(ilep.eq.1)) then
        paruse(152,ihrd)=0
      endif
c Then read in decay parameters specific to the imode=2 run:
 15   write(6,*) ' '
      write(6,*)'Input decay params to replace defaults: type and value'
      write(6,*)
     $   '(input ''print 1'' to display the list of parameter types and
     $ their current values)'
      write(6,*)
     $     '(input ''ctrl-D'' to terminate the input sequence)'
 20   read(5,*,end=200,err=100) chparam,chvalue
 100  if(chparam(1:3).eq.'eoi') goto 200
      if(chparam(1:1).eq.'*') goto 20
      itmp=aluisl(chparam)
      call alfdpa(chparam,chvalue,iflag)
      if(iflag.le.1) then
        goto 20
      elseif(iflag.eq.2) then
        goto 15
      elseif(iflag.eq.3) then
        write(6,*) 'param ',chparam(1:itmp),' not accepted in imode=2'
        goto 15
      endif
 200  continue
      return
      end

      
c-------------------------------------------------------------------
      subroutine Alfpar(chparam,chvalue,iflag)
c     deposit parameter values into common blocks and store in fname.par
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
CC      double precision chvalue
      double precision tmp
      character chvalue*20
      character chparam*8,chtmp*8
      character*50 tmpstr
      integer i,iflag,ipar,itmp,aluisl,istall
      data istall/0/
      iflag=0
      itmp=aluisl(chparam)
      if(chparam(1:5).eq.'print') then
CCCCCCCMODIFIED
        read(chvalue,*) tmp
        call alppar(int(tmp))
CCCCCCCMODIFIED
        iflag=2
        return
      elseif(chparam(1:3).eq.'***') then
        iflag=1
        return
      else
        do i=1,npar
          chtmp=chpar(i)
          if(chparam(1:itmp).eq.chtmp(1:itmp)) then
            ipar=i
            if(paruse(ipar,ihrd).eq.0) then
              iflag=3
              return
            endif
            goto 100
          endif
        enddo
      endif
      iflag=3
      return
 100  parval(ipar)=chvalue
CCCCCCCMODIFIED
      if(ipar.eq.4) then
         read(chvalue,*) tmp
         if(tmp.eq.-1) then
CCCCCCCMODIFIED
            if(istall.gt.1000) then
               write(6,*) 'PDF code unspecified, stop'
               stop
            endif
         call prntsf(6)
         write(6,*) ' enter ''ndns'' followed by value'
         iflag=1
         istall=istall+1
         return
         endif
      endif
      if(paruse(ipar,ihrd).eq.1) then
CCCCCCCMODIFIED
C        if(imode.eq.1) write(niopar,*) chpar(ipar),chvalue
        if(imode.eq.1) write(niopar,*) chpar(ipar),' ',chvalue
c
        call alustc(fname,'.err',tmpstr)
        call alugun(nioerr)
        open(unit=nioerr,file=tmpstr,status='unknown',access='append')
        write(nioerr,*) chpar(ipar),' ',chvalue
        close(nioerr)
c
        if(imode.eq.2) write(6,*) chpar(ipar),'=',chvalue
      endif
      end

c-------------------------------------------------------------------
      subroutine Alfdpa(chparam,chvalue,iflag)
c     deposit parameter values for imode=2 into common blocks and store
C     in fname.par 
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
CC      double precision chvalue
      character chvalue*20
      character chparam*8,chtmp*8
      integer i,iflag,ipar,itmp,aluisl,istall
      data istall/0/
      iflag=0
      itmp=aluisl(chparam)
      if(chparam(1:5).eq.'print') then
        call alppar(6)
        iflag=2
        return
      elseif(chparam(1:3).eq.'***') then
        iflag=1
        return
      else
        do i=1,npar
          chtmp=chpar(i)
          if(chparam(1:itmp).eq.chtmp(1:itmp)) then
            ipar=i
            if(paruse(ipar,ihrd).eq.0) then
              iflag=3
              return
            endif
            goto 100
          endif
        enddo
      endif
      iflag=3
      return
c modified_b
c 100  if(ipar.le.150) then
 100  if(ipar.le.150.or.ipar.gt.nparoff) then
c modified_e
        iflag=3
        return
      endif 
      parval(ipar)=chvalue
      if(paruse(ipar,ihrd).eq.1) write(6,*) chpar(ipar),'=',chvalue
      end
c-------------------------------------------------------------------
      subroutine Aldpar(n)
c     set list of parameters types and assign default values
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
CC      double precision chvalue
      character chvalue*20
      character chparam*8,nchar*10
      integer n,iunit,i,j,k,aluisl,itmp
      nchar='1234567890'
      do i=1,npar
        chpar(i)='***'
        chpdes(i)='parameter not assigned'
        do itmp=1,nprocs
          paruse(i,itmp)=0
        enddo
      enddo
c     beam parameters
c
      chpar(2)='ih2'
      chpdes(2)='Select pp (1) or ppbar (-1) collisions'
      partyp(2)=1
      parval(2)='1'
      do i=1,nactprc
        paruse(2,i)=1
      enddo
c
      chpar(3)='ebeam'
      chpdes(3)='beam energy in CM frame (default: 7000 for LHC14)'
      partyp(3)=0
      parval(3)='7000'
      do i=1,nactprc
        paruse(3,i)=1
      enddo
c
      chpar(4)='ndns'
      chpdes(4)='parton density set'
c currently available:
c ndns= 1      2      3       4      5      6       7      8
c pdf = cteq4m cteq4l cteq4hj cteq5m cteq5l cteq5hj cteq6m cteq6l
c ndns= 101    102        103        104        105       
c pdf = mrst99 mrst2002-1 mrst2002-2 mrst2002-3 mrst2002-4
      partyp(4)=1
      parval(4)='5'
      do i=1,nactprc
        paruse(4,i)=1
      enddo
c
      chpar(5)='iqopt'
      chpdes(5)='scale option (process dependent)'
      partyp(5)=1
      parval(5)='1'
      do i=1,nactprc
        paruse(5,i)=1
      enddo
c
      chpar(6)='qfac'
      chpdes(6)='Q scale rescaling factor'
      partyp(6)=0
      parval(6)='1d0'
      do i=1,nactprc
        paruse(6,i)=1
      enddo
c
      chpar(7)='ickkw'
      chpdes(7)
     $     ='CKKW scale option: set to 1 to enable jet-parton matching'
      partyp(7)=1
      parval(7)='0'
      do i=1,7
        paruse(7,i)=1
      enddo
      do i=9,12
        paruse(7,i)=1
      enddo
      do i=14,18
        paruse(7,i)=1
      enddo
      paruse(7,20)=1
c
      chpar(8)='ktfac'
      chpdes(8) ='scale factor for ckkw alphas scale'
      partyp(8)=0
      parval(8)='1d0'
      do i=1,7
        paruse(8,i)=1
      enddo
      do i=9,12
        paruse(8,i)=1
      enddo
      do i=14,18
        paruse(8,i)=1
      enddo
      paruse(8,20)=1
c
      chpar(10)='njets'
      chpdes(10)='number of light jets'
      partyp(10)=1
      if(ihrd.eq.3.or.
     +   ihrd.eq.4.or.
     +   ihrd.eq.11.or.
     +   ihrd.eq.13) then
        parval(10)='1'
      elseif(ihrd.eq.9) then
        parval(10)='2'
      else
        parval(10)='0'
      endif
      do i=1,nactprc
        paruse(10,i)=1
      enddo
c
      chpar(11)='ihvy'
      chpdes(11)='heavy flavour type for procs like WQQ, ZQQ, 2Q, etc'/
     $     /'(4=c, 5=b, 6=t)'
      partyp(11)=1
      parval(11)='5'
      paruse(11,1)=1
      paruse(11,2)=1
      paruse(11,6)=1
      paruse(11,7)=1
      paruse(11,8)=1
      paruse(11,15)=1
      paruse(11,16)=1
      paruse(11,18)=1
c
      chpar(12)='ihvy2'
      chpdes(12)='2nd heavy flavour type for procs like 4Q'
      partyp(12)=1
      parval(12)='5'
      paruse(12,7)=1
c
      chpar(13)='nw'
      chpdes(13)='number of W bosons'
      partyp(13)=1
      parval(13)='0'
      if(ihrd.eq.5) parval(13)='2'      
      paruse(13,5)=1
c
      chpar(14)='nz'
      chpdes(14)='number of Z bosons'
      partyp(14)=1
      parval(14)='0'      
      paruse(14,5)=1
c
      chpar(15)='nh'
      chpdes(15)='number of H bosons'
      partyp(15)=1
      parval(15)='0' 
      if(ihrd.eq.12) parval(15)='1'
      paruse(15,5)=1
c      paruse(15,12)=1
c
      chpar(16)='nph'
      chpdes(16)='number of photons'
      partyp(16)=1
      parval(16)='0'
      if(ihrd.eq.11.or.ihrd.eq.14.or.ihrd.eq.15.or.ihrd.eq.16.
     &   or.ihrd.eq.17.or.ihrd.eq.18) parval(16)='1'
      paruse(16,5)=1
      paruse(16,11)=1
      paruse(16,14)=1
      paruse(16,15)=1
      paruse(16,16)=1
      paruse(16,17)=1
      paruse(16,18)=1
c
      chpar(17)='zfstate'
      chpdes(17)='Final states for Z bosons in imode=1'
      partyp(17)=1
      parval(17)='0'
      paruse(17,5)=1
c
c     masses
      chpar(20)='mc'
      chpdes(20)='charm mass'
c     the charm quark is considered massless unless explicitly requested
C     (e.g. in processes like W c cbar or c cbar)
c     it is 0 in W c processes
      partyp(20)=0
      parval(20)='0'
      paruse(20,1)=1
      paruse(20,2)=1
      paruse(20,6)=1
      paruse(20,7)=1
      paruse(20,10)=1
      paruse(20,15)=1
      paruse(20,16)=1
      paruse(20,18)=1
      if(paruse(20,ihrd).eq.1) parval(20)='1.5'
c
      chpar(21)='mb'
      chpdes(21)='bottom mass'
      partyp(21)=0
      parval(21)='4.7d0'   
      paruse(21,1)=1
      paruse(21,2)=1
      paruse(21,6)=1
      paruse(21,7)=1
      paruse(21,8)=1
      paruse(21,13)=1
      paruse(21,15)=1
      paruse(21,16)=1
      paruse(21,18)=1
c
      chpar(22)='mt'
      chpdes(22)='top mass'
      partyp(22)=0
      parval(22)='173.2d0'
      paruse(22,1)=1
      paruse(22,2)=1
      paruse(22,6)=1
      paruse(22,7)=1
      paruse(22,8)=1
      paruse(22,13)=1
      paruse(22,15)=1
      paruse(22,16)=1
      paruse(22,18)=1
c
      chpar(23)='mh'
      chpdes(23)='higgs mass'
      partyp(23)=0
      parval(23)='125d0'
      paruse(23,5)=1
      paruse(23,8)=1
      paruse(23,12)=1
      paruse(23,20)=1
c
c     pt cuts
      chpar(30)='ptjmin'
      chpdes(30)='minimum pt for light jets'
      partyp(30)=0
      parval(30)='20d0'
      do i=1,nactprc
        paruse(30,i)=1
      enddo
c
      chpar(31)='ptbmin'
      chpdes(31)='ptmin for bottom quarks (in procs with explicit b)'
      partyp(31)=0
      parval(31)='20d0'
      paruse(31,1)=1
      paruse(31,2)=1
      paruse(31,6)=1
      paruse(31,7)=1
      paruse(31,8)=1
      paruse(31,13)=1
      paruse(31,15)=1
      paruse(31,16)=1
      paruse(31,18)=1
c
      chpar(32)='ptcmin'
      chpdes(32)='ptmin for charm quarks (in procs with explicit c)'
      partyp(32)=0
      parval(32)='20d0'
      paruse(32,1)=1
      paruse(32,2)=1
      paruse(32,6)=1
      paruse(32,7)=1
      paruse(32,10)=1
      paruse(32,15)=1
      paruse(32,16)=1
      paruse(32,18)=1
c
      chpar(33)='ptlmin'
      chpdes(33)='minimum pt for charged leptons'
      partyp(33)=0
      parval(33)='0d0'
      do i=1,4
        paruse(33,i)=1
      enddo
      paruse(33,10)=1
      paruse(33,12)=1
      paruse(33,14)=1
      paruse(33,15)=1
      paruse(33,17)=1
      paruse(33,18)=1
      paruse(33,20)=1
c
      chpar(34)='metmin'
      chpdes(34)='minimum missing et'
      partyp(34)=0
      parval(34)='0d0'
      do i=1,4
        paruse(34,i)=1
      enddo
      paruse(34,10)=1
      paruse(34,14)=1
      paruse(34,15)=1
      paruse(34,17)=1
      paruse(34,18)=1
      paruse(34,20)=1
c
      chpar(35)='ptphmin'
      chpdes(35)='minimum pt for photons'
      partyp(35)=0
      parval(35)='20d0'
      paruse(35,5)=1
      paruse(35,11)=1
      paruse(35,14)=1
      paruse(35,15)=1
      paruse(35,16)=1
      paruse(35,17)=1
      paruse(35,18)=1
c
      chpar(36)='ptcen'
      chpdes(36)='min pt for central jet in VBF 3-jet final states'/
     $ /' (used if irapgap = 1 and njets= 3)'
      partyp(36)=0
      parval(36)=parval(30)
      paruse(36,4)=1
      paruse(36,5)=1
      paruse(36,12)=1
      paruse(36,20)=1
c
c     eta cuts
      chpar(40)='etajmax'
      chpdes(40)='max|eta| for light jets'
      partyp(40)=0
      parval(40)='2.5'
      do i=1,nactprc
        paruse(40,i)=1
      enddo
c
      chpar(41)='etabmax'
      chpdes(41)='max|eta| for b quarks (in procs with explicit b)'
      partyp(41)=0
      parval(41)='2.5'
      paruse(41,1)=1
      paruse(41,2)=1
      paruse(41,6)=1
      paruse(41,7)=1
      paruse(41,8)=1
      paruse(41,13)=1
      paruse(41,15)=1
      paruse(41,16)=1
      paruse(41,18)=1
c
      chpar(42)='etacmax'
      chpdes(42)='max|eta| for c quarks (in procs with explicit c)'
      partyp(42)=0
      parval(42)='2.5'
      paruse(42,1)=1
      paruse(42,2)=1
      paruse(42,6)=1
      paruse(42,7)=1
      paruse(42,15)=1
      paruse(42,10)=1
      paruse(42,16)=1
      paruse(42,18)=1
c
      chpar(43)='etalmax'
      chpdes(43)='max abs(eta) for charged leptons'
      partyp(43)=0
      parval(43)='10d0'
      do i=1,4
        paruse(43,i)=1
      enddo
      paruse(43,10)=1
      paruse(43,12)=1
      paruse(43,14)=1
      paruse(43,15)=1
      paruse(43,17)=1
      paruse(43,18)=1
      paruse(43,20)=1
c
      chpar(44)='etaphmax'
      chpdes(44)='max abs(eta) for photons'
      partyp(44)=0
      parval(44)='2.5'
      paruse(44,5)=1
      paruse(44,11)=1
      paruse(44,14)=1
      paruse(44,15)=1
      paruse(44,16)=1
      paruse(44,17)=1
      paruse(44,18)=1
c
      chpar(45)='irapgap'
      chpdes(45)='enable central rap-gap in VBF >=2 jet events'
      partyp(45)=1
      parval(45)='0'
      paruse(45,4)=1
      paruse(45,5)=1
      paruse(45,12)=1
      paruse(45,20)=1
c
      chpar(46)='etagap'
      chpdes(46)='min rap for 2 "fwd" jets in VBF >=2 jet events'/
     $ /' (used if irapgap = 1)'
      partyp(46)=0
      parval(46)='2.5'
      paruse(46,4)=1
      paruse(46,5)=1
      paruse(46,12)=1
      paruse(46,20)=1
c
c     4 lepton final state specification
      chpar(47)='i4l'
      chpdes(47)='4 lepton final state according to the following: '/
     $ /'[0]  mu+ mu- e+ e-'/
     $ /'[1]  mu+ mu- mu+ mu-'/
     $ /'[2]  tau+ tau- e+ e-'/
     $ /'[3]  tau+ tau- tau+ tau-'/
     $ /'[4]  nu_mu bar nu_mu e+ e-'/
     $ /'[5]  nu_mu bar nu_mu mu+ mu-'/
     $ /'[6]  nu_mu bar nu_mu nu_e tau+ tau-'/
     $ /'[7]  nu_tau bar nu_tau nu_e tau+ tau-'/
     $ /'[8]  nu_mu bar nu_mu nu_e bar nu_e'/
     $ /'[9]  nu_mu bar nu_mu nu_mu bar nu_mu'/
     $ /'[10] nu_mu mu+ nu_e bar e-'
      partyp(47)=1
      parval(47)='0'
      paruse(47,20)=1
c
c     isolation cuts
      chpar(50)='drjmin'
      chpdes(50)='min deltaR(j-j), deltaR(Q-j) [j=light jet, Q=c/b]'
      partyp(50)=0
      parval(50)='0.7'
      do i=1,nactprc
        paruse(50,i)=1
      enddo
c
      chpar(51)='drbmin'
      chpdes(51)='min deltaR(b-b) (procs with explicit b)'
      partyp(51)=0
      parval(51)='0.7'
      paruse(51,1)=1
      paruse(51,2)=1
      paruse(51,6)=1
      paruse(51,7)=1
      paruse(51,8)=1
      paruse(51,13)=0
      paruse(51,15)=1
      paruse(51,16)=1
      paruse(51,18)=1
c
      chpar(52)='drcmin'
      chpdes(52)='min deltaR(c-c) (procs with explicit charm)'
      partyp(52)=0
      parval(52)='0.7'
      paruse(52,1)=1
      paruse(52,2)=1
      paruse(52,6)=1
      paruse(52,7)=1
      paruse(52,8)=1
      paruse(52,15)=1
      paruse(52,16)=1
      paruse(52,18)=1
c      paruse(52,10)=1
c
      chpar(55)='drlmin'
      chpdes(55)='min deltaR between charged lepton and light jets'
      partyp(55)=0
      parval(55)='0d0'
      do i=1,4
        paruse(55,i)=1
      enddo
      paruse(55,10)=1
      paruse(55,14)=1
      paruse(55,15)=1
      paruse(55,17)=1
      paruse(55,18)=1
      paruse(55,20)=1
c
      chpar(56)='drphjmin'
      chpdes(56)='min deltaR between photon and light jets'
      partyp(56)=0
      parval(56)='0.7'
      paruse(56,5)=1
      paruse(56,11)=1
      paruse(56,14)=1
      paruse(56,15)=1
      paruse(56,16)=1
      paruse(56,17)=1
      paruse(56,18)=1
c
      chpar(57)='drphlmin'
      chpdes(57)='min deltaR between photon and charged lepton'
      partyp(57)=0
      parval(57)='0.4'
      paruse(57,14)=1
      paruse(57,15)=1
      paruse(57,17)=1
      paruse(57,18)=1
c
      chpar(58)='drphmin'
      chpdes(58)='min deltaR between photons'
      partyp(58)=0
      parval(58)='0.7'
      paruse(58,5)=1
      paruse(58,11)=1
      paruse(58,14)=1
      paruse(58,15)=1
      paruse(58,16)=1
      paruse(58,17)=1
      paruse(58,18)=1
c
c     dilepton cuts
      chpar(60)='ilep'
      chpdes(60)='Z*/gamma fin state: 0=lept (1 family) 1=nu (3 fam)'
      partyp(60)=1
      parval(60)='0'
      paruse(60,2)=1
      paruse(60,4)=1
      paruse(60,17)=1
      paruse(60,18)=1
c
      chpar(61)='mllmin'
      chpdes(61)='min dilepton inv mass'
      partyp(61)=0
      parval(61)='40d0' 
      paruse(61,2)=1
      paruse(61,4)=1
      paruse(61,17)=1
      paruse(61,18)=1
      paruse(61,20)=1
c
      chpar(62)='mllmax'
      chpdes(62)='max dilepton inv mass'
      partyp(62)=0
      parval(62)='200d0'
      paruse(62,2)=1
      paruse(62,4)=1
      paruse(62,17)=1
      paruse(62,18)=1
      paruse(62,20)=1
c
c     pt cuts for leading jet/heavy quark
      chpar(65)='ptj1min'
      chpdes(65)='minimum pt of hardest light jet'
      partyp(65)=0
      parval(65)='0d0'
      do i=1,nactprc
        paruse(65,i)=1
      enddo
c
      chpar(66)='ptj1max'
      chpdes(66)='maximum pt of hardest light jet'
      partyp(66)=0
      parval(66)='1d8'
      do i=1,nactprc
        paruse(66,i)=1
      enddo
c
      chpar(67)='ptq1min'
      chpdes(67)='minimum pt for hardest heavy quark'
      partyp(67)=0
      parval(67)='0d0'
      paruse(67,1)=1
      paruse(67,2)=1
      paruse(67,6)=1
      paruse(67,7)=1
      paruse(67,8)=1
      paruse(67,13)=1
      paruse(67,15)=1
      paruse(67,16)=1
      paruse(67,18)=1
c
      chpar(68)='ptq1max'
      chpdes(68)='maximum pt for hardest heavy quark'
      partyp(68)=0
      parval(68)='1d8'
      paruse(68,1)=1
      paruse(68,2)=1
      paruse(68,6)=1
      paruse(68,7)=1
      paruse(68,10)=1
      paruse(68,15)=1
      paruse(68,16)=1
      paruse(68,18)=1
c
      chpar(69)='ptph1min'
      chpdes(69)='minimum pt for hardest photon'
      partyp(69)=0
      parval(69)='0d0'
      paruse(69,5)=1
      paruse(69,11)=1
      paruse(69,14)=1
      paruse(69,15)=1
      paruse(69,16)=1
      paruse(69,17)=1
      paruse(69,18)=1
c
      chpar(70)='ptph1max'
      chpdes(70)='maximum pt for hardest photon'
      partyp(70)=0
      parval(70)='1d8'
      paruse(70,5)=1
      paruse(70,11)=1
      paruse(70,14)=1
      paruse(70,15)=1
      paruse(70,16)=1
      paruse(70,17)=1
      paruse(70,18)=1
c
c     pt cuts for softest jet
      chpar(71)='ptjsmin'
      chpdes(71)='minimum pt of softest light jet'
      partyp(71)=0
      parval(71)='0d0'
      do i=1,nactprc
        paruse(71,i)=1
      enddo
c
      chpar(72)='ptjsmax'
      chpdes(72)='maximum pt of softest light jet'
      partyp(72)=0
      parval(72)='1d8'
      do i=1,nactprc
        paruse(72,i)=1
      enddo
c
c     seeds
      chpar(90)='iseed1'
      chpdes(90)='first random number seed (5-digit integer)'
      partyp(90)=1
      parval(90)='12345'
      do i=1,nactprc
        paruse(90,i)=1
      enddo

c
      chpar(91)='iseed2'
      chpdes(91)='second random number seed (5-digit integer)'
      partyp(91)=1
      parval(91)='67890'
      do i=1,nactprc
        paruse(91,i)=1
      enddo
c
c
c     anomalous couplings
*     changing cosvma between -1 and 1 gives a coupling of the form 
*     cosvma*(V-A) + sinvma*(V+A)   (sinvma= sqrt(1-cosvma^2))
      chpar(100)='cosvma'
      chpdes(100)='top-W coupling: cosvma*(V-A)+sinvma*(V+A)'
      partyp(100)=0
      parval(100)='1d0'
      paruse(100,6)=1
      paruse(100,13)=1
      paruse(100,16)=1
c
      chpar(101)='itdec'
      chpdes(101)='forces top decays, with spin-correlations'
      partyp(101)=1
c
      if(ihrd.eq.1.or.ihrd.eq.2.or.
     +   ihrd.eq.6.or.ihrd.eq.7.or.ihrd.eq.8.or.
     +   ihrd.eq.13.or.ihrd.eq.15.or.ihrd.eq.16.or.
     +   ihrd.eq.18) then
        parval(101)='1'
      else 
        parval(101)='0'
      endif 
      paruse(101,1)=1
      paruse(101,2)=1
      paruse(101,6)=1
      paruse(101,7)=1
      paruse(101,8)=1
      paruse(101,13)=1
      paruse(101,15)=1
      paruse(101,16)=1
      paruse(101,18)=1
c
      chpar(102)='itopprc'
      chpdes(102)='Selection of single-top process'
      partyp(102)=1
      parval(102)='1'
      paruse(102,13)=1
c
      chpar(103)='weqch'
      chpdes(103)='Select like-sign W pairs in vbjet (only for nw=2)'
      partyp(103)=1
      parval(103)='0'
      paruse(103,5)=1
c
      chpar(104)='hww_coup'
      chpdes(104)='hWW kfactor: the input value multiplies the SM one'
      partyp(104)=0
      parval(104)='1d0'
      paruse(104,5)=1
c
      chpar(105)='hzz_coup'
      chpdes(105)='hZZ kfactor: the input value multiplies the SM one'
      partyp(105)=0
      parval(105)='1d0'
      paruse(105,5)=1
c
      chpar(106)='runwidth'
      chpdes(106)='option for running width [runwidth=1]'
      partyp(106)=1
      parval(106)='0'
      paruse(106,1)=1
      paruse(106,2)=1
      paruse(106,3)=1
      paruse(106,4)=1
c
      chpar(110)='xlclu'
      chpdes(110)='lambda value for ckkw alpha (match shower alpha)'
      partyp(110)=0
c
c     default to -1. If it does not get replaced with a positive value
c     in the input file, it will be set equal to xlam in alinit
      parval(110)='-1d0'
      do i=1,6
        paruse(110,i)=1
      enddo
      do i=9,12
        paruse(110,i)=1
      enddo
      do i=14,18
        paruse(110,i)=1
      enddo
      paruse(110,20)=1
c
      chpar(111)='lpclu'
      chpdes(111)='loop order for ckkw alpha (match shower alpha)'
      partyp(111)=1
c     default to -1. If it does not get replaced with a positive value
c     in the input file, it will be set equal to nloop in alinit
      parval(111)='-1'
      do i=1,6
        paruse(111,i)=1
      enddo
      do i=9,12
        paruse(111,i)=1
      enddo
      do i=14,18
        paruse(111,i)=1
      enddo
      paruse(111,20)=1
c
      do k=1,10
        j=140+k
        chpar(j)='par'//nchar(k:k)
        chpdes(j)='auxiliary parameter (real*8)'
        partyp(j)=0
        parval(j)='0'
        do i=1,nactprc
          paruse(j,i)=1
        enddo
      enddo


c parameters for run with imode=2, npar>150
      chpar(151)='iwdecmode'
      chpdes(151)='W decay modes, in imode=2'
      partyp(151)=1
      parval(151)='1'
      if(ihrd.eq.5) parval(151)='11'
      paruse(151,1)=1
      paruse(151,3)=1
      paruse(151,5)=1
      paruse(151,10)=1
      paruse(151,13)=1
      paruse(151,14)=1
      paruse(151,15)=1
c
      chpar(152)='izdecmode'
      chpdes(152)='Z->l+l- (ilep=0) decay modes, in imode=2'
      partyp(152)=1
      parval(152)='1'
      paruse(152,2)=1
      paruse(152,4)=1
      paruse(152,17)=1
      paruse(152,18)=1
c
      chpar(153)='itdecmode'
      chpdes(153)='top decay modes, in imode=2'
      partyp(153)=1
      parval(153)='11'
      paruse(153,1)=1
      paruse(153,2)=1
      paruse(153,6)=1
      paruse(153,7)=1
      paruse(153,8)=1
      paruse(153,13)=1
      paruse(153,15)=1
      paruse(153,16)=1
      paruse(153,18)=1
c
      chpar(160)='cluopt'
      chpdes(160)
     $     ='kt scale option. 1:kt propto pt, 2:kt propto mt'
      partyp(160)=1
      parval(160)='1'
      do i=1,7
        paruse(160,i)=1
      enddo
      do i=9,12
        paruse(160,i)=1
      enddo
      paruse(160,14)=1
      paruse(160,15)=1
      paruse(160,16)=1
      paruse(160,17)=1
      paruse(160,18)=1
      paruse(160,20)=1
c     seeds for unweighting
      chpar(190)='iseed3'
      chpdes(190)
     $   ='first random number seed for unweighting (5-digit integer)'
      partyp(190)=1
      parval(190)='12345'
      do i=1,nactprc
        paruse(190,i)=1
      enddo

c
      chpar(191)='iseed4'
      chpdes(191)
     $   ='second random number seed for unweighting (5-digit integer)'
      partyp(191)=1
      parval(191)='67890'
      do i=1,nactprc
        paruse(191,i)=1
      enddo
c
c output formats
      chpar(192)='ilhe'
      chpdes(192)='output data format. unw=0, lhe=1'
      partyp(192)=1
      parval(192)='0'
      do i=1,nactprc
        paruse(192,i)=1
      enddo
c hidden parameters (paruse=2 => won't be printed out)
      chpar(195)='impunw'
      chpdes(195)
     $   ='impunw=2: maxwgt<0, impunw=1: ask user, impunw=0: default'
      partyp(195)=1
      parval(195)='0'
      do i=1,nactprc
        paruse(195,i)=2
      enddo
c
      call Aldpar_external
c
c if the process is unavailable paruse= 0
      do i=1,npar
        do itmp=1,nprocs
          if(chprc(itmp)(1:6).eq.'unavai') paruse(i,itmp)=0 
        enddo
      enddo
c modified_e
      do i=1,npar
        parlen(i)=aluisl(chpar(i))
      enddo
      return
c
      entry Alppar(n)
c print definitions and current values of parameters
      if(n.eq.2.or.n.eq.4.or.n.eq.5) then
        call alugun(iunit)
        if(n.eq.5) then
          open(iunit,file='prc.list',status='unknown')
        else
          open(iunit,file='par.list',status='unknown')
        endif
      else
        iunit=6
      endif
      if(n.lt.5) then
        write(iunit,*) '------'
        write(iunit,*) 'hard process code (not to be changed):'
        write(iunit,*) 'ihrd=',ihrd
        do i=1,npar
          if(chpar(i).ne.'***'.and.paruse(i,ihrd).eq.1) then
            itmp=aluisl(chpdes(i))
            write(iunit,*) '------'
            write(iunit,*) chpdes(i)(1:itmp),':'
            itmp=aluisl(chpar(i))
            if(partyp(i).eq.0) then
              write(iunit,*) chpar(i)(1:itmp),'=',parval(i)
            else
c              write(iunit,*) chpar(i)(1:itmp),'=',int(parval(i))
              write(iunit,*) chpar(i)(1:itmp),'=',parval(i)
            endif
            if(i.eq.4) call prntsf(iunit)
            if(i.eq.5) call alhsca(ihrd,iunit)
          endif
        enddo
      elseif(n.eq.5) then
c print list of all processes
        write(iunit,*) '======'
        write(iunit,*) 'List of processes'
        do i=1,nactprc
          if(chprc(i)(1:6).ne.'unavai') then
            write(iunit,'(I2,1x,a)') i,chprc(i)
          endif
        enddo
        write(iunit,*) '======'
        write(iunit,*) 'List of parameters'
        do i=1,npar
          do j=1,nactprc
            if(paruse(i,j).eq.2) goto 490
          enddo
          if(chpar(i).ne.'***') then
            write(iunit,*) '------'
            itmp=aluisl(chpar(i))
            if(partyp(i).eq.0) then
              write(iunit,*) i,' ',chpar(i),parval(i)
            else
c              write(iunit,*) i,' ',chpar(i),int(parval(i))
              write(iunit,*) i,' ',chpar(i),parval(i)
            endif
            itmp=aluisl(chpdes(i))
            write(iunit,*) chpdes(i)(1:itmp),':'
            write(iunit,'(20(I2,1x))') (paruse(i,j),j=1,nactprc)
          endif
 490      continue
        enddo
        write(iunit,*) '======'
        write(iunit,*) 'Scale setting choices for each process:'
        do j=1,nactprc
          if(chprc(j)(1:6).eq.'unavai') goto 500
          write(iunit,*) '------'
          itmp=aluisl(chprc(j))
          write(iunit,'(I2,1x,a)') j,chprc(j)(1:itmp)
          call alhsca(j,iunit)
 500      continue 
        enddo
        write(iunit,*) '======'
        write(iunit,*) 'imode=1 process-specific documentation'
        do j=1,nactprc
          if(chprc(j)(1:6).eq.'unavai') goto 600
          write(iunit,*) '------'
          itmp=aluisl(chprc(j))
          write(iunit,'(I2,1x,a)') j,chprc(j)(1:itmp)
          call alpdoc(j,iunit,1)
 600      continue 
        enddo
        write(iunit,*) '======'
        write(iunit,*) 'imode=2 process-specific documentation'
        do j=1,nactprc
          if(chprc(j)(1:6).eq.'unavai') goto 700
          write(iunit,*) '------'
          itmp=aluisl(chprc(j))
          write(iunit,'(I2,1x,a)') j,chprc(j)(1:itmp)
          call alpdoc(j,iunit,2)
 700      continue 
        enddo
        write(iunit,*) '======'
        write(iunit,*) 'Structure function menu'
        call prntsf(iunit)
      elseif(n.eq.6) then
        do i=151,npar
          if(chpar(i).ne.'***'.and.paruse(i,ihrd).eq.1) then
            itmp=aluisl(chpdes(i))
            write(iunit,*) '------'
            write(iunit,*) chpdes(i)(1:itmp),':'
            itmp=aluisl(chpar(i))
            if(partyp(i).eq.0) then
              write(iunit,*) chpar(i)(1:itmp),'=',parval(i)
            else
CCCCCCCMODIFIED
c              write(iunit,*) chpar(i)(1:itmp),'=',int(parval(i))
              write(iunit,*) chpar(i)(1:itmp),'=',parval(i)
CCCCCCCMODIFIED
            endif
          endif
        enddo
      endif
      if(imode.lt.5) then
        write(6,*) ' '
        call alpdoc(ihrd,iunit,imode)
      endif
      if(iunit.ge.10) close(iunit)
      end

c-------------------------------------------------------------------
      subroutine Alpdoc(ihrd,iunit,imode)
c     print information regarding hard process ihrd
c-------------------------------------------------------------------
      implicit none
      integer ihrd,iunit,imode
c modified_b 
      call Alpdoc_external(ihrd,iunit,imode)
c modified_e
      if(imode.le.1.or.imode.gt.2) then
C Z+jets and VBJET PROCESSES:
        if(ihrd.eq.4.or.ihrd.eq.5) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=0,1:'
          write(iunit,*) 
     & 'GENERATION OF RAPIDITY-GAP CONFIGURATIONS:'
          write(iunit,*) 
     & 'If njets>=2 the option exists to generate events with fwd jets.'
          write(iunit,*) 
     & 'Select irapgap=1 to force two jets to satisfy the constraints:'
          write(iunit,*) 
     & '|eta(j1)|>etagap   |eta(j2)|>etagap  eta(j1)*eta(j2)<0'
          write(iunit,*) 
     & 'In the case of njets=3, irapgap=1 assumes the third jet to be'
          write(iunit,*) 
     & 'central, with abs(eta)<etagap and with pt>ptcen'
          write(iunit,*) ' '
        endif 
        if(ihrd.eq.5) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=0,1:'
          write(iunit,*) 
     &         'SELECTION OF DECAY MODES:'
          write(iunit,*) 
     & 'the choice of decay modes for the Z bosons should be performed'
          write(iunit,*) 
     & 'already in imode=1 running (the EW couplings of the Z depend on'
          write(iunit,*) 
     & 'flavour). Input the integer string "zfstate" describing  the '
          write(iunit,*) 
     & 'decay modes of the individual Zs:'
          write(iunit,*) 
     & '1: Z->nu nubar (summed over all flavours)' 
          write(iunit,*) 
     & '2: Z->l+l- (summed over all charged leptons)'
          write(iunit,*) 
     & '3: Z->q qbar (summed over all quark flavours < top)'
          write(iunit,*) 
     & '4: Z->b bbar'
          write(iunit,*) 
     & '5: Z->all f fbar modes'
          write(iunit,*) 
     & 'Example: input "zfstate 24" for ZZ -> l+l- b bar'
          write(iunit,*) 
     & 'Example: input "zfstate 25" for ZZ -> l+l- + (Z->anything)'
          write(iunit,*) 
     & 'Example: input "zfstate 234" for ZZZ -> l+l- q qbar b bar'
          write(iunit,*) 
     & 'NB: The decay modes of the W boson are entered in imode=2'
        endif
        if(ihrd.eq.1.or.ihrd.eq.2.or.ihrd.eq.6.or.ihrd.eq.7.or.
     &     ihrd.eq.8.or.ihrd.eq.13.or.ihrd.eq.15.or.
     &     ihrd.eq.16.or.ihrd.eq.18) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=0,1:'
          write(iunit,*) 
     &         'SELECTION OF STABLE/UNSTABLE TOP:'
          write(iunit,*) 
     & 'Input the integer string "itdec":'
          write(iunit,*) 
     & '1: decaying top' 
          write(iunit,*) 
     & '2: stable top' 
          write(iunit,*) 
     & 'NB: The decay modes of the top quark are entered in imode=2'
        endif
        if(ihrd.eq.13) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=0,1:'
          write(iunit,*) 
     &         '4 single-top processes can be selected:'
          write(iunit,*)'itopprc=1: t+q (njets=0)'
          write(iunit,*)'itopprc=2: t+b (njets=0)'
          write(iunit,*)'itopprc=3: t+W(W->f fbar")+jets (njets=0,1)'
          write(iunit,*)"itopprc=4: t+b+W(W->f fbar')+jets (njets=0,1)"
        endif
        if(ihrd.eq.20) then
          write(iunit,*) ' '
          write(iunit,*) 'The following 4lepton final states available:'
          write(iunit,*) 'i4l= 0: mu+ mu- e+ e-'
          write(iunit,*) 'i4l= 1: mu+ mu- mu+ mu-'
          write(iunit,*) 'i4l= 2: tau+ tau- e+ e-'
          write(iunit,*) 'i4l= 3: tau+ tau- tau+ tau-'
          write(iunit,*) 'i4l= 4: nu_mu bar nu_mu e+ e-'
          write(iunit,*) 'i4l= 5: nu_mu bar nu_mu mu+ mu-'
          write(iunit,*) 'i4l= 6: nu_mu bar nu_mu nu_e tau+ tau-'
          write(iunit,*) 'i4l= 7: nu_tau bar nu_tau nu_e tau+ tau-'
          write(iunit,*) 'i4l= 8: nu_mu bar nu_mu nu_e bar nu_e'
          write(iunit,*) 'i4l= 9: nu_mu bar nu_mu nu_mu bar nu_mu'
          write(iunit,*) 'i4l= 10: nu_mu mu+ nu_e bar e-'
       endif
      elseif(imode.eq.2.or.imode.gt.2) then
        write(iunit,*) ' '
        if(ihrd.eq.1.or.ihrd.eq.3.or.ihrd.eq.10) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=2:'
          write(iunit,*) 'select W decay modes, set iwdecmode to:'
          write(iunit,*) '1: e nu'
          write(iunit,*) '2: mu nu'
          write(iunit,*) '3: tau nu'
          write(iunit,*) '4: e/mu/tau nu'
          write(iunit,*) '5: q q''bar'
          write(iunit,*) '6: fully inclusive'
        elseif(ihrd.eq.2.or.ihrd.eq.4.or.ihrd.eq.17.or.ihrd.eq.18) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=2 and ilep=0:'
          write(iunit,*)
     $         'select Z->l+l- decay modes, set izdecmode to:'
          write(iunit,*) '1: e e'
          write(iunit,*) '2: mu mu'
          write(iunit,*) '3: tau tau'
          write(iunit,*) '4: ell+ ell-'
        elseif(ihrd.eq.5) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=2:'
          write(iunit,*)
     $         'select decay mode for each W, set iwdecmode as:' 
          write(iunit,*) '1: e nu'
          write(iunit,*) '2: mu nu'
          write(iunit,*) '3: tau nu'
          write(iunit,*) '4: e/mu/tau nu '
          write(iunit,*) '5: q qbar'' '
          write(iunit,*) '6: fully inclusive'
          write(iunit,*)
     $         'E.g.: input iwdecmode 24 for WW-> mu l nu_mu nu_l' 
          write(iunit,*
     $         )'E.g.: input iwdecmode 45 for WW-> l nu q qbar'''  
          write(iunit,*)
     $     'E.g.: input iwdecmode 126 for WWW-> e nue mu nu(mu) ff''bar'
        elseif(ihrd.eq.13) then
          write(iunit,*) 'if itopprc.ge.3, then'
          write(iunit,*) 'select W decay modes, set iwdecmode to:'
          write(iunit,*) '1: e nu'
          write(iunit,*) '2: mu nu'
          write(iunit,*) '3: tau nu'
          write(iunit,*) '4: e/mu/tau nu'
          write(iunit,*) '5: 2 jets'
          write(iunit,*) '6: fully inclusive'
        elseif(ihrd.eq.14.or.ihrd.eq.15) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=2:'
          write(iunit,*) 'select W decay modes, set iwdecmode to:'
          write(iunit,*) '1: e nu'
          write(iunit,*) '2: mu nu'
          write(iunit,*) '3: tau nu'
          write(iunit,*) '4: e/mu/tau nu'
        endif
        if(ihrd.eq.1.or.ihrd.eq.2.or.ihrd.eq.6.or.ihrd.eq.7.or.
     &     ihrd.eq.8.or.ihrd.eq.13.or.ihrd.eq.15.or.
     &     ihrd.eq.16.or.ihrd.eq.18) then
          write(iunit,*) ' '
          write(iunit,*) 'For imode=2:'
          write(iunit,*) 'select top decay modes, set itdecmode to:'
          write(iunit,*) '1: e nu b'
          write(iunit,*) '2: mu nu b'
          write(iunit,*) '3: tau nu b'
          write(iunit,*) '4: e/mu/tau nu b'
          write(iunit,*) '5: b + 2 jets'
          write(iunit,*) '6: fully inclusive'
          write(iunit,*) 'One choice for each decaying (anti-)top'
          write(iunit,*) 'Example: input "itdecmode 23" '
          write(iunit,*) 'for t tbar -> mu nu b tau nu b'
       endif
      endif
      end

c-------------------------------------------------------------------
      subroutine Alspar
c     set list of parameters types and assign default values
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      integer i
c
      read(parval(2),*)   ih2
      read(parval(3),*)   ebeam
      read(parval(4),*)      ndns
      read(parval(5),*)      iqopt
      read(parval(6),*)       qfac
      read(parval(7),*)      ickkw
      read(parval(8),*)      ktfac
      read(parval(10),*)      njets
      read(parval(11),*)      ihvy
      read(parval(12),*)      ihvy2
      read(parval(13),*)      nw
      read(parval(14),*)      nz
      read(parval(15),*)      nh
      read(parval(16),*)      nph
      read(parval(17),*)      zfstate
      read(parval(20),*)      mc
      read(parval(21),*)      mb
      read(parval(22),*)      mt
      read(parval(23),*)      mh
      read(parval(30),*)      ptjmin
      read(parval(31),*)      ptbmin
      read(parval(32),*)      ptcmin
      read(parval(33),*)      ptlmin
      read(parval(34),*)      metmin
      read(parval(35),*)      ptphmin
      read(parval(36),*)      ptcen
      read(parval(40),*)      etajmax
      read(parval(41),*)      etabmax
      read(parval(42),*)      etacmax
      read(parval(43),*)      etalmax
      read(parval(44),*)      etaphmax
      read(parval(45),*)      irapgap
      read(parval(46),*)      etagap
      read(parval(47),*)      i4l
      read(parval(50),*)      drjmin
      read(parval(51),*)      drbmin
      read(parval(52),*)      drcmin
      read(parval(55),*)      drlmin
      read(parval(56),*)      drphjmin
      read(parval(57),*)      drphlmin
      read(parval(58),*)      drphmin
      read(parval(60),*)      ilep
      read(parval(61),*)      mllmin
      read(parval(62),*)      mllmax
      read(parval(65),*)      ptj1min
      read(parval(66),*)      ptj1max
      read(parval(67),*)      ptQ1min
      read(parval(68),*)      ptQ1max
      read(parval(69),*)      ptph1min
      read(parval(70),*)      ptph1max
      read(parval(71),*)      ptjsmin
      read(parval(72),*)      ptjsmax
      read(parval(90),*)      iseed(1)
      read(parval(91),*)      iseed(2)
      read(parval(100),*)      cosvma
      read(parval(101),*)      itdec
      read(parval(102),*)      itopprc
      read(parval(103),*)      weqch
      read(parval(104),*)      hww_coup
      read(parval(105),*)      hzz_coup
      read(parval(106),*)      runwidth
      read(parval(110),*)      xlclu
      read(parval(111),*)      lpclu
      do i=1,10
        read(parval(140+i),*)      xpar(i)
      enddo
c
c parameters for imode=2, n>150
      read(parval(151),*)      iWdecmode
      read(parval(152),*)      iZdecmode
      read(parval(153),*)      itdecmode
      read(parval(160),*)      cluopt
      read(parval(190),*)      iseed2(1)
      read(parval(191),*)      iseed2(2)
      read(parval(192),*)      ilhe
      read(parval(195),*)      impunw
c modified_b
      call Alspar_external
c modified_e
      end

c-------------------------------------------------------------------
      subroutine alinit
c     process the input parameters and fills the couplings constants
c     needed by ALPHA
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      character*5 nversion
      common/version/nversion
      integer iasclu
      common/iascl/iasclu
c local variables
      double precision aluhiw,alfas,alfas_clu
      integer i,iret,iunit,aluisl,itmp
      character*50 ewopt(0:4)
c
      if(iewopt.eq.0) then
c     input GF, sin-thetaW, alphaem
c     calculate the rest
         gw=sqrt(4.d0*pi*aem)/stw
         mw=gw/sqrt(4d0*sqrt(2d0)*gfermi)
         ctw=sqrt(1d0-stw**2)
         mz=mw/ctw
         gbar=gw/ctw 
         g1=gbar*stw
      elseif(iewopt.eq.1) then
c     input mW, GF and sin_thetaW
c     calculate the rest
         ctw=sqrt(1-stw**2)
         mZ=mW/ctw
         gw=sqrt ( 4d0*sqrt(2d0)*gfermi* mw**2 )
         gbar=gw/ctw 
         g1=gbar*stw
      elseif(iewopt.eq.2) then
c     input mZ, alpha(em) and sin_thetaW
c     calculate the rest
         ctw=sqrt(1-stw**2)
         mw=mz*ctw 
         gw=sqrt(4.d0*pi*aem)/stw
         gbar=gw/ctw 
         g1=gbar*stw
      elseif(iewopt.eq.3) then
c     input mZ, mW and GF
c     calculate the rest
         gw=sqrt ( 4d0*sqrt(2d0)*gfermi* mw**2 )
         ctw=mw/mz
         stw=sqrt(1-ctw**2)
         gbar=gw/ctw 
         g1=gbar*stw
         aem=gw**2/(4.d0*pi)*stw**2
      elseif(iewopt.eq.4) then
c     input alpha, GF and MZ
c     calculates the rest  
         stw= sqrt((1.d0-sqrt(1.d0-2.d0*
     +                        sqrt(2.d0)*pi*aem/gfermi/mz/mz))/2.d0)
         ctw= sqrt(1.d0-stw**2)
         mw= ctw*mz
         gw= sqrt(4.d0*pi*aem)/stw
         gbar=gw/ctw
         g1=gbar*stw
         aem=gw**2/(4.d0*pi)*stw**2
      else
        write(6,*)
     $       'option for selection of EW parameters not specified'
        stop
      endif
      wwid = 3d0/16d0/pi * gw**2 * mw
      zwid = 
c     leptons
     +     3*(1-2*stw**2+4*stw**4) +
c     u,d,c,s
     +     2*3*(1-2*stw**2+20d0/9d0*stw**4) +
c     b
     +     3*sqrt(1-(2*mb/mz)**2)*(1d0-(2d0*mb/mz)**2 + (-1d0+4d0/3d0
     $     *stw**2)**2*(1d0+2d0*(mb/mz)**2))/4d0
      zwid = 1/(48*pi) * mz * gbar**2 * zwid
      hwid = aluhiw(mh,mw,mz,mt)
c initialise QCD and PDF parameters
 1    call pdfconv(ndns,nmnr,pdftyp)
      call pdfpar(nmnr,ih1,xlam,sche,nloop,iret)
      if(iret.eq.1) then
         write(6,*) 'pdf set',ndns,' not available, re-input ndns:'
         call prntsf(6)
         read(5,*) ndns
         goto 1
       endif
c       write(*,*) 'lhapdf lamda=',xlam
       as=alfas(mt**2,xlam,nloop,-1)
c       write(*,*) 'lhapdf as,xlam=',as,xlam
c glu-glu-higgs coupling in mt-> infinity limit
      ggh= as/3.d0/246.d0/pi
c      ggh= as/3.d0/246.d0/pi*(1.d0+11.d0/4.d0*as/pi)
c lambda and loop order in alphas for reweighting after clustering
c     Set to xlam and nloop if default was not changed in input
      iasclu=0
      if(xlclu.gt.0.or.lpclu.gt.0) iasclu=1
      if(xlclu.lt.0) xlclu=xlam
      if(lpclu.lt.0) lpclu=nloop
c fill in  alpha common with mass and coupling parameters
      do i=1,100
         apar(i)=0d0
      enddo
c
      apar(1)=mz
      apar(2)=zwid
      apar(3)=mw
      apar(4)=wwid
      apar(5)=mh
      apar(6)=hwid
      apar(14)=mc
      apar(15)=mb
      apar(16)=mt
      apar(51)=gw
      apar(52)=g1
      apar(53)=ctw
      apar(54)=stw
      apar(55)=gbar
      apar(61)=ggh
c
      apar(71)=hww_coup          !anomalous hww coupling, default (SM) 1.d0 
      apar(72)=hzz_coup 
c
      do i=1,24
        amass(i)=0
      enddo
      amass(4)=mc
      amass(5)=mb
      amass(6)=mt
      amass(11)=mlep(1)
      amass(13)=mlep(2)
      amass(15)=mlep(3)
      amass(23)=mz
      amass(24)=mw
c
      roots=2d0*ebeam
      s=roots*roots
      ptjmax=0.8*ebeam
      ptbmax=ptjmax
      ptcmax=ptjmax
c
      ewopt(0)='input GF, sin-thetaW, alphaem, calculate the rest'
      ewopt(1)='input mW, GF, sin-thetaW, calculate the rest'
      ewopt(2)='input mZ sin-thetaW, alphaem, calculate the rest'
      ewopt(3)='input mW, mZ, GF calculate the rest'
      write(6,*) ' '
      write(6,*) '---------------------'
      do i=1,2
        if(i.eq.1) iunit=6
        if(imode.eq.2) goto 100
        if(i.eq.2) iunit=niosta
        if(i.eq.1) write(6,*) ' ALPGEN, Version: ',nversion
        write(iunit,*) ' RUNNING PARAMETERS'
        write(iunit,*) '.....................'
        write(iunit,*) ' '
        write(iunit,*) '        Electroweak parameters:'
        write(iunit,*) 'iewopt=',iewopt,':'
        write(iunit,*) ewopt(iewopt)
        write(iunit,*) 'M(W)=',mw,' Gamma(W)=',wwid
        write(iunit,*) 'M(Z)=',mz,' Gamma(Z)=',zwid
        write(iunit,*) 'M(H)=',mh,' Gamma(H)=',hwid
        write(iunit,*) ' gW=',gw,'; sin^2(thetaW)=',stw**2,
     $       '; 1/aem(mZ)=',1/aem
        write(iunit,*) ' '
        write(iunit,*) '        Quark masses:'
        write(iunit,*) 'm(top)=',mt,' m(b)=',mb
        write(iunit,*) ' '
        write(iunit,*) '       Beams'' parameters:'
        if(ih2.eq.1) write(iunit,*) 'beam1=proton, beam2=proton'
        if(ih2.eq.-1) write(iunit,*) 'beam1=proton, beam2=antiproton'
        write(iunit,*) 'Ebeam=',ebeam,' PDF set=',pdftyp
        as=alfas(mz**2,xlam,nloop,-1)
        write(iunit,*) 'as(MZ)[nloop=',nloop,'] = ',as
        write(iunit,*) ' '
      enddo
      return
c write out parameters for shower evolution
 100  continue
      write(niosta,'(a)') '************** run parameters '
      write(niosta,9991) ihrd,' ! hard process code'
      write(niosta,9992) mc,mb,mt,mw,mz,mh,' ! mc,mb,mt,mw,mz,mh'
 9991 FORMAT(I4,A20)
 9992 FORMAT(6F8.3,A20)  
      do i=1,npar
        if(chpar(i).ne.'***'.and.paruse(i,ihrd).eq.1) then
          itmp=aluisl(chpar(i))
          write(niosta,*) i,' ',parval(i),'  ! ',chpar(i)(1:itmp)
        endif
      enddo
      write(niosta,'(a)') '************** end parameters '
c     
      end


c-------------------------------------------------------------------
      subroutine aligrd
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      integer nvar,nv
      common/psopt/nvar,nv
      real*8 totalmax,tmpfac
      common/mxfact/totalmax,tmpfac
      character*50 tmpstr
c local variables
      integer n,nit,iunit
      double precision ntmp
*
c-    If  unweighting run:
*
      if(imode.eq.2) then
         avgrew=0d0
         colgen=.true.
         evopt=.false.
         evgen=.true.
         return
      endif
*
c     input grid, if so required
      if(igrid.ne.0) then
c first of all save old grids
c grid 1:
        call alugun(iunit)
        if(igrid.eq.1) then
          call alustc(fname,'.grid1',tmpstr)
          open(unit=iunit,file=tmpstr,status='unknown')
          call readgrid(iunit)
          close(iunit)
          call alustc(fname,'.grid1-old',tmpstr)
          open(unit=iunit,file=tmpstr,status='unknown')
          do n= 1,nvar+1
            call grid1W(iunit,n)
          enddo
          write (iunit,20) totalmax,tmpfac
          close(iunit)
          call alustc(fname,'.grid1',tmpstr)
        elseif(igrid.eq.2) then
c grid 2
          call alustc(fname,'.grid2',tmpstr)
          open(unit=iunit,file=tmpstr,status='unknown')
          call readgrid(iunit)
          close(iunit)
          call alustc(fname,'.grid2-old',tmpstr)
          open(unit=iunit,file=tmpstr,status='unknown')
          do n= 1,nvar+1
            call grid1W(iunit,n)
          enddo
          write (iunit,20) totalmax,tmpfac
          close(iunit)
          call alustc(fname,'.grid2',tmpstr)
        endif
c now open the grid required for the run:
        open(unit=iunit,file=tmpstr,status='unknown')
        call readgrid(iunit)
        close (iunit)
      endif
 20   format(/,'  totalmax and tmpfac',//,2(d20.9))
*     
c-    Optimize grid with first iterations, if required:
*
      ntmp=maxev
      if(nopt.ge.1) then
         evgen=.false.          ! do NOT store event info (e.g. histograms)
         evopt=.true.           ! perform optimization procedures
         colgen=.false.         ! do NOT generate colour flows
         maxev=nopt
         do nit=1,niter
            call alegen
         enddo
      endif
*
c-    Reset generation parameters, move to event generation
*
      maxev=ntmp
      evgen=.true.              ! store event info (e.g. histograms)
      evopt=.true.              ! perform optimization procedures
      colgen=.false.            ! do NOT generate colour flows
      end


c-------------------------------------------------------------------
      subroutine alegen
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
c     commons
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      real*8  psum(1000), p2sum(1000)
      common/subproc/psum,p2sum
      integer maxn,nct,nx1,nx2
      parameter (maxn= 500)
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      integer nvar,nv
      common/psopt/nvar,nv
c     
c     local variables
      integer keeps
      real *8 djpd,djg
      real *8 totalmax,tmpfac,totalmin,ncount,ncount1,extra
      integer nmon
      data nmon/100000/
      data extra/1.5d0/
      common/mxfact/totalmax,tmpfac
      integer navgusd,maxiter,navgin
      data avgwgt/-1.d0/
      data navgin/40/,maxiter/20/
      real *8 rn1,rn2,rn3,pswgt,xlum,total,alphap(4,maxpar)
      real *8 wgcol,dummy,ranget,ranset,xrn
      integer colfl(2*maxpar)
      integer i,n,iflag,ip(maxpar),ipinv(maxpar),
     $     afl(maxpar),alphafl(maxpar),instate(2),ndummy
c     local counters
      real *8 ipscount,nev
      real *8 weff,weff1
      real *8 ipscounti(1000),ispncounti(1000),ni(1000)
      save ipscounti,ispncounti,ni
      save ipscount
      real *8 gtotal,g2total,effps,effev
      real *8 sigerr,cumavg,cumerr,cumsig
      real*8  gsum(1000),g2sum(1000),gg2sum(1000)
      data cumsig/0d0/,cumerr/0d0/,cumavg/0d0/
c     
      real *8 newwmax(5)
      real *8 rspin,evtme,evtme0,evwgt,djproc,evtme1
      logical numericalcanc
c     local file naming for grids
      character*50 tmpstr
c
      logical ewk
c
      integer iunit
c     random numbers defining the event
      integer jseed(2)
      common/js/jseed
      data jseed/2*0/
c     local weight variables
      real rwgt,rx1
      double precision nevwgt
      logical unwgt
      real *8 avgwgt_in,g2total_in,maxwgt_in
      integer mxnavg
      parameter (mxnavg=100)
      real*8 refwgt,tmp,avgsig(mxnavg)
      double precision wusr
c     local spin variables
      integer hel(100)
      integer iavg,iavg_store,iavst(mxnavg)
c     clustering reweight
      double precision rewgt,shat,alfas,alfas_clu,rn,nrewgt
      integer color(2*maxpar)     !coulor string
      common/colore/color
c     local momentum variables
      integer maxdec
      parameter (maxdec=40)
      integer nwrt,iflwrt(maxdec),icuwrt(2,maxdec)
      double precision pwrt(5,maxdec),decbr
c     
c     for debugging:    
c     character*3 pname(-12:12)
c     data pname/'nb ','e+ ',5*'   ','bb ','cb ','sb ','ub ','db ','gl '
c     $     ,'d  ','u  ','s  ','c  ','b  ',5*'   ','e- ','nu '/
c     
c     reset counters
      real*8  aextra(0:maxpar-2)
      data aextra/1.d0,1.15d0,1.3d0,1.45d0,1.6d0,
     $     1.75d0,1.9d0,2.05d0,2.2d0,10*2.2d0/
*     
      real*8 qref(4)/221.d0,0.d0,0.d0,21.d0/,qout(4),msqref/48400.d0/
      integer j1,j2,ifail(3,2)
      logical flg90
      common/f90/flg90          ! to avoid spurious divergences in temporal gauge
*     
      integer jpr
      common/jp/jpr
      save qref,msqref
ccc modified fulvio
      double precision q2modified,as1,as2
ccc modified fulvio
*     
      if(imode.le.1.and.maxev.eq.0) return
c initialise averaging over multiple color/spin configurations
      if (imode.eq.2) nprtns= maxpar-2
      maxiter= imaxiter(nprtns)
      extra  = aextra(nprtns) 
      navgin = 1
      if (maxiter.ne.1) navgin = 2*maxiter
      if(navgin.gt.mxnavg) then
        write(6,*)
     $       'bookeeping of multiple colour/spin calls not properly',
     $       ' initialised: '
        write(6,*)
     $       'number of calls exceeds ',mxnavg 
        write(6,*) 'Redimension array ``avgsig'' in routine ``alegen'''
        stop
      endif
      gtotal=0d0
      g2total=0d0
c     reset counters for individual processes to 0
      do i=1,1000
        gsum(i)=0d0
        g2sum(i)=0d0
      enddo
      ipscount=0d0
      do i= 1,jprocmax
        ipscounti(i) = 0d0
        ispncounti(i)= 0d0
        ni(i)        = 0d0
      enddo
      maxwgt=-1d0
      unwev=0
      nevwgt=0d0
      nrewgt=0d0
c     if writing events to a file for later unweighting,
c     dump the grids first
      if(imode.eq.1.and.evgen) then
        do n= 1,nvar+1
          call grid1W(niopar,n)
        enddo
        write (niopar,20) totalmax,tmpfac
      endif
c
      call alustc(fname,'.err',tmpstr)
      call alugun(nioerr)
      open(unit=nioerr,file=tmpstr,status='unknown',access='append')
        do n= 1,nvar+1
          call grid1W(nioerr,n)
        enddo
        write (nioerr,20) totalmax,tmpfac
      close(nioerr)
c     if reading events from a file for unweighting,
c     input the grids first
      if(imode.eq.2) then
        call readgrid(niopar)
        read(niopar,*)
        read(niopar,10) maxev,avgwgt_in,g2total_in,maxwgt_in
        read(niopar,*)
        read(niopar,'(e12.4)') tmp
c user sets the maximum weight:
        read(niopar,'(5(e12.6,1x))') newwmax
        do i=1,5
          if(newwmax(i).gt.maxwgt_in) newwmax(i)=maxwgt_in
        enddo
        if(impunw.eq.1) then
          write(6,*) ' '
          write(6,*)
     $         'Expected approx # of unweighted events:',int(avgwgt_in
     $         /maxwgt_in*tmp)
 8        FORMAT(1X,A,100(/,1X,A))
          write(6,8)
     $         'To improve unweigthing efficiency, you can rescale',
     $         'the true maximum weight.',
     $         'This may slightly bias a fraction of the events.',
     $  'The following factors lead to the quoted eff improvements',
     $  'affecting at most the quoted fraction of events:'
          write(6,'(a)')
     $'bias:      1%          2%           3%           4%           5%'  
          write(6,'(a,5(e12.6,1x))') 'rescale: ',(newwmax(i)/maxwgt_in,i
     $         =1,5)
          write(6,'(a,5(e12.6,1x))') 'eff impr:',(maxwgt_in
     $         /newwmax(i),i=1,5)
 9        write(6,*) ' '
          write(6,*) 'input desired increase in efficiency',
     $         ' (1 for the default max wgt)'
          read(5,*) maxwgt
c     
c     keep fraction of affected events <5% :
          if(maxwgt.gt.maxwgt_in/newwmax(5)) maxwgt=maxwgt_in/newwmax(5)
c     keep efficiency boost < 100, not to be biased by pre-unweighting:
          maxwgt=min(maxwgt,1d2)
          write(6,*) 'implemented increase in efficiency =',maxwgt
          maxwgt=maxwgt_in/maxwgt
        elseif(impunw.eq.0) then
c          write(6,'(a,5(e12.6,1x))') 'rescale: ',(newwmax(i)/maxwgt_in,i
c     $         =1,5)
c          write(6,'(a,5(e12.6,1x))') 'eff impr:',(maxwgt_in
c     $         /newwmax(i),i=1,5)
          maxwgt=maxwgt_in/newwmax(5)          
          maxwgt=min(maxwgt,1d2)
c          write(6,*) 'implemented increase in efficiency =',maxwgt
          maxwgt=maxwgt_in/maxwgt
        elseif(impunw.eq.2) then
          maxwgt=-1
        else
          write(6,*) 'unspecified choice of maximum weight setting'
          stop
        endif
      endif 
 10   format(f15.1,1x,3(e12.6,1x))
 11   format(e12.6)
 20   format(/,'  totalmax and tmpfac',//,2(d20.9))
c     generate events
      write(6,*) ' '
      if(imode.le.1) then
        write(6,*) 'starting generation of',maxev,' events'
        write(niosta,*) ' '
        write(niosta,*) 'starting generation of',maxev,' events'
      else
        write(6,*) 'starting scan/unweighting of',maxev,' events'
      endif 
      iflag= 0
      n=0
      ncount= 0.d0
      ncount1= 0.d0
c     counters for PDF lumi . LE . 0 events
      ifail(3,1)=0
      ifail(3,2)=3
c     start loop over events, between 99 and 100
      nev=1d0
 99   continue
      n=n+1
c     reset local integer counter, to avoid overruns, Used only for
C     monitoring
      if(n.gt.1000000000) n=n-1000000000
c     generate kinematical configuration, return phase-space weight

 1    if(imode.eq.0) then
        dummy= ranget(jseed)
        continue
      elseif(imode.eq.1) then
        dummy= ranget(jseed)
      elseif(imode.eq.2) then 
        read(niowgt,15,end=1000) jseed,iavg_store,rwgt,rx1
        evwgt=dble(rwgt)
        if(.not.unwgt(evwgt,maxwgt)) goto 100
        dummy= ranset(jseed)
        unwev=unwev+1
      endif
 15   format(i12,1x,i12,1x,i4,2(1x,e12.6))
c
      if(mod(n,1000).eq.0) then
        call alustc(fname,'.tmp',tmpstr)
        call alugun(nioerr)
        open(unit=nioerr,file=tmpstr,status='unknown')
          write(nioerr,*)'n', n
          write(nioerr,*)jseed,'   jseed'
        close(nioerr)
      endif
c
c     MC over channels
c     generate jproc in range 1-jprocmax, 
c     return inverse jacobian djproc

      call onedimbin(1,jproc,djproc,1,ndummy,dummy)
      ni(jproc)= ni(jproc)+1d0
      jpr= jproc

      total= 0.d0
c     counters to avoid running into infinite loops in the cut-checking
c     routines
      do i=1,2
        ifail(i,1)=0
        ifail(i,2)=1
      enddo
 16   call phspace(iflag,pswgt,djpd,djg)
      if(iflag.eq.0) call momcheck(iflag,jseed)
      if(iflag.eq.1) then
        ipscount=ipscount+1d0
        ipscounti(jproc)= ipscounti(jproc)+1d0
        ifail(1,1)=ifail(1,1)+1
        if(ifail(1,1).gt.1e6) then
          write(*,*) 'Generation fails kin cuts',ifail(1,2),'M times'
          ifail(1,1)=0
          ifail(1,2)=ifail(1,2)+1
          if(ifail(1,2).eq.101) then
            write(*,*) '1E8 kin cut failures: too low eff, terminate'
            goto 1000
          endif
        endif
        goto 16
      endif
c     evaluate parton densities
      call setpdf
c     evaluate parton luminosities, alpha_s factor,  select flavours and
c     apply possible flavour-dependent cuts
      call selflav(jproc,xlum,afl)
      if(xlum.lt.0d0) then
        ipscount=ipscount+1d0
        ipscounti(jproc)= ipscounti(jproc)+1d0
        goto 16
      endif
c      if(xlum.eq.0) then
c        ifail(3,1)=ifail(3,1)+1
c        if(log10(real(ifail(3,1))).eq.real(ifail(3,2))) then
c          write(*,*) 'Warning: PDFlumi=0 more than',ifail(3,1)
c     $         ,' times'
cc     $         ,' xlum set to 0'
c          if(ifail(3,2).lt.9) ifail(3,2)=ifail(3,2)+1
c        endif
c      endif
c     Apply user-defined selection cuts;
c     wusr allows the user to reweight the event as a function of the
c     phase-space cuts, flavour assignment, etc
      call usrcut(iflag,wusr)
      pswgt=pswgt*wusr
      if(iflag.eq.1) then
        ipscount=ipscount+1d0
        ipscounti(jproc)= ipscounti(jproc)+1d0
        ifail(2,1)=ifail(2,1)+1
        if(ifail(2,1).gt.1e6) then
          write(*,*) 'event didnt pass user cuts',ifail(2,2),'M times'
          ifail(2,1)=0
          ifail(2,2)=ifail(2,2)+1
          if(ifail(2,2).eq.101) then
            write(*,*) '1E8 user cut failures: too low eff, terminate'
            goto 1000
          endif
        endif
        goto 16
      endif
c     events passed all kinematical cuts
      weff= ni(jproc)/(ni(jproc)+ipscounti(jproc))
c     check that weighted-events file was not contaminated
      if(imode.eq.2) then
        if(abs(rx1-real(x1))/(rx1+real(x1)).gt.1e-5) then
          write(6,*) 'Weighted-events file has problems:'
          write(6,*) 'x1 from phase-space=',real(x1)
          write(6,*) 'x1 from file =',rx1
          write(6,*) 'skip event',unwev
          unwev=unwev-1
          goto 100
        endif
      endif
c     setup input for alpha computation:
      call setalp(npart,afl,p,alphafl,alphap,instate,ip,ipinv)
c     generate colour      
      tmp=0d0
      navg    = navgin
      totalmin=-1.d0
      if (avgwgt.lt.0.d0.and.igrid.eq.0) then
        navg    = 1
        totalmax= 1.d25
        tmpfac  = 100.d0
      else
        if (.not.evgen) totalmax= avgwgt*tmpfac
      endif
      if(imode.eq.2) navg=navgin
      navgusd = 0
      keeps= 0
      numericalcanc=.false.
      do 111 iavg=1,navg
 2      call randa(rn1)
        call randa(rn2)
        if(npart.gt.10) call randa(rn3)
        call selcol(alphafl,instate,rn1,rn2,rn3,iflag)
        if(iflag.eq.1)   goto 2
c     generate spin
 17     call randa(rspin)
        call selspin(rspin,hel)
        call fltspn(alphafl,hel,instate,iflag)
        if(iflag.eq.1) then
          if(keeps.eq.0) ispncounti(jproc)= ispncounti(jproc)+1d0
          goto 17
        endif
        do j1= 1,npart
c mauro 30/05/2014
           spnwrt(j1)= hel(ipinv(j1))
           if((ifl(j1).ge.1.and.ifl(j1).le.6).or.
     .          (ifl(j1).ge.11.and.ifl(j1).le.16)) 
     .          spnwrt(j1)=-spnwrt(j1)
c mauro 30/05/2014
        enddo
c     non-zero spin configuration selected, efficiency updated
        weff1= weff*ni(jproc)/(ni(jproc)+ispncounti(jproc))     
c     
c     reevaluate alphas using ktmin as a scale
        if(ickkw.eq.1) then
c          write(6,*) ' '
c          write(6,*) (alphafl(i),i=1,npart)
c          write(6,*) (color(2*i-1),i=1,npart)
c          write(6,*) (color(2*i),i=1,npart)
          call setcol(npart,ip,ipinv,color,ifl,icu)
c          write(6,*) 'qsq=',qsq
          call cktmin
c          write(6,*) 'kres=',kres
          asmax=alfas_clu(kres,xlclu,lpclu,-1)
c     multiply alphas by a rescaling factor, to allow separation of qcd
C     and ew
c     contributions to jet production: 
c     
c     sig=1/(resc**N) * sig(as->as*resc)
c     
c     for processes where maximum power of alphas is N.
c     resc>>1 will suppress ew processes
c     resc<<1 will suppress qcd processes
          apar(56)=sqrt(4*pi*as*resc)
        endif
        
        if(imode.eq.2.and.iavg.eq.iavg_store) then 
          total=evwgt
c debug
c          write(2,*) ' '
c          write(2,*) (alphafl(i),i=1,npart)
c          write(2,*) (color(2*i-1),i=1,npart)
c          write(2,*) (color(2*i),i=1,npart)
c          write(2,*) 'ktmin=',sqrt(kres),' Q=',sqrt(qsq)
c          write(2,*) 'as(Q)/as(kres)=',alfas(qsq,xlam,nloop,-1)
c     $         /alfas(kres,xlam,nloop,-1) 
c end debug
          goto 112
        elseif(imode.le.1) then
          resonance= 'n'
          flg90=.false.
          call matrix(alphafl,instate,alphap,hel,evtme,0,wgcol
     $         ,colfl)
          if (resonance.eq.'y') then
c            print*,'taglio',evtme
            evtme= 0.d0
          endif
          if(flg90.and.resonance.ne.'y') then
            do j1=1,npart
              call boost(1,msqref,qref,qout,alphap(1,j1))
              do j2=1,4
                alphap(j2,j1)=qout(j2)
              enddo
            enddo
            flg90=.false.
            call matrix(alphafl,instate,alphap,hel,evtme,0,wgcol
     +           ,colfl)
          endif
          ncount= ncount+1
          ncount1= ncount1+1
c     
c     include MC sum over jprocs, gev->pb
          total = dble(jprocmax)/djproc * pbfac * evtme
c     combine with parton luminosities, phase-space weight and
c     efficiency
          total = total*xlum*pswgt*weff1
c
c     rescale alphas for the light jets
          if(ickkw.eq.1) total=total*(asmax/as)**naspow
          if ((iavg.eq.1).and.
     .         (total.le.totalmax).and. 
     .         (total.ge.totalmin)) then
            navgusd= navgusd+1
            avgsig(navgusd)=total
            iavst(navgusd)=iavg
            tmp=tmp+total
            goto 115
          else
            keeps= 1
            if ((total.ge.totalmax).or.
     .           (total.le.totalmin)) then
              navgusd= navgusd+1
              avgsig(navgusd)=total
              iavst(navgusd)=iavg
              tmp=tmp+total
              if(total.ge.1.d2*totalmax) then
                call setgauge(0.d0)
                call matrix(alphafl,instate,alphap,hel,evtme1,0,wgcol
     $               ,colfl)
                call setgauge(1.d0)
                if(abs(evtme1-evtme)/evtme.gt.1.d-4) then
                   call alustc(fname,'.err',tmpstr)
                   call alugun(nioerr)
                   open(unit=nioerr,file=tmpstr,status='unknown',
     +                  access='append')
                     if(.not.numericalcanc) write(nioerr,*)'  *** '
                     if(.not.numericalcanc) write(nioerr,*)'  '
                     if(.not.numericalcanc) write(nioerr,*)'jseed',jseed
                     if(.not.numericalcanc) write(nioerr,*)'  *** '
                     write(nioerr,*)'evtme', evtme
                     write(nioerr,*)'evtme xi=0', evtme1
                     write(nioerr,*)
                   close(nioerr)
                   numericalcanc=.true.
                   evtme= 0.d0
                endif
              endif
              if (navgusd.ge.maxiter) goto 115 
c     else
c     cycle
            endif
          endif
        endif     
 111  continue
c     omment  navg -> navgusd
 115  total=tmp/dble(navgusd)
      if(imode.eq.1) then
        call randa(rn1)
        rn1=rn1*tmp
        tmp=0d0
        do iavg=1,navgusd
          tmp=tmp+avgsig(iavg)
          if(tmp.ge.rn1) then
            iavg_store=iavst(iavg)
            goto 112
          endif
        enddo
      endif
 112  continue
c
c     generate color flow, if required
      if(colgen) then
        call randa(wgcol) 
        flg90=.false.
        ewk=.false.
        if(ihrd.eq.5.and.ickkw.eq.1.and.njets.ge.2) then
!          resc0=resc
!          resc=1.d0
          apar(56)=sqrt(4*pi*as)
          call matrix(alphafl,instate,alphap,hel,evtme0,0
     >       ,wgcol,colfl)
          apar(56)=sqrt(4*pi*as*1.d6)
!          resc=1.d6
          call matrix(alphafl,instate,alphap,hel,evtme,0
     >       ,wgcol,colfl)
          evtme0=abs(evtme0-evtme*1.d-6**(njets))/evtme0
          if(abs(evtme0).gt.0.5d0) then
            ewk=.true.
          else
            ewk=.false.
          endif 
          apar(56)=sqrt(4*pi*as*resc)
        endif
        call matrix(alphafl,instate,alphap,hel,evtme0,1
     >       ,wgcol,colfl)
        if(flg90) then
          call selcol(alphafl,instate,rn1,rn2,rn3,iflag)
          do j1=1,npart
            call boost(1,msqref,qref,qout,alphap(1,j1))
            do j2=1,4
              alphap(j2,j1)=qout(j2)
            enddo
          enddo
          flg90=.false.
          call matrix(alphafl,instate,alphap,hel,evtme,0,wgcol
     +         ,colfl)
        endif
        call setcol(npart,ip,ipinv,colfl,ifl,icu)
        call colchk(npart,icu,ifl,nev,total)
      endif
      gtotal=gtotal+total                     
      g2total=g2total+total**2
      gsum(jproc) =gsum(jproc)+total
      g2sum(jproc)=g2sum(jproc)+total**2
      if(imode.le.1) maxwgt=max(maxwgt,total)
      if(evgen) then
        if(imode.eq.0) then 
          call aleana(jproc,total)
        elseif(imode.eq.1) then 

c     write events after a pre-unweighting, using as maximum weight 
c     refwgt=1d-2*maxwgt 
          refwgt=1d-2*maxwgt
          if(unwgt(total,refwgt)) then
            write(niowgt,15) jseed,iavg_store,real(max(total,refwgt))
     $           ,real(x1)
c     if(total.gt.1d-25) write(niowgt,15) jseed,jproc,real(evtme)
c     ,real(x1)
            nevwgt=nevwgt+1d0
          endif
          call aleana(jproc,total)
          call alhbkk(total)
        elseif(imode.eq.2) then ! dump event for later Herwig processing
c     First cluster and evaluate reweighting factor
          if(ickkw.eq.1) then
            rewgt=1d0
            call alpclu(ewk,rewgt)
            avgrew=avgrew+rewgt
            nrewgt=nrewgt+1d0
            call mfill(190,sngl(rewgt),1e0)
c     next calls for debugging only
c            if(rewgt.gt.1d0) then
c              write(2,*) ' '
c              write(2,*) 'nev=',nev,' rewgt:',rewgt
c              write(2,*) (ifl(i),i=1,npart)
c              write(2,*) (icu(1,i),i=1,npart)
c              write(2,*) (icu(2,i),i=1,npart)
c              write(2,*) (sqrt(clkt(2,i)),i=1,ncltot)
c            endif
c            write(2,*)
c     $         '   fl col1  col2  cmot cdau1 cdau2  p(3)   ktin  ktout' 
c            write(2,*) 'BEFORE CLUSTERING'
c            do i=1,nclext
c             write(2,"(6(I5,1x),2(f10.2,1x))") icfl(i),iccol(1,i)
c     $             ,iccol(2,i),icmot(i),icdau(1,i),icdau(2,i)
c     $             ,sqrt(clkt(1,i)),sqrt(clkt(2,i))
c            enddo
c            write(2,*) ' '
c            write(2,*) 'CLUSTERS'
c            do i=nclext+1,ncltot
c              write(2,"(6(I5,1x),2(f10.2,1x))") icfl(i),iccol(1,i)
c     $             ,iccol(2,i),icmot(i),icdau(1,i),icdau(2,i)
c     $             ,sqrt(clkt(1,i)),sqrt(clkt(2,i))
c            enddo
c            write(2,*) ' '

c     end debug
c     3: unweight w.r.t. to reweighting factor:
            if(.not.unwgt(rewgt,1d0)) then
              unwev=unwev-1
              goto 100
            endif
          endif
c     First add to the output event record decay particles, when present
          decbr=1d0
          call setdec(nwrt,iflwrt,icuwrt,pwrt,decbr)
          call aleana(jproc,avgwgt_in*decbr)
          if(ilhe.eq.0) then
            call evdump(niounw,nwrt,jproc,iflwrt,icuwrt,pwrt,qsq,1d0)
          elseif(ilhe.eq.1) then
            call evdlhe
          else
            write(*,*) 'undefined event output format, stop'
            stop
          endif
        endif
      endif
c     weight bookkeeping and various optimizations

      if (evopt) then
        call wgtopt(total)
      endif

c     end generation of n-th event
c     Produce monitor file, to report partial results during run
      if(imode.ne.2.and.mod(n,nmon).eq.0) then
        if (.not.evgen) then
          tmpfac= tmpfac*ncount1/(extra*nmon)
          ncount1= 0.d0
        endif
        call alustc(fname,'.mon',tmpstr)
        call alugun(iunit)
        open(unit=iunit,file=tmpstr,status='unknown')
        write (iunit,51) nev
 51     format(3x,' processed=',f15.1,' events') 
        write (iunit,*) '       '
        totwgt= gtotal
        effps=nev/(nev+ipscount)
        effev=nev
        sigerr=sqrt(g2total-gtotal**2/(effev))/(effev)
        avgwgt=totwgt/effev
        write(iunit,*) 'average ph-space eff=',effps
        write(iunit,*) 'avgwgt(pb)=',avgwgt,'+-',sigerr
     $       ,' maxwgt=',maxwgt
        if(maxwgt.gt.0) write(iunit,*) 'unwgt eff=',totwgt/nev
     $       /maxwgt
        write(iunit,*) '                 '
        write(iunit,*) ' sub-processes:  '
        write(iunit,*) '                 '
        do i= 1,jprocmax
          gg2sum(i)=sqrt(g2sum(i)-gsum(i)**2/(effev))/(effev)
          write(iunit,*) 'jproc=',i,' total: ',gsum(i)/(effev),'+-'
     $         ,gg2sum(i)
          xsprc(i)=gsum(i)/effev
          errprc(i)=gg2sum(i)
        enddo
c     write(iunit,*) '               '
c     write(iunit,*) 'n. of ME calls=',ncount,avgwgt,tmpfac,totalmax
c     write(iunit,*) '               '
        close (iunit)
        call monitor(n,tmpstr)
      endif
 100  continue
      nev=nev+1d0
      if(nev.le.maxev) goto 99
c     
c     event generation completed
c     
 1000 continue
c     
c     print out final statistics
      totwgt=gtotal
      if(imode.ne.2) then
        effps=maxev/(maxev+ipscount)
        effev=maxev
        sigerr=sqrt(g2total-gtotal**2/(effev))/(effev)
        avgwgt=totwgt/effev
        if (avgwgt.eq.0.d0) then
          write(6,*) 'cross section = 0, stop'
          stop
        endif
        write(6,*) 'average ph-space eff=',effps
        write(6,*) 'avgwgt(pb)=',avgwgt,'+-',sigerr,' maxwgt=',maxwgt
        write(6,*) 'unwgt eff =',totwgt/maxev/maxwgt
        write(6,*) '                 '
        write(6,*) ' sub-processes:  '
        write(6,*) '                 '
        do i= 1,jprocmax
          gg2sum(i)=sqrt(g2sum(i)-gsum(i)**2/(effev))/(effev)
          write(6,*) 'jproc=',i,' total(pb): ',gsum(i)/(effev),'+-'
     $         ,gg2sum(i) 
          xsprc(i)=gsum(i)/effev
          errprc(i)=gg2sum(i)
        enddo
c     
c     same info, written to file
        write(niosta,*) 'average ph-space eff=',effps
        write(niosta,*) 'avgwgt(pb)=',avgwgt,'+-',sigerr,' maxwgt='
     $       ,maxwgt
        write(niosta,*) 'unwgt eff =',totwgt/maxev/maxwgt
        write(niosta,*) '                 '
        write(niosta,*) ' sub-processes:  '
        write(niosta,*) '                 '
        do i= 1,jprocmax
          gg2sum(i)=sqrt(g2sum(i)-gsum(i)**2/(effev))/(effev)
          write(niosta,*) 'jproc=',i,' total(pb): ',gsum(i)/(effev)
     $         ,'+-',gg2sum(i) 
        enddo
        if(cumerr.eq.0d0) then
          cumerr=1d0/sigerr**2
        else
          cumerr=1d0/cumerr**2+1d0/sigerr**2
        endif
        cumsig=cumsig+avgwgt/sigerr**2
        cumavg=cumsig/cumerr
        cumerr=1d0/sqrt(cumerr)
        write(niosta,*) '                 '
        write(niosta,*) 'cumulated cross-section:'
        write(niosta,*) 'avgwgt(pb)=',cumavg,'+-',cumerr
        write(6,*) '                 '
        write(6,*) 'cumulated cross-section:'
        write(6,*) 'avgwgt(pb)=',cumavg,'+-',cumerr
      elseif(imode.eq.2) then
        if(ickkw.eq.1) then 
c rescale by average reweighting factor with CKKW-like scales
          avgrew=avgrew/nrewgt
        else
          avgrew=1d0
        endif
        avgrew=avgrew*decbr
        sigerr=g2total_in*avgrew
        avgwgt=avgwgt_in*avgrew
        write(6,*)'Crosssection(pb)=',avgwgt,'+-',sigerr
        write(6,*) 'Generated ',unwev,' unweighted events, lum=',unwev
     $       /avgwgt,'pb-1' 
        write(niosta,9995) avgwgt,sigerr,' ! sigma +- error (pb)'
        write(niosta,*) unwev,unwev/avgwgt,' ! unwtd events, lum (pb-1)' 
        if(ilhe.eq.1) call lhehea(NIOSTA,ebeam,ih2,avgwgt,sigerr)
        close(niosta)
c merge files
        close(niounw)
        call mrgunw(ilhe,fname)
      endif
 9995 FORMAT(2F25.8,A30)     
c     write out updated grid, if igridw properly initialised
      if(imode.ne.2) then
        if (evgen) tmpfac= tmpfac*ncount1/(extra*maxev)
        totalmax= avgwgt*tmpfac
        if(evopt) then
          if(evgen) then
            call alustc(fname,'.grid2',tmpstr)
          else
            call alustc(fname,'.grid1',tmpstr)
          endif
          call alugun(iunit)
          open(unit=iunit,file=tmpstr,status='unknown')
          call dumpgrid(iunit)
          close (iunit)
        endif
      endif
      if(imode.eq.1.and.evgen) then
        write(niopar,*)
     $'number wgted evts in the file,  sigma(pb), error(pb),  maxwgt'
         write(niopar,10) nevwgt,cumavg,cumerr,maxwgt
         write(niopar,*)
     $'number weights'
         write(niopar,'(e12.4)') real(maxev)
      endif
c next call for debugging only
c      if(imode.ne.2) call prfinal(jprocmax,effev)
c
      end

      subroutine mrgunw(ilhe,fname)
      implicit none
      integer ilhe,status
      character*40 fname
      character*80 tmpf1,tmpf2,tmpf3
      if(ilhe.eq.0) then
        return
c        call alustc(fname,'_unw.par',tmpf1)
c        call alustc(fname,'.unw',tmpf2)
      elseif(ilhe.eq.1) then
        call alustc(fname,'.lhh',tmpf1)
        call alustc(fname,'.lhu',tmpf2)
        call alustc(fname,'.lhe',tmpf3)
      endif
      call system("cat "//tmpf1 // tmpf2 // ">" // tmpf3 , status)
      if(status.eq.0) then
        call system("rm "//tmpf1 // tmpf2 )
      else
        write(*,*)"Could not merge para
     $m's and event files at end of run"
      endif
      end

cc-------------------------------------------------------------------
c      subroutine prstat(n,jseed,jprocmax,jproc,evtme,xlum,pswgt
c     $     ,djproc,wgtmax)
cc     routine for the debugging of individual large-weight events
cc     Not functional to the running of the code
cc-------------------------------------------------------------------
c
c      implicit none
c      real*8  evtme,xlum,pswgt,djproc,wgtmax,total,effev
c      integer n,jprocmax,jproc,i     
c      real*8 psum(1000), p2sum(1000)
c      data psum/1000*0/,p2sum/1000*0/
c      integer jseed(2)
c      save 
c      total=evtme*xlum*pswgt
c      if(mod(n,100000).eq.0) then
c          open (unit= 31,file= 'monitor',status='unknown')
c          write (31,51) n
c   51     format('  event=',i10) 
c          close (31)
c      endif
c      if(total.gt.1d5) then
c         write(6,*) ' '
c         write(6,*) ' '
c         write(6,*) 'large weight: event=',n,' jproc=',jproc
c         write(6,*) 'total wgt=',total,' pswgt=',pswgt,' jproc jac='
c     $        ,1/djproc,'xlum=',xlum,' ME^2=',evtme*djproc
c 8       format(5(a,e12.4))
c         if(wgtmax.eq.total) write(niosta,*) 'NEW MAX WEIGHT'
c         write(niosta,*) 'large weight: event=',n,' jseeds=',jseed(1)
c     $        ,jseed(2),' jproc=',jproc
c         write(niosta,8) 'total wgt=',total,' pswgt=',pswgt,' jproc jac='
c     $        ,1/djproc,'xlum=',xlum,' ME^2=',evtme*djproc
c         write(6,*) ' '
c         call dumpwgt
c      endif
cc bookkeeping of the grandtotals
c      psum(jproc)=psum(jproc)+total
cc     bookkeeping of the phase-space weight
c      psum(jprocmax+jproc)=psum(jprocmax+jproc)+pswgt
cc bookkeeping of the jacobian of jproc selection
c      psum(2*jprocmax+jproc)=psum(2*jprocmax+jproc)+1/djproc
cc bookkeeping of the ME^2
c      psum(3*jprocmax+jproc)=psum(3*jprocmax+jproc)+evtme*djproc
c      p2sum(jproc)=p2sum(jproc)+total**2
c      p2sum(jprocmax+jproc)=p2sum(jprocmax+jproc)+pswgt**2
c      p2sum(2*jprocmax+jproc)=p2sum(2*jprocmax+jproc)+(1/djproc)**2
c      p2sum(3*jprocmax+jproc)=p2sum(3*jprocmax+jproc)+(evtme*djproc)
c     $     **2
c      return
c      entry prfinal(jprocmax,effev)
cc      write(6,*) 'bookkeeping of the grandtotals'
c      write(niosta,*) 'bookkeeping of the grandtotals'
c      do i=1,jprocmax
c         p2sum(i)=sqrt(p2sum(i)-psum(i)**2/(effev))/(effev)
cc         write(6,*) 'jproc=',i,' total: ',psum(i)/(effev),'+-'
cc     $        ,p2sum(i) 
c         write(niosta,*) 'jproc=',i,' total: ',psum(i)/(effev),'+-'
c     $        ,p2sum(i) 
c      enddo
cc      write(6,*) 'bookkeeping of the  phspace wgt'
c      write(niosta,*) 'bookkeeping of the  phspace wgt'
c      do i=jprocmax+1,2*jprocmax
c         p2sum(i)=sqrt(p2sum(i)-psum(i)**2/(effev))/(effev)
cc         write(6,*) 'jproc=',i-jprocmax,' pswgt: ',psum(i)/(effev),'+-'
cc     $        ,p2sum(i) 
c         write(niosta,*) 'jproc=',i-jprocmax,' pswgt: ',psum(i)/(effev),'+-'
c     $        ,p2sum(i) 
c      enddo
cc      write(6,*) 'bookkeeping of the jproc jacobian'
c      write(niosta,*) 'bookkeeping of the jproc jacobian'
c      do i=2*jprocmax+1,3*jprocmax
c         p2sum(i)=sqrt(p2sum(i)-psum(i)**2/(effev))/(effev)
cc         write(6,*) 'jproc=',i-2*jprocmax,' prcwgt: ',psum(i)/(effev)
cc     $        ,'+-',p2sum(i) 
c         write(niosta,*) 'jproc=',i-2*jprocmax,' prcwgt: ',psum(i)/(effev)
c     $        ,'+-',p2sum(i) 
c      enddo
cc      write(6,*) 'bookkeeping of the ME^2'
c      write(niosta,*) 'bookkeeping of the ME^2'
c      do i=3*jprocmax+1,4*jprocmax
c         p2sum(i)=sqrt(p2sum(i)-psum(i)**2/(effev))/(effev)
cc         write(6,*) 'jproc=',i-3*jprocmax,' ME^2: ',psum(i)/(effev)
cc     $        ,'+-',p2sum(i) 
c         write(niosta,*) 'jproc=',i-3*jprocmax,' ME^2: ',psum(i)/(effev)
c     $        ,'+-',p2sum(i) 
c      enddo
c      do i=1,100 
c         psum(i)=0d0
c         p2sum(i)=0d0
c      enddo
c      end
c


      subroutine alccut(lnot,pt,eta,dr,nfirst)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Applies kinematical cuts on light jets, common to all processes
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'alpgen.inc'
      integer lnot,nfirst,i,j
      real*8 pt(maxpar),eta(maxpar),dr(maxpar,maxpar)
      real*8 ptmin,ptmax
      lnot= 0
c
c      write(11,*) (pt(i),i=nfirst,nfirst+njets-1),dr(nfirst,nfirst+1)
      if(njets.gt.0) then
        ptmax=-1d0
        ptmin= 1d10
        do i=nfirst,nfirst+njets-1
          if (pt(i).lt.ptjmin)           goto 10
          if (abs(eta(i)).gt.etajmax)    goto 10
          ptmax=max(ptmax,pt(i))
          ptmin=min(ptmin,pt(i))
        enddo 
        if(ptmax.lt.ptj1min.or.ptmax.gt.ptj1max) goto 10
        if(ptmin.lt.ptjsmin.or.ptmin.gt.ptjsmax) goto 10
        do i=nfirst,nfirst+njets-1
          do j=i+1,nfirst+njets-1
            if(dr(i,j).lt.drjmin)  goto 10
          enddo
        enddo
      endif

 5    return
      
 10   lnot= 1
      return
      end



      subroutine wgtopt(total)
      implicit none
c inputs
      real *8 total
c common
      integer nvar,nv
      common/psopt/nvar,nv
      integer maxn
      parameter (maxn= 500)  
      integer mask,mmask
      real*8 peropt
      common/psopt1/mask(maxn),mmask(maxn),peropt(maxn)
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
c locals
      real *8 dummy,wgtt
      integer n,ndummy,nn
      data ndummy/-1/
      wgtt= total

      call onedimbin(2,mmask(1),wgtt,1,ndummy,dummy)
      nn=jgrid(jproc)
      do n= (nn-1)*nv+2,nn*nv+1
        call onedimbin(2,mmask(n),wgtt,n,ndummy,dummy)
      enddo
      end

c-------------------------------------------------------------------
      subroutine setpdf
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      real *8 alfas
      real xx1,xx2,q2,tmp
      integer i
c
      q2=real(qsq)                   
      xx1=real(x1)
      xx2=real(x2)
      call genericpdf(nmnr,ih1,q2,xx1,f1,5)
      call genericpdf(nmnr,ih2,q2,xx2,f2,5)      
c     enforce positivity of pdfs
      do i=-5,5
        if(f1(i).lt.0) f1(i)=0
        if(f2(i).lt.0) f2(i)=0
      enddo
c     adopt pdg numbering scheme for up and down quarks
      tmp=f1(1)
      f1(1)=f1(2)
      f1(2)=tmp
      tmp=f2(1)
      f2(1)=f2(2)
      f2(2)=tmp
      tmp=f1(-1)
      f1(-1)=f1(-2)
      f1(-2)=tmp
      tmp=f2(-1)
      f2(-1)=f2(-2)
      f2(-2)=tmp
c
      as=alfas(qsq,xlam,nloop,-1)
c     multiply alphas by a rescaling factor, to allow separation of qcd and ew
c     contributions to jet production: 
c
c                 sig=1/(resc**N) * sig(as->as*resc)
c
c     for processes where maximum power of alphas is N.
c     resc>>1 will suppress ew processes
c     resc<<1 will suppress qcd processes
      apar(56)=sqrt(4*pi*as*resc)
      end                                    

      subroutine grans(n1,nx1,x1)      
c----------------------------------------------------------------c
c                                                                c       
c     x1 is uniformely genetared in                              c
c     x1a < x1 < x1b.                                            c
c                                                                c
c     n1 is the bin, nx1 is the total number of bins             c
c                                                                c       
c----------------------------------------------------------------c

      implicit none
      real*8 x1,ran1,x1a,x1b
      integer n1,nx1   
      call rans(ran1)

      x1a= dfloat(n1-1)/dfloat(nx1)
      x1b= x1a+1.d0/dfloat(nx1)
      x1= x1b*ran1+(1.d0-ran1)*x1a

      return
      end

      subroutine photsm(lflag,cn,cxm,cxp,s,dj,ran1)
c----------------------------------------------------------------c
c                                                                c
c   Massless particle propagator:                                c
c                                                                c
c   the invariant mass squared has the distribution              c
c                                                                c
c                       1/s^cn                                   c
c                       cxm < s < cxp                            c
c                                                                c
c              INPUT                           OUTPUT            c
c                                                                c
c  lflag= 0:   cn, cxm, cxp                    s, dj             c
c                                                                c
c  lflag= 1:   cn, cxm, cxp , s                dj                c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 cn,cxm,cxp,s,dj,ran1,tjm,hj
      integer lflag
      dj= 0.d0
      if (lflag.eq.0) then
        s= tjm(0.d0,cn,cxm,cxp,1,ran1)
      else
        if (s.le.cxm.or.s.ge.cxp) return
      endif
      dj= 1.d0/(s**cn*hj(0.d0,cn,cxm,cxp))
      return
      end

      subroutine resonm(lflag,rm,ga,lim,cxm,cxp,s,dj,ran1)
c----------------------------------------------------------------c
c                                                                c
c   Resonant propagator:                                         c
c                                                                c
c   the invariant mass squared has the distribution              c
c                                                                c
c                    1/[(s-rm^2)^2+rm^2*ga^2]                    c
c                                                                c
c   if lim= 0 -inf < s < +inf                                    c
c   if lim= 1  cxm < s < cxp                                     c
c                                                                c
c              INPUT                           OUTPUT            c
c                                                                c
c  lflag= 0:   rm, ga, lim, (cxm, cxp)         s, dj             c
c                                                                c
c  lflag= 1:   rm, ga, lim, s, (cxm, cxp)      dj                c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 rm,ga,cxm,cxp,s,dj,ran1,pi,rms,ard,arp,arm,ym,ypmym
      integer lflag,lim    
      parameter(pi= 3.14159265358979324d0)
      dj= 0.d0
      rms= rm*rm
      if (lim.eq.0) then
        if (lflag.eq.0) then
          s= rms+rm*ga*tan(pi*(ran1-0.5d0))
        endif
        dj= 1.d0/(pi/rm/ga*((s-rms)*(s-rms)+rms*ga*ga))
      else
        if (lflag.eq.1) then
          if (s.lt.cxm.or.s.gt.cxp) return
        endif
        ard= cxp-cxm
        arp= (cxp-rms)/rm/ga
        arm= (cxm-rms)/rm/ga
        ym= atan(arm)
        ypmym= atan(ard/rm/ga/(1.d0+arp*arm))
        if (arp*arm.lt.-1.d0) then
          if (arp.gt.0.d0) ypmym= ypmym+pi
          if (arp.lt.0.d0) ypmym= ypmym-pi
        endif
        if (lflag.eq.0) then
          s= rms+rm*ga*tan(ran1*(ypmym)+ym)
        endif
        dj= 1.d0/(ypmym/rm/ga*((s-rms)*(s-rms)+rms*ga*ga))
      endif
      return
      end

      subroutine dec2fm(lflag,s,p,s1,s2,p1,p2,dj,ran)
c----------------------------------------------------------------c
c                                                                c
c   Isotropic 2-body decay:                                      c
c                                                                c
c                       p ----> p1 + p2                          c
c                       s1= p1^2                                 c
c                       s2= p2^2                                 c
c                                                                c
c              INPUT                          OUTPUT             c
c                                                                c
c  lflag= 0:   s, p, s1, s2                   p1, p2, dj         c
c                                                                c
c  lflag= 1:   s, s1, s2                      dj                 c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 pi,s,s1,s2,dj,sqlam,rs,p1m,ran1,ran2,ct,st,phi,cp,sp
      integer lflag,k    
      parameter (pi= 3.14159265358979324d0)
      real*8 p(0:3),p1(0:3),p2(0:3),p1h(0:3)
      real*8 ran(1:2)      
      dj= 2.d0/pi/sqlam(s,s1,s2)
      if (lflag.eq.0) then
        rs= sqrt(dabs(s))
        p1h(0)= (s+s1-s2)/rs/2.d0
        p1m= rs/dj/pi
        ran1= ran(1)
        ran2= ran(2)
        ct= 2.d0*ran1-1.d0
        st= sqrt(1.d0-ct*ct)
        phi= 2.d0*pi*ran2
        cp= cos(phi)
        sp= sin(phi)
        call qvec(p1h(0),p1m,st,ct,sp,cp,p1h)
        call boost(0,s,p,p1h,p1)
        do k= 0,3
           p2(k)= p(k)-p1(k)
        enddo
      endif
      return
      end

      function tjm(a,cn,cxm,cxp,k,ran1)
c----------------------------------------------------------------c
c                                                                c
c   tjm is produced according to the distribution                c
c                                                                c
c                    (1/tjm)^cn     where:                       c
c                                                                c
c                    tjm= a+k*ct                                 c  
c                    cxm < ct < cxp                              c
c                    k= +1 or -1                                 c  
c                                                                c  
c----------------------------------------------------------------c

      implicit none
      real*8 a,cn,cxm,cxp,ran1,ce,tjm
      integer k 
      ce= 1.d0-cn
      if (abs(ce).gt.1.d-8) then
        if (k.eq.1) then
           tjm= (ran1*(a+cxp)**ce+(1.d0-ran1)*
     *         (a+cxm)**ce)**(1.d0/ce)
        else
           tjm= (ran1*(a-cxm)**ce+(1.d0-ran1)*
     *         (a-cxp)**ce)**(1.d0/ce)
        endif
      else
        if (k.eq.1) then
          if(cxp.lt.-a) then
            tjm= -exp(ran1*log(-a-cxp)+
     *          (1.d0-ran1)*log(-a-cxm))
          else
            tjm= exp(ran1*log(a+cxp)+
     *          (1.d0-ran1)*log(a+cxm))
          endif
        else
          if(cxp.lt.a) then
            tjm= exp(ran1*log(a-cxm)+
     *          (1.d0-ran1)*log(a-cxp))
          else
            tjm= -exp(ran1*log(-a+cxm)+
     *          (1.d0-ran1)*log(-a+cxp))
          endif
        endif
      endif
      return
      end

      function rantjm(a,cn,cxm,cxp,k,tjm)
c----------------------------------------------------------------c
c                                                                c
c     given tjm produced according to the distribution           c
c                                                                c
c                    (1/tjm)^cn     where:                       c
c                                                                c
c                    tjm= a+k*ct                                 c  
c                    cxm < ct < cxp                              c
c                    k= +1 or -1                                 c  
c                                                                c  
c     rantjm is the corresp. uniform random number in [0,1]      c 
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 a,cn,cxm,cxp,rantjm,ce,tjm
      integer k 
      ce= 1.d0-cn
      if (abs(ce).gt.1.d-8) then
        if (k.eq.1) then
          rantjm= (tjm**ce-(a+cxm)**ce)/((a+cxp)**ce-(a+cxm)**ce)
        else
          rantjm= (tjm**ce-(a-cxp)**ce)/((a-cxm)**ce-(a-cxp)**ce)
        endif
      else
        if (k.eq.1) then
          if(cxp.lt.-a) then
            rantjm= (log(-tjm)-log(-a-cxm))/
     *              (log(-a-cxp)-log(-a-cxm))
          else
            rantjm= (log(tjm)-log(a+cxm))/
     *              (log(a+cxp)-log(a+cxm))
          endif
        else
          if(cxp.lt.a) then
            rantjm= (log(tjm)-log(a-cxp))/
     *              (log(a-cxm)-log(a-cxp))
          else
            rantjm= (log(-tjm)-log(-a+cxp))/
     *              (log(-a+cxm)-log(-a+cxp))
          endif
        endif
      endif
      return
      end

      subroutine boost(lflag,sq,q,ph,p)
c----------------------------------------------------------------c
c                                        _                       c
c   Boost of a 4-vector ( relative speed q/q(0) ):               c
c                                                                c
c   ph is the 4-vector in the rest frame of q                    c
c   p is the corresponding 4-vector in the lab frame             c
c                                                                c  
c              INPUT                               OUTPUT        c
c                                                                c
c  lflag= 0:  sq, q, ph                            p             c
c                                                                c
c  lflag= 1:  sq, q, p                             ph            c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 sq,c1,rsq
      integer lflag,i   
      real*8 q(0:3),ph(0:3),p(0:3)
      rsq= sqrt(sq)
      if (lflag.eq.0) then
        p(0)= (q(0)*ph(0)+q(1)*ph(1)+q(2)*ph(2)+q(3)*ph(3))/rsq
        c1= (ph(0)+p(0))/(rsq+q(0))
        do i= 1,3
          p(i)= ph(i)+q(i)*c1
        enddo
      else
        ph(0)= (q(0)*p(0)-q(1)*p(1)-q(2)*p(2)-q(3)*p(3))/rsq
        c1= (p(0)+ph(0))/(rsq+q(0))
        do i= 1,3
          ph(i)= p(i)-q(i)*c1
        enddo
      endif
      return
      end

      subroutine qvec(q0,qm,st,ct,sp,cp,q)
c----------------------------------------------------------------c
c                                                                c
c   A 4-vector is built                                          c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 q0,qm,st,ct,sp,cp
      real*8 q(0:3)
      q(0)= q0
      q(1)= qm*st*sp        
      q(2)= qm*st*cp       
      q(3)= qm*ct
      return
      end 

      function sqlam(s,s1,s2)
      implicit none
      real*8 s,s1,s2,x1,x2,arg1,sqlam,arg
      x1= s1/s
      x2= s2/s
      arg= 1.d0-x1-x2
      arg1= arg*arg-4.d0*x1*x2
      if (arg1.gt.0.d0) then
        sqlam= sqrt(arg1)
      else
        sqlam= 1.d-30
      endif
      return
      end

      function hj(a,cn,cxm,cxp)
c----------------------------------------------------------------c
c                                                                c
c   Normalization factor for tj                                  c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 a,cn,cxm,cxp,ce,hj
      ce= 1.d0-cn
      if (abs(ce).gt.1.d-8) then
        hj= ((a+cxp)**ce-(a+cxm)**ce)/ce
      else
        hj= log((a+cxp)/(a+cxm))
      endif
      return
      end

      subroutine onedimbin(lflag,n1,w,n,igridw,peropt)
c-
c-    Self optimizing binning
c- 
c-    n labels the variable for which one wishes to
c-    have a self optimized binning
c-
c-    if lflag= 1:
c-        n1 is generated in [1:nx1] and w= djtot1 (the inverse jacobian)
c-        is an output.
c-
c-    if lflag= 2:
c-        w is the input weight for the optimization. The optimization
c-        is performed only if n1= 1.
c-
c-    if lflag= 3:
c-        given n1 and n w= djtot1 (the inverse jacobian) is computed
c-
c-    if lflag= 0:
c-        the new set of a-priori weights is computed and
c-        written to unit= igridw if igridw.ne.0.
c-
      implicit none
      real*8 peropt,w,beta,rcho,wt1,ni
      real*8 al1,bet1,v1,aptot1,alsum1,alsum2
      integer lflag,n1,n,ncmax,maxn,modopt,igridw
      integer nx1,nx2,nct1,nct,init,i,m,j,k,lrun
      parameter (beta= 0.5d0, ncmax= 5000, maxn= 500)
c-
c-    if modopt= 0 [1] the grid optimization is done according to the
c-    variance [cross section]
c-
      parameter (modopt= 1)
c-
c-    peropt is the percentage [0:1] of allowed optimization
c-
c-    WARNING: ncmax must be .le. 6000
c-
      real*8 s1(ncmax)
      integer iopt(maxn),jad1(maxn)
      integer ntmp  
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      common/rdandwrt/al1(ncmax,maxn),nct1(maxn)
      common/ausil/init(ncmax),bet1(0:ncmax,maxn)
      common/book/v1(ncmax),ni(maxn)
comment
      save 
*
c-    self-optimization of the a-priori weights   
*                
      if (lflag.eq.0) then
*
c-      only if there was optimization the new set is computed
*
        if (iopt(n).eq.1) then
           do i= nct(n)+1,nct(n)+nx1(n)
             s1(i)= v1(i)/ni(n)
           enddo
           aptot1= 0.d0                                             
           alsum1= 0.d0                                             
           alsum2= 0.d0                                             
           do m= 1,nx1(n)                                             
             if (modopt.eq.0) then
               aptot1= aptot1+al1(m,n)*s1(nct(n)+m)**beta                  
             else 
               if (s1(nct(n)+m).lt.1.d-10) then
                  aptot1= aptot1+1.d-10 
               else
                  aptot1= aptot1+s1(nct(n)+m)                 
               endif
             endif
           enddo                                                    
           do m= 1,nx1(n)                                             
             if (modopt.eq.0) then
               al1(m,n)= al1(m,n)*(s1(nct(n)+m)**beta)/aptot1       
             else
               if (s1(nct(n)+m).lt.1.d-10) then
                 al1(m,n)=  1.d-10/aptot1
               else
                 al1(m,n)=  s1(nct(n)+m)/aptot1
               endif 
             endif
           enddo                                                    
*
c-         rescaling  according to the chosen value of peropt
*
           do m= 1,nx1(n)                                             
             al1(m,n)= al1(m,n)*peropt
             al1(m,n)= al1(m,n)+(1.d0-peropt)/dfloat(nx1(n))
           enddo                                                    
        endif
*
        if (igridw.ne.0) then
*
c-        if igridw= 0 the a-priori weights are not written 
c-        if igridw.ne.0 they are written in unit= igridw
*
          write (igridw,20) n,(m,al1(m,n),m= 1,nx1(n))
        endif
*
c-      computes the new bet1(j) 
*
        do j= 1,nx1(n)                                             
          bet1(j,n)= 0                                             
          do k= 1,j                                              
            bet1(j,n)= bet1(j,n)+al1(k,n)                              
          enddo
        enddo
        return 
*
c-    generation of n1 and computation of w= djtot1
*
      elseif (lflag.eq.1) then
         call rans(rcho)                                            
         ntmp = nx1(n) 
         do lrun= 1,ntmp                                            
            if (rcho.le.bet1(lrun,n)) then
               jad1(n)= lrun
               goto 10                        
            endif
         enddo                                                      
 10      n1= jad1(n)
         w = ntmp*al1(n1,n)                             
         return
*
c-    bookkeeping for the self-optimization, w is an input 
*
      elseif(lflag.eq.2) then  
         if (n1.eq.1) then
            iopt(n)= 1
            if (modopt.eq.0) then
              wt1= w*w/al1(jad1(n),n)
            else
              wt1= w
            endif
            ntmp       = nct(n)+jad1(n)
            v1(ntmp)   = v1(ntmp)+wt1
            ni(n)      = ni(n)+1
         endif
*
c-    given n1 and n w= djtot1 (the inverse jacobian) is computed
*
      elseif(lflag.eq.3) then  
         w = nx1(n)*al1(n1,n) 
      endif
 20   format(/,'  a-priori weights for variable n.',i4,
     +       //,(i3,1x,d20.9))
 21   format(///,(4x,d20.9))
      return
      end
*
      subroutine grid1W(nwrt,n)
      implicit none
      real*8 al1
      integer nwrt,n,ncmax,maxn
      integer nct,nx1,nx2,nct1,m     	
      parameter (ncmax= 5000, maxn= 500)
c-
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      common/rdandwrt/al1(ncmax,maxn),nct1(maxn)
      save 

      write (nwrt,20) n,(m,al1(m,n),m= 1,nx1(n))
 20   format(/,'  a-priori weights for variable n.',i4,
     +       //,(i3,1x,d20.9))
      end

      subroutine grid1R(nrd,n)
      implicit none 
      real*8 al1,bet1
      integer nct,nx1,nx2,nct1,init,j,k,nrd,i,n
      integer ncmax,maxn
      parameter (ncmax= 5000, maxn= 500)
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      common/rdandwrt/al1(ncmax,maxn),nct1(maxn)
      common/ausil/init(ncmax),bet1(0:ncmax,maxn)
      save 

c-    if nrd.ne.0
c-    the a-priori weights are read from unit= nrd
      if (nrd.ne.0) read (nrd,21) (al1(i,n),i= 1,nx1(n))
         
      do j= 1,nx1(n)                                         
        bet1(j,n)= 0                                         
        do k= 1,j                                          
          bet1(j,n)= bet1(j,n)+al1(k,n)                          
        enddo                                  
      enddo
cmlm
      if(abs(bet1(nx1(n),n)-1d0).lt.1d-6) then
         bet1(nx1(n),n)=1d0
      else
         write(6,*) 'grid entries for variable',n,' don''t sum to 1:'
         write(6,*) 'beta1=',bet1(nx1(n),n)
         write(6,*) 'error initialising grid, stop'
         stop
      endif

 21   format(///,(4x,d20.9))
      end

c-------------------------------------------------------------------
      subroutine dumpgrid(igridw)
c-------------------------------------------------------------------
      implicit none
c input variables
      integer igridw
      double precision peropt,totalmax,tmpfac
c commons
      integer nvar,nv
      common/psopt/nvar,nv
      integer maxn
      parameter (maxn= 500)  
      integer mask,mmask 
      common/psopt1/mask(maxn),mmask(maxn),peropt(maxn)
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      common/mxfact/totalmax,tmpfac
c locals
      integer ndummy,n
      real*8 dummy
c
      do n= 1,nvar+1
        call onedimbin(0,ndummy,dummy,n,igridw,peropt(n))
      enddo
      write (igridw,20) totalmax,tmpfac
 20   format(/,'  totalmax and tmpfac',//,2(d20.9))
      end

c-------------------------------------------------------------------
      subroutine readgrid(nrd)
c-------------------------------------------------------------------
      implicit none
c inputs
      integer nrd
c common
      integer nvar,nv
      common/psopt/nvar,nv
      real*8 totalmax,tmpfac
      common/mxfact/totalmax,tmpfac
c locals
      integer  n
c
c      if (nrd.ne.niopar) open (unit= nrd,status='old')
      do n= 1,nvar+1
        call grid1r(nrd,n)
      enddo
      read(nrd,20) totalmax,tmpfac
 20   format(///,2(d20.9))
c      if (nrd.ne.niopar) close(nrd)

      end
*
      subroutine spk(lflag,a,cn,pt2,pt3,phi2,phi3,xmr2,xmr3,eta2,
     .               etc0,etc1,eta3,wjac,ran0,lw)
c----------------------------------------------------------------c
c                                                                c
c     It returns eta3, generated within etc0 and etc1,           c
c     such that the invariant mass s23 has the distribution      c
c                                                                c
c     1/(a+s_23)^cn                                              c
c                                                                c
c     If lflag.eq.0:                                             c
c                                                                c
c     Input: a,cn,pt2,pt3,phi2,phi3,xmr2,xmr3,eta2,              c
c            etc0,etc1,ran0                                      c
c                                                                c
c     Ouput: eta3,wjac,lw                                        c
c                                                                c
c     If lflag.eq.1:                                             c
c                                                                c
c     Input: a,cn,pt2,pt3,phi2,phi3,xmr2,xmr3,eta2,eta3          c
c            etc0,etc1                                           c
c                                                                c
c     Ouput: ran0,wjac,lw                                        c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 pt(2:3),eta(2:3),etc0,etc1,xmr(2:3)
      real*8 al2,al3,xp,xm,x2p,x2m,zc0,zc1,x3p0,x3p1,x3mi,sm
      real*8 x3m,x3p,pt23s,yp,wjac,sp,spp,smm,x,al,fla,cn,wj,s
      real*8 rnn,ran0,wf,sig,sc,a
      real*8 pt2,pt3,phi2,phi3,xmr2,xmr3,eta2,eta3
      integer lw,lflag
*
c-    Set local variables 
*      
      pt(2)  = pt2
      pt(3)  = pt3
      xmr(2) = xmr2
      xmr(3) = xmr3
      eta(2) = eta2
      if (lflag.eq.1) eta(3)= eta3
      pt23s  = pt(2)**2+pt(3)**2+2.d0*pt(2)*pt(3)*cos(phi2-phi3)
*
      al2= pt(2)**2+xmr(2)**2
      al3= pt(3)**2+xmr(3)**2
      if (eta(2).ge.0.d0) then
        xp = sqrt(pt(2)**2*(cosh(eta(2)))**2+xmr(2)**2)
     .           +pt(2)*sinh(eta(2))
        x2p= xp
        x2m= al2/x2p
      else
        xm = sqrt(pt(2)**2*(cosh(eta(2)))**2+xmr(2)**2)
     .           -pt(2)*sinh(eta(2))
        x2m= xm
        x2p= al2/x2m
      endif
      if (lflag.eq.1) then
        if (eta(3).ge.0.d0) then
          xp = sqrt(pt(3)**2*(cosh(eta(3)))**2+xmr(3)**2)
     .             +pt(3)*sinh(eta(3))
          x3p= xp
          x3m= al3/x3p
        else
          xm = sqrt(pt(3)**2*(cosh(eta(3)))**2+xmr(3)**2)
     .             -pt(3)*sinh(eta(3))
          x3m= xm
          x3p= al3/x3m
        endif
        s=-pt23s+al2+al3+x3p*x2m+x2p*x3m 
      endif
      zc0 = 2.d0*pt(3)*sinh(etc0)
      zc1 = 2.d0*pt(3)*sinh(etc1)
      if (zc0.ge.0.d0) then
        x3p0= 0.5d0*(zc0+sqrt(zc0*zc0+4.d0*al3))
      else
        x3p0=-al3/(0.5d0*(zc0-sqrt(zc0*zc0+4.d0*al3)))
      endif
      if (zc1.ge.0.d0) then
        x3p1= 0.5d0*(zc1+sqrt(zc1*zc1+4.d0*al3))
      else
        x3p1=-al3/(0.5d0*(zc1-sqrt(zc1*zc1+4.d0*al3)))
      endif
      x3mi= sqrt(al2*al3)/x2m
      sm  = -pt23s+al2+al3+x2m*x3p0+al2*al3/x2m/x3p0
      sp  = -pt23s+al2+al3+x2m*x3p1+al2*al3/x2m/x3p1
      if     (x3p0.ge.x3mi) then
        smm= sm
        spp= sp
        wf = 1.d0
        if (lflag.eq.0) then
          sig= 1.d0
          rnn= ran0
        else
          sig= 0.d0
        endif
      elseif (x3p1.le.x3mi) then
        smm= sp
        spp= sm
        wf = 1.d0
        if (lflag.eq.0) then
          sig=-1.d0
          rnn= ran0
        else
          sig= 0.d0
        endif
      elseif (x3p0.lt.x3mi.and.x3p1.gt.x3mi) then
        sc= -pt23s+al2+al3+2.d0*sqrt(al2*al3)
        wf = 2.d0
        if (lflag.eq.0) then
          if (ran0.le.0.5d0) then
            smm= sc
            spp= sm
            sig=-1.d0
            rnn= ran0*2.d0
          else
            smm= sc
            spp= sp
            sig= 1.d0
            rnn= ran0*2.d0-1.d0
          endif
        else
          if (x3p.lt.x3mi) then
            smm= sc
            spp= sm
            sig= 0.d0
          else
            smm= sc
            spp= sp
            sig= 1.d0
          endif
        endif
      else
        print*,'error'
        stop
      endif
      if (smm.le.0.d0) goto 100
      if (spp.le.0.d0) goto 100
      if (smm.ge.spp)  goto 100
      call ppeaka(lflag,a,cn,smm,spp,s,wj,rnn,lw)
      if (lw.eq.1) goto 100
      if (lflag.eq.1) then
        if (pt(3).lt.1.d-20) goto 100
        fla = dabs(x3p*x2m-x3m*x2p)
        if (fla.lt.1.d-20)   goto 100
        wjac= (x3p+x3m)/(2.d0*pt(3)*cosh(eta(3)))/
     .       fla*wf
        wjac= wjac*wj
        ran0= (rnn+sig)/wf 
        return
      endif
      al  = pt23s+s
      fla = al**2+al2**2+al3**2-2.d0*(al*al2+al*al3+al2*al3)
      fla = sqrt(dabs(fla))
      yp  = 0.5d0*((al-al2-al3)+fla)
      if (dabs(sig-1.d0).lt.1.d-5) x3p= yp/x2m
      if (dabs(sig+1.d0).lt.1.d-5) x3p= al2*al3/yp/x2m
      x3m = al3/x3p
      x   = 0.5d0/pt(3)*(x3p-al3/x3p)
      if (x.gt.0.d0) then
        eta(3)= log(sqrt(x*x+1.d0)+x)
      else
        eta(3)=-log(sqrt(x*x+1.d0)-x)
      endif
      wjac=dabs(x3p*x2m-x3m*x2p)
      if(wjac.lt.1.d-20) goto 100
      wjac= (x3p+x3m)/(2.d0*pt(3)*cosh(eta(3)))/
     .      wjac*wf
      wjac= wjac*wj
      eta3= eta(3)
*
      return
 100  wjac= 0.d0
      lw= 1
      end

*
      subroutine etagen(lflag,ph1,ph2,eta1,eta2m,eta2p,eta2,wjac,rnd,lw)
c---------------------------------------------------------------c
c
c     Given 2 massless 4-momenta p1 and p2, eta2 is generated 
c     in such a way that  s12 has the distribution 1/s12:
c
c             1
c        --------------
c        cosh(eta2-a)-b
c
c     When p1 is along z use a=b=0 --> eta1=0, ph1=0, ph2=pi/2.
c
c     When lflag.eq.0
c
c       input:    ph1,ph2,eta1,eta2m,eta2p,rnd
c       
c       output:   eta2,wjac
c
c     When lflag.eq.1
c
c       input:    ph1,ph2,eta1,eta2,eta2m,eta2p
c       
c       output:   wjac,rnd
c 
c---------------------------------------------------------------c

      implicit none
      real*8 ph1,ph2,eta1,eta2,eta2m,eta2p,rnd,wjac
      real*8 a,b,c,xm,xp,bb,aa,den
      integer lflag,lw
      lw= 0
*
      if (lflag.eq.1) then
        if (eta2.lt.eta2m) goto 100
        if (eta2.gt.eta2p) goto 100
      endif
*
c-    local variables:
*      
      a = eta1
      b = cos(ph1-ph2)
*
c-    Protection:
*
      if (b.ge.0.d0) then
        b= min(b,0.99d0)
      else
        b= max(b,-0.99d0)
      endif
*
      c = sqrt(dabs(1.d0-b*b))
      xm= eta2m
      xp= eta2p
*
      bb = atan((exp(xm-a)-b)/c)
      den= atan((exp(xp-a)-b)/c)-bb
      if (dabs(den).lt.1.d-12) goto 100
      aa = c/2.d0/den
      if (lflag.eq.0) then
        den= b+c*tan(0.5d0*c*rnd/aa+bb)
        if (den.le.0.d0) goto 100
        eta2= a+log(den)
      endif
      wjac= (cosh(eta2-a)-b)/aa
      if (lflag.eq.1) rnd= 2.d0*aa/c*(atan((exp(eta2-a)-b)/c)-bb)
*
      return
 100  wjac= 0.d0
      lw= 1
      return
      end

      SUBROUTINE RAMBO(LFLAG,N,ET,XM,P,DJ)
C------------------------------------------------------
C
C                       RAMBO
C
C    RA(NDOM)  M(OMENTA)  B(EAUTIFULLY)  O(RGANIZED)
C
C    A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR
C    AUTHORS:  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING
C    THIS IS VERSION 1.0 -  WRITTEN BY R. KLEISS
C    (MODIFIED BY R. PITTAU)
C
C                INPUT                 OUTPUT
C
C    LFLAG= 0:   N, ET, XM             P, DJ
C    LFLAG= 1:   N, ET, XM, P          DJ
C
C    N  = NUMBER OF PARTICLES (>1, IN THIS VERSION <101)
C    ET = TOTAL CENTRE-OF-MASS ENERGY
C    XM = PARTICLE MASSES ( DIM=100 )
C    P  = PARTICLE MOMENTA ( DIM=(4,100) )
C    DJ = 1/(WEIGHT OF THE EVENT)
C
C------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XM(100),P(4,100),Q(4,100),Z(100),R(4),
     .   B(3),P2(100),XM2(100),E(100),V(100),IWARN(5)
      SAVE ACC,ITMAX,IBEGIN,IWARN,Z,TWOPI,PO2LOG
      DATA ACC/1.D-14/,ITMAX/10/,IBEGIN/0/,IWARN/5*0/
C
C INITIALIZATION STEP: FACTORIALS FOR THE PHASE SPACE WEIGHT
      IF(IBEGIN.NE.0) GOTO 103
      IBEGIN=1
      TWOPI=8.*DATAN(1.D0)
      PO2LOG=LOG(TWOPI/4.)
      Z(2)=PO2LOG
      DO 101 K=3,100
  101 Z(K)=Z(K-1)+PO2LOG-2.*LOG(DFLOAT(K-2))
      DO 102 K=3,100
  102 Z(K)=(Z(K)-LOG(DFLOAT(K-1)))
C
C CHECK ON THE NUMBER OF PARTICLES
  103 IF(N.GT.1.AND.N.LT.101) GOTO 104
      PRINT 1001,N
      STOP
C
C CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES
  104 XMT=0.
      NM=0
      DO 105 I=1,N
      IF(XM(I).NE.0.D0) NM=NM+1
  105 XMT=XMT+ABS(XM(I))
      IF(XMT.LE.ET) GOTO 201
      PRINT 1002,XMT,ET
      STOP

  201 CONTINUE 
      if (lflag.eq.1) then
        w0= exp((2.*N-4.)*LOG(ET)+Z(N))
        do j= 1,N
          v(j)= sqrt(p(1,j)**2+p(2,j)**2+p(3,j)**2)
        enddo

        a1= 0.d0
        a3= 0.d0
        a2= 1.d0
        do j= 1,N
          a1= a1+v(j)/ET
          a2= a2*v(j)/p(4,j)
          a3= a3+v(j)*v(j)/p(4,j)/ET
        enddo
        wm= a1**(2*N-3)*a2/a3
        dj= 1.d0/w0/wm
        return
      endif
C
C THE PARAMETER VALUES ARE NOW ACCEPTED
C
C GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE

      DO 202 I=1,N
      call rans(RAN1)
      call rans(RAN2)
      call rans(RAN3)
      call rans(RAN4)
      C=2.*RAN1-1.
      S=SQRT(1.-C*C)
      F=TWOPI*RAN2
      Q(4,I)=-LOG(RAN3*RAN4)
      Q(3,I)=Q(4,I)*C
      Q(2,I)=Q(4,I)*S*COS(F)
  202 Q(1,I)=Q(4,I)*S*SIN(F)
C
C CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION
      DO 203 I=1,4
  203 R(I)=0.
      DO 204 I=1,N
      DO 204 K=1,4
  204 R(K)=R(K)+Q(K,I)
      RMAS=SQRT(R(4)**2-R(3)**2-R(2)**2-R(1)**2)
      DO 205 K=1,3
  205 B(K)=-R(K)/RMAS
      G=R(4)/RMAS
      A=1./(1.+G)
      X=ET/RMAS
C
C TRANSFORM THE Q'S CONFORMALLY INTO THE P'S
      DO 207 I=1,N
      BQ=B(1)*Q(1,I)+B(2)*Q(2,I)+B(3)*Q(3,I)
      DO 206 K=1,3
  206 P(K,I)=X*(Q(K,I)+B(K)*(Q(4,I)+A*BQ))
  207 P(4,I)=X*(G*Q(4,I)+BQ)
C
C CALCULATE WEIGHT AND POSSIBLE WARNINGS
      WT=PO2LOG
      IF(N.NE.2) WT=(2.*N-4.)*LOG(ET)+Z(N)
      IF(WT.GE.-180.D0) GOTO 208
      IF(IWARN(1).LE.5) PRINT 1004,WT
      IWARN(1)=IWARN(1)+1
  208 IF(WT.LE. 174.D0) GOTO 209
      IF(IWARN(2).LE.5) PRINT 1005,WT
      IWARN(2)=IWARN(2)+1
C
C RETURN FOR WEIGHTED MASSLESS MOMENTA
  209 IF(NM.NE.0) GOTO 210
      WT=EXP(WT)
      DJ= 1.d0/WT
      RETURN
C
C MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X
  210 XMAX=SQRT(1.-(XMT/ET)**2)
      DO 301 I=1,N
      XM2(I)=XM(I)**2
  301 P2(I)=P(4,I)**2
      ITER=0
      X=XMAX
      ACCU=ET*ACC
  302 F0=-ET
      G0=0.
      X2=X*X
      DO 303 I=1,N
      E(I)=SQRT(XM2(I)+X2*P2(I))
      F0=F0+E(I)
  303 G0=G0+P2(I)/E(I)
      IF(ABS(F0).LE.ACCU) GOTO 305
      ITER=ITER+1
      IF(ITER.LE.ITMAX) GOTO 304
      PRINT 1006,ITMAX
      GOTO 305
  304 X=X-F0/(X*G0)
      GOTO 302
  305 DO 307 I=1,N
      V(I)=X*P(4,I)
      DO 306 K=1,3
  306 P(K,I)=X*P(K,I)
  307 P(4,I)=E(I)
C
C CALCULATE THE MASS-EFFECT WEIGHT FACTOR
      WT2=1.
      WT3=0.
      DO 308 I=1,N
      WT2=WT2*V(I)/E(I)
  308 WT3=WT3+V(I)**2/E(I)
      WTM=(2.*N-3.)*LOG(X)+LOG(WT2/WT3*ET)
C
C RETURN FOR  WEIGHTED MASSIVE MOMENTA
      WT=WT+WTM
      IF(WT.GE.-180.D0) GOTO 309
      IF(IWARN(3).LE.5) PRINT 1004,WT
      IWARN(3)=IWARN(3)+1
  309 IF(WT.LE. 174.D0) GOTO 310
      IF(IWARN(4).LE.5) PRINT 1005,WT
      IWARN(4)=IWARN(4)+1
  310 WT=EXP(WT)
      DJ= 1.d0/WT
      RETURN
C
 1001 FORMAT(' RAMBO FAILS: # OF PARTICLES =',I5,' IS NOT ALLOWED')
 1002 FORMAT(' RAMBO FAILS: TOTAL MASS =',D15.6,' IS NOT',
     . ' SMALLER THAN TOTAL ENERGY =',D15.6)
 1004 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY UNDERFLOW')
 1005 FORMAT(' RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY  OVERFLOW')
 1006 FORMAT(' RAMBO WARNS:',I3,' ITERATIONS DID NOT GIVE THE',
     . ' DESIRED ACCURACY =',D15.6)
      END

c RP: I have added subroutine flat
      subroutine flat(x,x0,x1,w,ran0)
      implicit none
      real*8 x,x0,x1,w,ran0
*
      x= x1*ran0+(1.d0-ran0)*x0
      w= (x1-x0)
      return
      end
*
      subroutine fflat(lflag,x,x0,x1,w,ran0,lw)
      implicit none
      real*8 x,x0,x1,w,ran0
      integer lflag,lw
      w= (x1-x0)
      if (lflag.eq.1) then
        if (x.lt.x0.or.x.gt.x1) goto 100
        ran0= (x-x0)/(x1-x0)
        return
      endif
*
      x= x1*ran0+(1.d0-ran0)*x0
      return
 100  w= 0.d0
      lw= 1
      return
      end
*

c RP: I have added function asinh
      function asinh(x)
      implicit none
      real*8 x,asinh
      if (x.ge.0.d0) then
        asinh= log(sqrt(x*x+1.d0)+x)
      else
        asinh=-log(sqrt(x*x+1.d0)-x)
      endif
      return
      end     
*
c RP: I have added subroutine peaka
      subroutine peaka(a,cn,cxm,cxp,x,wj,ran1)
c----------------------------------------------------------------c
c                                                                c
c   Peaked distribution:                                         c
c   The variable x has the distribution                          c
c                                                                c
c                   1/(a+x)^cn                                   c
c                   cxm < x < cxp                                c
c                                                                c
c              INPUT                      OUTPUT                 c 
c                                                                c
c              a, cn, cxm, cxp            x, wj                  c
c                                                                c
c----------------------------------------------------------------c
*
      implicit none
      real*8 a,cn,cxm,cxp,x,wj,ran1,tjm,apx,hj
      apx= tjm(a,cn,cxm,cxp,1,ran1)
      x= -a+apx
      if (abs(1.d0-cn).gt.1.d-8) then
         wj= (apx)**cn*hj(a,cn,cxm,cxp)
      else
         wj= (apx*hj(a,cn,cxm,cxp))
      endif
      return
      end
* 
      subroutine triangle(x0,xm,xp,x,wj,ran0)
c----------------------------------------------------------------c
c                                                                c
c     Triangular distribution                                    c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      real*8 x0,xm,xp,x,wj,ran0,root,ran1,al1,alsum
      if (xm.lt.0.d0.and.xp.lt.0.d0) then
        root= (x0+xm)**2+ran0*(xp-xm)*(xp+xm+2.d0*x0)
        x   =-x0+sqrt(root)
        wj  = 0.5d0*(xp-xm)/(x+x0)*(xp+xm+2.d0*x0)
      elseif(xm.gt.0.d0.and.xp.gt.0.d0) then
        root= (-x0+xm)**2+ran0*(xp-xm)*(xp+xm-2.d0*x0)
        x   = x0-sqrt(root)
        wj  = 0.5d0*(xp-xm)/(x-x0)*(xp+xm-2.d0*x0)
      elseif(xm.lt.0.d0.and.xp.gt.0.d0) then
        alsum=-xm+xp
        al1  = -xm/alsum
        if (ran0.le.al1) then
          ran1= ran0/al1
          root= (x0+xm)**2-ran1*xm*(xm+2.d0*x0)
          x   =-x0+sqrt(root)
          wj  = 0.5d0*(-xm)/(x+x0)*(xm+2.d0*x0)/al1
        else
          ran1= (ran0-al1)/(1.d0-al1)
          root= x0**2+ran1*xp*(xp-2.d0*x0)
          x   = x0-sqrt(root)
          wj  = 0.5d0*xp/(x-x0)*(xp-2.d0*x0)/(1.d0-al1)
        endif
      else
        print*,'ERROR in subroutine triangle'
        stop 
      endif
      return
      end
*
      subroutine ppeaka(lflag,a,cn,cxm,cxp,x,wj,ran1,lw)
c----------------------------------------------------------------c
c                                                                c
c   Peaked distribution:                                         c
c   The variable x has the distribution                          c
c                                                                c
c                   1/(a+x)^cn                                   c
c                   cxm < x < cxp                                c
c                                                                c
c              INPUT                      OUTPUT                 c 
c                                                                c
c  lflag= 0:   a, cn, cxm, cxp, ran1      x, wj                  c
c                                                                c
c  lflag= 1:   a, cn, cxm, cxp, x         wj, ran1               c
c                                                                c
c----------------------------------------------------------------c
*
      implicit none
      real*8 a,cn,cxm,cxp,x,wj,ran1,tjm,apx,hj,rantjm
      integer lflag,lw
      if (lflag.eq.0) then
        apx= tjm(a,cn,cxm,cxp,1,ran1)
        x= -a+apx
      else
        apx= a+x
        if (x.lt.cxm.or.x.gt.cxp) goto 100
        ran1= rantjm(a,cn,cxm,cxp,1,apx)
      endif
      if (abs(1.d0-cn).gt.1.d-8) then
         wj= (apx)**cn*hj(a,cn,cxm,cxp)
      else
         wj= (apx*hj(a,cn,cxm,cxp))
      endif
      return
 100  wj= 0.d0
      lw= 1
      return
      end

      subroutine ttriangle(lflag,x0m,x0p,xm,xp,x,wj,ran0,lw)
c----------------------------------------------------------------c
c                                                                c
c     Triangular distribution                                    c
c                                                                c
c----------------------------------------------------------------c

      implicit none
      integer lflag,lw
      real*8 x0,x0m,x0p,xm,xp,x,wj,ran0,root,ran1,al1,alsum,prot
      real*8 a,b
      data prot/1.d-20/
*
      if (lflag.eq.1) then
        if (x.lt.xm.or.x.gt.xp) goto 100
      endif
*
      if (xm.lt.0.d0.and.xp.lt.0.d0) then
        x0= dabs(x0m)
        if (lflag.eq.0) then
          root= (x0+xm)**2+ran0*(xp-xm)*(xp+xm+2.d0*x0)
          x   =-x0+sqrt(dabs(root))
        else
          a   =  0.5d0*(xp-xm)*(xp+xm+2.d0*x0)
          b   =  0.5d0*xm*(xm+2.d0*x0)
          ran0= (x*(0.5d0*x+x0)-b)/a
        endif
        if (dabs(x+x0).lt.prot) goto 100
        wj  = 0.5d0*(xp-xm)/(x+x0)*(xp+xm+2.d0*x0)
      elseif(xm.gt.0.d0.and.xp.gt.0.d0) then
        x0= dabs(x0p)
        if (lflag.eq.0) then
          root= (-x0+xm)**2+ran0*(xp-xm)*(xp+xm-2.d0*x0)
          x   = x0-sqrt(dabs(root))
        else
          a   =  0.5d0*(xp-xm)*(xp+xm-2.d0*x0)
          b   =  0.5d0*xm*(xm-2.d0*x0)
          ran0= (x*(0.5d0*x-x0)-b)/a
        endif
        if (dabs(x-x0).lt.prot) goto 100
        wj  = 0.5d0*(xp-xm)/(x-x0)*(xp+xm-2.d0*x0)
      elseif(xm.lt.0.d0.and.xp.gt.0.d0) then
        alsum=-xm+xp
        if (dabs(alsum).lt.prot) goto 100
        al1  = -xm/alsum
        if (lflag.eq.0) then
          if (ran0.le.al1) then
            x0= dabs(x0m)
            if (dabs(al1).lt.prot) goto 100
            ran1= ran0/al1
            root= (x0+xm)**2-ran1*xm*(xm+2.d0*x0)
            x   =-x0+sqrt(dabs(root))
            if (dabs(x+x0).lt.prot) goto 100
            wj  = 0.5d0*(-xm)/(x+x0)*(xm+2.d0*x0)/al1
          else
            x0= dabs(x0p)
            if (dabs(1.d0-al1).lt.prot)    goto 100
            ran1= (ran0-al1)/(1.d0-al1)
            root= x0**2+ran1*xp*(xp-2.d0*x0)
            x   = x0-sqrt(dabs(root))
            if (dabs(x-x0).lt.prot) goto 100
            wj  = 0.5d0*xp/(x-x0)*(xp-2.d0*x0)/(1.d0-al1)
          endif
        elseif (lflag.eq.1) then
          if (x.lt.0.d0) then
            x0= dabs(x0m)
            if (dabs(al1).lt.prot)  goto 100
            if (dabs(x+x0).lt.prot) goto 100
            wj  = 0.5d0*(-xm)/(x+x0)*(xm+2.d0*x0)/al1
            a   =  0.5d0*(-xm)*(xm+2.d0*x0)
            b   =  0.5d0*xm*(xm+2.d0*x0)
            ran1= (x*(0.5d0*x+x0)-b)/a
            ran0= ran1*al1
          else
            x0= dabs(x0p)
            if (dabs(1.d0-al1).lt.prot)  goto 100
            if (dabs(x-x0).lt.prot)      goto 100
            wj  = 0.5d0*xp/(x-x0)*(xp-2.d0*x0)/(1.d0-al1)
            a   =  0.5d0*(xp)*(xp-2.d0*x0)
            b   =  0.d0
            ran1= (x*(0.5d0*x-x0)-b)/a
            ran0= al1+ran1*(1.d0-al1)
          endif 
        else
          goto 101
        endif
      else
        goto 100
      endif
*
      return
 100  wj= 0.d0
      lw= 1
      return
 101  print*,'ERROR IN SUBROUTINE TTRIANGLE'
      stop
      end
*
      function aluhiw(mh,mw,mz,mt)
      implicit none
      double precision aluhiw,mh,mw,mz,mt,tmp
      double precision gtow,gtoz,gtot
      double precision b,x
      double precision pi,gfermi
      parameter (pi=3.14159265d0,gfermi=1.16639d-5)
c
      tmp=0d0
c     H-> WW
      if(mh.gt.2d0*mw) then
         x=mw**2/mh**2
         b=sqrt(1d0-4d0*x)
         gtow=2d0*sqrt(2d0)*gfermi/32d0/pi * mh**3*(1d0-4d0*x+12d0*x**2)
     $        *b
         tmp=tmp+gtow
      endif
c     H-> ZZ
      if(mh.gt.2d0*mz) then
         x=mz**2/mh**2
         b=sqrt(1d0-4d0*x)
         gtoz=sqrt(2d0)*gfermi/32d0/pi * mh**3*(1d0-4d0*x+12d0*x**2)*b
         tmp=tmp+gtoz
      endif
c     H-> t tbar
      if(mh.gt.2d0*mt) then
         x=mt**2/mh**2
         b=sqrt(1d0-4d0*x)
         gtot= gfermi /4.d0/pi/dsqrt(2.d0)*mh*3.d0*mt**2 * b
         tmp=tmp+gtot
      endif
      aluhiw=tmp
      end


      subroutine rans(ran)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     Random number generator                                    c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 ran,rangen
      ran= rangen(1)
      end

      subroutine randa(rn)
      implicit none
      real *8 rn,rangen
      rn=rangen(1)
      end

*-- Author :    F. James, modified by Mike Seymour
C-----------------------------------------------------------------------
      FUNCTION RANGEN(I)
C-----------------------------------------------------------------------
C     MAIN RANDOM NUMBER GENERATOR
C     USES METHOD OF l'Ecuyer, (VIA F.JAMES, COMP PHYS COMM 60(1990)329)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION RANGEN,RANSET,RANGET
      INTEGER I,ISEED(2),K,IZ,JSEED(2)
      SAVE ISEED
      DATA ISEED/12345,67890/
      K=ISEED(1)/53668
      ISEED(1)=40014*(ISEED(1)-K*53668)-K*12211
      IF (ISEED(1).LT.0) ISEED(1)=ISEED(1)+2147483563
      K=ISEED(2)/52774
      ISEED(2)=40692*(ISEED(2)-K*52774)-K*3791
      IF (ISEED(2).LT.0) ISEED(2)=ISEED(2)+2147483399
      IZ=ISEED(1)-ISEED(2)
      IF (IZ.LT.1) IZ=IZ+2147483562
      RANGEN=DBLE(IZ)*4.656613001013252D-10
C--->                (4.656613001013252D-10 = 1.D0/2147483589)
      RETURN
C-----------------------------------------------------------------------
      ENTRY RANSET(JSEED)
C-----------------------------------------------------------------------
      IF (JSEED(1).EQ.0.OR.JSEED(2).EQ.0) then
         write(6,*) 'Jseeds=0, wrong settings for RANSET'
         stop
      endif
      ISEED(1)=JSEED(1)
      ISEED(2)=JSEED(2)
 999  RETURN
C-----------------------------------------------------------------------
      ENTRY RANGET(JSEED)
C-----------------------------------------------------------------------
      JSEED(1)=ISEED(1)
      JSEED(2)=ISEED(2)
      RETURN
      END

*-- Author :    F. James, modified by Mike Seymour
C-----------------------------------------------------------------------
      FUNCTION RANGEN2(I)
C-----------------------------------------------------------------------
C     MAIN RANDOM NUMBER GENERATOR
C     USES METHOD OF l'Ecuyer, (VIA F.JAMES, COMP PHYS COMM 60(1990)329)
C-----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION RANGEN2,RANSET2,RANGET2
      INTEGER I,ISEED(2),K,IZ,JSEED(2)
      SAVE ISEED
      DATA ISEED/12347,67892/
      K=ISEED(1)/53668
      ISEED(1)=40014*(ISEED(1)-K*53668)-K*12211
      IF (ISEED(1).LT.0) ISEED(1)=ISEED(1)+2147483563
      K=ISEED(2)/52774
      ISEED(2)=40692*(ISEED(2)-K*52774)-K*3791
      IF (ISEED(2).LT.0) ISEED(2)=ISEED(2)+2147483399
      IZ=ISEED(1)-ISEED(2)
      IF (IZ.LT.1) IZ=IZ+2147483562
      RANGEN2=DBLE(IZ)*4.656613001013252D-10
C--->                (4.656613001013252D-10 = 1.D0/2147483589)
      RETURN
C-----------------------------------------------------------------------
      ENTRY RANSET2(JSEED)
C-----------------------------------------------------------------------
      IF (JSEED(1).EQ.0.OR.JSEED(2).EQ.0) then
         write(6,*) 'Jseeds=0, wrong settings for RANSET'
         stop
      endif
      ISEED(1)=JSEED(1)
      ISEED(2)=JSEED(2)
 999  RETURN
C-----------------------------------------------------------------------
      ENTRY RANGET2(JSEED)
C-----------------------------------------------------------------------
      JSEED(1)=ISEED(1)
      JSEED(2)=ISEED(2)
      RETURN
      END


C-----------------------------------------------------------------------
      subroutine aluend(iunit)
c     go to end of file
C-----------------------------------------------------------------------
      ios = 0    
      dowhile(ios.eq.0)
         read(unit=iunit,fmt='(1x)',iostat=ios)
      enddo                        
      end
C-----------------------------------------------------------------------
      subroutine alugun(n)
c     look for next available unit
C-----------------------------------------------------------------------
      implicit none
      integer n,i
      logical yes
      do i=10,100
        inquire(unit=i,opened=yes)
        if(.not.yes) goto 10
      enddo
      write(6,*) 'no free units to write to available, stop'
      stop
 10   n=i
      end
C---------------- -------------------------------------------------------
      subroutine alustn(string,num)
c- writes the number num on the string string starting at the blank
c- following the last non-blank character
      character * (*) string
      character * 20 tmp
      integer aluisl
      l = len(string)
      write(tmp,'(i15)')num
      j=1
      dowhile(tmp(j:j).eq.' ')
        j=j+1
      enddo
      ipos = aluisl(string)
      ito = ipos+1+(15-j)
      if(ito.gt.l) then
         write(6,*)'error, string too short'
         write(6,*) string
         stop
      endif
      string(ipos+1:ito)=tmp(j:)
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

      subroutine alustc(str1,str2,str)
c concatenates str1 and str2 into str. Ignores trailing blanks of str1,str2
      integer aluisl
      character *(*) str1,str2,str
      l1=aluisl(str1)
      l2=aluisl(str2)
      l =len(str)
      if(l.lt.l1+l2) then
          write(6,*) 'error: l1+l2>l in alustc'
          write(6,*) 'l1=',l1,' str1=',str1
          write(6,*) 'l2=',l2,' str2=',str2
          write(6,*) 'l=',l
          stop
      endif
      if(l1.ne.0) str(1:l1)=str1(1:l1)
      if(l2.ne.0) str(l1+1:l1+l2)=str2(1:l2)
      if(l1+l2+1.le.l) str(l1+l2+1:l)= ' '
      end

c-------------------------------------------------------------------
      function unwgt(w,wmax)
c     routine for the unweighting of weighted events.
c     unwgt=.true.   ==> event was selected
c-------------------------------------------------------------------
      implicit none
      logical unwgt
      double precision w,wmax,rn,rangen2
      unwgt=.false.
      rn=rangen2(1)
      if(w.ge.rn*wmax) unwgt=.true.
      end



      subroutine alsinp
      implicit none
      include 'alpgen.inc'
      character*50 tmpstr
      integer l,aluisl
c
      call alustc(fname,'.err',tmpstr)
      call alugun(nioerr)
      open(unit=nioerr,file=tmpstr,status='unknown')
        write(nioerr,*)'Log file for raised exceptions'
        write(nioerr,*)
      close(nioerr)
c

      if(imode.eq.1) then
         call alustc(fname,'.par',tmpstr)
         call alugun(niopar)
         open(unit=niopar,file=tmpstr,status='unknown')
         l=aluisl(tmpstr)
         write(6,*) 'Run parameters and diagnostics written to '
     $        ,tmpstr(1:l)
         call alustc(fname,'.wgt',tmpstr)
         call alugun(niowgt)
         open(unit=niowgt,file=tmpstr,status='unknown')
         l=aluisl(tmpstr)
         write(6,*) 'Wgted events are written to ',tmpstr(1:l)
       elseif(imode.eq.2) then
        call alustc(fname,'.par',tmpstr)
        call alugun(niopar)
        open(unit=niopar,file=tmpstr,status='old')
        call alustc(fname,'.wgt',tmpstr)
        call alugun(niowgt)
        open(unit=niowgt,file=tmpstr,status='old')
        l=aluisl(tmpstr)
        write(6,*) 'Wgted events are read from ',tmpstr(1:l)
      endif
      end

      subroutine alsout
      implicit none
      include 'alpgen.inc'
      character*5 nversion
      common/version/nversion
      character*50 tmpstr
      integer l,aluisl
      if(imode.lt.2) then
         call alustc(fname,'.stat',tmpstr)
         l=aluisl(tmpstr)
         write(6,*) 'Statistics and results written to ',tmpstr(1:l)
         call alugun(niosta)
         open(unit=niosta,file=tmpstr,status='unknown')
         write(niosta,*) 'ALPGEN, Version: ',nversion
         call alustc(fname,'.top',tmpstr)
         l=aluisl(tmpstr)
         write(6,*) 'Topdrawer plots (if any) written to ',tmpstr(1:l)
         topfile=tmpstr
      elseif(imode.eq.2) then
        call alustc(fname,'_unw.top',tmpstr)
        l=aluisl(tmpstr)
        write(6,*) 'Topdrawer plots (if any) written to ',tmpstr(1:l)
        topfile=tmpstr
c open output files
        if(ilhe.eq.0) then
          call alustc(fname,'_unw.par',tmpstr)
          call alugun(niosta)
          open(unit=niosta,file=tmpstr,status='unknown')
        elseif(ilhe.eq.1) then
          call alustc(fname,'.lhh',tmpstr)
          call alugun(niosta)
          open(unit=niosta,file=tmpstr,status='unknown')
          write(niosta,'(a)') '<LesHouchesEvents version="1.0">'
          write(niosta,'(a)') '<!--'
        else
          write(*,*) 'undefined out format, stop'
          stop
        endif
        l=aluisl(tmpstr)
        write(6,*) 'Parameters and results written to ',tmpstr(1:l)
        call alugun(niounw)
        if(ilhe.eq.0) then
          call alustc(fname,'.unw',tmpstr)
          open(unit=niounw,file=tmpstr,status='unknown')
        elseif(ilhe.eq.1) then
          call alustc(fname,'.lhu',tmpstr)
          open(unit=niounw,file=tmpstr,status='unknown')
        endif
        l=aluisl(tmpstr)
        write(6,*) 'Unwgted events are written to ',tmpstr(1:l)
        if(ilhe.eq.1) then
          call alustc(fname,'.lhe',tmpstr)
          l=aluisl(tmpstr)
          write(6,*) 'Params and event files merged at the end',
     $    ' of the run into  ',tmpstr(1:l)
        endif
      endif
      end
      
      
      subroutine refinemom(n,p,iflag)                           
*                                                                    
c-    Refine the momenta in such a way they are                 
c-    on-shell and energy-momentum conserving up to             
c-    the computer precision.                                   
c-    n is the total number of momenta (initial+final state)    
*                                                                    
      implicit none                                             
      integer maxpar,n,i,iflag                                  
      parameter (maxpar= 20)                                    
      real*8 p(5,maxpar),ptotx,ptoty,ptotz,etot                 
*                                                               
      iflag= 0                                                  
      ptotx= 0.d0                                               
      ptoty= 0.d0                                               
      ptotz= 0.d0                                               
      do i= 3,n                                                 
         ptotx= ptotx+p(1,i)                                    
         ptoty= ptoty+p(2,i)                                    
         ptotz= ptotz+p(3,i)                                    
      enddo                                                     
*                                                                    
      etot= 0.d0                                                
      do i= 3,n                                                 
        p(1,i)= p(1,i)-ptotx/dfloat(n-2)                        
        p(2,i)= p(2,i)-ptoty/dfloat(n-2)                        
        p(4,i)= sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2+p(5,i)**2)   
        etot  = etot+p(4,i)                                     
      enddo                                                     
*                                                                    
      if (dabs(ptotz).gt.etot) then                             
        iflag= 1                                                
      else                                                      
        p(4,1)= 0.5d0*(etot+ptotz)                              
        p(3,1)= p(4,1)                                          
        p(4,2)= 0.5d0*(etot-ptotz)                              
        p(3,2)=-p(4,2)                                          
      endif                                                     
      return                                                    
      end                                                       




c-------------------------------------------------------------------
      SUBROUTINE USRTIME(NOW)
c-------------------------------------------------------------------
c     this subroutine should return the current time, in the format
c     "Thu Apr  4 12:10:45 2002". 
c     The present version is based on use of the Intrinsic Time and
c     CTime functions, supported by most Unix and GNU compilers.  The
c     user may need to change them should they not be supported by her
c     compiler. The sole purpose of this routine is to provide the time
c     information when outputting the topdrawer histograms, so if you
c     don't want to bother finding a Time routine supported by your
c     compiler, you can uncomment the dummy value assigned below to NOW
c     ,and comment away the call to ctime
c     
      IMPLICIT NONE
      integer time
      character*24 CTIME,now
c      now='Day Mon XX hh:mm:ss yyyy'
      now=''
c      now=ctime(time())
      end


*CMZ :-        -05/11/95  19.33.42  by  Mike Seymour
*-- Author :    Adapted by Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE ALULB4(PS,PI,PF)
C-----------------------------------------------------------------------
C     TRANSFORMS PI (GIVEN IN REST FRAME OF PS) INTO PF (IN LAB)
C     N.B. P(1,2,3,4) = (PX,PY,PZ,E); PS(5)=M
C-----------------------------------------------------------------------
      DOUBLE PRECISION PF4,FN,PS(5),PI(4),PF(4)
      IF (PS(4).EQ.PS(5)) THEN
        PF(1)= PI(1)
        PF(2)= PI(2)
        PF(3)= PI(3)
        PF(4)= PI(4)
      ELSE
        PF4  = (PI(1)*PS(1)+PI(2)*PS(2)
     &         +PI(3)*PS(3)+PI(4)*PS(4))/PS(5)
        FN   = (PF4+PI(4)) / (PS(4)+PS(5))
        PF(1)= PI(1) + FN*PS(1)
        PF(2)= PI(2) + FN*PS(2)
        PF(3)= PI(3) + FN*PS(3)
        PF(4)= PF4
      END IF
      END
*CMZ :-        -05/11/95  19.33.42  by  Mike Seymour
*-- Author :    Adapted by Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE ALULF4(PS,PI,PF)
C-----------------------------------------------------------------------
C     TRANSFORMS PI (GIVEN IN LAB) INTO PF (IN REST FRAME OF PS)
C     N.B. P(1,2,3,4) = (PX,PY,PZ,E); PS(5)=M
C-----------------------------------------------------------------------
      DOUBLE PRECISION PF4,FN,PS(5),PI(4),PF(4)
      IF (PS(4).EQ.PS(5)) THEN
        PF(1)= PI(1)
        PF(2)= PI(2)
        PF(3)= PI(3)
        PF(4)= PI(4)
      ELSE
        PF4  = (PI(4)*PS(4)-PI(3)*PS(3)
     &         -PI(2)*PS(2)-PI(1)*PS(1))/PS(5)
        FN   = (PF4+PI(4)) / (PS(4)+PS(5))
        PF(1)= PI(1) - FN*PS(1)
        PF(2)= PI(2) - FN*PS(2)
        PF(3)= PI(3) - FN*PS(3)
        PF(4)= PF4
      END IF
      END
*CMZ :-        -05/11/95  19.33.42  by  Mike Seymour
*-- Author :    Adapted by Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE ALULOB(PS,PI,PF)
C-----------------------------------------------------------------------
C     TRANSFORMS PI (GIVEN IN REST FRAME OF PS) INTO PF (IN LAB)
C     N.B. P(1,2,3,4,5) = (PX,PY,PZ,E,M)
C-----------------------------------------------------------------------
      DOUBLE PRECISION PS(5),PI(5),PF(5)
      CALL ALULB4(PS,PI,PF)
      PF(5)= PI(5)
      END
CDECK  ID>, ALULOF.
*CMZ :-        -05/11/95  19.33.42  by  Mike Seymour
*-- Author :    Adapted by Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE ALULOF(PS,PI,PF)
C-----------------------------------------------------------------------
C     TRANSFORMS PI (GIVEN IN LAB) INTO PF (IN REST FRAME OF PS)
C     N.B. P(1,2,3,4,5) = (PX,PY,PZ,E,M)
C-----------------------------------------------------------------------
      DOUBLE PRECISION PS(5),PI(5),PF(5)
      CALL ALULF4(PS,PI,PF)
      PF(5)= PI(5)
      END
CDECK  ID>, ALULOR.
*CMZ :-        -26/04/91  11.11.56  by  Bryan Webber
*-- Author :    Giovanni Abbiendi & Luca Stanco
C-----------------------------------------------------------------------
      SUBROUTINE ALULOR (TRANSF,PI,PF)
C-----------------------------------------------------------------------
C     Makes the ALULOR transformation specified by TRANSF on the
C     quadrivector PI(5), giving PF(5).
C-----------------------------------------------------------------------
      DOUBLE PRECISION TRANSF(4,4),PI(5),PF(5)
      INTEGER I,J
      DO 1 I=1,5
        PF(I)=0.D0
    1 CONTINUE
      DO 3 I=1,4
       DO 2 J=1,4
         PF(I) = PF(I) + TRANSF(I,J) * PI(J)
    2  CONTINUE
    3 CONTINUE
      PF(5) = PI(5)
      RETURN
      END
CDECK  ID>, ALUMAS.
*CMZ :-        -26/04/91  11.11.56  by  Bryan Webber
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      SUBROUTINE ALUMAS(P)
C-----------------------------------------------------------------------
C     PUTS INVARIANT MASS IN 5TH COMPONENT OF VECTOR
C     (NEGATIVE SIGN IF SPACELIKE)
C-----------------------------------------------------------------------
      DOUBLE PRECISION ALUSQR,P(5)
      P(5)=SQRT((P(4)+P(3))*(P(4)-P(3))-P(1)**2-P(2)**2)
      END
*CMZ :-        -26/04/91  11.11.56  by  Bryan Webber
*-- Author :    Bryan Webber
C-----------------------------------------------------------------------
      FUNCTION ALUPCM(EM0,EM1,EM2)
C-----------------------------------------------------------------------
C     C.M. MOMENTUM FOR DECAY MASSES EM0 -> EM1 + EM2
C     SET TO -1 BELOW THRESHOLD
C-----------------------------------------------------------------------
      DOUBLE PRECISION ALUPCM,EM0,EM1,EM2,EMS,EMD
      EMS=ABS(EM1+EM2)
      EMD=ABS(EM1-EM2)
      IF (EM0.LT.EMS.OR.EM0.LT.EMD) THEN
        ALUPCM=-1.
      ELSEIF (EM0.EQ.EMS.OR.EM0.EQ.EMD) THEN
        ALUPCM=0.
      ELSE
        ALUPCM=SQRT((EM0+EMD)*(EM0-EMD)*
     &              (EM0+EMS)*(EM0-EMS))*.5/EM0
      ENDIF
      END

      subroutine rescms(p,p1,p2,m1,m2)
      implicit none
      integer il
      double precision p(5),p1(5),p2(5),m1,m2
      double precision po1(5),po2(5),m,mo1,mo2,pcm,pcmo
      double precision alupcm
      m=p(5)
      mo1=p1(5)
      mo2=p2(5)
      call alulof(p,p1,po1)
      call alulof(p,p2,po2)
      if(max(mo1,mo2).lt.1d-6) then
        pcmo=m/2d0
      else 
        pcmo=alupcm(m,mo1,mo2)
      endif 
      pcm=alupcm(m,m1,m2)
      do il=1,3
        po1(il)=pcm/pcmo*po1(il)
        po2(il)=pcm/pcmo*po2(il)
      enddo
      po1(4)=sqrt(pcm**2+m1**2)
      po2(4)=sqrt(pcm**2+m2**2)
      po1(5)=m1
      po2(5)=m2
      call alulob(p,po1,p1)
      call alulob(p,po2,p2)
      end


c     start clustering routines
c
      subroutine alsclu
c     setup variables for clustering:
c     iopt=1: quarks and gluons only
c     iopt=2: flavour clustering, quarks/gluons and EW particles
      implicit none
      include 'alpgen.inc'
c     locals
      integer i,j,if1,ini,nfcol,nweak
c
      integer mode
c
      data ini/0/
c
c      if(ihrd.gt.4) return
c     count coloured particles (quarks and gluons)
      if(ini.eq.0) then
        nclext=0
        do i=1,npart
          if(abs(ifl(i)).le.6.or.ifl(i).eq.21) nclext=nclext+1
        enddo
c non strongly interacting final state particles
c        nweak=nw+nh+nz+nph
c strongly interacting final state particles
c        nfcol=nclext-2
c        ncltot=2*nfcol+min(2,nweak)
c
c determine total number of clusters 
c     This is equal to the number of colore particles for the LO tree
c     process, plus twice the number of jets beyond the tree: 
c     Ncol(LO tree)+2njets
c
c w/z qq+jets
        if(ihrd.le.2) then
          ncltot=4+2*njets
        elseif(ihrd.le.4.or.ihrd.eq.10) then
          ncltot=2+2*njets
        elseif(ihrd.eq.5) then
c          if(nh.eq.0) then
c the LO tree process is always qqbar->nw W + nh H + nz Z + nph Photons
            ncltot=2+2*njets
c     for EW processes, e.g. qq -> H qq, ncltot=4+2(njets-2)=2njets =>
c     ncltot-2, which is implemented by the ewk='true' option in alcclu
c
c old code:
c          else
c            if(nh.eq.nweak) then
c     qq->Hqq => purely EW
c              ncltot=2*(nfcol-2)+4
c            else
c              ncltot=2+2*nfcol
c            endif
c          endif
c     2Q
        elseif(ihrd.eq.6) then
          ncltot=4+2*njets
c     4Q
        elseif(ihrd.eq.7) then
          if(max(ihvy,ihvy2).eq.6.and.min(ihvy,ihvy2).ne.6) then
            ncltot=4+2*(njets+2)
          else
            ncltot=6+2*njets
          endif
c     QQh
        elseif(ihrd.eq.8) then
          ncltot=4+2*njets
c     Njet
        elseif(ihrd.eq.9) then
          ncltot=4+2*(njets-2)
c     gamma jets
        elseif(ihrd.eq.11) then
          ncltot=3+2*(njets-1)
c     H jets
        elseif(ihrd.eq.12) then
          ncltot=2+2*njets
c     W gamma jets
        elseif(ihrd.eq.14) then
          ncltot=2+2*njets
c     W gamma QQ jets
        elseif(ihrd.eq.15) then
          ncltot=4+2*njets
        endif
c
        naspow=ncltot-nclext
        ini=1
      endif
      if(nclext.eq.ncltot) return
c
      kres=qsq
c
      mode=0
c
      do i=1,nclext
        clkt(1,i)=kres
        clkt(2,i)=kres
        icst(i)=0
        icmot(i)=0
        if1=ifl(i)
        if(if1.eq.21) if1=0
        icfl(i)=if1
        icfs(i)=1
        if(i.le.2) icfs(i)=-1
        do j=1,4
          clp(j,i)=p(j,i)
        enddo
        do j=1,2
          icdau(j,i)=i
          iccol(j,i)=icu(j,i)
        enddo
      enddo
      do i=nclext+1,ncltot
        clkt(1,i)=kres
        clkt(2,i)=kres
        icst(i)=0
        icmot(i)=0
        icfl(i)=0
        icfs(i)=1
        do j=1,4
          clp(j,i)=0
        enddo
        do j=1,2
          icdau(j,i)=i
          iccol(j,i)=0
        enddo
      enddo
      end
c********************************************************************
      subroutine antennas(ewk,icu,ncl,idxant)
c********************************************************************
      implicit none
c
      logical ewk
      integer icu(2,20),ncl,idxant(20)
c
      integer j1,nq,nqb,iq(20),iqb(20),i1,i2,j2,icu0(2,20)
c
      do j1=1,ncl
        if(j1.le.2) then
          icu0(1,j1)=icu(2,j1)
          icu0(2,j1)=icu(1,j1)
        else
          icu0(1,j1)=icu(1,j1)
          icu0(2,j1)=icu(2,j1)
        endif
      enddo
c
      nq=0
      nqb=0
      do j1=1,ncl
        if(icu0(1,j1).eq.0) then
          if(icu0(2,j1).eq.0) then
            write(*,*)'something going wrong with clustering'
            stop
          endif
          idxant(j1)=-1
          nq=nq+1
          iq(nq)=j1
        elseif(icu0(2,j1).eq.0) then
          idxant(j1)=-1
          nqb=nqb+1
          iqb(nqb)=j1
        else
          idxant(j1)=1
        endif
      enddo
      if(nq.ne.nqb) then
        write(*,*)'something wrong with clustering'
        stop
      endif
c
      if(nq.le.1) then
        ewk=.false.
        return
      elseif(nq.eq.2) then
c
        i1=min(icu0(2,iq(1)),icu0(2,iq(2)))
        i2=min(icu0(1,iqb(1)),icu0(1,iqb(2)))
        j1=min(i1,i2)
        if(i2.eq.j1) then
          i2=i1
          i1=j1
        endif
        do j1=1,ncl
          j2=max(icu0(1,j1),icu0(2,j1))
          if(j2.le.i2.and.j2.ge.i1) then
          else
            idxant(j1)=idxant(j1)*2
          endif
        enddo
c
      else
        write(*,*)'ew clustering with six quarks not yet implemented'
        stop
      endif
c
      return
      end
c      
      
      subroutine clumrg(p1,p2,is1,is2,pclu,mrgopt)
c     merge pseudoparticle momenta p1&p2 into new cluster momentum
      implicit none
c     inputs
      integer is1,is2,mrgopt
      double precision p1(4),p2(4),pclu(4)
c     locals
      integer k,is
      is=is1*is2
      do k=1,4
        pclu(k)=is1*p1(k)+is2*p2(k)
        if(is.lt.0) pclu(k)=-pclu(k)
      enddo
      if(mrgopt.eq.2) then
c keep massless cluster
        pclu(4)=sqrt((is*p1(1)+p2(1))**2+(is*p1(2)+p2(2))**2+(is*p1(3)
     $       +p2(3))**2)
      endif
      end


      subroutine alcclu(ewk,iflag)
c colour driven clustering
      implicit none
      include 'alpgen.inc'
c inputs
      integer iflag
      logical ewk
c local
      integer i,j,nmax,nclmax,is,imin,jmin,istart,if1,if2,ic1(2)
     $     ,ic2(2),icol1,icol2,icmin1,icmin2
      double precision ktmin,kperp,tmp,tiny
      parameter (tiny=1d-2)
c
      integer mrgall(20),ncltot0
c debugging
c      integer nev
c      data nev/0/
c
c initialise
c      nev=nev+1
      iflag=0
c
      if(ewk) then
        call antennas(ewk,icu,nclext,mrgall)
      else
        do i=1,nclext
          mrgall(i)=1
        enddo
      endif
      ncltot0=ncltot
      if(ewk) ncltot0=ncltot0-2
      if(nclext.gt.ncltot0) then
        write(*,*)'problems with clustering'
        stop
      elseif(nclext.eq.ncltot0) then
        return
      endif
c
      nclmax=ncltot0
c      write(6,*) 'event=',nev
      nmax=nclext
      if(nmax.eq.nclmax) return
 1    continue
      ktmin=1d12
c
      do i=1,nmax-1
        if(icst(i).eq.1) goto 10
        istart=i
        if1=icfl(i)*icfs(i)
        ic1(1)=iccol(1,i)
        ic1(2)=iccol(2,i)
c inhibit clustering of initial-initial states
        if(i.eq.1) istart=2
        do j=istart+1,nmax
          if(icst(j).eq.1) goto 9
          if(icfs(i).lt.0. and .icfs(j).lt.0) goto 9
          if2=icfl(j)*icfs(j)
          ic2(1)=iccol(1,j)
          ic2(2)=iccol(2,j)
          is=icfs(i)*icfs(j)
c     check clusterability of i+j and assign colours to potential new
C     cluster
c
c     inhibit clustering among partons belonging to different colour singlets
          if(abs(mrgall(i)).ne.abs(mrgall(j))) goto 9
c     inhibit q qbbr clustering
          if(mrgall(i)+mrgall(j).lt.0) goto 9
c     do not allow W-g clusters
          if( (if1.eq.0.and.if2.eq.-1000) . or.
     $        (if2.eq.0.and.if1.eq.-1000))  goto 9
          if(iccol(1,i).eq.-1) then
c     i is an EW object, cluster inherits j colours
            icol1=iccol(1,j)
            icol2=iccol(2,j)
          elseif(iccol(1,j).eq.-1) then
c     j is an EW object, cluster inherits i colours
            icol1=iccol(1,i)
            icol2=iccol(2,i)
          elseif(is.lt.0) then
c initial-final cluster first
            if((if1+if2)*if1*if2.ne.0) goto 9
            if(ic1(1).ne.0.and.ic1(1).eq.ic2(1)) then
              icol1=ic1(2)*(1+icfs(i))/2+ic2(2)*(1+icfs(j))/2
              icol2=ic1(2)*(1-icfs(i))/2+ic2(2)*(1-icfs(j))/2
c              write(6,*) 'passed 1'
            elseif(ic1(2).ne.0.and.ic1(2).eq.ic2(2)) then
              icol1=ic1(1)*(1-icfs(i))/2+ic2(1)*(1-icfs(j))/2
              icol2=ic1(1)*(1+icfs(i))/2+ic2(1)*(1+icfs(j))/2
c              write(6,*) 'passed 2'
            elseif(ic1(1).eq.ic1(2)) then
              icol1=ic2(1)*(1-icfs(j))/2+ic2(2)*(1+icfs(j))/2
              icol2=ic2(2)*(1-icfs(j))/2+ic2(1)*(1+icfs(j))/2
c              write(6,*) 'passed 3'
            elseif(ic2(1).eq.ic2(2)) then
              icol1=ic1(1)*(1-icfs(i))/2+ic1(2)*(1+icfs(i))/2
              icol2=ic1(2)*(1-icfs(i))/2+ic1(1)*(1+icfs(i))/2
c              write(6,*) 'passed 4'
c q -> q g(cluster)
            elseif(if1+if2.eq.0.and.if1.ne.0) then
              if(ic1(1).ne.0) then
                icol1=ic1(1)
                icol2=ic2(1)
              elseif(ic1(2).ne.0) then
                icol1=ic2(2)
                icol2=ic1(2)
              endif
c              write(6,*) 'passed 5'
            else
              goto 9
            endif
          else
c final-final cluster then
            if((if1+if2)*if1*if2.ne.0) goto 9
            if(ic1(1).ne.0.and.ic1(1).eq.ic2(2)) then
              icol1=ic2(1)
              icol2=ic1(2)
c              write(6,*) 'passed 11'
            elseif(ic1(2).ne.0.and.ic1(2).eq.ic2(1)) then
              icol1=ic1(1)
              icol2=ic2(2)
c              write(6,*) 'passed 12'
            elseif(ic1(1).eq.ic1(2)) then
              icol1=ic2(1)
              icol2=ic2(2)
c              write(6,*) 'passed 13'
            elseif(ic2(1).eq.ic2(2)) then
              icol1=ic1(1)
              icol2=ic1(2)
c              write(6,*) 'passed 14'
c g-> q qbar clusters
            elseif(if1+if2.eq.0.and.if1.ne.0) then
              if(ic1(1).ne.0) then
                icol1=ic1(1)
                icol2=ic2(2)
              elseif(ic2(1).ne.0) then
                icol1=ic2(1)
                icol2=ic1(2)
              endif
c              write(6,*) 'passed 15'
            else
              goto 9
            endif
          endif
c
 8        tmp=kperp(clp(1,i),clp(1,j),icfs(i),icfs(j),cluopt,ktfac)
c          write(6,*) (clp(imin,i),imin=1,3)
c          write(6,*) (clp(imin,j),imin=1,3)
c          write(6,*) 'i,j:',i,j,' kt=',sqrt(tmp)
c          if(tmp.lt.kres) then
c below resolution, reject event
c            iflag=1
c            return
c          endif
c     resolve the symmetry between merging with either beam-like
C     particle. Make more probable the configuration with jet and beam
C     traveling in the same z direction
c          if(is.lt.0) then
c            tmp=tmp-sign(tiny,clp(3,i)*clp(3,j))
c          endif
          if(ktmin.gt.tmp) then
            ktmin=tmp
            imin=i
            jmin=j
            if(icol1.eq.icol2) then
              icmin1=0
              icmin2=0
            else
              icmin1=icol1
              icmin2=icol2
            endif
c            write(6,*) 'new min kt:',ktmin,' pair:',i,j
          endif
 9        continue
        enddo
 10     continue 
      enddo
c create new cluster
c      write(6,*) 'new cluster:',imin,jmin,' ktmin=',ktmin
c no new cluster found
      if(ktmin.eq.1d12) goto 100
      nmax=nmax+1
c
      mrgall(nmax)=min(mrgall(imin),mrgall(jmin))
c
      icdau(1,nmax)=imin
      icdau(2,nmax)=jmin
      clkt(2,nmax)=ktmin
      is=icfs(imin)*icfs(jmin)
      icfl(nmax)=is*(icfs(imin)*icfl(imin)+icfs(jmin)*icfl(jmin))
      icfs(nmax)=is
c     merge pesudoparticles into cluster
      call clumrg(clp(1,imin),clp(1,jmin),icfs(imin),icfs(jmin),clp(1
     $     ,nmax),mrgopt)
c update i,j status
      icst(imin)=1
      icst(jmin)=1
      clkt(1,imin)=ktmin
      clkt(1,jmin)=ktmin
      icmot(imin)=nmax
      icmot(jmin)=nmax
c cluster colour
      iccol(1,nmax)=icmin1
      iccol(2,nmax)=icmin2
      if(nmax.lt.nclmax) goto 1
      return
 100  continue
c leave incomplete clusters with scale at ktres
c      write(6,*) 'incomplete cluster'
c completed clustering
      end

      subroutine cktmin
c determine ktmin for the event
      implicit none
      include 'alpgen.inc'
c inputs
      double precision ktmin
c local
      integer i,j,nmax,nclmax,is,imin,jmin,istart,ifs1,ifs2,ic1(2)
     $     ,ic2(2),ifl1,ifl2
      double precision kperp,tmp
c
      call alsclu
      nmax=nclext
      ktmin=qsq
      do i=1,nmax-1
        istart=i
        ifs1=icfs(i)
        ifl1=icfl(i)*ifs1
        ic1(1)=iccol(1,i)
        ic1(2)=iccol(2,i)
        if(ic1(1).eq.0.and.ic1(2).eq.0) goto 10
        do j=istart+1,nmax
          ifs2=icfs(j)
          ifl2=icfl(j)*ifs2
          ic2(1)=iccol(1,j)
          ic2(2)=iccol(2,j)
          if(ic2(1).eq.0.and.ic2(2).eq.0) goto 9
c          is=ifs1+ifs2
c     check clusterability of i+j:
c accept all  qg or gg clusters, check flavour for  qqbar clusters
          if(ifl1*ifl2*(ifl1+ifl2).ne.0) goto 9
c            if(ifl1+ifl2.eq.0) goto 8
c          else
c            goto 8
c            if(is.eq.0) then
c initial-final
c              if(ic1(1).eq.ic2(1).or.ic1(2).eq.ic2(2)) goto 8
c            else
c final-final or initial-initial
c              if(ic1(1).eq.ic2(2).or.ic1(2).eq.ic2(1)) goto 8
c            endif
c          endif
c          goto 9
c
 8        tmp=kperp(clp(1,i),clp(1,j),ifs1,ifs2,cluopt,ktfac)
c          write(6,*) 'i,j:',i,j,' kt=',tmp
          ktmin=min(ktmin,tmp)
 9        continue
        enddo
 10     continue 
      enddo
      kres=ktmin
      end



      function kperp(p1,p2,is1,is2,iopt,qscale)
      implicit none
      double precision kperp,p1(4),p2(4),qscale
      integer is1,is2,iopt
c local variables
      double precision pt1,pt2,y1,y2,dphi,mt1,mt2,q1,q2,m1,m2,tmp
c commons
c
      pt1=p1(1)**2+p1(2)**2
      pt2=p2(1)**2+p2(2)**2
      mt1=p1(4)**2-p1(3)**2
      mt2=p2(4)**2-p2(3)**2
      m1=mt1-pt1
      m2=mt2-pt2
      if(iopt.eq.1) then
        q1=pt1
        q2=pt2
      elseif(iopt.eq.2) then
        q1=mt1
        q2=mt2
      endif
      if(is1.gt.0.and.is2.gt.0) then
c     final-final
        y1=-0.5*log((p1(4)-p1(3))/(p1(4)+p1(3)))
        y2=-0.5*log((p2(4)-p2(3))/(p2(4)+p2(3)))
        dphi=acos( (p1(1)*p2(1)+p1(2)*p2(2) )/sqrt(pt1*pt2+1d-3))
        kperp=((y2-y1)**2+dphi**2)*min(q1,q2)
c     initial-final
      elseif(is1.lt.0.and.is2.gt.0) then
c     kperp=mt2
        kperp=q2+abs(q1-pt1)
c     final-initial
      elseif(is2.lt.0.and.is1.gt.0) then
c     kperp=mt1
        kperp=q1+abs(q2-pt2)
c     initial-initial
      else
        kperp=(p1(4)+p2(4))**2-(p1(1)+p2(1))**2-(p1(2)+p2(2))**2-(p1(3
     $       )+p2(3))**2
      endif
c     protect against sqrt(kperp) falling below the mass of final state
C     and below 2 gev
c
      kperp=qscale**2*max(kperp,abs(m1)+abs(m2))
      kperp=max(kperp,4d0)
      end

      subroutine alasrs(ewk,rewgt)
      implicit none
      include 'alpgen.inc'
c arguments
      double precision rewgt
      logical ewk
c locals
      integer i,ncltot0
      double precision alfas_clu
c
      ncltot0=ncltot
      if(ewk) ncltot0=ncltot-2
      rewgt=1d0
      do i=nclext+1,ncltot0
c        if(clkt(2,i).ne.kres) write(2,*) 'rewgt',i,alfas(clkt(2,i),xlclu
c     $       ,lpclu,-1)/asmax
        if(abs(clkt(2,i)-kres).gt.0.1) rewgt=rewgt*alfas_clu(clkt(2,i),
     $       xlclu,lpclu,-1)/asmax
      enddo
c
      if(ewk) then
        rewgt=rewgt*(as/asmax)**2
      endif
c
      end


      subroutine rstcol(icol)
c set colour flags for and from clustering
      implicit none
      include 'alpgen.inc'
      integer i,j,icol
      if(icol.eq.1) then
        do j=1,nclext
          do i=1,2
            icu(i,j)=0
          enddo
        enddo
      elseif(icol.eq.2) then
        do j=1,nclext
          do i=1,2
            icu(i,j)=iccol(i,j)
          enddo
        enddo
      endif
      end

      subroutine alpclu(ewk,rewgt)
      implicit none
      include 'alpgen.inc'
      integer iflag
      logical ewk
      double precision rewgt,w
      rewgt=1d0
      call alsclu
      call alcclu(ewk,iflag)
      if(iflag.eq.1) then
        rewgt=0d0
        return
      endif
c rescale alphas only during the weighted event generation
c      if(imode.le.1) then
      call alasrs(ewk,w)
      rewgt=rewgt*w
c      endif
      end
c
      subroutine rndtd
      implicit none
      include 'alpgen.inc'
      real*8 rnd
      integer ptr(4,24),tdec(4),j
      data ptr/1,2,3,4,
     +         1,2,4,3,
     +         1,3,2,4,
     +         1,3,4,2,
     +         1,4,2,3,
     +         1,4,3,2,
     +         2,1,3,4,
     +         2,1,4,3,
     +         2,3,1,4,
     +         2,3,4,1,
     +         2,4,1,3,
     +         2,4,3,1,
     +         3,1,2,4,
     +         3,1,4,2,
     +         3,2,1,4,
     +         3,2,4,1,
     +         3,4,1,2,
     +         3,4,2,1,
     +         4,1,2,3,
     +         4,1,3,2,
     +         4,2,1,3,
     +         4,2,3,1,
     +         4,3,1,2,
     +         4,3,2,1/
      if(itdecmode.lt.10) then
         return
      elseif(itdecmode.lt.100) then
         call randa(rnd)
         if(rnd.gt.0.5) then
            tdec(2)= itdecmode/10
            tdec(1)= itdecmode-tdec(2)*10
            itdecmode= tdec(1)*10+tdec(2)
         endif
      elseif(itdecmode.lt.10000) then
         call randa(rnd)
         j= rnd*24
         j= min(j,23)+1
         tdec(1)= itdecmode/1000
         tdec(2)= (itdecmode-tdec(1)*1000)/100
         tdec(3)= (itdecmode-tdec(1)*1000-tdec(2)*100)/10
         tdec(4)=  itdecmode-tdec(1)*1000-tdec(2)*100-tdec(3)*10
         itdecmode= tdec(ptr(1,j))*1000
     +            + tdec(ptr(2,j))*100
     +            + tdec(ptr(3,j))*10
     +            + tdec(ptr(4,j))
      else
         write(*,*) 'itdecmode value not allowed'
         stop
      endif
      return
      end
c***********************************************************************
      subroutine momcheck(iflag,jseed)
c***********************************************************************
      implicit none
c
      integer iflag,jseed(2)
c
      include 'alpgen.inc'
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      integer j1,j2
      character*50 tmpstr
      real*8 tmp
c
      do j1=1,4
         tmp=p(j1,1)+p(j1,2)
         do j2=3,npart
            tmp=tmp-p(j1,j2)
         enddo
         if(abs(tmp)/(p(4,1)+p(4,2)).gt.1.d-6) then
           call alustc(fname,'.err',tmpstr)
           call alugun(nioerr)
           open(unit=nioerr,file=tmpstr,status='unknown',
     +          access='append')
             write(nioerr,*)'momentum conservation not fulfilled'
             write(nioerr,*)'p  ',p
             write(nioerr,*)'jseed',jseed
             write(nioerr,*)
           close(nioerr)
           iflag=1
           return
         endif
      enddo
C
      do j1=1,npart
        tmp=p(4,j1)**2-p(1,j1)**2-p(2,j1)**2-p(3,j1)**2-p(5,j1)**2
        if(abs(tmp)/p(4,j1)**2.gt.1.d-6) then
           call alustc(fname,'.err',tmpstr)
           call alugun(nioerr)
           open(unit=nioerr,file=tmpstr,status='unknown',
     +          access='append')
             print*,' '
             print*,'momentum onshellness not fulfilled'
             print*,'jproc= ',jproc
             print*,'j1= ',j1
             print*,'mj^2= ',p(4,j1)**2-p(1,j1)**2-p(2,j1)**2-p(3,j1)**2
             print*,'abs(tmp)/p4**2',abs(tmp)/p(4,j1)**2
             print*,'tmp= ',tmp
             print*,'p(4,j1)= ',p(4,j1)
             print*,' '
             write(nioerr,*)'momentum onshellness not fulfilled'
             write(nioerr,*)'p  ',p
             write(nioerr,*)'jseed',jseed
             write(nioerr,*)
           close(nioerr)
           iflag=1
           return
        endif
      enddo
C
      return
      end
c**************************************************************************************
      subroutine raiseexception1(md,nprt,poscol,posfld,col,flv)
c**************************************************************************************
      implicit none
      include 'alpgen.inc'
c
      character*4 md
      integer nprt,poscol,posfld,col(200),flv(100)
c
      integer j1
      character*50 tmpstr
c
      if (nprt.gt.100) then
        write(*,*)'ill dimensioned array in raiseexception1'
        stop
      endif
c
      call alustc(fname,'.err',tmpstr)
      call alugun(nioerr)
      open(unit=nioerr,file=tmpstr,status='unknown',access='append')
c
      write(nioerr,*)
      write(nioerr,*)'ill dim array OFFSH and/or OFFSHCOL in AMP'
      write(nioerr,*)md,' poscol,posfld',poscol,posfld
      write(nioerr,*)'color string', (col(j1),j1=1,2*nprt)
      write(nioerr,*)'flv string', (flv(j1),j1=1,nprt)
c
      close(nioerr)
c
      return
      end
c**************************************************************************************
      subroutine raisenumexcep(label)
c**************************************************************************************
      implicit none
      include 'alpgen.inc'
c
      integer label
c
      integer jseed(2)
      common/js/jseed
      character*50 tmpstr
c
      call alustc(fname,'.err',tmpstr)
      call alugun(nioerr)
      open(unit=nioerr,file=tmpstr,status='unknown',access='append')
c
      write(nioerr,*)
      write(nioerr,*)'raised exception number', label
      write(nioerr,*)'jseed', jseed
c
      close(nioerr)
c
      return
      end
     
