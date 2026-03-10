c process code
      integer ihrd
      integer itopprc
c data format
      integer ilhe
      common/lhefmt/ilhe
c--event file data
      integer Nunit,NunitOut,NunitIni
      character*100 filename,tmpstr
c total number of input events
      common/AHio/Nunit,NunitOut,NunitIni,filename
c inputs for the matching algorithm
      integer iexc,npfst,nplst,nljets,njstart,njlast
     $     ,ickkw
      double precision etclus,rclus,etaclmax
      common/AHopts/etclus,rclus,etaclmax,iexc,npfst
     $     ,nplst,nljets,njstart,njlast,ickkw
c process and particles parameters
      double precision mc,mb,mt,mw,mz,mh
      double precision ebeam
      integer ndns,ih1,ih2       
      integer nw,nz,nh,nph
      integer ihvy,ihvy2,itdec
c pdf set type
      character pdftyp*25
c total number of partons
      integer npart
      common/AHppara/mc,mb,mt,mw,mz,mh,
     & ebeam,ih1,ih2,
     & ihrd,itopprc,itdec,
     & nw,nz,nh,nph,
     & ihvy,ihvy2,
     & npart,ndns,pdftyp
c weight information
      real *8 maxwgt,avgwgt,errwgt,totlum
      integer unwev
      common/AHwgts/maxwgt,avgwgt,errwgt,totlum,unwev
c parton level input parameters
      integer nparam
      parameter (nparam=200)
      integer parlen,partyp
      character chpar*8,chpdes*70
      double precision parval
      common/AHPLpr/parval(nparam),chpar(nparam),chpdes(nparam)
     $     ,parlen(nparam),partyp(nparam)
c shower-related user parameters
      integer nSHparam
      parameter (nSHparam=100)
      integer SHparlen,SHpartyp,SHparuse
      character chSHpar*8,chSHpdes*70
      double precision SHparval
      common/AHSHpr/SHparval(nSHparam),chSHpar(nSHparam)
     $     ,chSHpdes(nSHparam),SHparlen(nSHparam),SHpartyp(nSHparam)
     $     ,SHparuse(nSHparam,0:1)
      integer nevmax,ihvyv,ncmax,nbmax,nvetoc,nvetob,nvetomch,nvetoexc,
     $     nvetohvy,nvetonoj
      common/AHVpar/nevmax,ihvyv,ncmax,nbmax,nvetonoj,nvetoc(100)
     $     ,nvetob(100),nvetomch(0:10),nvetoexc(10),nvetohvy(10)
c other unspecified input params available to user to drive analysis, run, etc
      double precision xpar(0:9)
      common/AHApar/xpar
c global event cuts
      double precision ptjmin,ptjmax,etajmax,drjmin,
     +     ptbmin,ptbmax,etabmax,drbmin,
     +     ptcmin,ptcmax,etacmax,drcmin,
     +     ptphmin,etaphmax,drphjmin,drphmin,drphlmin,
     +     ptlmin,etalmax,drlmin,metmin,mllmin,mllmax
      common/AHcuts/ptjmin,ptjmax,etajmax,drjmin,
     +     ptbmin,ptbmax,etabmax,drbmin,
     +        ptcmin,ptcmax,etacmax,drcmin,
     +        ptphmin,etaphmax,drphjmin,drphmin,drphlmin,
     +        ptlmin,etalmax,drlmin,metmin,mllmin,mllmax
