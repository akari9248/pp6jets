      subroutine Aldpar_external(n)
      include 'alpgen.inc'
      integer n,i,j
c
c THINGS TO DO:
C 1. check the values of the paruse(ivar,ihrd) set by default in Aldpar
C
C     EXAMPLE:  
C      paruse(10,nexternal)= 0
C      paruse(13,nexternal)= 1
C
C 2. define the new proc-specific parameters (only use ipar > nparoff)
C
C     EXAMPLE:  
C      chpar(nparoff+1)='mn'
C      chpdes(nparoff+1)='mass of the heavy neutrino N'
C      partyp(nparoff+1)=0
C      parval(nparoff+1)=100.d0
C      paruse(nparoff+1,nexternal)=1
      return
      end 

      subroutine Alspar_external
C     EXAMPLE:  
c      include 'alpgen.inc'
C      include 'hvyN.inc'
C      mn =parval(nparoff+1) 
      return
      end 

      subroutine alhsca_external(ihrd,iunit)
      implicit none
      integer ihrd,iunit
c$$$      if(ihrd.eq.1) then
c$$$        write(iunit,*) 'iqopt=1 => Q=qfac*sqrt{m_W^2+ sum_jets(m_tr^2)}'
c$$$        write(iunit,*) 'iqopt=2 => Q=qfac*mW'
c$$$        write(iunit,*) 'iqopt=3 => Q=qfac*sqrt{m_W^2+ pt_W^2}'
c$$$        write(iunit,*) 'iqopt=4 => Q=qfac*sqrt{sum_jets(m_tr^2)}'
c$$$        write(iunit,*) 'where:'
c$$$        write(iunit,*) 
c$$$     $ '- m_tr^2=m^2+pt^2, summed over heavy quarks and light jets'
c$$$      endif
      return
      end 

      subroutine Alpdoc_external(ihrd,iunit,imode)
      implicit none
      integer ihrd,iunit,imode
c$$$  EXAMPLE
c$$$      if(imode.le.1.or.imode.gt.2) then
c$$$      elseif(imode.eq.2.or.imode.gt.2) then
c$$$        write(iunit,*) ' '
c$$$        write(iunit,*) ' imode= ',imode
c$$$        write(iunit,*) ' '
c$$$        write(iunit,*) 'For imode=2 and indec= 1:'
c$$$        write(iunit,*) 'select W decay modes, set iwdecmode to:'
c$$$        write(iunit,*) '1: e nu'
c$$$        write(iunit,*) '2: mu nu'
c$$$        write(iunit,*) '3: tau nu'
c$$$        write(iunit,*) '4: e/mu/tau nu'
c$$$        write(iunit,*) '5: q q''bar'
c$$$        write(iunit,*) '6: fully inclusive'
c$$$      endif
      return
      end 
c modified_e
c-------------------------------------------------------------------
      subroutine setalp(npart,flin,pin,flout,pout,instate,ip,ipinv)
c
c     First reorder particles according to alpha ordering
c     and then reorder the respective momenta, etc.
c     Alpha ordering:
c     i= 1    2    3     4  5    6    7    8     9 10 11  12 13 14 15 16
c        ubar dbar nubar e+ cbar sbar tbar bbar  u d  nu  e- c  s  t  b  
c        17     18   19  20  21     22
c        Higgs  Z0   W-  W+  photon gluon
c     (all particles outgoing)
c     ip(i): position of the "i-th" particle of the Alpha array,
c            in the array of momenta defined by routine alegen
c     ipinv(i): inverse of ip, i.e. where to find the "i"-th particle
c            of the alegen array in the alpha momentum array 
c
c     flin is the sequence of flavours required by alpha to fully
c     define a specific process, with the generation ordering.
c
c     flout is the sequence of flavours ordered according to the alpha
c     ordering: flout(i)=flin(ip(i))
c     
c     pout is the sequence of momenta, ordered according to the alpha
c     ordering. All energies are positive, and the pointers to the two
c     momenta corresponding to the intial state partons are instate(i),
c     i=1,2
c     
c-------------------------------------------------------------------
      implicit none
      integer i,j,k,ix,jx,kp,iff,npart,maxpar
      parameter (maxpar=20)
      integer flin(maxpar),flout(maxpar),ip(maxpar),ipinv(maxpar)
      integer instate(2)
      real *8 pin(5,maxpar),pout(4,maxpar)
c     Alpha ordering of outgoing momenta
      integer aref,aord(22)
c     Labels of the particles according to the PDG
c
      data aord/-2,-1,-12,-11,-4,-3,-6,-5,2,1,12,11,4,3,6,5,
     .          25,23,-24,24,22,0/
c
      do i=1,maxpar
       flout(i)=1001
      enddo
c
      if(npart.gt.maxpar) then
         write(6,*) 'error in setper:'
         write(6,*) 'number of particles larger than ',maxpar
         stop
      endif
c     First reorder particles according to alpha ordering
c
c     scan particle types in the Alpha order
      k=1
      do i=1,22
         aref=aord(i)
c     search for all particles of flavour aref in the alegen array, and
c     order them in the "ip" array
         do j=1,npart
            iff=flin(j)
c     turn incoming to outgoing flavours
            if(j.le.2.and.abs(iff).lt.20) iff=-iff
c
            if(iff.eq.aref) then
               ip(k)=j
               ipinv(j)=k
               flout(k)=flin(j)
               if(j.le.2) instate(j)=k
               do ix=1,4
c     alegen momentum conventions:     p(1)=Px p(2)=Py  p(3)=Pz  p(4)=E   
c     alpha momentum conventions:      p(1)=E  p(2)=Px  p(3)=Py  p(4)=Pz
                  jx=ix+1
                  if(jx.eq.5) jx=1
                  pout(jx,k)=pin(ix,j)
               enddo
               if(k.eq.npart) goto 99
               k=k+1
            endif
         enddo
      enddo
 99   continue
c reassign momentum labels to agree with alpha conventions
c                             
      do k=1,npart
        kp=ip(k)
        do ix=1,4
c     alegen momentum conventions:
c     p(1)=Px   p(2)=Py    p(3)=Pz p(4)=E   
c     alpha momentum conventions:
c     p(1)=E  p(2)=Px   p(3)=Py    p(4)=Pz
           jx=ix+1
          if(jx.eq.5) jx=1
          pout(jx,k)=pin(ix,kp)
        enddo
      enddo
c
c
      do k=1,npart
         if(abs(flout(k)).eq.24) flout(k)= -flout(k)
      enddo
      end
C****************************************************************************
      subroutine selcol(flvmlm,posi,rnd1,rnd2,rnd3,allowed)
c     select colours
C****************************************************************************
C
      implicit none
C     
      integer nmax
      parameter (nmax=20)       !maximum number of external particles. 
      integer flvmlm(nmax),posi(2)
      integer color(2*nmax)     !coulor string
      common/colore/color
      integer nlb
      parameter (nlb=4)            
      integer nrep
      parameter (nrep=100)
      integer allowed           !0 if coulor flow possible 1 otherwise
      integer j1,m1,m2,count,col(nmax),cola(nmax)
      integer gell(2,8),spin
      integer spin1(8),spin2(8),nqq
      integer rep(nmax),nqrk,nprt,nglu,nlep,ngb,nphot,nh
      real*8 rnd1,rnd2,rnd3          ! random numbers to select coulor
      real*8 n1,n2,qq(3,3)
      common/process/rep,nqrk,nprt,nglu,nlep,ngb,nphot,nh
      data  gell/1,1,1,2,1,3,2,1,2,2,2,3,3,1,3,2/
      data  spin1/0,1,2,-1,0,1,-2,-1/,spin2/0,-1,0,1,0,1,0,-1/
      data qq/1,4,7,2,5,8,3,6,1/
c     data repa/-1,0,-1,1,0,1,2,2,2,2,   !d nu b ub e- bb g g g g   incoming
c     >         -1,0,-1,-1,1,0,1,1,2,2,  !d nu c b ub e- cb bb g g  incoming
c     >         -1,-1,0,-1,1,1,0,1,2,2   !d d nu b ub db e- bb g g  incoming
c     >             ,970*10  /
C
         integer flvax(nmax)
         common/flv/flvax
C
         do j1=1,nmax
           flvax(j1)=flvmlm(j1)
         enddo      
C
      call flvtocol(flvmlm,posi)
C     
      n1=3.d0**nqrk
      nqq=nqrk/2
      if (nqrk.eq.0.) then
      elseif (nqrk.eq.2.) then
      elseif (nqrk.eq.4) then
      elseif (nqrk.eq.6) then
      else
         write(6,*)'wrong flavour assignment',flvmlm
         stop
      endif
C     
      count=0
      m1=int(rnd1*n1)+1
      do j1=1,nqrk
         count=count+1
         cola(count)=mod(m1,3)+1
         m1=m1/3
      enddo
      do j1=1,nqq
         col(j1)=qq(cola(j1),cola(j1+nqq))
      enddo
      if(nglu.le.8) then
        n2=8.d0**nglu
        m2=int(rnd2*n2)+1
        do j1=nqrk+1,nqrk+nglu
          count=count+1
          cola(count)=mod(m2,8)+1
          m2=m2/8
        enddo
        do j1=nqq+1,nqq+nglu
          col(j1)=cola(nqq+j1)
        enddo
      else
        if(nglu.gt.16) then
          write(*,*)'too many gluons in selcol'
          stop
        endif
        n2=8.d0**8
        m2=int(rnd2*n2)+1
        do j1=nqrk+1,nqrk+8
          count=count+1
          cola(count)=mod(m2,8)+1
          m2=m2/8
        enddo
        do j1=nqq+1,nqq+8
          col(j1)=cola(nqq+j1)
        enddo
        n2=8.d0**(nglu-8)
        m2=int(rnd3*n2)+1
        do j1=nqrk+9,nqrk+nglu
          count=count+1
          cola(count)=mod(m2,8)+1
          m2=m2/8
        enddo
        do j1=nqq+9,nqq+nglu
          col(j1)=cola(nqq+j1)
        enddo
      endif
C     
      allowed=0
      spin=0
      do j1=1,nqq+nglu
         spin=spin+spin1(col(j1))
      enddo
      if(spin.ne.0) then
         allowed=1
         return
      endif
      spin=0
      do j1=1,nqq+nglu
         spin=spin+spin2(col(j1))
      enddo
      if(spin.ne.0) then
         allowed=1
         return
      endif
C     
      count=0
      do j1=1,nprt
         if (rep(j1).eq.0) then
            color(2*j1-1)=0
            color(2*j1)=0
         elseif (rep(j1).eq.-1) then
            count=count+1
            color(2*j1-1)=cola(count)
            color(2*j1)=0
         elseif (rep(j1).eq.1) then
            count=count+1
            color(2*j1-1)=0
            color(2*j1)=cola(count)
         elseif (rep(j1).eq.2) then
            count=count+1
            color(2*j1-1)=gell(1,cola(count))
            color(2*j1)=gell(2,cola(count))
         else
            write(6,*)'wrong coulor representation in INCOLBBW'
            stop
         endif
      enddo
C     
      return
      end                                                                 
C
C**************************************************************************
         subroutine selspin(rspin,hel)
C**************************************************************************
C
C Input parameters:   RSPIN random number to select spin configuration
C Output parameters:  HEL helicity configuration 
C
         implicit none
C
         integer nmax
         parameter (nmax=20)
C
         real*8 rspin
         integer hel(nmax)
C
         integer rep(nmax),nqrk,nprt,nglu,nlep,ngb,nphot,nh
         common/process/rep,nqrk,nprt,nglu,nlep,ngb,nphot,nh
         integer j1,j2
C
         nh=nprt-nglu-nqrk-ngb-nlep-nphot
C
         j2=rspin*2**(nglu+nqrk+nphot+nlep)*3**ngb
         do j1=1,nqrk+nlep
           hel(j1)=mod(j2,2)
           if(hel(j1).eq.0)hel(j1)=-1
           j2=j2/2       
         enddo
         do j1=nqrk+nlep+1,nqrk+nlep+nh
            hel(j1)=0
         enddo
         do j1=nqrk+nlep+nh+1,nqrk+nlep+nh+ngb
           hel(j1)=mod(j2,3)-1
           j2=j2/3       
         enddo
         do j1=nqrk+nlep+nh+ngb+1,nqrk+nlep+nh+ngb+nglu+nphot
            hel(j1)=mod(j2,2)
            if(hel(j1).eq.0)hel(j1)=-1
            j2=j2/2       
         enddo
C
         return
         end
C*********************************************************************
        subroutine fltspn(flvmlm,hel,posi,allowed)
C*********************************************************************
C
        implicit none
C
        integer nmax
        parameter (nmax=20)
        integer flvmlm(nmax), hel(nmax),allowed,posi(2)
C
        integer j1,hs(-12:12),j2,nf(-12:12),flv(nmax)
C
        do j1=1,nmax
         flv(j1)=flvmlm(j1)
        enddo
        flv(posi(1))=-flv(posi(1))
        flv(posi(2))=-flv(posi(2))
C
        do j1=-12,12
         hs(j1)=0
         nf(j1)=0
        enddo
        j1=1
        do j1=1,nmax
           j2=flv(j1)
           if(j2.eq.1001) goto 111
           if(abs(j2).le.12) then
              hs(j2)=hs(j2)+hel(j1)
              nf(j2)=nf(j2)+1
           endif
           
        enddo
C
 111    allowed=1
        j2=abs(abs(nf(1)-nf(-1))+abs(nf(2)-nf(-2)))
        if(mod(j2,2).ne.0) then
         write(*,*)'wrong flavour assignment in fltspn'
         stop
        endif
        j2=j2/2
        if(nf(1).gt.nf(-1)) then
         if(hs(1)-hs(-1).ne.j2) return
        else
         if(hs(-1)-hs(1).ne.j2) return
        endif
        if(nf(2).gt.nf(-2)) then
         if(hs(2)-hs(-2).ne.j2) return
        else
         if(hs(-2)-hs(2).ne.j2) return
        endif
c
c        j2=abs(abs(nf(3)-nf(-3))+abs(nf(4)-nf(-4)))
c        if(mod(j2,2).ne.0) then
c         write(*,*)'wrong flavour assignment in fltspn'
c         stop
c        endif
c        j2=j2/2
c        if(nf(3).gt.nf(-3)) then
c         if(hs(3)-hs(-3).ne.j2) return
c        else
c         if(hs(-3)-hs(3).ne.j2) return
c        endif
c        if(nf(4).gt.nf(-4)) then
c         if(hs(4)-hs(-4).ne.j2) return
c        else
c         if(hs(-4)-hs(4).ne.j2) return
c        endif
c
        j2=abs(abs(nf(7)-nf(-7))+abs(nf(8)-nf(-8)))
        if(mod(j2,2).ne.0) then
         write(*,*)'wrong flavour assignment in fltspn'
         stop
        endif
        j2=j2/2
        if(nf(7).gt.nf(-7)) then
         if(hs(7)-hs(-7).ne.j2) return
        else
         if(hs(-7)-hs(7).ne.j2) return
        endif
        if(nf(8).gt.nf(-8)) then
         if(hs(8)-hs(-8).ne.j2) return
        else
         if(hs(-8)-hs(8).ne.j2) return
        endif
c
        j2=abs(abs(nf(9)-nf(-9))+abs(nf(10)-nf(-10)))
        if(mod(j2,2).ne.0) then
         write(*,*)'wrong flavour assignment in fltspn'
         stop
        endif
        j2=j2/2
        if(nf(9).gt.nf(-9)) then
         if(hs(9)-hs(-9).ne.j2) return
        else
         if(hs(-9)-hs(9).ne.j2) return
        endif
        if(nf(10).gt.nf(-10)) then
         if(hs(10)-hs(-10).ne.j2) return
        else
         if(hs(-10)-hs(10).ne.j2) return
        endif
c
        j2=abs(abs(nf(11)-nf(-11))+abs(nf(12)-nf(-12)))
        if(mod(j2,2).ne.0) then
         write(*,*)'wrong flavour assignment in fltspn'
         stop
        endif
        j2=j2/2
        if(nf(12).ne.hs(12)) return
        if(nf(-12).ne.hs(-12)) return
        if(nf(11).gt.nf(-11)) then
         if(hs(11)-hs(-11).ne.j2) return
        else
         if(hs(-11)-hs(11).ne.j2) return
        endif
C
        allowed=0
C
        return
        end
c-------------------------------------------------------------------
      subroutine setcol(npart,ip,ipinv,colfl,ifl,icu)
c     input colour flow from alpha (colfl) , output colour flow string
c     for herwig (icu)
c-------------------------------------------------------------
      implicit none
c inputs
      integer maxpar
      parameter(maxpar=20)
      integer npart,ip(maxpar),ipinv(maxpar),colfl(2*maxpar),icu(2
     $     ,maxpar),ifl(maxpar)
c locals
      integer nfree,nfree0,iswap,iend(maxpar),ist(maxpar)
      integer i,j,k,kb,it,nev, idbg,iq(maxpar), nq,ng
c
      data nev/0/,idbg/0/
c
c      nev=nev+1
c      if(idbg.eq.1.and.nev.gt.30) stop
      nq=0
      ng=0
      do i=1,maxpar
         icu(1,i)=0
         icu(2,i)=0
         iend(i)=1
         iq(i)=0
         if(i.le.2) then
c initial state
            ist(i)=1
         else
c final state
            ist(i)=-1
         endif
      enddo
      do i=1,npart
         icu(1,i)=colfl(2*ipinv(i)-1)
         icu(2,i)=colfl(2*ipinv(i))
      enddo
      nfree=npart
      do i=1,npart
        if(ifl(i).lt.0.and.ifl(i).gt.-7) then         
c     antiquark
          iq(i)=-1
          nq=nq+1
          iend(i)=0
          nfree=nfree-1
          if(icu(2,i).eq.0) then
            icu(2,i)=icu(1,i)
            icu(1,i)=0
          endif
        elseif(ifl(i).gt.0.and.ifl(i).lt.7) then         
c     quark
          iq(i)=1
          nq=nq+1
          iend(i)=0
          nfree=nfree-1
          if(icu(1,i).eq.0) then
            icu(1,i)=icu(2,i)
            icu(2,i)=0
          endif
        elseif(ifl(i).eq.21) then
          ng=ng+1
        endif
      enddo
c
c special for gluon-only events, flip colours of initial state gluons
c      if(nq.eq.0) then
c        do 20 i=1,2
c          if(ifl(i).eq.21) then
c            k=icu(2,i)
c            icu(2,i)=icu(1,i)
c            icu(1,i)=k
c          endif
c 20     enddo
c        return
c      endif
c     

c
c flip colours for initial state gluons
      do 20 i=1,2
        if(ifl(i).eq.21) then
          k=icu(2,i)
          icu(2,i)=icu(1,i)
          icu(1,i)=k
        endif
 20   enddo
c wrap it up if no quarks
      if(nq.eq.0) return


      if(idbg.eq.1) then
         write(6,*) 'Quark status: nfree= ',nfree
         do i=1,npart
            write(6,*) 'ifl=',ifl(i),' colfl(1)=',colfl(2*ipinv(i)-1)
     $           ,' colfl(2)=',colfl(2*ipinv(i)),' icup(1)=',icu(1,i)
     $           ,' icup(2)=',icu(2,i),' iend=',iend(i)
         enddo
      endif
c     
c
c swap colour-anticolour for gluons whose colour line ends on antiquarks
      do i=1,npart
         if(abs(iq(i)).eq.1) then
            k=icu(1,i)
            kb=icu(2,i)
c     found a quark/antiquark
            j=0
 100        j=j+1
            if(j.gt.npart) goto 120
            if(ifl(j).ne.21.or.iend(j).eq.0) goto 100
            if(k.ne.icu(1,j).and.k.ne.icu(2,j).and.kb.ne.icu(1,j).and.kb
     $           .ne.icu(2,j)) goto 100
            iend(j)=0
            nfree=nfree-1
            iswap=0
c     colour line between init-init or final-final states: swap the gluon
            if(k.eq.icu(1,j).and.ist(i)*ist(j).gt.0) iswap=1
c     colour-anticol line between init-final states: swap the gluon
            if(k.eq.icu(2,j).and.ist(i)*ist(j).lt.0) iswap=1
c     anticol line between init-init or final-final states: swap the gluon
            if(kb.eq.icu(2,j).and.ist(i)*ist(j).gt.0) iswap=1
c     colour-anticol line between init-final states: swap the gluon
            if(kb.eq.icu(1,j).and.ist(i)*ist(j).lt.0) iswap=1
            if(iswap.eq.1) then
               k=icu(2,j)
               icu(2,j)=icu(1,j)
               icu(1,j)=k
            endif
         endif
 120  continue
      enddo
c
c     go through the remaining gluons
      it=0
 180  continue
      nfree0=nfree
      if(idbg.eq.1) then
         it=it+1
         write(6,*) 'Gluon status, IT=',it,' nfree= ',nfree
         do i=1,npart
            write(6,*) 'ifl=',ifl(i),' colfl(1)=',colfl(2*ipinv(i)-1)
     $           ,' colfl(2)=',colfl(2*ipinv(i)),' icup(1)=',icu(1,i)
     $           ,' icup(2)=',icu(2,i),' iend=',iend(i)
         enddo
      endif
      do i=1,npart
         if(ifl(i).eq.21.and.iend(i).eq.0) then
            k=icu(1,i)
            kb=icu(2,i)
            j=0
 200        j=j+1
            if(j.gt.npart) goto 220
            if(ifl(j).ne.21.or.iend(j).eq.0) goto 200
            if(k.ne.icu(1,j).and.k.ne.icu(2,j).and.kb.ne.icu(1,j).and.kb
     $           .ne.icu(2,j)) goto 200
            iend(j)=0
            nfree=nfree-1
            iswap=0
            if(k.eq.icu(1,j).and.ist(i)*ist(j).gt.0) iswap=1
            if(k.eq.icu(2,j).and.ist(i)*ist(j).lt.0) iswap=1
            if(kb.eq.icu(2,j).and.ist(i)*ist(j).gt.0) iswap=1
            if(kb.eq.icu(1,j).and.ist(i)*ist(j).lt.0) iswap=1
            if(iswap.eq.1) then
               k=icu(2,j)
               icu(2,j)=icu(1,j)
               icu(1,j)=k
            endif
         endif
 220  continue
      enddo
      if(nfree.ne.nfree0) goto 180
c check for close colour clusters

c     
c
      if(idbg.eq.1) then
         write(6,*) 'Final, IT=',it,' nfree= ',nfree
         do i=1,npart
            write(6,*) 'ifl=',ifl(i),' colfl(1)=',colfl(2*ipinv(i)-1)
     $           ,' colfl(2)=',colfl(2*ipinv(i)),' icup(1)=',icu(1,i)
     $           ,' icup(2)=',icu(2,i),' iend=',iend(i)
         enddo
         write(6,*) ' '
      endif
      end

      subroutine colchk(npart,icu,ifl,nev,w)
      implicit none
c arguments
      integer nmax
      parameter (nmax=20)
      integer npart,icu(2,nmax),ifl(nmax)
      double precision nev,w
c locals
      integer i,j,icount(0:nmax),itmp
c
      do i=1,nmax
        icount(i)=0
      enddo
      do i=1,2
        do j=1,npart
          itmp=icu(i,j)
          icount(itmp)=icount(itmp)+1
        enddo
      enddo
      itmp=0
      do i=1,nmax
        if(icount(i).ne.2.and.icount(i).ne.0) itmp=1
      enddo
      if(itmp.ne.0) then
        write(6,*) 'Problem with colour flow, event=',nev,' ME=',w
        write(6,*) (ifl(j),j=1,npart)
        write(6,*) (icu(1,j),j=1,npart)
        write(6,*) (icu(2,j),j=1,npart)
      endif
      end

c-------------------------------------------------------------------
      subroutine evdump(Nunit,npart,iproc,ifl,icolup,p,qsq,wgt)
c     dump event details to a file, for future reading by herwig
c-------------------------------------------------------------------
      implicit none
      integer maxdec
      parameter (maxdec=40)
c inputs
      integer Nunit,npart,iproc,ifl(maxdec),icolup(2,maxdec)
      real*8 qsq,p(5,maxdec),wgt
c locals
      real Sq,Sp(3,maxdec),Sm(maxdec),Swgt
      integer nev,i,j
      data nev/0/
      save
c
      if(npart.gt.maxdec) then
         write(6,*) 'Sp(4,maxdec) in evdump underdimensioned'
         stop
      endif
      nev=nev+1
      Sq=real(sqrt(qsq))
      Swgt=real(wgt)
      do i=1,npart
         do j=1,3
            Sp(j,i)=real(p(j,i))
         enddo
         Sm(i)=real(p(5,i))
      enddo
c     Nevent, iproc, npart, wgt, Q scale 
      write(Nunit,2) nev,iproc,npart,Swgt,Sq
 2    format(i8,1x,i4,1x,i2,2(1x,e12.6))
c     flavour, colour and z-momentum of incoming partons
      write(Nunit,8) ifl(1),icolup(1,1),icolup(2,1),Sp(3,1)
      write(Nunit,8) ifl(2),icolup(1,2),icolup(2,2),Sp(3,2)
c     flavour, colour, 3-momentum and mass of outgoing partons
      do i=3,npart
      write(Nunit,9) ifl(i),icolup(1,i),icolup(2,i),Sp(1,i),Sp(2,i)
     $     ,Sp(3,i),Sm(i)
      enddo
 8    format(i8,1x,2(i4,1x),f10.3)
 9    format(i8,1x,2(i4,1x),4(1x,f10.3))
      end


      SUBROUTINE LHEHEA(NIOSTA,ebeam,ih2,avgwgt,sigerr)
      IMPLICIT NONE
      double precision ebeam,avgwgt,sigerr
      integer NIOSTA,ih2,iproc,ipr
      INTEGER MAXPUP
      PARAMETER(MAXPUP=100)
      INTEGER IDBMUP(2),PDFGUP(2),PDFSUP(2),IDWTUP,NPRUP,LPRUP(MAXPUP)
      DOUBLE PRECISION EBMUP(2),XSECUP(MAXPUP),XERRUP(MAXPUP)
     $     ,XMAXUP(MAXPUP)
C WRITE THE HEADERS OF THE LHE EVENT OUTPUT
      WRITE(NIOSTA,'(A)') '-->'
C
      WRITE(NIOSTA,'(A)') '<header>'
      WRITE(NIOSTA,'(A)') '</header>'
C
      WRITE(NIOSTA,'(A)') '<init>'
C
      IDBMUP(1) = 2212
      IF(IH2.EQ.1) THEN
        IDBMUP(2) = 2212
      ELSEIF(IH2.EQ.-1) THEN
        IDBMUP(2) =-2212
      ENDIF
      EBMUP(1) = ABS(EBEAM)
      EBMUP(2) = ABS(EBEAM)
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
      XERRUP(1) = sigerr
C--   MAXIMUM WEIGHT
      XMAXUP(1) = avgwgt
C-- THE FOLLOWING IS HERWIG/PYTHIA DEPENDENT, AND NEEDS TO BE SETUP IN ALPSHO
      LPRUP(1)=1111
C
      WRITE(NIOSTA,200) IDBMUP(1),IDBMUP(2),EBMUP(1),
     &EBMUP(2),PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),IDWTUP,NPRUP
      DO 120 IPR=1,NPRUP
        WRITE(NIOSTA,201) XSECUP(IPR),XERRUP(IPR),XMAXUP(IPR),LPRUP(IPR)
  120 CONTINUE
      WRITE(NIOSTA,'(A)') '</init>'
 200  FORMAT(2(I5,1X),2(E12.6,1X),4(I6,1X),2(I2,1X))
 201  FORMAT(3(E12.6,1X),I8)
      END

c-------------------------------------------------------------------
      subroutine evdlhe
c     dump event details directly in LHE format
c-------------------------------------------------------------------
      implicit none
      INTEGER I,J,IPR,NEV,iproc
      DATA NEV/0/
      DOUBLE PRECISION EWGT
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &     IDUP(MAXNUP),ISTUP(MAXNUP),MOTHUP(2,MAXNUP),
     &     ICOLUP(2,MAXNUP),PUP(5,MAXNUP),VTIMUP(MAXNUP),
     &     SPINUP(MAXNUP)
      CHARACTER*2 VTI,SPI
      INCLUDE 'alpgen.inc'
      save
c
      NEV=NEV+1
      EWGT=1D0
ccc modified fulvio 09/01/2014 - begin
c      IPROC=1
      IDPRUP=1111
ccc modified fulvio 09/01/2014 - end
      AQEDUP=-1
      AQCDUP=-1
ccc modified fulvio 19/05/2014 - begin
c      VTI='0'
c      SPI='9'
ccc modified fulvio 19/05/2014 - end
      WRITE(NIOUNW,*) "<event>" 
ccc modified fulvio 09/01/2014 - begin
c      WRITE(NIOUNW,1) NUP,IPROC,EWGT,SQRT(QSQ),AQEDUP,AQCDUP
      WRITE(NIOUNW,1) NUP,IDPRUP,EWGT,SQRT(QSQ),AQEDUP,AQCDUP
ccc modified fulvio 09/01/2014 - end
 1    FORMAT(I2,1X,I4,1X,4(E12.6,1X))
      DO 100 I=1,NUP
ccc modified fulvio 19/05/2014 - begin
        VTIMUP(I)=0.D0
c        SPINUP(I)=SPNWRT(I)
        IF(ISTUP(I).eq.2) SPINUP(I)=9.D0
ccc modified fulvio 19/05/2014 - end
        WRITE(NIOUNW,101) 
     &   IDUP(I),ISTUP(I),
     &       MOTHUP(1,I),MOTHUP(2,I),ICOLUP(1,I),ICOLUP(2,I),
ccc modified fulvio 19/05/2014 - begin
c     &       (PUP(J,I),J=1,5),VTI,SPINUP(I)
     &       (PUP(J,I),J=1,5),VTIMUP(I),SPINUP(I)
ccc modified fulvio 19/05/2014 - end
 100  CONTINUE
      WRITE(NIOUNW,*) "</event>" 
ccc modified fulvio 19/01/2014 - begin
cccc modified fulvio 09/01/2014 - begin
 101  FORMAT(1P,I8,5(1X,I5),5(1X,E16.9),1X,E12.5,1X,E10.3)
c 101  FORMAT(I8,1X,(I4,1X),4(I3,1X),7(E16.9,1X))
c 101  FORMAT(I8,1X,(I4,1X),4(I3,1X),5(E16.9,1X),A,1x,E10.3)
cc 101  FORMAT(I8,1X,(I4,1X),4(I3,1X),5(F10.3,1X),A,1x,F10.3)
cccc modified fulvio 09/01/2014 - end
ccc modified fulvio 19/01/2014 - end
C 101  FORMAT(I8,1X,(I4,1X),4(I3,1X),5(1PE11.4,1X),2(A,1x))
      end

C***********************************************************************
       subroutine matrix(flvmlm,posi,impul,hel,rst,labcol
     >                    ,wgcol,colfl)
C***********************************************************************
C
C Input parameters: FLVMLM(NMAX) particles according
C                   to michelangelo convention, POSI(NMAX) positions of
C                   incoming particles IMPUL(4,8) particles four momenta,
C                   HEL, array of particle helicities 
C                   LABCOL to choose whether su(3) mode only (LABCOL=0) or
C                   color flow is required (LABCOL=1) or dual mode only
C                   (LABCOL=2), WGCOL random number for color flow unweighting
C Output parameters:   RESULT squared matrix element (modulus), COLFLOW 
C                       string representing the selected coulor flow 
C
C
C At present:
C
C     the color flow is returned in the following way: each particle is
C     represented by a pair of integers (n1,n2); two particle (n1,n2), (n3,n4)
C     are colour connected if either n2=n3 or n4=n1 (in this case the coulor
C     flows from 2 to 1). COLFLOW will contain a string of pairs of integers 
C     ordered as follows: d (or dbar) nu_ebar bbar ubar (or u) e^- b glu glu
C     (with gluons ordered according to the ordering of momenta).
C     (COLFLOW=(n1,....,nn,....) only the first nn=2*particles_number
C       elements to be used)   
C
C
C   IMPUL(J,NPART) ===   J=1 Energy of the NPART-th particle
C   IMPUL(J,NPART) ===   J=2,3,4 x,y,z components  of the NPART-th particle 
C                                                         three momentum
C    the order of the momenta is: 
C        dbar nubar bbar u e- b ( glu glu glu)  (processes 1,2,3,4)
C        dbar nubar cbar bbar u e- c b (glu)  (processes > 5 )
C    where all particles are assumed outcoming and incoming gluons
C    must be before outcoming ones
C
         implicit none
C
         integer nmax        !maximum number of external particles, 
         parameter (nmax=20)
         integer nlb
         parameter (nlb=4)    
         real*8 impul(4,nmax),wgcol    !the contribution to the integral
         real*8 rst
c         complex*16 rst
         integer labcol,colfl(2*nmax),flvmlm(nmax)
c
         integer posi(2)
         integer hel(nmax)             !particles helicities
C
         integer flgdual          !dual (0) or su3 (1) amplitudes
         common/dual/flgdual
         integer color(2*nmax)    !coulor string
         integer colst(2*nmax)
         common/colore/color
         integer ncls,ant4(4,2),ant6a(6,6),ant6b(6,6),ant2(2,1)
         integer nant,nantl,j3
         parameter (ncls=2)       !different class of processes
         integer nx,proc(nmax),antns(6),ndual,j1,j2
         parameter (ndual=4040000)
         real*8 dualamp(ndual),colfac,damp,dampref,avgspin
         integer colaux(ndual,2*nmax),ndl,ndla(0:12,0:3),class
         integer rep(nmax),nqrk,nprt,nglu,nlep,ngb,nphot,nh
         common/process/rep,nqrk,nprt,nglu,nlep,ngb,nphot,nh
         integer fq(nmax),pq(nmax)
         data ndla/0,0,0,1,5,23,119,719,5039,40319,362879,3628799,
     >                                                    39916799,
     >          0,0,1,5,23,119,719,5039,40319,362879,3628799,0,0,
     >          0,1,5,23,119,719,5039,40319,362879,0,0,0,0,
     >          0,2,11,59,359,8*0/         
c         data ndla/0,0,0,1,5,23,119,719,5039,40319,362879,
c     >          0,0,1,5,23,119,719,5039,40319,0,0,
c     >          0,1,5,23,119,719,5039,0,0,0,0,
c     >          0,2,11,59,359,6*0/

         save ndla,ant4,ant6a,ant2
         complex*16 result
C
         data ant2  /1,2/
         data ant4  /1,4,2,3,         ! 2Q 2Qph, 2t2b, hjet, Njet, QQh, phjet, top,vbjet old
     >               1,3,2,4/
c         data ant4  /1,3,2,4,
c     >               1,4,2,3/          ! 4t 4b old
c         data ant6a  /1,5,2,4,3,6,     ! 4t 4b old
c     >               1,5,2,6,3,4,
c     >               1,6,2,4,3,5,
c     >               1,6,2,5,3,4,
c     >               1,4,2,5,3,6,
c     >               1,4,2,6,3,5/    
         data ant6a  /1,5,2,6,3,4,     ! 2Q 2Qph, 2t2b, hjet, Njet, phjet, QQh, top,vbjet  old
     >               1,6,2,4,3,5,
     >               1,6,2,5,3,4,
     >               1,4,2,6,3,5,
     >               1,4,2,5,3,6,
     >               1,5,2,4,3,6/    
C
         call setinteraction77
C
         avgspin=0.d0
C
         if (labcol.eq.0) then
C
          flgdual=1
          call matrix0(flvmlm,posi,impul,hel,result)
c
         elseif (labcol.eq.1) then
c
c     
          j1=0
          do j3=1,nmax
           if (abs(flvmlm(j3)).gt.0.and.abs(flvmlm(j3)).le.6) then
            j1=j1+1
            fq(j1)=flvmlm(j3)
            pq(j1)=j3
           endif
          enddo
          if(j1.ne.nqrk) then
           write(*,*)'something wrong in matrix, NQRK not ok',nqrk,j3
           stop
          endif
c
          flgdual=0
          do j1=1,2*nmax
           colst(j1)=color(j1)
          enddo 
c
          do j1=1,nmax
           proc(j1)=abs(rep(j1))
          enddo
c
          ndl=0
          if(nqrk.eq.0) then
           nantl=1
          elseif(nqrk.eq.2) then
           nantl=1
          elseif(nqrk.eq.4) then
           nantl=2
          elseif(nqrk.eq.6) then
           nantl=6
          else
           write(*,*)'Cannot deal with',nqrk,'quarks'
           stop
          endif
          do j3=1,nantl
           do j1=1,nqrk
            if(nqrk.eq.2) then
             antns(j1)=pq(ant2(j1,j3))
            elseif(nqrk.eq.4) then
             antns(j1)=pq(ant4(j1,j3))
            elseif(nqrk.eq.6) then
             antns(j1)=pq(ant6a(j1,j3))
            endif
           enddo
           j2=nqrk/2
           do j1=0,ndla(nglu,j2)
            nx=j1
            if(nqrk.eq.0) then
               call gendual0q(proc,nglu,colst,nx,colfac) 
            elseif(nqrk.eq.2) then
               call gendual2q(proc,antns,nglu,colst,nx,colfac) 
            elseif(nqrk.eq.4) then
               call gendual4q(proc,antns,nglu,colst,nx,colfac) 
            elseif(nqrk.eq.6) then
               call gendual6q(proc,antns,nglu,colst,nx,colfac) 
            endif 
            if(abs(colfac).gt.1.d-20) then
             ndl=ndl+1
             call matrix0(flvmlm,posi,impul,hel,result)
             dualamp(ndl)=abs(result*colfac)**2
             do j2=1,2*nmax
              colaux(ndl,j2)=color(j2)
             enddo
            endif
           enddo
          enddo
c
          damp=0
          do j1=1,ndl
           damp=damp+dualamp(j1)
          enddo
          dampref=damp*wgcol
          j1=0
          damp=0.
          do while(damp.lt.dampref.and.j1.lt.ndl)
           j1=j1+1
           damp=damp+dualamp(j1)
          enddo
          ndl=j1
          if(ndl.eq.0) then
             do j1=1,2*nmax
                colfl(j1)=colst(j1)
             enddo
          else 
             do j1=1,2*nmax
                colfl(j1)=colaux(ndl,j1)
             enddo
          endif
C
         else
          write(6,*)'wrong LABCOL in matrix',labcol
          stop
         endif
C
         rst=abs(result)**2
c         rst=result
C
         return
         end
c


      SUBROUTINE LHEHEA_old(NIOSTA,ebeam,ih2,avgwgt,sigerr)
      IMPLICIT NONE
      double precision ebeam,avgwgt,sigerr
      integer NIOSTA,ih2,iproc,ipr
      INTEGER MAXPUP
      PARAMETER(MAXPUP=100)
      INTEGER IDBMUP(2),PDFGUP(2),PDFSUP(2),IDWTUP,NPRUP,LPRUP(MAXPUP)
      DOUBLE PRECISION EBMUP(2),XSECUP(MAXPUP),XERRUP(MAXPUP)
     $     ,XMAXUP(MAXPUP)
C WRITE THE HEADERS OF THE LHE EVENT OUTPUT
      WRITE(NIOSTA,'(A)') '-->'
C
      WRITE(NIOSTA,'(A)') '<header>'
      WRITE(NIOSTA,'(A)') '</header>'
C
      WRITE(NIOSTA,'(A)') '<init>'
C
      IDBMUP(1) = 2212
      IF(IH2.EQ.1) THEN
        IDBMUP(2) = 2212
      ELSEIF(IH2.EQ.-1) THEN
        IDBMUP(2) =-2212
      ENDIF
      EBMUP(1) = ABS(EBEAM)
      EBMUP(2) = ABS(EBEAM)
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
      XERRUP(1) = sigerr
C--   MAXIMUM WEIGHT
      XMAXUP(1) = avgwgt
C-- THE FOLLOWING IS HERWIG/PYTHIA DEPENDENT, AND NEEDS TO BE SETUP IN ALPSHO
      LPRUP(1)=1111
C
      WRITE(NIOSTA,200) IDBMUP(1),IDBMUP(2),EBMUP(1),
     &EBMUP(2),PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),IDWTUP,NPRUP
      DO 120 IPR=1,NPRUP
        WRITE(NIOSTA,201) XSECUP(IPR),XERRUP(IPR),XMAXUP(IPR),LPRUP(IPR)
  120 CONTINUE
      WRITE(NIOSTA,'(A)') '</init>'
 200  FORMAT(2(I5,1X),2(E12.6,1X),4(I6,1X),2(I2,1X))
 201  FORMAT(3(E12.6,1X),I8)
      END

c-------------------------------------------------------------------
      subroutine evdlhe_old
c     dump event details directly in LHE format
c-------------------------------------------------------------------
      implicit none
      INTEGER I,J,IPR,NEV,iproc
      DATA NEV/0/
      DOUBLE PRECISION EWGT
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &     IDUP(MAXNUP),ISTUP(MAXNUP),MOTHUP(2,MAXNUP),
     &     ICOLUP(2,MAXNUP),PUP(5,MAXNUP),VTIMUP(MAXNUP),
     &     SPINUP(MAXNUP)
      CHARACTER*2 VTI,SPI
      INCLUDE 'alpgen.inc'
      save
c
      NEV=NEV+1
      EWGT=1D0
      IPROC=1
      AQEDUP=-1
      AQCDUP=-1
      VTI='0'
      SPI='9'
      WRITE(NIOUNW,*) "<event>" 
      WRITE(NIOUNW,1) NUP,IPROC,EWGT,SQRT(QSQ),AQEDUP,AQCDUP
 1    FORMAT(I2,1X,I4,1X,4(E12.6,1X))
      DO 100 I=1,NUP
        WRITE(NIOUNW,101) 
     &   IDUP(I),ISTUP(I),
     &       MOTHUP(1,I),MOTHUP(2,I),ICOLUP(1,I),ICOLUP(2,I),
     &       (PUP(J,I),J=1,5),VTI,SPI
 100  CONTINUE
      WRITE(NIOUNW,*) "</event>" 
 101  FORMAT(I8,1X,(I4,1X),4(I3,1X),5(F10.3,1X),2(F6.3,1x))
c 101  FORMAT(I8,1X,(I4,1X),4(I3,1X),5(F10.3,1X),2(A,1x))
      end

