c-------------------------------------------------------------------
      subroutine alsprc
c     assigns the hard process code
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      ihrd=9
      end

c-------------------------------------------------------------------
      subroutine alhset
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'Njet.inc'
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
c     ngrid is the total number of grids allowd for P.S. variables.
c     jgrid(jproc) labels the grid associated to a given jproc.
c       jgrid= 1   -->   g g
c       jgrid= 2   -->   q qbar and qbar q
c       jgrid= 3   -->   g q    and g qbar
c       jgrid= 4   -->   q g    and qbar g
c       jgrid= 5   -->   q q     and qbar qbar
      data jgrid/1,2,3,4,1,5,5,2,2,2,3,4,3,4,1,1,
     +           2,2,2,2,5,5,5,5,5,2,2,2,
     +           3,4,3,4,3,4,3,4,1,1,1,61*0/
      integer i
c     parameters for the gauge invariance prescription:
      winsize  = 2.d0/pi
      resonance= 'n'
      wmode    = 'yy'
c     process input parameters
      npart=njets+2
      nprtns=njets
      if(njets.eq.2) then
         jprocmax= 10
         ngrid   = 5
      elseif(njets.eq.3) then
         jprocmax= 14
         ngrid   = 5
      elseif(njets.eq.4) then
         jprocmax= 28
         ngrid   = 5
      elseif(njets.eq.5) then
         jprocmax= 36
         ngrid   = 5
      elseif(njets.eq.6) then
         jprocmax= 39
         ngrid   = 5
      elseif(njets.le.16) then
         jprocmax= 39
         ngrid   = 5
      else 
         print*,'njets=',njets,' not yet available'
         stop
      endif 
c masses
      do i=1,njets+2
         p(5,i)=0
      enddo
c run parameters:
      if(njets.gt.1) then
        write(niosta,*) njets,' jets'
      endif 
      write(niosta,*)
     $     '======================================='
      write(niosta,*)
     $     'Generation cuts for the partonic event sample:' 
      write(niosta,*) '     Light jets:'
      write(niosta,*) 'ptmin=',ptjmin,' |etamax|=',etajmax
     $     ,' dR(j-j)>',drjmin 
      end


c-------------------------------------------------------------------
      subroutine selflav(jproc,xlum,afl)
c     evaluates parton luminosities, and assigns PDG-code flavours
c     gluon=21 bbar=-5 cbar=-4 sbar=-3 ubar=-2 dbar=-1 d=1 u=2 s=3 c=4 b=5
c     jproc
c
c---n gluon final states
c
c  1  g g    -> ng
c  2  q qbar -> ng
c     qbar q -> ng
c
c---1 q jet + (n-1 gluons) final state
c
c  3  g q    -> q    (+ gluons)
c     g qbar -> qbar (+ gluons)
c  4  q g    -> q    (+ gluons)
c     qbar g -> qbar (+ gluons)
c
c---2 q jets + (n-2 gluons) final state
c
c  5  g g    -> q qbar   (+ gluons)
c  6  q q    -> q q      (+ gluons)
c  7  q q'   -> q q'     (+ gluons)
c  8  q qbar'-> q qbar'  (+ gluons)
c  9  q qbar -> q qbar   (+ gluons)
c 10  q qbar -> q' qbar' (+ gluons)
c
c---3 q jets + (n-3 gluons) final state
c
c 11  g q    -> q q qbar   (+ gluons)
c     qbar g -> qbar q qbar(+ gluons)
c 12  q g    -> q q qbar   (+ gluons)
c     qbar g -> qbar q qbar (+ gluons)
c 13  g q    -> q q' qbar'   (+ gluons)
c     qbar g -> qbar q' qbar'(+ gluons)
c 14  q g    -> q q' qbar'   (+ gluons)
c     qbar g -> qbar q' qbar' (+ gluons)
c
c---4 q jets + (n-4 gluons) final state
c
c 15  g g    -> q qbar q qbar (+gluons)
c 16  g g    -> q qbar q' qbar' (+gluons)
c
c --- processes with 6 quarks
c 17  q qbar -> q' qbar' q'' qbar'' 
c 18  q qbar -> q' qbar' q' qbar' 
c 19  q qbar -> q qbar q' qbar' 
c 20  q qbar -> q qbar q qbar
c 21  q q    -> q q q qbar 
c 22  q q    -> q q q' qbar'
c 23  q q'   -> q q' q qbar
c 24  q q'   -> q q' q' qbar'
c 25  q q'   -> q q' q'' qbar''
c 26  q qbar'-> q qbar' q qbar 
c 27  q qbar'-> q qbar' q' qbar' 
c 28  q qbar'-> q qbar' q'' qbar''
c
c---5 q jets + (n-5 gluons) final state
c
c 29  g q    -> q q qbar q qbar 
c 30  q g    -> q q qbar q qbar
c 31  g q    -> q q' qbar' q' qbar'
c 32  q g    -> q q' qbar' q' qbar'
c 33  g q    -> q q qbar q' qbar'
c 34  q g    -> q q qbar q' qbar'
c 35  g q    -> q q' qbar' q'' qbar''
c 36  q g    -> q q' qbar' q'' qbar''
c
c---6 q jets + (n-6 gluons) final state
c
c 37  g g    -> q qbar q qbar q qbar 
c 38  g g    -> q qbar q' qbar' q'' qbar''
c 39  g g    -> q qbar q qbar q' qbar'
c-----------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'Njet.inc'
      integer afl(maxpar)
      real tmp(1000),slum,cwgt,swgt,rn,tmptot
c
      double precision xlum,xrn                                 
      integer i,k,itmp,icount,init,jproc,ng,j
      real cfac(0:6),ifact(0:maxpar)
      data cfac/8e0,6*3e0/,init/0/
      integer imap(-5:5)
      data imap/-1,-2,-1,-2,-1,0,1,2,1,2,1/
      integer nlqp,kj
      double precision effco
      save init,cfac,ifact
c
      nlqp= 5
c
      if(init.eq.0) then
         njets=npart-2
         init=1
         ifact(0)=1e0
         do i=1,maxpar
            ifact(i)=ifact(i-1)/real(i)
         enddo
      endif
c
      do i=1,npart
         ifl(i)=0
      enddo 
 1    call randa(xrn)
      rn=real(xrn)
      if(1e0-rn.lt.1e-7) goto 1

      icount=0
      slum=0e0
      tmptot=0e0
c
c      print*,'jproc = ',jproc
c
c     jproc=1   njets>=2
c     g g -> (njets) g
      if(jproc.eq.1) then
         slum = f1(0)*f2(0)
         ifl(1) = 0
         ifl(2) = 0
         do k=3,njets+2
            ifl(k)=0
         enddo
         itmp= 1
         effco = ccoef(0,njets+2)
         goto 100
c
c     jproc=2   njets>=2
c  2  q qbar -> (njets) g
c     and qbar q -> (njets) g   q=u,d,c,s,b
      else if(jproc.eq.2) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(i)*f2(-i)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=i
                  ifl(2)=-i
                  do k=3,njets+2
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(1))
                  effco = ccoef(1,njets)
                  goto 100
               endif 
            endif
         enddo 
c
c     jproc=3   njets>=2
c  3  g q    -> q (njets-1) g 
c     and g qbar  -> qbar (njets-1) g
      else if(jproc.eq.3) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(0)*f2(i)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=0
                  ifl(2)=i
                  ifl(3)=0
                  ifl(4)=i
                  do k=5,njets+2
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(2))
                  effco = ccoef(1,njets)
                  goto 100
               endif 
            endif
         enddo 
c
c     jproc=4   njets>=2
c  4  q g    -> q (njets-1) g
c     and qbar g  -> qbar (njets-1) g
      else if(jproc.eq.4) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(i)*f2(0)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=i
                  ifl(2)=0
                  ifl(3)=0
                  ifl(4)=i
                  do k=5,njets+2
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(1))
                  effco = ccoef(1,njets)
                  goto 100
               endif 
            endif
         enddo 
c
c     jproc=5   njets>=2
c  5  g g -> q qbar (njets-2 g) q=u,d,c,s,b
      else if(jproc.eq.5) then
         slum=f1(0)*f2(0)
c
c  slum x number of possible f.s. q qbar pairs
c
         slum=slum*dble(nlqp)
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=1.
            endif
         enddo 
         icount=0
         rn=rn*2.*dble(nlqp)
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=0
                  ifl(2)=0
                  ifl(3)=i
                  ifl(4)=-i
                  do k=5,njets+2
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(3))
                  effco = ccoef(1,njets)
                  goto 100
               endif 
            endif
         enddo 
c     jproc=6   njets>=2
c  6  q q -> q q (njets-2) g
c     and qbar qbar -> qbar qbar (njets-2) g   q=u,d,c,s,b
      else if(jproc.eq.6) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(i)*f2(i)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=i
                  ifl(2)=i
                  ifl(3)=i
                  ifl(4)=i
                  do k=5,njets+2
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(1))
                  effco = ccoef(2,njets-2)
c statistical factor to account for indistinguishable quarks
                  slum= slum/2.d0
                  goto 100
               endif 
            endif
         enddo 
c     jproc=7   njets>=2
c  7  q q' -> q q' (njets-2) g
c     and qbar qbar' -> qbar qbar' (njets-2) g   q=u,d,c,s,b
      else if(jproc.eq.7) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(j)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=j
                        ifl(3)=i
                        ifl(4)=j
                        do k=5,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(1))
                        effco = ccoef(2,njets-2)
                        goto 100
                     endif
                  endif
               enddo 
            endif
         enddo 
c     jproc=8   njets>=2
c  8  q qbar' -> q qbar' (njets-2) g
c     and qbar q' -> qbar q' (njets-2) g   q=u,d,c,s,b
      else if(jproc.eq.8) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.lt.0) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(j)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.lt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=j
                        ifl(3)=i
                        ifl(4)=j
                        do k=5,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(1))
                        effco = ccoef(2,njets-2)
                        goto 100
                     endif
                  endif
               enddo 
            endif
         enddo 
c
c     jproc=9   njets>=2
c  9  q qbar -> q qbar (njets-2) g
c     and qbar q -> qbar q (njets-2) g   q=u,d,c,s,b
      else if(jproc.eq.9) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(i)*f2(-i)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=i
                  ifl(2)=-i
                  ifl(3)=i
                  ifl(4)=-i
                  do k=5,njets+2
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(1))
                  effco = ccoef(2,njets-2)
                  goto 100
               endif
            endif
         enddo 
c
c     jproc=10   njets>=2
c 10  q qbar -> q' qbar' (njets-2) g
c     and qbar q -> qbar' q' (njets-2) g   q=u,d,c,s,b
      else if(jproc.eq.10) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j= -nlqp,nlqp
                  if(j.ne.0.and.abs(j).ne.abs(i)) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(-i)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j= -nlqp,nlqp
                  if(j.ne.0.and.abs(j).ne.abs(i)) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=-i
                        ifl(3)=j
                        ifl(4)=-j
                        do k=5,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(1))
                        effco = ccoef(2,njets-2)
c to remove double counting (j <-> -j)
                        slum= slum/2.d0
                        goto 100
                     endif
                  endif
               enddo
            endif
         enddo 
c
c     jproc=11   njets>=3
c 11  g q    -> q q qbar (njets-3) g 
c     and g qbar  -> qbar q qbar (njets-3) g
      else if(jproc.eq.11) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(0)*f2(i)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=0
                  ifl(2)=i
                  ifl(3)=i
                  ifl(4)=i
                  ifl(5)=-i
                  do k=6,njets+2
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(2))
                  effco = ccoef(2,njets-2)
c statistical factor to account for indistinguishable quarks
                  slum= slum/2.d0
                  goto 100
               endif 
            endif
         enddo 
c
c     jproc=12   njets>=3
c 12  q g    -> q q qbar (njets-3) g
c     and qbar g -> qbar q qbar (njets-3) g
      else if(jproc.eq.12) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(i)*f2(0)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=i
                  ifl(2)=0
                  ifl(3)=i
                  ifl(4)=i
                  ifl(5)=-i
                  do k=6,njets+2
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(1))
                  effco = ccoef(2,njets-2)
c statistical factor to account for indistinguishable quarks
                  slum= slum/2.d0
                  goto 100
               endif 
            endif
         enddo 
c
c     jproc=13   njets>=3
c 13  g q    -> q q' qbar' (njets-3) g 
c     and g qbar'  -> qbar q' qbar' (njets-3) g
      else if(jproc.eq.13) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.(abs(i).ne.abs(j))) then
                     icount=icount+1
                     tmp(icount)=f1(0)*f2(i)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.(abs(i).ne.abs(j))) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=0
                        ifl(2)=i
                        ifl(3)=i
                        ifl(4)=j
                        ifl(5)=-j
                        do k=6,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(2))
                        effco = ccoef(2,njets-2)
c to remove double counting (j <-> -j)
                        slum= slum/2.d0
                        goto 100
                     endif
                  endif
               enddo
            endif
         enddo 
c
c     jproc=14   njets>=3
c 14  q g    -> q q' qbar' (njets-3) g
c     and qbar g -> qbar q' qbar' (njets-3) g
      else if(jproc.eq.14) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.(abs(i).ne.abs(j))) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(0)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.(abs(i).ne.abs(j))) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=0
                        ifl(3)=i
                        ifl(4)=j
                        ifl(5)=-j
                        do k=6,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(1))
                        effco = ccoef(2,njets-2)
c to remove double counting (j <-> -j)
                        slum= slum/2.d0
                        goto 100
                     endif
                  endif
               enddo
            endif
         enddo 
c
c     jproc=15   njets>=4
c  15  g g -> q qbar q qbar (njets-4 g) q=u,d,c,s,b
      else if(jproc.eq.15) then
         slum=f1(0)*f2(0)
c
c  slum x number of possible f.s. q qbar pairs
c
         slum=slum*dble(nlqp)
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=1.
            endif
         enddo 
         icount=0
         rn=rn*nlqp*2.d0       !to count plus and minus signs of q
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=0
                  ifl(2)=0
                  ifl(3)=i
                  ifl(4)=-i
                  ifl(5)=i
                  ifl(6)=-i
                  do k=7,njets+2
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(3))
                  effco = ccoef(2,njets-2)
c     statistical factor to account for indistinguishable quarks
                  slum= slum/4.d0
                  goto 100
               endif 
            endif
         enddo 
c
c     jproc=16   njets>=4
c  16  g g -> q qbar q' qbar' (njets-4 g) q=u,d,c,s,b
      else if(jproc.eq.16) then
         slum=f1(0)*f2(0)
c
c  slum x number of possible f.s. q qbar pairs
c
         slum=slum*dble(nlqp)*dble(nlqp-1)/2.d0
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.(abs(i).ne.abs(j))) then
                     icount=icount+1
                     tmp(icount)=1.
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*dble(nlqp)*dble(nlqp-1) * 4.d0  
c     to count plus and minus signs of q and q'
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.(abs(i).ne.abs(j))) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=0
                        ifl(2)=0
                        ifl(3)=i
                        ifl(4)=-i
                        ifl(5)=j
                        ifl(6)=-j
                        do k=7,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(3))
                        effco = ccoef(2,njets-2)
                        goto 100
                     endif
                  endif 
               enddo
            endif
         enddo 
c
c     jproc=17   njets>=4
c 17  q qbar -> q' qbar' q'' qbar'' (njets-4) g
c     and qbar q -> qbar' q' q'' qbar'' (njets-4) g   q=u,d,c,s,b
      else if(jproc.eq.17) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j= -nlqp,nlqp
                  if(j.ne.0.and.abs(j).ne.abs(i)) then
                     do kj= -nlqp,nlqp
                        if((kj.ne.0).and.(abs(kj).ne.abs(j)).and.
     +                     (abs(kj).ne.abs(i))) then
                           icount=icount+1
                           tmp(icount)=f1(i)*f2(-i)
                           slum=slum+tmp(icount)
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j= -nlqp,nlqp
                  if(j.ne.0.and.abs(j).ne.abs(i)) then
                     do kj= -nlqp,nlqp
                        if((kj.ne.0).and.(abs(kj).ne.abs(j)).and.
     +                     (abs(kj).ne.abs(i))) then
                           icount=icount+1
                           tmptot=tmptot+tmp(icount)
                           if(tmptot.ge.rn) then 
                              ifl(1)=i
                              ifl(2)=-i
                              ifl(3)=j
                              ifl(4)=-j
                              ifl(5)=kj
                              ifl(6)=-kj
                              do k=7,njets+2
                                 ifl(k)=0
                              enddo
                              itmp=imap(ifl(1))
                              effco = ccoef(3,npart-6)
c to remove double counting (j <-> -j, kj <-> kj) and label exchange 
c                of q' and q''
                              slum= slum/2.d0/2.d0 /2.d0
                              goto 100
                           endif
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo 
c     jproc=18   njets>=4
c 18  q qbar -> q' qbar' q' qbar' (njets-4) g
c     and qbar q -> qbar' q' q' qbar' (njets-4) g   q=u,d,c,s,b
      else if(jproc.eq.18) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j= -nlqp,nlqp
                  if(j.ne.0.and.abs(j).ne.abs(i)) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(-i)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j= -nlqp,nlqp
                  if(j.ne.0.and.abs(j).ne.abs(i)) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=-i
                        ifl(3)=j
                        ifl(4)=-j
                        ifl(5)=j
                        ifl(6)=-j
                        do k=7,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(1))
                        effco = ccoef(3,npart-6)
c to remove double counting (j <-> -j)
                        slum= slum/2.d0
                        slum= slum/4.d0 !to account for identical f.s. particles
                        goto 100
                     endif
                  endif
               enddo
            endif
         enddo 
c
c     jproc=19   njets>=4
c 19  q qbar -> q qbar q' qbar' (njets-4) g
c     and qbar q -> qbar q q' qbar' (njets-4) g   q=u,d,c,s,b
c
      else if(jproc.eq.19) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j= -nlqp,nlqp
                  if(j.ne.0.and.abs(j).ne.abs(i)) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(-i)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j= -nlqp,nlqp
                  if(j.ne.0.and.abs(j).ne.abs(i)) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=-i
                        ifl(3)=i
                        ifl(4)=-i
                        ifl(5)=j
                        ifl(6)=-j
                        do k=7,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(1))
                        effco = ccoef(3,npart-6)
c to remove double counting (j <-> -j)
                        slum= slum/2.d0
                        goto 100
                     endif
                  endif
               enddo
            endif
         enddo 
c
c     jproc=20   njets>=4
c 20  q qbar -> q qbar q qbar (njets-4) g
c     and qbar q -> qbar q q qbar (njets-4) g   q=u,d,c,s,b
c
      else if(jproc.eq.20) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(i)*f2(-i)
               slum=slum+tmp(icount)
            endif
         enddo
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=i
                  ifl(2)=-i
                  ifl(3)=i
                  ifl(4)=-i
                  ifl(5)=i
                  ifl(6)=-i
                  do k=7,njets+2
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(1))
                  effco = ccoef(3,npart-6)
                  slum= slum/4.d0 !to account for identical f.s. particles
                  goto 100
               endif
            endif
         enddo 
c     jproc=21   njets>=4
c 21  q q -> q q q qbar (njets-4) g
c     and qbar qbar -> qbar qbar qbar q (njets-4) g   q=u,d,c,s,b
      else if(jproc.eq.21) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(i)*f2(i)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=i
                  ifl(2)=i
                  ifl(3)=i
                  ifl(4)=i
                  ifl(5)=i
                  ifl(6)=-i
                  do k=7,njets+2
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(1))
                  effco = ccoef(3,npart-6)
c statistical factor to account for indistinguishable quarks
                  slum= slum/3.d0/2.d0 
                  goto 100
               endif 
            endif
         enddo 
c     jproc=22   njets>=4
c 22  q q -> q q q' qbar' (njets-4) g
c     and qbar qbar -> qbar qbar qbar' q' (njets-4) g   q=u,d,c,s,b
      else if(jproc.eq.22) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j)) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(i)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j)) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=i
                        ifl(3)=i
                        ifl(4)=i
                        ifl(5)=j
                        ifl(6)=-j
                        do k=7,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(1))
                        effco = ccoef(3,npart-6)
c     statistical factor to account for indistinguishable quarks
                        slum= slum/2.d0 
c     to remove double counting on j (sign + and -)
                        slum= slum/2.d0
                        goto 100
                     endif
                  endif
               enddo 
            endif
         enddo 
c     jproc=23   njets>=4
c 23  q q' -> q q' q qbar (njets-4) g
c     and qbar qbar' -> qbar qbar' q qbar (njets-4) g   q=u,d,c,s,b
      else if(jproc.eq.23) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(j)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=j
                        ifl(3)=i
                        ifl(4)=j
                        ifl(5)=i
                        ifl(6)=-i
                        do k=7,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(1))
                        effco = ccoef(3,npart-6)
c   to account for indistingushable particles
                        slum= slum/2.d0
                        goto 100
                     endif
                  endif
               enddo 
            endif
         enddo 
c     jproc=24   njets>=4
c 24  q q' -> q q' q' qbar' (njets-4) g
c     and qbar qbar' -> qbar qbar' q' qbar' (njets-4) g   q=u,d,c,s,b
      else if(jproc.eq.24) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(j)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=j
                        ifl(3)=i
                        ifl(4)=j
                        ifl(5)=j
                        ifl(6)=-j
                        do k=7,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(1))
                        effco = ccoef(3,npart-6)
c   to account for indistingushable particles
                        slum= slum/2.d0
                        goto 100
                     endif
                  endif
               enddo 
            endif
         enddo 
c     jproc=25   njets>=4
c 25  q q' -> q q' q'' qbar'' (njets-4) g
c     and qbar qbar' -> qbar qbar' q'' qbar'' (njets-4) g   q=u,d,c,s,b
      else if(jproc.eq.25) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     do kj=-nlqp,nlqp
                        if(kj.ne.0.and.abs(kj).ne.abs(j).and.
     +                     abs(kj).ne.abs(i)) then
                           icount=icount+1
                           tmp(icount)=f1(i)*f2(j)
                           slum=slum+tmp(icount)
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.gt.0) then
                     do kj=-nlqp,nlqp
                        if(kj.ne.0.and.abs(kj).ne.abs(j).and.
     +                     abs(kj).ne.abs(i)) then
                           icount=icount+1
                           tmptot=tmptot+tmp(icount)
                           if(tmptot.ge.rn) then 
                              ifl(1)=i
                              ifl(2)=j
                              ifl(3)=i
                              ifl(4)=j
                              ifl(5)=kj
                              ifl(6)=-kj
                              do k=7,njets+2
                                 ifl(k)=0
                              enddo
                              itmp=imap(ifl(1))
                              effco = ccoef(3,npart-6)
c to account for sign of q''
                              slum= slum/2.d0
                              goto 100
                           endif
                        endif
                     enddo
                  endif
               enddo 
            endif
         enddo 
c     jproc=26   njets>=4
c 26  q qbar' -> q qbar' q qbar (njets-4) g
c     and qbar q' -> qbar q' qbar q (njets-4) g   q=u,d,c,s,b
      else if(jproc.eq.26) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.lt.0) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(j)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.lt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=j
                        ifl(3)=i
                        ifl(4)=j
                        ifl(5)=i
                        ifl(6)=-i
                        do k=7,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(1))
                        effco = ccoef(3,npart-6)
c to account for indistinguishable final state particles
                        slum= slum/2.d0
                        goto 100
                     endif
                  endif
               enddo 
            endif
         enddo 
c     jproc=27   njets>=4
c 27  q qbar' -> q qbar' q' qbar' (njets-4) g
c     and qbar q' -> qbar q' qbar' q' (njets-4) g   q=u,d,c,s,b
      else if(jproc.eq.27) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.lt.0) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(j)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.lt.0) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=j
                        ifl(3)=i
                        ifl(4)=j
                        ifl(5)=j
                        ifl(6)=-j
                        do k=7,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(1))
                        effco = ccoef(3,npart-6)
                        goto 100
                     endif
                  endif
               enddo 
            endif
         enddo 
c     jproc=28   njets>=4
c 28  q qbar' -> q qbar' q'' qbar'' (njets-4) g
c     and qbar q' -> qbar q' qbar'' q'' (njets-4) g   q=u,d,c,s,b
      else if(jproc.eq.28) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.lt.0) then
                     do kj=-nlqp,nlqp
                        if(kj.ne.0.and.abs(kj).ne.abs(i).and.
     +                     abs(kj).ne.abs(j)) then
                           icount=icount+1
                           tmp(icount)=f1(i)*f2(j)
                           slum=slum+tmp(icount)
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j).and.i*j.lt.0) then
                     do kj=-nlqp,nlqp
                        if(kj.ne.0.and.abs(kj).ne.abs(i).and.
     +                     abs(kj).ne.abs(j)) then
                           icount=icount+1
                           tmptot=tmptot+tmp(icount)
                           if(tmptot.ge.rn) then 
                              ifl(1)=i
                              ifl(2)=j
                              ifl(3)=i
                              ifl(4)=j
                              ifl(5)=kj
                              ifl(6)=-kj
                              do k=7,njets+2
                                 ifl(k)=0
                              enddo
                              itmp=imap(ifl(1))
                              effco = ccoef(3,npart-6)
c     to account for the sign of q''
                              slum= slum/2.d0
                              goto 100
                           endif
                        endif
                     enddo 
                  endif
               enddo
            endif
         enddo 
c     jproc=29   njets>=5
c 29  g q    -> q q qbar q qbar (njets-5) g 
c     and g qbar  -> qbar q qbar q qbar (njets-5) g
      else if(jproc.eq.29) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(0)*f2(i)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=0
                  ifl(2)=i
                  ifl(3)=i
                  ifl(4)=i
                  ifl(5)=-i
                  ifl(6)=i
                  ifl(7)=-i
                  do k=8,njets+2
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(2))
                  effco = ccoef(3,npart-6)
c statistical factor to account for indistinguishable quarks
                  slum= slum/3.d0/2.d0 /2.d0   ! /3! /2|
                  goto 100
               endif 
            endif
         enddo 
c     jproc=30   njets>=5
c 30  q g    -> q q qbar q qbar (njets-5) g 
c     and qbar g  -> qbar q qbar q qbar (njets-5) g
      else if(jproc.eq.30) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=f1(i)*f2(0)
               slum=slum+tmp(icount)
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=i
                  ifl(2)=0
                  ifl(3)=i
                  ifl(4)=i
                  ifl(5)=-i
                  ifl(6)=i
                  ifl(7)=-i
                  do k=8,njets+2
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(2))
                  effco = ccoef(3,npart-6)
c statistical factor to account for indistinguishable quarks
                  slum= slum/3.d0/2.d0 /2.d0   ! /3! /2|
                  goto 100
               endif 
            endif
         enddo 
c     jproc=31   njets>=5
c 31  g q    -> q q' qbar' q' qbar' (njets-5) g 
c     and g qbar  -> qbar q' qbar' q' qbar' (njets-5) g
      else if(jproc.eq.31) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j)) then
                     icount=icount+1
                     tmp(icount)=f1(0)*f2(i)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j)) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=0
                        ifl(2)=i
                        ifl(3)=i
                        ifl(4)=j
                        ifl(5)=-j
                        ifl(6)=j
                        ifl(7)=-j
                        do k=8,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(2))
                        effco = ccoef(3,npart-6)
c statistical factor to account for indistinguishable quarks
                        slum= slum/2.d0 /2.d0 ! /2! /2|
c to account for the sign of q'
                        slum= slum/2.d0
                        goto 100
                     endif 
                  endif
               enddo
            endif
         enddo 
c     jproc=32   njets>=5
c 32  q g    -> q q' qbar' q' qbar' (njets-5) g 
c     and qbar g  -> qbar q' qbar' q' qbar' (njets-5) g
      else if(jproc.eq.32) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j)) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(0)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j)) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=0
                        ifl(3)=i
                        ifl(4)=j
                        ifl(5)=-j
                        ifl(6)=j
                        ifl(7)=-j
                        do k=8,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(2))
                        effco = ccoef(3,npart-6)
c statistical factor to account for indistinguishable quarks
                        slum= slum/2.d0 /2.d0 ! /2! /2|
c to account for the sign of q'
                        slum= slum/2.d0
                        goto 100
                     endif 
                  endif
               enddo
            endif
         enddo 
c     jproc=33   njets>=5
c 33  g q    -> q q qbar q' qbar' (njets-5) g 
c     and g qbar  -> qbar q qbar q' qbar' (njets-5) g
      else if(jproc.eq.33) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j)) then
                     icount=icount+1
                     tmp(icount)=f1(0)*f2(i)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j)) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=0
                        ifl(2)=i
                        ifl(3)=i
                        ifl(4)=i
                        ifl(5)=-i
                        ifl(6)=j
                        ifl(7)=-j
                        do k=8,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(2))
                        effco = ccoef(3,npart-6)
c statistical factor to account for indistinguishable quarks
                        slum= slum/2.d0    ! /2!
c to account for the sign of q'
                        slum= slum/2.d0
                        goto 100
                     endif 
                  endif
               enddo
            endif
         enddo 
c     jproc=34   njets>=5
c 34  q g    -> q q qbar q' qbar' (njets-5) g 
c     and qbar g  -> qbar q qbar q' qbar' (njets-5) g
      else if(jproc.eq.34) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j)) then
                     icount=icount+1
                     tmp(icount)=f1(i)*f2(0)
                     slum=slum+tmp(icount)
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j)) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=i
                        ifl(2)=0
                        ifl(3)=i
                        ifl(4)=i
                        ifl(5)=-i
                        ifl(6)=j
                        ifl(7)=-j
                        do k=8,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(2))
                        effco = ccoef(3,npart-6)
c statistical factor to account for indistinguishable quarks
                        slum= slum/2.d0 ! /2|
c to account for the sign of q'
                        slum= slum/2.d0
                        goto 100
                     endif 
                  endif
               enddo
            endif
         enddo 
c     jproc=35   njets>=5
c 35  g q    -> q q' qbar' q'' qbar'' (njets-5) g 
c     and g qbar  -> qbar q' qbar' q'' qbar'' (njets-5) g
      else if(jproc.eq.35) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j)) then
                     do kj=-nlqp,nlqp
                        if(kj.ne.0.and.abs(kj).ne.abs(i).and.
     +                     abs(kj).ne.abs(j)) then
                           icount=icount+1
                           tmp(icount)=f1(0)*f2(i)
                           slum=slum+tmp(icount)
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j)) then
                     do kj=-nlqp,nlqp
                        if(kj.ne.0.and.abs(kj).ne.abs(i).and.
     +                     abs(kj).ne.abs(j)) then
                           icount=icount+1
                           tmptot=tmptot+tmp(icount)
                           if(tmptot.ge.rn) then 
                              ifl(1)=0
                              ifl(2)=i
                              ifl(3)=i
                              ifl(4)=j
                              ifl(5)=-j
                              ifl(6)=kj
                              ifl(7)=-kj
                              do k=8,njets+2
                                 ifl(k)=0
                              enddo
                              itmp=imap(ifl(2))
                              effco = ccoef(3,npart-6)
c to account for sign of q' and q''
                              slum= slum/2.d0/2.d0
c to account for exchange of q' and q'' current
                              slum= slum/2.d0
                              goto 100
                           endif 
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo 
c     jproc=36   njets>=5
c 36  q g    -> q q' qbar' q'' qbar'' (njets-5) g 
c     and qbar g  -> qbar q' qbar' q'' qbar'' (njets-5) g
      else if(jproc.eq.36) then
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j)) then
                     do kj=-nlqp,nlqp
                        if(kj.ne.0.and.abs(kj).ne.abs(i).and.
     +                     abs(kj).ne.abs(j)) then
                           icount=icount+1
                           tmp(icount)=f1(i)*f2(0)
                           slum=slum+tmp(icount)
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*slum
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.abs(i).ne.abs(j)) then
                     do kj=-nlqp,nlqp
                        if(kj.ne.0.and.abs(kj).ne.abs(i).and.
     +                     abs(kj).ne.abs(j)) then
                           icount=icount+1
                           tmptot=tmptot+tmp(icount)
                           if(tmptot.ge.rn) then 
                              ifl(1)=i
                              ifl(2)=0
                              ifl(3)=i
                              ifl(4)=j
                              ifl(5)=-j
                              ifl(6)=kj
                              ifl(7)=-kj
                              do k=8,njets+2
                                 ifl(k)=0
                              enddo
                              itmp=imap(ifl(2))
                              effco = ccoef(3,npart-6)
c to account for sign of q' and q''
                              slum= slum/2.d0/2.d0
c to account for exchange of q' and q'' current
                              slum= slum/2.d0
                              goto 100
                           endif 
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo 
c     jproc=37   njets>=6
c  37  g g -> q qbar q qbar q qbar (njets-6 g) q=u,d,c,s,b
      else if(jproc.eq.37) then
         slum=f1(0)*f2(0)
c
c  slum x number of possible f.s. q qbar pairs
c
         slum=slum*dble(nlqp)
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmp(icount)=1.
            endif
         enddo 
         icount=0
         rn=rn*nlqp*2.d0       !to count plus and minus signs of q
         do i=-nlqp,nlqp
            if(i.ne.0) then
               icount=icount+1
               tmptot=tmptot+tmp(icount)
               if(tmptot.ge.rn) then 
                  ifl(1)=0
                  ifl(2)=0
                  ifl(3)=i
                  ifl(4)=-i
                  ifl(5)=i
                  ifl(6)=-i
                  ifl(7)=i
                  ifl(8)=-i
                  do k=9,njets+2
                     ifl(k)=0
                  enddo
                  itmp=imap(ifl(3))
                  effco = ccoef(3,npart-6)
c     statistical factor to account for indistinguishable quarks
                  slum= slum/3.d0/2.d0/3.d0/2.d0  ! /3! /3!
                  goto 100
               endif 
            endif
         enddo 
c     jproc=38   njets>=6
c  38  g g -> q qbar q' qbar' q'' qbar'' (njets-6 g) q=u,d,c,s,b
      else if(jproc.eq.38) then
         slum=f1(0)*f2(0)
c
c  slum x number of possible f.s. q qbar pairs
c
         slum=slum*dble(nlqp)*dble(nlqp-1)*dble(nlqp-2)/6.d0
c 6 is the number of permutations of the three quark currents
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.(abs(i).ne.abs(j))) then
                     do kj=-nlqp,nlqp
                        if(kj.ne.0.and.(abs(kj).ne.abs(j)).and.
     +                     (abs(kj).ne.abs(i))) then
                           icount=icount+1
                           tmp(icount)=1.
                        endif
                     enddo
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*dble(nlqp)*dble(nlqp-1)*dble(nlqp-2) * 8.d0  !2^3 to count plus and minus signs of q, q' and q''
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.(abs(i).ne.abs(j))) then
                     do kj=-nlqp,nlqp
                        if(kj.ne.0.and.(abs(kj).ne.abs(j)).and.
     +                     abs(kj).ne.abs(i)) then
                           icount=icount+1
                           tmptot=tmptot+tmp(icount)
                           if(tmptot.ge.rn) then 
                              ifl(1)=0
                              ifl(2)=0
                              ifl(3)=i
                              ifl(4)=-i
                              ifl(5)=j
                              ifl(6)=-j
                              ifl(7)=kj
                              ifl(8)=-kj
                              do k=9,njets+2
                                 ifl(k)=0
                              enddo
                              itmp=imap(ifl(3))
                              effco = ccoef(3,npart-6)
                              goto 100
                           endif
                        endif 
                     enddo
                  endif
               enddo
            endif
         enddo 
c     jproc=39   njets>=6
c  39  g g -> q qbar q qbar q' qbar' (njets-6 g) q=u,d,c,s,b
      else if(jproc.eq.39) then
         slum=f1(0)*f2(0)
c
c  slum x number of possible f.s. q qbar pairs
c
         slum=slum*dble(nlqp)*dble(nlqp-1)/2.d0 
c  2 to count the number of permutations of q and q' currents
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.(abs(i).ne.abs(j))) then
                     icount=icount+1
                     tmp(icount)=1.
                  endif
               enddo
            endif
         enddo 
         icount=0
         rn=rn*dble(nlqp)*dble(nlqp-1) * 4.d0  
c  4 to count plus and minus signs of q and q'
         do i=-nlqp,nlqp
            if(i.ne.0) then
               do j=-nlqp,nlqp
                  if(j.ne.0.and.(abs(i).ne.abs(j))) then
                     icount=icount+1
                     tmptot=tmptot+tmp(icount)
                     if(tmptot.ge.rn) then 
                        ifl(1)=0
                        ifl(2)=0
                        ifl(3)=i
                        ifl(4)=-i
                        ifl(5)=i
                        ifl(6)=-i
                        ifl(7)=j
                        ifl(8)=-j
                        do k=9,njets+2
                           ifl(k)=0
                        enddo
                        itmp=imap(ifl(3))
                        effco = ccoef(3,npart-6)
c to account for 2 indistinguishable quark pairs
                        slum= slum/2.d0/2.d0
                        goto 100
                     endif
                  endif 
               enddo
            endif
         enddo 
c
      else
         write(*,*) 'jproc not defined, slum=0'
         slum=0
         stop
      endif
c
      xlum=-1d0
      return
 100  continue
c
      do i=1,maxpar
         afl(i)=ifl(i)
      enddo
c
c     evaluate colour weight factors
      do i=1,2
         if(ifl(i).eq.0) ifl(i)=21
      enddo
      ng=0
      cwgt=1e0
      do i=3,njets+2
         cwgt=cwgt*cfac(abs(ifl(i)))
         if(ifl(i).eq.0) then
            ifl(i)=21
            ng=ng+1
         endif
      enddo
c     evaluate spin weight factors
      swgt=2e0
      swgt=swgt**(njets)
*
      xlum=dble(slum*cwgt*swgt*ifact(ng)) /resc**(njets)
      xlum=xlum*effco
c
c     complete flavour-dependent momenta assignements in the user common block:
      call usrfll()
      end                                                
c
c
      subroutine usrfll()
      implicit none
      include 'alpgen.inc'
      include 'Njet.inc'
      integer i,j
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
c
      do i=1,4
         do j=1,2
            pin(i,j)=p(i,j)
         enddo
         do j=1,npart-2
            pout(i,j)=p(i,j+2)
         enddo
      enddo
      do i=1,4
         do j=1,njets
            pjet(i,j)=p(i,j+2)
            ptj(j)=pt(j+2)
            etaj(j)=eta(j+2)
         enddo
      enddo
      do j=1,njets
         do i=j+1,njets
            drjj(i,j)=dr(i+2,j+2)
            drjj(j,i)=drjj(i,j)
         enddo
      enddo
      end
c
      subroutine setdec(nwrt,iflwrt,icuwrt,pwrt,decbr)
      implicit none
      include 'alpgen.inc'
      include 'Njet.inc'
C--   GUP Event common block
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,
     &     IDUP(MAXNUP),ISTUP(MAXNUP),MOTHUP(2,MAXNUP),
     &     ICOLUP(2,MAXNUP),PUP(5,MAXNUP),VTIMUP(MAXNUP),
     &     SPINUP(MAXNUP)
      integer iup,i,j
c debug
      double precision rate(16)
      common/dbg/rate
      data rate/16*0d0/
c locals
      integer maxdec
      parameter (maxdec=40)
      integer ip, ic, il
c arguments
      integer nwrt,iflwrt(maxdec),icuwrt(2,maxdec)
      double precision pwrt(5,maxdec),decbr
c
      do ip=1,npart
        iflwrt(ip)=ifl(ip)
        do ic=1,2
          icuwrt(ic,ip)=icu(ic,ip)
        enddo
        do il=1,5
          pwrt(il,ip)=p(il,ip)
        enddo
      enddo
      nwrt=npart
      decbr=1d0
 100  if(ilhe.eq.0) return
C prepare LHE format
      nup=nwrt
      do iup=1,nup
        idup(iup)=iflwrt(iup)
        do j=1,5
          pup(j,iup)=pwrt(j,iup)
        enddo
        vtimup(iup)=0.d0
        spinup(iup)=spnwrt(iup)
      enddo
c---  initial state lines
      do iup=1,2
        istup(iup)=-1
        mothup(1,iup)=0
        mothup(2,iup)=0
      enddo
c---  outgoing lines
      do iup=3,nup
        istup(iup)=1
        mothup(1,iup)=1
        mothup(2,iup)=2
      enddo
c---  set up colour structure labels
      Do iup=1,nup
        do i=1,2
          icolup(i,iup)=0
        enddo
        if(icuwrt(1,iup).ne.0) icolup(1,iup)=icuwrt(1,iup)+500
        if(icuwrt(2,iup).ne.0) icolup(2,iup)=icuwrt(2,iup)+500
      Enddo
      end
c
      subroutine phspace(lnot,pswgt,djpd,djg)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                               c
c  (njets)-body phase space for the process:                    c
c                                                               c
c  h(1)  h(2) -> h(3) ... h(njets)                              c
c                                                               c
c  where h is any light-flavoured quark or gluon                c
c                                                               c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'alpgen.inc'
      include 'Njet.inc'
      real *8 dummy,djg
      real *8 pswgt
      real *8 djpd,factor
      real *8 cutkin(10)
      real *8 wgt
      common/loccut/cutkin
      real *8 pl(maxpar),y(maxpar)
      real *8 pcm(0:3,maxpar)
*
c-    debugging variables
*
      real *8 totpt
      integer nx,ninit,i,l,m,lnot
      parameter (nx= 20)
      data ninit/0/
      save ninit
      if(ninit.eq.0) then
*
c-       setup local generation cuts
*
         cutkin(1)=ptjmin
         cutkin(3)=etajmax
         cutkin(5)=drjmin
         cutkin(7)=0.d0
         cutkin(8)=0.d0
         cutkin(9)=ptj1min
         cutkin(10)=ptj1max
c-
         ninit=1
      endif
*
c-    The generation starts
*
      pswgt=0.d0
      call momgen(njets,roots,x1,x2,pcm,wgt,lnot)
      djg= 1.d0 ! dummy variable
      if (lnot.eq.1) then
         pswgt= 0.d0
         goto 100
      endif
*
c-    will write factor=factor0/(x1*x2), with factor0 function of njets etc.
*
      factor= 1d0/(2.d0*pi)**(3*(njets)-4)/2.d0/s/x1/x2
*
c-    initial state momenta in the LAB frame:
*
      pcm(0,1)= roots/2.d0*x1   
      pcm(1,1)= 0.d0   
      pcm(2,1)= 0.d0   
      pcm(3,1)= roots/2.d0*x1   

      pcm(0,2)= roots/2.d0*x2   
      pcm(1,2)= 0.d0   
      pcm(2,2)= 0.d0   
      pcm(3,2)=-roots/2.d0*x2   

c-    Rapidities and pseudo-rapidities (in the LAB system), pt's, deltar's:

      do l= 3,njets+2
          pt(l) = sqrt(pcm(1,l)**2+pcm(2,l)**2)
          pl(l) = pcm(3,l)
          eta(l)= -log(tan(0.5d0*atan2(pt(l),pl(l))))
          y(l)=(pcm(0,l)+pcm(3,l))/(pcm(0,l)-pcm(3,l))
          y(l)=dabs(y(l))
          y(l)=0.5d0*log(y(l))
      enddo
*
c-    Calculates jet-jet distances:
*
      do l= 3,njets+1
         do m= l+1,njets+2
            if(min(pt(l),pt(m)).gt.0d0) then
               dphi(l,m)= 
     .              (pcm(1,l)*pcm(1,m)+pcm(2,l)*pcm(2,m))
     .              /pt(l)/pt(m)
               if(dabs(dphi(l,m)).gt.1.d0) then
c                  write(*,*) 'partons',l,m,', cos(Dphi)=', dphi(l,m)
c     .                 ,', set to +-1'
c                  write(*,*) 'pt(',l,')=',pt(l),'p(',l,')=',  (pcm(i,l)
c     .                 ,i=0,3)
c                  write(*,*) 'pt(',m,')=',pt(m),'p(',m,')=',  (pcm(i,m)
c     .                 ,i=0,3)
                  if (dphi(l,m).gt.0.d0) dphi(l,m)= 1.d0
                  if (dphi(l,m).lt.0.d0) dphi(l,m)=-1.d0
               endif
               dphi(l,m)=acos(dphi(l,m))
               dphi(m,l)=dphi(l,m)
            else
c***fix***
c                dphi(l,m)=pi
c                dphi(m,l)=pi
               goto 100
            endif
         enddo
      enddo
*
      do l= 3,njets+1
         do m= l+1,njets+2
            dr(l,m)= sqrt(dphi(l,m)**2+(eta(l)-eta(m))**2)
            dr(m,l)=dr(l,m)
         enddo
      enddo
*
c-    Redefine the momenta:           
*                                     
      do l= 1,njets+2            
         do m=1,2                     
            p(m,l)=pcm(m,l)           
         enddo                        
         p(4,l)= pcm(0,l)             
         p(3,l)= pcm(3,l)             
      enddo                           
comment
c      do l= 1,npart
c        print*,'osh',l,'=',p(4,l)**2-p(1,l)**2-p(2,l)**2-p(3,l)**2
c      enddo     
c
c      do l= 1,4
c        dummy= -p(l,1)-p(l,2)
c        do m= 3,npart
c          dummy= dummy+p(l,m)
c        enddo
c        print*,'summ',l,'=',dummy
c      enddo
c      print*,'   '
comment
*                                     
c-    Initial state pt, eta, and y:         
*                                     
      do l= 1,2                       
         pt(l) = 0.d0                 
         y(l)  = 1.d6                 
         eta(l)= y(l)                 
      enddo                           
*
c-    Call to the the cut routine:
*
c     first gobal cuts on light jets, common to all processes (routine
c     alccut in alpgen.f. Penultimate argument: position of first jet
      call alccut(lnot,pt,eta,dr,3)
      if (lnot.eq.1) then
         pswgt= 0.d0
         goto 100
      endif
c     then cuts specific to this hard process
      call chkcut(lnot,pt,p,eta,dr,njets)
      if (lnot.eq.1) then
         pswgt= 0.d0
         goto 100
      endif
*
      pswgt = factor*wgt
*
c-    Evaluate q2: will include several possible options:
*
      totpt= 0.d0
*
c-    Total et^2 of jets
*
      do i=3,njets+2
         totpt=totpt+pt(i)**2
      enddo 
      if(iqopt.eq.0) then 
         qsq=1d0
      elseif(iqopt.eq.1) then
         qsq=totpt
ccc modified fulvio
         qsq= qsq/float(njets)
ccc modifief fulvio
      elseif(iqopt.eq.2) then
         qsq=roots*roots*x1*x2
      endif
      qsq=qfac**2*qsq

 100  continue
      end
*
      subroutine chkcut(lnot,pt,p,eta,dr,njets)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                c
c     Applies kinematical cuts to the final state during the phase
c     -space generation                                          c
c                                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      real*8 ptjmin,drjmin,etajmax,ptj1min,ptj1max,ptmax
      integer maxpar,ninit,njets,j,lnot,i
      real*8 cutkin(10)
      common/loccut/cutkin
      data ninit/0/
      parameter (maxpar=20)
      real*8 pt(maxpar),eta(maxpar),dr(maxpar,maxpar),p(5,maxpar)
      save ninit,ptjmin,etajmax,drjmin
      if(ninit.eq.0) then
         ninit=1
         ptjmin=cutkin(1)
         etajmax=cutkin(3)
         drjmin=cutkin(5)
         ptj1min=cutkin(9)
         ptj1max=cutkin(10)
      endif

      lnot= 0
c
 5    return

 10   lnot= 1
      return
      end
*
      subroutine momgen(njets,roots,x1,x2,p,wgt,lw)
*
c-    Generator of np particles in the LAB frame.
*
      implicit none
      real*8 roots,x1,x2,wgt,zero,etacut
      real*8 ptjmin,etajmax,drjmin,cutkin,ptj1min,ptj1max
      double precision s,xmsum0,pt0lmax,pt0lsum,eta0min
      double precision wtau,wjr,en,pz
      double precision ag,tau0,tau
      real*8 rootsh,ranram,sq,y0,yr
      real*8 rcnmin,rcxml,cxml
      real*8 xmsum,ptlim
      real*8 wgt1,wgt2
c
c global event cuts
c      double precision ptjmin,ptjmax,etajmin,etajmax,drjmin
c      common/glocut/ptjmin,ptjmax,etajmin,etajmax,drjmin
      integer lw,npar,mpar,np,njets,j,k,lw1
      parameter (npar= 20)
      double precision pm(0:4,npar),pt0(npar),pt1(npar),xm(npar),
     .       xmr(npar),eta0(npar),pt0_l(npar),pt1_l(npar),p(0:3,npar)
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      common/loccut/cutkin(10)
      dimension ranram(111)
      integer nvar,nv,m,l,nbin,ndummy,mask,mmask,md
      integer nct,nx1,nx2,maxn,lmin,lmax
      parameter (maxn= 500)
      real*8 dj,djbintot,djbin,dummy,apw(1:2)
      real*8 djb0,djb1,ran0
      common/psopt/nvar,nv
      real*8 peropt
      common/psopt1/mask(maxn),mmask(maxn),peropt(maxn)
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      data mpar/0/
      save
*
      if ((njets.gt.npar-2).or.(njets.lt.0)) then
        write (6,*) 'WRONG NJETS VALUE'
        stop
      endif
      if (mpar.eq.0) then
        mpar   = 1
        np     = njets
        zero   = 0.d0
        etacut = 40.d0
        s      = roots*roots
*
        ptjmin = cutkin(1)
        etajmax= cutkin(3)
        drjmin = cutkin(5)
        ptj1min = cutkin(9)
        ptj1max = cutkin(10)
*
        do j= 1,np
           pt0_l(j) = ptjmin
           pt1_l(j) = min(ptj1max,roots/2.d0)
           eta0(j)  = etajmax
           xm(j)    = zero        
        enddo
*
c-      Parameters for the x1,x2 integration:
*
        ag  = 0.98d0
*
        rcnmin= 2.d0*ptjmin*sin(drjmin/2.d0)
*
c-      l is the number of light jets
*
        l    = njets                               
        rcxml= int(l/2)*rcnmin
        cxml = rcxml*rcxml
*
        xmsum0 = 0.d0
        pt0lmax= 0.d0
        pt0lsum= 0.d0
        eta0min= 100.d0
        do j= 1,np
           xmsum0 = xmsum0+xm(j)
           pt0lmax= max(pt0lmax,pt0_l(j))
           pt0lsum= pt0lsum+pt0_l(j)
           eta0min= min(eta0min,eta0(j))
        enddo
        xmsum= xmsum0
*
c-      sets the ratios of calls mom/momr
*
        apw(1)= 0.9d0
        apw(2)= 1.d0-apw(1)
      endif
*
      lw=   0
      lmin= (jgrid(jproc)-1)*nv+2
      lmax= jgrid(jproc)*nv+1
*
      dj= 0.d0
      call rans(ran0)
      if (ran0.le.apw(1)) then
         md= 1
      else
         md= 2 
      endif
*
      djb1= 1.d0
      djb0= 1.d0
      m= 0
      do l= lmin,lmax
        m= m+1
        if (md.eq.1) then
          mmask(l)= 1
        else
          mmask(l)= mask(l)
        endif
        if (mmask(l).eq.1) then
          call onedimbin(1,nbin,djbin,l,ndummy,dummy)
          call grans(nbin,nx1(l),ranram(m))
          if (mask(l).eq.1) then
            djb1= djb1*djbin
          elseif (mask(l).eq.0) then
            djb0= djb0*djbin
          endif
        else
          ranram(m)= 0.d0
        endif
      enddo
      djbintot= djb1
      if (md.eq.1) djbintot= djbintot*djb0 
*
c-    tau0 is the lower cut on tau= x1*x2
*
      tau0   = 1.d0/s*(xmsum)**2
      tau0   = max(tau0,pt0lsum**2/s)
      tau0   = max(tau0,1.d0/s*
     +        (float(njets)*ptjmin)**2)
*
      tau0   = max(tau0,cxml/s)

      if(tau0.gt.1.) goto 100
*
c-    Generates x1 and x2: 
*
      call ppeaka(0,tau0,ag,tau0,1.d0,tau,wtau,ranram(1),lw1)
      sq = sqrt(tau)
      y0 =-0.5d0*log(tau) 
      call ttriangle(0,-y0,y0,-y0,y0,yr,wjr,ranram(2),lw1)
*
      x1 = sq*exp(yr) 
      x2 = sq*exp(-yr) 
      rootsh= roots*sq
*
c-    Protection:
*
      if (rootsh.lt.(xmsum+pt0lsum)) goto 100
*
c-    Rescalings and transformations to feed mom:
*
      do j= 1,np
        ptlim = s*s
     .           +(xm(j)**2-(xmsum-xm(j))**2)**2
     .        -2.d0*s*(xm(j)**2+(xmsum-xm(j))**2)   
        if (ptlim.le.0.d0) goto 100
        ptlim =  0.5d0/roots*sqrt(ptlim)
        pt0(j)= pt0_l(j)/roots
        pt1(j)= min(pt1_l(j),ptlim)/roots
*
c-      Protection:
*
        if (pt0(j).gt.pt1(j)) goto 100
        xmr(j)= xm(j)/roots
      enddo
*
c-    momenta:
*
      en    = 0.5d0*(x1+x2)      ! Rescaled initial total energy
      pz    = 0.5d0*(x1-x2)      ! Rescaled initial longitudinal momentum
*
      if (md.eq.1) then
        call mom(0,np,pt0,pt1,eta0,xmr,en,pz,pm,wgt1,ranram,lw)
        if (lw.eq.0) then
          call momr(1,np,xmr,en,pz,pm,wgt2,lw)
          if (lw.eq.0) then
            if (wgt1.eq.0.d0.and.wgt2.eq.0.d0) goto 100
            wgt= wgt1*wgt2/
     +          (wgt2*apw(1)*djbintot+wgt1*apw(2)*djb1)
          else
            wgt= wgt1/(apw(1)*djbintot)
          endif
        else
          goto 100
        endif
      else
        call momr(0,np,xmr,en,pz,pm,wgt2,lw)
        if (lw.eq.0) then
          call mom(1,np,pt0,pt1,eta0,xmr,en,pz,pm,wgt1,ranram,lw)
          if (lw.eq.0) then
            m= 0
            do l= lmin,lmax
              m= m+1
              if (mask(l).eq.0) then
                nbin= min(nx1(l),int(1.d0+dfloat(nx1(l))*ranram(m)))
                call onedimbin(3,nbin,djbin,l,ndummy,dummy)
                djb0= djb0*djbin
              endif
            enddo
            if (wgt1.eq.0.d0.and.wgt2.eq.0.d0) goto 100
            wgt= wgt1*wgt2/
     +          (wgt1*apw(2)*djbintot+wgt2*apw(1)*djb0*djb1)
          else
            wgt= wgt2/(apw(2)*djbintot)
          endif
        else
          goto 100
        endif
      endif
*
c-    from the x1,x2 integration
*
      wgt= wgt*wtau*wjr
*
c-    Rescaling of momenta and weight:
*
      do j= 1,np
        do k= 0,3
           p(k,j+2)= pm(k,j)*roots
        enddo
      enddo
*
      wgt= wgt*s**(np-2)
*
      return
 100  lw= 1
      wgt= 0.d0
      return
      end
*
c-------------------------------------------------------------------
      subroutine alsgrd
c-------------------------------------------------------------------
      implicit none
      include 'alpgen.inc'
      include 'Njet.inc'
      integer nct1
      integer nvar,nch,n,nv
      integer init,j,k
      real*8 v1,ni
      real*8 al1,bet1 
      integer nct,nx1,nx2
      integer maxn,ncmax
      parameter (maxn= 500,ncmax= 5000)  
      common/printout/nct(maxn),nx1(maxn),nx2(maxn)
      common/rdandwrt/al1(ncmax,maxn),nct1(maxn)
      common/ausil/init(ncmax),bet1(0:ncmax,maxn)
      integer jprocmax,jproc,ngrid,jgrid,maxsubproc
      parameter (maxsubproc= 100)
      common/procct/jprocmax,jproc,ngrid,jgrid(maxsubproc)
      common/psopt/nvar,nv
      common/book/v1(ncmax),ni(maxn)
      data v1/ncmax*0d0/,ni/maxn*0d0/
      integer mask,mmask 
      real*8 peropt
      common/psopt1/mask(maxn),mmask(maxn),peropt(maxn)
      data mask/maxn*0/,mmask/maxn*1/,peropt/maxn*0.9d0/
* 
c-    initialise size of grids:
*
c-    for MC over jproc:

      nct(1)   = 0              ! first bin of variable 1
      nx1(1)   = jprocmax       ! number of jprocs
*
c-    for the Phase-space reweighting:
*
      nv  = 2*(njets)-1 
      nvar= nv*ngrid
      do n= 2,nvar+1
        nch= 10
*
        nct(n)   =  nct(n-1)+nx1(n-1)   
        nx1(n)   =  nch 
*
c-      set the channels shared by mom and momr: 
*
        if (mod(n-1,nv).eq.1)    mask(n)= 1
        if (mod(n-1,nv).eq.2)    mask(n)= 1
*
c-      set the percentages of optimization: 
*
        if (mask(n).eq.1) peropt(n)= 0.9d0 
        if (mask(n).eq.0) peropt(n)= 1.d0-dfloat(nv-1)/200.d0
      enddo
c-
      do n= 1,nvar+1
        nct1(n)= nx1(n)
        do j= 1,nct1(n)                                         

          init(nct(n)+j)= 0                                     

c-        the weights are initialized to be all equal 

          al1(j,n)= 1.d0/nct1(n)                                  
        enddo

        do j= 1,nct1(n)                                         
          bet1(j,n)= 0                                         
          do k= 1,j                                          
            bet1(j,n)= bet1(j,n)+al1(k,n)                          
          enddo                                              
        enddo                                                
      enddo
*
c-    protection:
*
      if (nct(nvar+1)+nct1(nvar+1).gt.ncmax) then
        print*,'INCREASE NCMAX'
        stop
      endif
*
      if ((nvar+1).gt.maxn) then
        print*,'INCREASE MAXN'
        stop
      endif
      end

c-------------------------------------------------------------------
      subroutine dumpwgt
c-------------------------------------------------------------------
c     routine for the debugging of individual large-weight events
c     Not functional to the rnning of the code
      implicit none
      include 'alpgen.inc'
      include 'Njet.inc'
      real *8 tmpmtt,ptmp(4)
      integer i
      write(niosta,*) 'x1,x2=',x1,x2
      write(niosta,*) 'p(Q1)=',(p(i,3),i=1,3)
      write(niosta,*) 'p(Q2)=',(p(i,4),i=1,3)
      write(niosta,*) 'p(j1)=',(p(i,5),i=1,3)
      write(niosta,*) 'p(j2)=',(p(i,6),i=1,3)
      write(niosta,*) 'p(j3)=',(p(i,7),i=1,3)
      write(niosta,*) 'p(j4)=',(p(i,8),i=1,3)
c
c b-bbar inv mass
      ptmp(4)=p(4,3)+p(4,4)
      tmpmtt=ptmp(4)**2
      do i=1,3
         ptmp(i)=p(i,3)+p(i,4)
         tmpmtt=tmpmtt-ptmp(i)**2
      enddo
      tmpmtt=sqrt(tmpmtt)
      write(niosta,*) 'm(b-bbar)=',tmpmtt
c
c j-j inv mass
      ptmp(4)=p(4,5)+p(4,6)
      tmpmtt=ptmp(4)**2
      do i=1,3
         ptmp(i)=p(i,5)+p(i,6)
         tmpmtt=tmpmtt-ptmp(i)**2
      enddo
      tmpmtt=sqrt(tmpmtt)
      write(niosta,*) 'm(jj)=',tmpmtt
      end
c
C***********************************************************************
       subroutine setinteraction77
C***********************************************************************
       implicit none
       integer inpint(1000)
       common/initinter/inpint  
c
         data inpint/ 6,
***     >               11,1,1,1,1,  11,1,2,1,2,  11,3,2,3,2,   2,1,1,1,2,  ! guu, gdd, gbb, w+ud
     >           11,1,1,1,1,    11,3,1,3,1,  11,1,2,1,2,  11,3,2,3,2,    ! guu, gtt, gdd, gbb
     >           11,2,1,2,1,  11,2,2,2,2,                                !gcc,gss
***     >                3,1,4,1,3,  11,2,1,2,1,                            ! w-en, gcc
     >                0,                                                 ! N of yukawa
     >                2,                                                 ! N self-gauge
     >                11,11,11,  12,11,11,                               ! ggg Auxgg
     >                961*-100/
***     >                971*-100/
c
       return
       end
c
c***** alpha-specific routines
*
      subroutine mom(lflag,np,pt0,pt1,eta0,xmr,en0,pz0,
     .               pm,wgt,ranram,lw)
      implicit none
      real*8 pi,wgt,rpr,phr,ran0,wj,en0,pz0,prx,pry
      real*8 etai1,etai2,etaj1,etaj2,xmri2,bxmri2,wjeta
      real*8 cn,wjaci,wjacj,phi,ranj,ref,alimp,alimm
      real*8 p0m,p0p 
      real*8 phj,alim,phip,phim,etap,etam
      integer npar,n,nri,nr,iter,i,j,k,lw,icont,lflag
      parameter (pi= 3.14159265358979323846264338327950d0)
      parameter (npar= 20)
      real*8 pt0(npar),pt1(npar),eta0(npar),
     .       xmr(npar),bxmr(npar)
      real*8 pm(0:4,npar),pt(npar),eta(npar)
      integer jp(npar),init,np
      real*8 p0r(npar),p1r(npar),qcut(npar),pcut(npar),bpt0(npar)
      real*8 bet0(npar),ph(npar)
      real*8 ranram(111)
      real*8 en,pz,ga,de,p,q,a,b,cut,ausp,ausm,v,vmr,vpr,v2
      real*8 det,rx,x1m,x2m,pmod,qmod,den
      real*8 asinh,aus,betp,alpp
      real*8 etalim,csi
      data init/0/
      save
*
      if (init.eq.0) then
        init= 1
        do iter= 1,np
          jp(iter)= iter
        enddo
        n  = np-1
        nri= 3
        etalim= 100.d0
*
        do iter= 1,n
          bpt0(jp(iter))= 0.d0
          pcut(jp(iter))= 0.d0
          qcut(jp(iter))= 0.d0
          bet0(jp(iter))= etalim
          do k= iter+1,np
              bpt0(jp(iter))= bpt0(jp(iter))+pt0(jp(k))
              pcut(jp(iter))= max(pcut(jp(iter)),pt0(jp(k)))
              bet0(jp(iter))= min(bet0(jp(iter)),eta0(jp(k)))
          enddo
        enddo
        qcut(jp(n))= pt0(jp(np))
*
      endif
*
      do iter= 1,n
        bxmr(jp(iter))= 0.d0
        do k= iter+1,np
            bxmr(jp(iter))= bxmr(jp(iter))+xmr(jp(k))
        enddo
      enddo
*
      lw  = 0
      nr  = nri
      rpr = 0.d0
      wgt = 1.d0
      en  = en0
      pz  = pz0
      if (lflag.eq.0) then
        phr = 0.d0
      elseif (lflag.eq.1) then
        prx = 0.d0
        pry = 0.d0
      else
        goto 101
      endif
*
      do iter= 1,n
        i= jp(iter)
        j= jp(iter+1)
*
c-      Generation of transverse momenta when lflag= 0
*
        ga = en+pz
        de = en-pz 
*
        if (ga.lt.0.d0) goto 100
        if (de.lt.0.d0) goto 100
        v  = sqrt(ga*de)
*
c-      Protection:
*
        vmr= v-rpr
        if (vmr.lt.0.d0) then
          vmr= dabs(vmr)
          v  = vmr+2.d0*rpr
          de = ga/v/v
        endif
        vpr= vmr+2.d0*rpr
        v2 = v*v
*
        if (ga.ge.de) then
          alpp= 2.d0*ga*exp(-eta0(i))
        else
          alpp= 2.d0*de*exp(-eta0(i))
        endif
        betp= ((ga+de)*cosh(eta0(i))-dabs(ga-de)*sinh(eta0(i)))
        csi = min(alpp,betp)
*
        if (lflag.eq.1) then
          pt(i)= sqrt(pm(1,i)**2+pm(2,i)**2)
          if (iter.ne.1) then
            phi= ((pm(1,i)*prx+pm(2,i)*pry)/pt(i)/rpr)
            if (phi.gt. 1.d0) phi= 1.d0
            if (phi.lt.-1.d0) phi=-1.d0
            phi= acos(phi)
          endif
          prx  = prx-pm(1,i)
          pry  = pry-pm(2,i)
          pt(j)= sqrt(prx*prx+pry*pry)
*
          eta(i)= -log(tan(0.5d0*atan2(pt(i),pm(3,i))))
        endif
*
        xmri2 = xmr(i)**2
        bxmri2= bxmr(i)**2
*
        aus   =  vpr*vmr 
        if(abs(aus).lt.1.d-30) then
           write(*,*)'fix in MOM, aus too small', aus
           goto 100
        endif
        ausp  = (aus+xmri2-bxmri2)
        det   = (ausp**2-4.d0*xmri2*aus)
        if (det.lt.0.d0) then
          det= abs(det)
          if (det.gt.1.d-12) goto 100
        endif
        det   = sqrt(det)
        if (ausp.gt.0.d0) then
           p0p= (rpr*ausp+v*det)/2.d0/aus
           p0m=-(ausp**2-4.d0*v2*xmri2)/4.d0/aus/p0p
        else
           p0m= (rpr*ausp-v*det)/2.d0/aus
           p0p=-(ausp**2-4.d0*v2*xmri2)/4.d0/aus/p0m
        endif
        p0r(i)= max(pt0(i),p0m,qcut(i)-rpr)
        p1r(i)= min(p0p,v-bpt0(i))
        if (csi.gt.2.d0*v) then
           p1r(i)= min(p1r(i),(vpr*vmr)/(csi-2.d0*rpr))
        endif
        cn= 1.d0
        if (p0r(i).eq.0.d0) cn= 0.98d0
        if (p0r(i).gt.p1r(i)) goto 100
        call ppeaka(lflag,max(pt0(i),xmr(i),1.d-2),
     .              cn,p0r(i),p1r(i),pt(i),wjaci,ranram(nr),lw)
        nr= nr+1
        if (lw.eq.1) goto 100
*
        if (lflag.eq.0) call rans(ran0)
*
        if (iter.eq.1) then
          if (lflag.eq.0) then
            phi= 2.d0*pi*ran0
            pt(j)= pt(i)
            phj= phi+pi
          endif
          wgt= wgt*2.d0*pi*pt(i)*wjaci
        else
          if (lflag.eq.0) then
            if (ran0.le.0.5d0) then
              ref = 1.d0
            else
              ref =-1.d0
            endif
            call rans(ranj)
          endif
*
          p0r(j)= max(dabs(pt(i)-rpr),qcut(i),pt(i)-v+2.d0*pcut(i))
          det= (v-sqrt(pt(i)**2+xmri2))**2-bxmri2
          if (det.lt.0.d0) then
             det= abs(det)
             if (det.gt.1.d-12) goto 100
          endif
          det   = sqrt(det)
          p1r(j)= min(pt(i)+rpr,det)
          if (csi.gt.2.d0*v) then
            p1r(j)= min(p1r(j),sqrt(pt(i)*(pt(i)-csi)+v2))
          endif
          if (p0r(j).gt.p1r(j)) goto 100
*
          alimp= (rpr**2+pt(i)**2-p1r(j)**2)/2.d0/rpr/pt(i)
          alimm= (rpr**2+pt(i)**2-p0r(j)**2)/2.d0/rpr/pt(i)
*
c-        Protections:
*
          if (alimp.gt. 1.d0) alimp= 1.d0
          if (alimp.lt.-1.d0) alimp=-1.d0 
          if (alimm.gt. 1.d0) alimm= 1.d0
          if (alimm.lt.-1.d0) alimm=-1.d0
*          
          phip= acos(alimp)
          phim= acos(alimm)
          if (phim.gt.phip) goto 100
*
          call fflat(lflag,phi,phim,phip,wjacj,ranj,lw)
          if (lw.eq.1) goto 100
          wgt= wgt*wjaci*wjacj*pt(i)*2.d0
*
          if (lflag.eq.0) then
            phi= phi*ref
            pt(j)= dabs(rpr**2+pt(i)**2-2.d0*rpr*pt(i)*cos(phi))
            pt(j)= sqrt(pt(j))
*
c-          Protections:
*
            if (rpr.eq.0.d0)   rpr  = 1.d-25
            if (pt(j).eq.0.d0) pt(j)= 1.d-25
            alim= (rpr**2+pt(j)**2-pt(i)**2)/2.d0/rpr/pt(j)
*
c-          Protections:
*
            if (alim.gt. 1.d0) alim= 1.d0
            if (alim.lt.-1.d0) alim=-1.d0
*
            phj=-acos(alim)*ref
*
          endif
        endif  
*
        if (lflag.eq.0) then
          phi= phi+phr
          phj= phj+phr
          pm(1,i)= pt(i)*sin(phi)
          pm(2,i)= pt(i)*cos(phi)
          pm(1,j)= pt(j)*sin(phj)
          pm(2,j)= pt(j)*cos(phj)
          ph(i)  = phi
          if (iter.eq.n) ph(j)  = phj
*
          phr= phj
        else
          ph(i)= pm(2,i)/pt(i)
          if (ph(i).ge.0.d0) then
            ph(i)= min(1.d0,ph(i))
          else
            ph(i)= max(-1.d0,ph(i))
          endif
          ph(i)= acos(ph(i))
          if (pm(1,i).lt.0.d0) ph(i)= 2.d0*pi-ph(i)
          if (iter.eq.n) then
             ph(j)= pm(2,j)/pt(j)
            if (ph(j).ge.0.d0) then
              ph(j)= min(1.d0,ph(j))
            else
              ph(j)= max(-1.d0,ph(j))
            endif
            ph(j)= acos(ph(j))
            if (pm(1,j).lt.0.d0) ph(j)= 2.d0*pi-ph(j)
          endif
        endif
*
        rpr= pt(j)
*
        p  = pt(i)
        q  = pt(j)
        a  = 1.d0+xmri2/p/p
        b  = max(1.d0+bxmri2/q/q,1.d0+4.d0/q/q*pcut(i)*
     .          (pcut(i)-q),(bpt0(i)/q)**2)
*
        cut= (sqrt(a)*p+sqrt(b)*q)**2-v2
        if (cut.gt.1.d-12) goto 100
        ausm= a*p*p-b*q*q
        det = (v2+ausm)**2-4.d0*v2*a*p*p
        if (det.lt.0.d0) then
           det= abs(det)
           if (det.gt.1.d-12) goto 100
        endif
        det = sqrt(det)
        rx  = v2+ausm
        x1m= (rx+det)/2.d0/p/de  
        x2m= ga*a/de/x1m                  ! x2m x1m
        etai1= asinh(0.5d0*(x1m-a/x1m))
        etai2= asinh(0.5d0*(x2m-a/x2m))
        if (iter.le.n-1) then
*
c-        Generation of etas when lflag= 0
*
          etam= max(etai2,-eta0(i))        
          etap= min(etai1, eta0(i))        
          if (etam.ge.etap) goto 100
*
comment
c          call ttriangle(lflag,etai2,etai1,etam,etap,eta(i),
c     .                  wjeta,ranram(nr),lw)
          if (iter.ge.2) then
            k    = jp(iter-1)
            alimp= 4.d0*(xmr(k)+xmr(i))**2   
            call spk(lflag,alimp,0.98d0,pt(k),pt(i),ph(k),ph(i),
     .               xmr(k),xmr(i),eta(k),
     .               etam,etap,eta(i),wjeta,ranram(nr),lw)
c            call etagen(lflag,ph(k),ph(i),eta(k),
c     .                  etam,etap,eta(i),wjeta,ranram(nr),lw)
          else
            call ttriangle(lflag,etai2,etai1,etam,etap,eta(i),
     .                     wjeta,ranram(nr),lw)
          endif
comment
          nr= nr+1
          if (lw.eq.1) goto 100
          wgt   = wgt*wjeta
          if (lflag.eq.0) then
            pm(0,i)= sqrt(pt(i)**2*cosh(eta(i))**2+xmri2)   
            pm(3,i)= pt(i)*sinh(eta(i))                         
            pm(0,j)= en-pm(0,i)
            pm(0,j)= dabs(pm(0,j))
            pm(3,j)= pz-pm(3,i)
          endif
          en     = en-pm(0,i)
          pz     = pz-pm(3,i)
          wgt    = wgt*sqrt(1.d0-(xmr(i)/pm(0,i))**2)
        else
*
c-        Cuts on the last 2 etas
*
          etaj1= asinh((pz-pt(i)*sinh(etai1))/q)
          etaj2= asinh((pz-pt(i)*sinh(etai2))/q)
*
c-        Look whether 0,1 or 2 solutions may contribute
*
          icont= 0                        
          if ((dabs(etai1).lt.eta0(i)).and.
     .        (dabs(etaj1).lt.eta0(j))) icont= icont+1
          if ((dabs(etai2).lt.eta0(i)).and.
     .        (dabs(etaj2).lt.eta0(j))) icont= icont+2
*
          if (icont.eq.0) then
             goto 100
          elseif (icont.eq.1) then
             if (lflag.eq.0) then
               eta(i)= etai1
               eta(j)= etaj1
             endif
             wj= 1.d0
          elseif (icont.eq.2) then
             if (lflag.eq.0) then
               eta(i)= etai2
               eta(j)= etaj2
             endif
             wj= 1.d0
          else
             if (lflag.eq.0) then
               call rans(ran0)
               if (ran0.lt.0.5d0) then
                 eta(i)= etai1
                 eta(j)= etaj1
               else
                 eta(i)= etai2
                 eta(j)= etaj2
               endif
             endif
             wj= 2.d0
          endif
          if (lflag.eq.0) then
             pm(0,i)= sqrt(pt(i)**2*cosh(eta(i))**2+xmri2)   
             pm(3,i)= pt(i)*sinh(eta(i))                         
             pm(0,j)= en-pm(0,i)
             pm(0,j)= dabs(pm(0,j))
             pm(3,j)= pz-pm(3,i)                                  
          elseif (lflag.eq.1) then
             en     = en-pm(0,i)-pm(0,j)
             pz     = pz-pm(3,i)-pm(3,j)
             eta(j) =-log(tan(0.5d0*atan2(pt(j),pm(3,j))))
          else
             goto 101
          endif
          pmod= pm(0,i)**2-xmri2
          qmod= pm(0,j)**2-xmr(j)*xmr(j)
          if (pmod.lt.0.d0) goto 100    
          if (qmod.lt.0.d0) goto 100
          pmod= sqrt(pmod)
          qmod= sqrt(qmod)
          den = dabs(pm(0,j)*pm(3,i)-pm(0,i)*pm(3,j))
          if(abs(den).lt.1.d-20) then
             call raisenumexcep(1)
             goto 100
          endif
          wgt = wgt*wj*
     .          pmod*qmod/p/q/cosh(eta(j))/cosh(eta(i))/den
        endif
      enddo     
      wgt= wgt/2.d0**np      
*
      return
 100  wgt= 0.d0
      lw= 1
      return
 101  print*,'ERROR IN SUBROUTINE MOM'
      stop
      end
*
      subroutine momr(lflag,np,xmr,en0,pz0,
     .                pm,wgt,lw)
      implicit none
      integer lflag,np,npar,lw,j
      parameter (npar= 20)
      real*8 xmr(npar),en0,pz0,pm(0:4,npar),wgt,dj
      real*8 x1,x2,bvel,gvel,rootshr,pr(4,npar)
*
      lw= 0
      if (lflag.eq.1) then
        en0= 0.d0
        pz0= 0.d0
        do j= 1,np
          en0= en0+pm(0,j)          
          pz0= pz0+pm(3,j)
        enddo
      endif
      x1  = dabs(en0+pz0)
      x2  = dabs(en0-pz0)
      bvel= (x1-x2)/(x1+x2)
      gvel= 1.d0/sqrt(1.d0-bvel*bvel)
      rootshr= sqrt(x1*x2)
*
      if (lflag.eq.1) then
        do j= 1,np
          pr(4,j)= gvel*(pm(0,j)-bvel*pm(3,j))
          pr(1,j)= pm(1,j)
          pr(2,j)= pm(2,j)
          pr(3,j)= gvel*(pm(3,j)-bvel*pm(0,j))
        enddo 
      endif
      call rambo(lflag,np,rootshr,xmr,pr,dj)
      wgt= 1.d0/dj
      if (lflag.eq.1) return
*
      do j= 1,np
        pm(0,j)= gvel*(pr(4,j)+bvel*pr(3,j))
        pm(1,j)= pr(1,j)
        pm(2,j)= pr(2,j)
        pm(3,j)= gvel*(pr(3,j)+bvel*pr(4,j))
      enddo
*
      return
      end
*

