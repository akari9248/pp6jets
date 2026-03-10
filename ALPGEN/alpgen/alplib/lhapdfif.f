      subroutine lhapdfif(ndns,ih,xmu2,x,fx,nf)
c Interface to lhapdf package.
      implicit none
      integer ndns,ih,nf
      real xmu2,x,fx(-5:5)
      real * 8 fxp(-6:6)
      integer j,jinit
      real * 8 tmp
      character * 20 parm(20)
      double precision val(20)
      data jinit/0/
      if(jinit.eq.0) then
c        call Initpdf(mod(ndns,100))
        parm(1)='DEFAULT'
        val(1)=ndns
        call pdfset(parm,val)
        jinit=1
      endif
      call evolvePDF(dble(x),dble(sqrt(xmu2)),fxp)
c     pftopdg returns density times x
      do j=-nf,nf
         fx(j)=real(fxp(j))
         fx(j)=fx(j)/x
      enddo
      if (ih.eq.-1) then
         do j=1,5
            tmp=fx(j)
            fx(j)=fx(-j)
            fx(-j)=tmp
         enddo
      endif
c to mantain alppdf standards (u <-> d)
      tmp=fx(1)
      fx(1)=fx(2)
      fx(2)=tmp
      tmp=fx(-1)
      fx(-1)=fx(-2)
      fx(-2)=tmp
      end

      subroutine lhapdfpar(ndns,ih,xlam,sche,nl,iret)
      implicit none 
      integer ndns,ih,iret
      real * 8 xlam
      character * 2 sche
      character * 20 parm(20)
      double precision val(20)
      integer nl
      real * 8 qcdl4,qcdl5
      double precision alfas,as
      common/w50512/qcdl4,qcdl5
      character *64 name
      
      parm(1)='DEFAULT'
      val(1)=ndns
c      name='CT14nlo.LHgrid'
c      name='cteq66.LHgrid'
      call pdfset(parm,val)
c      call InitPDFsetByName(name)
      call Initpdf(mod(ndns,100))
c      call setlhaparm('EXTRAPOLATE')
c lhapdf order of as is 1 less then alpgen
      call GetOrderAs(nl)
c     the following call requires LHAPDF v6.1.6 or higher
      call GetLam5(0,qcdl5)
      xlam = qcdl5
      sche='MS'
      iret=0
      end


      subroutine lhapdfconv(nin,nout,type)
      implicit none
      integer nin, nout
      character*25 type
      nout=nin
      write (type, "(A7,I10)") "LHAPDF ", nin 
      end

      subroutine lhaprntsf(iunit)
      end


