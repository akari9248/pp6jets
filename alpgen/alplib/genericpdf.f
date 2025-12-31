
      subroutine genericpdf(ndns,ih,xmu2,x,fx,nf)
c Interface to lhapdf and alppdf package.
      implicit none
      integer ndns,ih,nf
      real xmu2,x,fx(-5:5)

      character*6 whichpdf
      common/mypdf/whichpdf

      if (whichpdf.eq.'lhapdf') then
          call lhapdfif(ndns,ih,xmu2,x,fx,nf)
      elseif (whichpdf.eq.'alppdf') then
          call mlmpdf(ndns,ih,xmu2,x,fx,nf)
      else
          print*, " ******** NO PDFs ******** "
          stop
      endif

      end

*********************************************

      subroutine pdfpar(ndns,ih,xlam,sche,nl,iret)
      implicit none
      integer ndns,ih
      character * 2 sche
      real * 8 xlam
      integer nl,iret

      character*6 whichpdf
      common/mypdf/whichpdf

      if (whichpdf.eq.'lhapdf') then
          call lhapdfpar(ndns,ih,xlam,sche,nl,iret)
      elseif (whichpdf.eq.'alppdf') then
          call alppdfpar(ndns,ih,xlam,sche,nl,iret)
      else
          print*, " ******** NO PDFs ******** "
          stop
      endif

      end

*********************************************

      subroutine pdfconv(nin,nout,type)
      implicit none
      character*25 type
      integer nin,nout

      character*6 whichpdf
      common/mypdf/whichpdf

      if (nin.gt.1000) then
          whichpdf = 'lhapdf'
      else
          whichpdf = 'alppdf'
      endif

      if (whichpdf.eq.'lhapdf') then
          call lhapdfconv(nin,nout,type)
      elseif (whichpdf.eq.'alppdf') then
          call alppdfconv(nin,nout,type)
      else
          print*, " ******** NO PDFs ******** "
          stop
      endif

      end

*********************************************

      subroutine prntsf(iunit)
      implicit none 
      integer iunit
 
      character*6 whichpdf
      common/mypdf/whichpdf
     
      if (whichpdf.eq.'lhapdf') then
          call lhaprntsf(iunit)
      elseif (whichpdf.eq.'alppdf') then
          call alpprntsf(iunit)
c      else
c          print*, " ******** NO PDFs ******** "
c          stop
      endif

      end


*********************************************

      function alfas(q2,xlam,nloop,inf)
      implicit none
      real*8 alfas
      real*8 q2,xlam
      integer inf,nloop

      real*8 alphasPDF
      external alphasPDF

      real*8 alpalfas

      character*6 whichpdf
      common/mypdf/whichpdf

      if (whichpdf.eq.'lhapdf') then
          alfas=alphasPDF(sqrt(q2))
      elseif (whichpdf.eq.'alppdf') then
          alfas=alpalfas(q2,xlam,nloop,inf)
      else
          print*, " ******** NO PDFs ******** "
          stop
      endif

      end

*********************************************

      function alfas_clu(q2,xlam,nloop,inf)
      implicit none
      real*8 alfas_clu
      real*8 q2,xlam
      integer inf,nloop
      real*8 alphasPDF
      external alphasPDF
      real*8 alpalfas_clu
      character*6 whichpdf
      common/mypdf/whichpdf
      integer iasclu,iopt
      common/iascl/iasclu
      data iopt/0/
c
      if(iopt.eq.1) then
        alfas_clu=alpalfas_clu(q2,xlam,nloop,inf)
      elseif(iopt.eq.2) then
        alfas_clu=alphasPDF(sqrt(q2))
      elseif(iopt.eq.0) then
        if(iasclu.eq.1.or.whichpdf.eq.'alppdf') then
          iopt=1
          alfas_clu=alpalfas_clu(q2,xlam,nloop,inf)
        elseif(whichpdf.eq.'lhapdf') then
          iopt=2
          alfas_clu=alphasPDF(sqrt(q2))
        else
          print*, " ******** alfas_clu not defined ******** "
          stop
        endif
      endif
      end


