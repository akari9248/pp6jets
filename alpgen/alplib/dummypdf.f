      subroutine lhapdfpar(ndns,ih,xlam,sche,nl,iret)
      implicit none
      integer ndns,ih,nl,iret
      CHARACTER * 2 sche
      real * 8 xlam

      print *, "Compiled without LHAPDF"
      stop

      end

*********************************************

      subroutine lhapdfif(ndns,ih,xmu2,x,fx,nf)
      implicit none
      integer ndns,ih,nf
      real xmu2,x,fx(-5:5)

      print *, "Compiled without LHAPDF"
      stop

      end

*********************************************

      subroutine lhapdfconv(nin,nout,type)
      implicit none
      integer nin, nout
      character*25 type

      print *, "Compiled without LHAPDF"
      stop

      end

*********************************************

      subroutine lhaprntsf(iunit)

      print *, "Compiled without LHAPDF"
      stop

      end

*********************************************

      subroutine alphaspdf(q2)
      implicit none
      real * 8 q2

      print *, "Compiled without LHAPDF"
      stop

      end


