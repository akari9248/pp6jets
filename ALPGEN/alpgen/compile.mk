# DEFINE FORTRAN COMPILATION DEFAULTS
#FFF = gfortran  -fno-automatic -g -fcheck=bounds -Wall
#FFF = gfortran  -fno-automatic -g -fbacktrace -ffpe-trap=zero,overflow,underflow
FFF = gfortran  -fno-automatic -ffast-math
FF90 = gfortran  -fno-automatic
#FFF = gfortran  -fno-automatic -ffast-math -fbounds-check -march=x86-64 -O2 
#FF90 = gfortran  -fno-automatic -ffast-math  -fbounds-check -arch=x86-64 -O2 
#FFF = gfortran -fno-automatic -g -ffpe-trap=zero,invalid,overflow
#FF90 = gfortran -fno-automatic -g -fbounds-check -ffpe-trap=zero,invalid
ifeq ($(shell uname),AIX)
 FFF = xlf -O4
 FF90 = xlf90 -O4 -qsuffix=f=f90:o=o90
endif

ifeq ($(shell uname),Linux)
FFF = gfortran  -fno-automatic -ffast-math -O2 
FF90 = gfortran  -fno-automatic -ffast-math -O2 
endif

ifeq ($(shell uname),Darwin)
  ifeq  ($(shell uname -m),i386) 
#     FFF = ifort -noautomatic  -warn  
#     FF90 = ifort -noautomatic  -warn
     FFF = gfortran  -fno-automatic
    FF90 = gfortran -fno-automatic
  endif
endif
# Unix-ALPHA fortran
# FFF = f90 -fast
# FF90 = f90 -fast  

