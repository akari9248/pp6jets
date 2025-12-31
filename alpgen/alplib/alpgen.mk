# DO NOT EDIT FROM HERE ON:
#
# DEFINE DIRECTORY AND FILE ALIASES
alp= ../alplib
her= ../herlib
prclib= ../$(prc)lib
prcusr=.
prcfile=$(prclib)/$(prc)
ifeq ($(strip $(usrfile)),)
  execfile=$(prc)gen
  usrfile=$(prcusr)/$(prc)usr
else
  execfile=$(usrfile)
endif

## Choosing PDFlib: native,lhapdf
## LHAPDF package has to be installed separately
##
## Assign the directory where the file 'lhapdf-config' is located
#CONFIG_FILE_DIR=/usr/local/bin
#CONFIG_FILE_DIR=/Users/mlm/herwig/bin
#CONFIG_FILE_DIR=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.5.1-389ee/x86_64-centos7-gcc11-opt/bin
#CONFIG_FILE_DIR=/cvmfs/sft.cern.ch/lcg/releases/MCGenerators/lhapdf/6.5.3-3fa11/x86_64-centos7-gcc11-opt/bin/
CONFIG_FILE_DIR=/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-pafccj3/bin
LHAPDF_CONFIG=$(CONFIG_FILE_DIR)/lhapdf-config

ifeq ("$(shell which $(LHAPDF_CONFIG))", "")
  $(info ******* LHAPDF NOT FOUND IN ********)
  $(info $(CONFIG_FILE_DIR))
  $(info ******* WILL INSTALL NATIVE PDFS ********)
  PDFlib = native
else
  $(info ******* COMPILING WITH LHAPDF *******)
  LHAPDFDIR=$(shell $(LHAPDF_CONFIG) --prefix)
  LHAPATH= $(shell $(LHAPDF_CONFIG) --datadir)  
  LIBSLHAPDF= $(shell $(LHAPDF_CONFIG) --libdir)  
  $(info LHAPDF directory is $(LHAPDFDIR))
  $(info LHAPATH is $(LHAPATH))
  $(info LHAPDF Libraries path is $(LIBSLHAPDF))
  $(info C++ Libraries path is $(lstdc++))
  $(info For my ref: List of PDF codes in /opt/homebrew/Cellar/lhapdf/6.5.3/share/LHAPDF/PDFcodes.txt)
  $(info *******)
  PDFlib = lhapdf
endif

ifeq ("$(PDFlib)","lhapdf")
  PDFPACK=$(alp)/alppdf.f $(alp)/lhapdfif.f $(alp)/genericpdf.f 
  PDFLIB= -L$(LIBSLHAPDF) -lLHAPDF
#  LIBSLHAPDF= -Wl,-rpath, /usr/local/lib  -L/usr/local/lib -lLHAPDF
  LIBS+=$(PDFLIB)
else
  PDFPACK=$(alp)/genericpdf.f $(alp)/alppdf.f $(alp)/dummypdf.f
endif

# DEFINE FILE GROUPS:
# Files used for the parton-level event generation:
ALPGEN= $(alp)/alpgen.f $(alp)/alputi.f \
	$(alp)/Aint90.f $(PDFPACK)
ALP77=  $(alp)/Aint.f $(alp)/Asu3.f
FJCORE= $(alp)/fjfort.cc $(alp)/fjcore.cc $(alp)/alpfj.f

#  $(alp)/Aint.f $(alp)/Asu3.f 
PARTON= $(prcfile).f $(usrfile).f $(prclib)/Aproc.f $(ALPGEN) $(ALP77)

# Include files
INC=  $(prclib)/$(prc).inc $(alp)/alpgen.inc

# include files' dependencies
$(PARTON): $(INC)
$(FJCORE): $(alp)/fjcore.hh
# object files
OBJ=$(PARTON:.f=.o) $(PARTON90:.f90=.o) $(FJCORE: .cc=.o)


# compilation
%.o: %.f $(INC)
	$(FFF) -c -o $*.o $*.f 
CF=gfortran
%.o: %.cc
	$(CF) -c -o $*.o $*.cc 

#$(alp)/A90.o90: $(alp)/A90.f90 $(INC)
#	$(FFF) -c -o $(alp)/A90.o90 $(alp)/A90.f90
#	mv *.mod $(prclib)
$(prclib)/ini_$(prc).o90:  $(alp)/A90.f90 $(prclib)/ini_$(prc).f90 $(INC)
	$(FFF) -c -o $(alp)/A90.o90 $(alp)/A90.f90
	mv *.mod $(prclib)
	$(FFF) -c -o $(prclib)/ini_$(prc).o90 $(prclib)/ini_$(prc).f90

gen90: 
	make gen

genfj: 
	make genfj90

# fortran90 version: now default, can be called with either gen or gen90
gen: $(usrfile).o $(prcfile).o $(prclib)/$(prc).inc\
	$(alp)/alpgen.o $(alp)/alputi.o $(PDFPACK:.f=.o) \
	$(alp)/Aint90.o $(prclib)/Aproc.o  \
        $(prclib)/ini_$(prc).o90 $(alp)/alpgen.inc 
	$(FF90) -o $(execfile) $(usrfile).o $(prcfile).o \
	$(alp)/alpgen.o $(alp)/alputi.o $(PDFPACK:.f=.o) \
	$(alp)/Aint90.o $(prclib)/Aproc.o $(alp)/A90.o90 \
	$(prclib)/ini_$(prc).o90 $(LIBS)
# -lstdc++
lhapdfex: fexample1.o 
	$(FF90) -o fexample fexample1.o $(LIBS) -lstdc++

# Link to fastjet core routines
genfj90: $(usrfile).o $(prcfile).o $(prclib)/$(prc).inc\
	$(alp)/alpgen.o $(alp)/alputi.o $(PDFPACK:.f=.o) \
	$(alp)/Aint90.o $(prclib)/Aproc.o  \
        $(prclib)/ini_$(prc).o90 $(alp)/alpgen.inc \
	$(alp)/fjfort.o $(alp)/fjcore.o $(alp)/alpfj.o 
	$(FF90) -o $(execfile) $(usrfile).o $(prcfile).o \
	$(alp)/alpgen.o $(alp)/alputi.o $(PDFPACK:.f=.o) \
	$(alp)/Aint90.o $(prclib)/Aproc.o $(alp)/A90.o90 \
	$(prclib)/ini_$(prc).o90 $(LIBS) \
	$(alp)/fjfort.o $(alp)/fjcore.o $(alp)/alpfj.o -lstdc++

# fortran77 version -- for debugging/validation ony
gen77: $(OBJ)
	$(FFF) -o $(execfile) $(usrfile).o $(prcfile).o \
	$(alp)/alpgen.o $(alp)/alputi.o $(PDFPACK:.f=.o) \
	$(alp)/Asu3.o $(alp)/Aint.o $(prclib)/Aproc.o 

# DIRECTORY CLEANUP UTILITIES:
#
# remove object files only
cleanobj:
	-rm $(PARTON:.f=.o) $(PARTON90:.f90=.o) $(prcusr)/../*/*.o90*

# remove object files, etc
cleanall:
	-rm $(OBJ) $(prcusr)/fort.* $(prcusr)/*.top $(prcusr)/*.par \
	$(prcusr)/*.wgt $(prcusr)/*.unw $(prcusr)/*.mon \
	$(prcusr)/*.stat $(prcusr)/../*/*.o90*

