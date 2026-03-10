# This makefile allows to pack and unpack the alpgen libraries
# 
targz = tar -zcvf
# define directory tree structure:
#
# libraries:
alp= alplib
her= herlib
pyt= pylib
val= validation

PROCLIST= wqq zqq wjet zjet vbjet Njet 2Q 4Q QQh wcjet phjet hjet top \
wphjet wphqq 2Qph zphjet zphqq hvyN 

LIBS= $(addsuffix lib,  $(PROCLIST)) 
WORK= $(addsuffix work, $(PROCLIST)) 

find_77=$(wildcard $(addsuffix lib,  $(dir))/$(dir).f) \
	$(wildcard $(addsuffix lib,  $(dir))/Aproc.f) \
	$(wildcard $(addsuffix work,  $(dir))/$(dir)usr.f) 
F77files=$(foreach dir, $(PROCLIST), $(find_77))

find_inc=$(wildcard $(addsuffix lib,  $(dir))/$(dir).inc) \
	$(wildcard $(addsuffix work,  $(dir))/$(dir).inc)  \
	$(wildcard $(addsuffix lib,  $(dir))/alpgen.inc)  \
	$(wildcard $(addsuffix work,  $(dir))/alpgen.inc) 
INCfiles=$(foreach dir, $(PROCLIST), $(find_inc))

find_make=$(wildcard $(addsuffix work,  $(dir))/Makefile)
MAKEfiles= $(foreach dir, $(PROCLIST), $(find_make)) 

EXTRAS= Makefile compile.mk $(alp)/alpgen.mk compare 

find_input=$(wildcard $(addsuffix work,  $(dir))/input)
INPUTfiles= $(foreach dir, $(PROCLIST), $(find_input))

find_pdf=$(wildcard $(addsuffix work,  $(dir))/pdflnk)
PDFlnk=$(foreach dir, $(PROCLIST), $(find_pdf))

find_lis=$(wildcard $(addsuffix lib,  $(dir))/par.list)
LISfiles=$(foreach dir, $(PROCLIST), $(find_lis)) prc.list

find_90=$(wildcard $(addsuffix lib,  $(dir))/ini_$(dir).f90)
F90files= $(foreach dir, $(PROCLIST), $(find_90))

# HARD PROCESSES' FILES
PROCS77=$(F77files) $(INCfiles) $(MAKEfiles) $(INPUTfiles) \
        $(PDFlnk) $(LISfiles)
# F90 COMPONENTS FOR THE HARD PROCESSES
PROCS90=$(F90files)

# ALPGEN F77 SOURCE
ALP77=$(alp)/Aproc.f $(alp)/alpgen.f \
     $(alp)/alppdf.f $(alp)/genericpdf.f $(alp)/dummypdf.f $(alp)/lhapdfif.f\
	$(alp)/alputi.f $(alp)/alpsho.f \
     $(alp)/alpgen.inc $(alp)/alpsho.inc
# ALPGEN F90 COMPONENTS
ALP90= $(alp)/Aint90.f $(alp)/A90.f90 
# FASTJET COMPONENTS
FJCORE= $(alp)/fjfort.cc $(alp)/fjcore.cc $(alp)/alpfj.f $(alp)/fjcore.hh

# PDF DATA SETS
PDF=  $(alp)/pdfdat/ctq45 $(alp)/pdfdat/ctq61 $(alp)/pdfdat/ctq66 \
$(alp)/pdfdat/mrst $(alp)/pdfdat/mstw  $(alp)/pdfdat/pdflnk

# HERWIG files
HWIG= $(her)/HERWIG65.INC $(her)/Makefile $(her)/atoher.f \
    $(her)/getjet.f $(her)/herwig6521.inc $(her)/hwuser.f \
    $(her)/herwig6521.f $(her)/pdfdummy.f $(her)/alpsho.inc
# PYTHIA files
PYT=$(pyt)/Makefile $(pyt)/pythia-6.4.24.f \
    $(pyt)/pyuser.f $(pyt)/atopyt.f $(pyt)/alpsho.inc

# DOCUMENTATION FILES
DOCS= DOCS/alpdoc.tex DOCS/alpdoc.pdf DOCS/JHEP3.cls DOCS/History.txt
#PS=$(patsubst%.tex,%.ps,$(TEXFILE))
#%.ps: %.tex
#	cp DOCS/JHEP3.cls .
#	latex $*.tex 
#	latex $*.tex 
#	dvips $*.dvi

VALfiles= $(val)/REF/REF*


pack: 
	cd DOCS;  pdflatex alpdoc ; pdflatex alpdoc ;
	$(targz) alpgen.tar.gz $(ALP77) $(ALP90) $(PROCS77) $(PROCS90)	\
	$(HWIG) $(PYT) $(PDF) $(val)/validate \
	$(EXTRAS) $(DOCS) $(FJCORE) $(val)/REF
# $(VALfiles) 

v2last: 
	cd DOCS;  pdflatex alpdoc ; pdflatex alpdoc ;
	$(targz) v2last.tar.gz $(ALP77) $(ALP90) $(PROCS77) $(PROCS90)	\
	$(HWIG) $(PYT) $(PDF) $(alp)/alpsho_fj.f $(val)/validate \
	$(EXTRAS) $(DOCS) $(FJCORE)

packfull: 
	cd DOCS;  pdflatex alpdoc ; pdflatex alpdoc ;
	$(targz) alpgen.tar.gz $(ALP77) $(ALP90) $(PROCS77) $(PROCS90)	\
	$(HWIG) $(PYT) $(PDF) $(val)/validate $(val)/REF $(VALfiles) \
	$(EXTRAS) $(DOCS)  $(FJCORE) \
	Transition.txt Validation.txt flowchart.key \
	flowchart.pdf
package:
	cd validation; ./validate "ref" ;
	cd DOCS; pdflatex alpdoc ; pdflatex alpdoc ;
	$(targz) alpgen.tar.gz $(ALP77) $(ALP90) $(PROCS77) $(PROCS90) \
	$(HWIG) $(PYT) $(PDF) $(val)/validate $(VALfiles) \
	$(EXTRAS) $(DOCS) $(FJCORE) ;

validate:
	cd validation; ./validate "val"
