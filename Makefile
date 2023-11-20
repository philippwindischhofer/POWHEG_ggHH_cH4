#-*- Makefile -*-
## Choose compiler: gfortran,ifort ( g77 no longer compatible )
COMPILER=gfortran
## Choose PDF: native,lhapdf5,lhapdf6
## LHAPDF package has to be installed separately
PDF = lhapdf6
#Choose Analysis: none, default
## default analysis may require FASTJET package, that has to be installed separately (see below)
ANALYSIS = default
## For static linking uncomment the following
#STATIC= -static
#

OBJ=obj-$(COMPILER)
OBJDIR:=$(OBJ)

MODDIR=mod-$(COMPILER)

ifeq ("$(COMPILER)","gfortran")
FORTCOMP=gfortran
FC=$(FORTCOMP)
FF=$(FORTCOMP) -cpp -fno-automatic -ffixed-line-length-132
# Additional flags for Mac OS X
UNAME_S := $(shell uname -s)
ifeq ("$(UNAME_S)","Darwin")
FF+= -lc++
endif
## -fbounds-check sometimes causes a weird error due to non-lazy evaluation
CC=gcc
CXX=g++
## of boolean in gfortran.
#FFLAGS= -Wall -Wimplicit-interface -fbounds-check -g -O2
FFLAGS= -g -O2 -std=legacy -fno-lto
## For floating point exception trapping  uncomment the following
#FPE=-ffpe-trap=invalid,zero,overflow,underflow
## Use -O2 optimization
OPT=-O2
## For debugging uncomment the following
#DEBUG= -ggdb -pg
MOD= -J $(MODDIR)
endif

ifeq ("$(COMPILER)","ifort")
FORTCOMP=ifort
FC=$(FORTCOMP)
FF=$(FORTCOMP) -save -cpp -extend_source -assume bscc
CC=icpc
CXX=icpc
LIBS=-limf
#FFLAGS =  -check
## For floating point exception trapping  uncomment the following
#FPE = -fpe0
OPT = -O3 #-fast
## For debugging uncomment the following
#DEBUG= -debug -g
MOD = -module $(MODDIR)
endif

GOSAMLIBS=$(shell cd GoSamlib ; echo *.f *.f90 ' ' | sed 's/qlonshellcutoff.f// ; s/qlconstants.f// ; s/.f /.o /g ; s/.f90 /.o /g ; ')
GOSAMLIBS+=$(shell cd GoSamlib ; echo *.cpp *.cc *.F90 ' ' | sed 's/.cc /.o /g ; s/.F90 /.o /g ; ')
GOSAMLIBS+=$(shell cd GoSamlib/quadsources ; echo *.cc *.F90 ' ' | sed 's/.cc /.o /g ; s/.F90 /.o /g ; ')

PWD=$(shell pwd)
WDNAME=$(shell basename $(PWD))
VPATH= ./:../:$(OBJDIR):./GoSamlib/:./Gosamlib/quadsources/:./Virtual/

INCLUDE0=$(PWD)
INCLUDE1=$(shell dirname $(PWD))/include
INCLUDE2=$(PWD)/GoSamlib
INCLUDE3=$(PWD)/Virtual
FF+= $(FFLAGS) $(FPE) $(OPT) $(DEBUG) -I$(INCLUDE0) -I$(INCLUDE1) -I$(INCLUDE2) -I$(INCLUDE3) 
LIBS+=-lz -lquadmath

INCLUDE =$(wildcard ../include/*.h *.h include/*.h)

ifeq ("$(PDF)","lhapdf6")
LHAPDF_CONFIG=lhapdf-config
FJCXXFLAGS+= $(shell $(LHAPDF_CONFIG) --cxxflags)
PDFPACK=lhapdf6if.o lhapdf6ifcc.o
LIBSLHAPDF= -Wl,-rpath,$(shell $(LHAPDF_CONFIG) --libdir)  -L$(shell $(LHAPDF_CONFIG) --libdir) -lLHAPDF -lstdc++
ifeq  ("$(STATIC)","-static")
## If LHAPDF has been compiled with gfortran and you want to link it statically, you have to include
## libgfortran as well. The same holds for libstdc++.
## One possible solution is to use fastjet, since $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
## does perform this inclusion. The path has to be set by the user.
# LIBGFORTRANPATH= #/usr/lib/gcc/x86_64-redhat-linux/4.1.2
# LIBSTDCPP=/lib64
# LIBSLHAPDF+=  -L$(LIBGFORTRANPATH)  -lgfortranbegin -lgfortran -L$(LIBSTDCPP) -lstdc++
endif
LIBS+=$(LIBSLHAPDF)
else
PDFPACK=mlmpdfif.o hvqpdfpho.o
endif

ifeq ("$(ANALYSIS)","default")
##To include Fastjet configuration uncomment the following lines.
FASTJET_CONFIG=$(shell which fastjet-config)
LIBSFASTJET += $(shell $(FASTJET_CONFIG) --libs --plugins ) -lstdc++
FJCXXFLAGS+= $(shell $(FASTJET_CONFIG) --cxxflags)
PWHGANAL=pwhg_analysis.o
## Also add required Fastjet drivers to PWHGANAL (examples are reported)
PWHGANAL+= fastjetfortran.o
else
PWHGANAL=pwhg_analysis-dummy.o
endif

all: pwhg_main

NINJAFLAGS=-DHAVE_CONFIG_H -I./GoSamlib  -O2
NINJAFLAGSQP=-DHAVE_CONFIG_H -I./GoSamlib -I./Gosamlib/quadsources -O2

# VirtualGrids
LIBS+= $(shell python3-config --ldflags)
CCFLAGS+=$(shell python3-config --cflags)

pwhg_cpHTO_reweight.o: pwhg_cpHTO_reweight.f $(INCLUDE)
	$(FF) -c $(MOD) -o $(OBJDIR)/$@ $<

%.o: %.F $(INCLUDE)  | $(OBJDIR) $(MODDIR)
	$(FF) -c $(MOD) -o $(OBJDIR)/$@ $<

%.o: %.f $(INCLUDE)  | $(OBJDIR) $(MODDIR)
	$(FF) -c $(MOD) -o $(OBJDIR)/$@ $<

%.o: %.f90 $(INCLUDE) | $(OBJDIR)
	$(FF) -c $(MOD) -o $(OBJDIR)/$@ $<

%_qp.o: %_qp.f90 $(INCLUDE) quadninjago.o | $(OBJDIR)
	$(FF) -c $(MOD) -o $(OBJDIR)/$@ $<

%.o: %.F90 $(INCLUDE) | $(OBJDIR)
	$(FF) -DHAVE_CONFIG_H -c $(MOD) -o $(OBJDIR)/$@ $<

%.o: quadsources/%.F90 $(INCLUDE) | $(OBJDIR)
	$(FF) -DHAVE_CONFIG_H -c $(MOD) -o $(OBJDIR)/$@ $<

%.o: %.c
	$(CC) $(DEBUG) -c -o $(OBJDIR)/$@ $^ $(CCFLAGS)

%.o: %.cpp
	$(CXX) $(DEBUG) -c -o $(OBJDIR)/$@ $^

%.o: %.cc
	$(CXX) $(DEBUG) -c -o $(OBJDIR)/$@ $^ $(NINJAFLAGS) $(FJCXXFLAGS) $(LIBHEPMC)

%.o: quadsources/%.cc
	$(CXX) $(DEBUG) -c -o $(OBJDIR)/$@ $^ $(NINJAFLAGSQP) $(FJCXXFLAGS)

# PYTHIA 8
#PYTHIA8LOCATION=$(shell pythia8-config --prefix)
PYTHIA8LOCATION=/usr/
FJCXXFLAGS+=-I$(PYTHIA8LOCATION)/include -I$(PYTHIA8LOCATION)/include/Pythia8 -I$(PYTHIA8LOCATION)/include/Pythia8Plugins
LIBPYTHIA8=-L$(PYTHIA8LOCATION)/lib64/ -lpythia8  -lstdc++ -ldl -lz
LIBHEPMC=$(shell pythia8-config --hepmc)

USER=init_couplings.o init_processes.o Born_phsp.o Born.o virtual.o MEborn.o grids.o D0functions.o \
        real.o btilde_gghh.o pt2maxreg.o LesHouchesreg.o pwhg_cpHTO_reweight.o hhgrid.o \
	$(PWHGANAL)

PWHG=pwhg_main.o pwhg_init.o bbinit.o lhefwrite.o LesHouches.o phi3m.o  \
	gen_Born_phsp.o find_regions.o test_Sudakov.o sigborn.o	boost.o	\
	gen_real_phsp.o maxrat.o gen_index.o gen_radiation.o		\
	Bornzerodamp.o sigremnants_gghh.o random.o boostrot.o		\
	bra_ket_subroutines.o cernroutines.o init_phys.o powheginput.o	\
	pdfcalls.o sigreal.o sigcollremn.o pwhg_bookhist-multi-new.o	\
	pwhg_analysis_driver.o checkmomzero.o setstrongcoupl.o		\
	integrator.o newunit.o mwarn.o sigsoftvirt.o sigcollsoft.o	\
	sigvirtual.o setlocalscales.o validflav.o pwhgreweight.o        \
	opencount.o ubprojections.o utils.o mint_upb.o lhefread.o       \
	pwhg_io_interface.o rwl_weightlists.o rwl_setup_param_weights.o \
	refine_momenta.o $(PDFPACK) $(USER) $(FPEOBJ)

# target to generate LHEF output
pwhg_main:gosamlib.a $(PWHG)
	$(FF) $(patsubst %,$(OBJDIR)/%,$(PWHG)) $(LIBS) $(LIBSFASTJET) $(STATIC) -static-libgfortran  $(OBJ)/gosamlib.a -o $@ -lstdc++


gosamlib.a: $(GOSAMLIBS)
	\rm $(OBJ)/gosamlib.a ; $(PWD)/GoSamlib/compile_gosamlib.sh $(PWD)/$(OBJ) $(PWD)/$(MODDIR)

LHEF=lhef_analysis.o boostrot.o random.o cernroutines.o		     \
     opencount.o powheginput.o pwhg_bookhist-multi-new.o $(PWHGANAL) \
     lhefread.o newunit.o pwhg_analysis_driver.o                     \
     pwhg_io_interface.o rwl_weightlists.o $(FPEOBJ)

# target to analyze LHEF output
lhef_analysis:$(LHEF)
	$(FF) $(patsubst %,$(OBJDIR)/%,$(LHEF)) $(LIBS) $(LIBSFASTJET) $(STATIC)  -o $@


# target to read event file, shower events with PYTHIA8 + analysis
PYTHIA8=main-PYTHIA8.o pythia8F77.o boostrot.o pwhg_bookhist-multi-new.o \
	powheginput.o $(PWHGANAL) opencount.o lhefread.o newunit.o	 \
	pdfdummies.o reshufflemoms.o random.o cernroutines.o 		 \
	Check_LesHouches.o pwhg_io_interface.o rwl_weightlists.o         \
	zlibdummy.o bra_ket_subroutines.o $(FPEOBJ)

main-PYTHIA8-lhef: $(PYTHIA8)
	$(FF) $(LIBPYTHIA8) $(patsubst %,$(OBJ)/%,$(PYTHIA8)) $(LIBSFASTJET) $(LIBPYTHIA8) $(LIBHEPMC) $(STATIC) $(LIBSLHAPDF) -o $@

# target to read event file, shower events with PYTHIA + analysis
PYTHIA6=main-PYTHIA.o setup-PYTHIA-lhef.o pythia.o boostrot.o powheginput.o    \
	$(PWHGANAL) lhefread.o pwhg_io_interface.o rwl_weightlists.o newunit.o \
	pwhg_analysis_driver.o random.o cernroutines.o opencount.o             \
	pwhg_bookhist-multi-new.o $(FPEOBJ)

main-PYTHIA6-lhef: $(PYTHIA6)
	$(FF) $(patsubst %,$(OBJ)/%,$(PYTHIA6)) $(LIBS) $(LIBSFASTJET)  $(STATIC) -o $@


# target to generate object directory if it does not exist
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(MODDIR):
	mkdir -p $(MODDIR)

# target to cleanup
.PHONY: clean veryclean

clean:
	rm -f $(patsubst %,$(OBJ)/%,$(USER) $(PWHG) $(LHEF) $(HERWIG) $(PYTHIA8) $(PYTHIA6)) pwhg_main lhef_analysis	\
	main-PYTHIA8-lhef main-PYTHIA6-lhef
	rm -f btilde_gghh.f sigremnants_gghh.f

veryclean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.* $(OBJDIR)/*.a \
	rm -f btilde_gghh.f sigremnants_gghh.f          \
	rm -f main-PYTHIA8-lhef main-PYTHIA6-lhef	\
	pwhg_main lhef_analysis main-HERWIG-lhef main-PYTHIA-lhef


##########################################################################

include GoSamlib/Makefile.virt.dep

MEborn.o:  precision_golem.o parametre.o matrice_s.o form_factor_type.o form_factor_4p.o form_factor_3p.o cache.o constante.o array.o spinor.o D0functions.o

Born.o: MEborn.o

btilde_gghh.f: btilde.f
	bash generate_files.sh btilde

sigremnants_gghh.f: sigremnants.f
	bash generate_files.sh sigremnants
