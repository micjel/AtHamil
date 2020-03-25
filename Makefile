#--------------------------------------------------
# Make file for NuHamil code
#
# When we use ifort,
# -heap-arrays option is needed for built in transpose function with
#  large dimension matrix
#--------------------------------------------------
TARGET=AtHamil
INSTLDIR=$(HOME)/bin
EXEDIR=$(PWD)/exe/$(shell hostname)
EXEDIR=$(PWD)/exe
MODDIR = mod
Host= $(shell if hostname|grep -q apt1; then echo apt; \
  elif hostname|grep -q oak; then echo oak; \
  elif hostname|grep -q cedar; then echo cedar; \
  else echo other; fi)
HOST=$(strip $(Host))
DEBUG_MODE=off
Gauss_Laguerre=off

OS = Linux
ifneq (,$(findstring arwin,$(shell uname)))
  OS = OSX
endif
$(info HOST is $(HOST).)
$(info OS is $(OS).)
$(info Debug mode is $(DEBUG_MODE).)
$(info Memory save mode is $(MemorySave).)
BRANCH=$(shell git branch -v | grep '*' | awk '{printf "%s",$$2}')
COMMIT=$(shell git log -1 | awk 'NR==1 {printf "%s",$$2}')
ifneq (,$(findstring HEAD,$(BRANCH)))
  BRANCH=HEAD
endif
VERSION=$(BRANCH):$(COMMIT)

FDEP=
FC=
LFLAGS=  # library
FFLAGS=  # option
DFLAGS=  # debug
LINT=    # 8-byte integer

#--------------------------------------------------
# Default Parameters
#--------------------------------------------------
ifeq ($(strip $(HOST)),other)
  FDEP=makedepf90
  FC=gfortran
  LFLAGS+= -I/usr/local/include -L/usr/local/lib
  LFLAGS+= -lblas -llapack -lgsl -lz
  FFLAGS= -O3
  CFLAGS= -O3
  FFLAGS+= -fopenmp
  FFLAGS+= -Dsingle_precision
  FFLAGS+= -DVERSION=\"$(VERSION)\"
  ifeq ($(Gauss_Laguerre),on)
    FFLAGS+= -Dgauss_laguerre
  endif
  ifeq ($(DEBUG_MODE),on)
    DFLAGS+=-Wall -pedantic -fbounds-check -O -Wuninitialized -fbacktrace
    #FDFLAGS+=-ffpe-trap=invalid,zero,overflow # Note: gsl larguerre signal
    ifneq ($(OS), OSX)
      DFLAGS+= -pg -g
    endif
  endif
endif

ifeq ($(strip $(HOST)),oak)
  FDEP=makedepf90
  FC=ifort
  MKL=-L$(MKLROOT)/lib/ -L$(MKLROOT)/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -Wl,-rpath,$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/../compiler/lib/
  LFLAGS+= -I/usr/local/include -L/usr/local/lib
  LFLAGS+= $(MKL) -lgsl -lz
  FFLAGS= -O3
  CFLAGS= -O3
  FFLAGS+= -fopenmp
  FFLAGS+= -Dsingle_precision
  FFLAGS+= -DVERSION=\"$(VERSION)\"
  ifeq ($(Gauss_Laguerre),on)
    FFLAGS+= -Dgauss_laguerre
  endif
  ifeq ($(DEBUG_MODE),on)
    DFLAGS+=-Wall -pedantic -fbounds-check -O -Wuninitialized -fbacktrace
    #FDFLAGS+=-ffpe-trap=invalid,zero,overflow # Note: gsl larguerre signal
    ifneq ($(OS), OSX)
      DFLAGS+= -pg -g
    endif
  endif
endif


ifeq ($(DEBUG_MODE),on)
endif



#--------------------------------------------------
# Source Files
#--------------------------------------------------

SRCDIR = src
SRCDIR_Main = main
SRCDIR_LinAlg = submodules/LinAlgf90/src
DEPDIR = .
OBJDIR = obj

SRCS=
OBJS=
MODS=

SRCF90:=$(wildcard $(SRCDIR)/*.f90)
SRCF95:=$(wildcard $(SRCDIR)/*.F90)
OBJF90:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90))))
OBJF95:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95))))
SRCS= $(SRCF90) $(SRCF95)
OBJS= $(OBJF90) $(OBJF95)

SRCF90_Main:=$(wildcard $(SRCDIR_Main)/*.f90)
SRCF95_Main:=$(wildcard $(SRCDIR_Main)/*.F90)
OBJF90_Main:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90_Main))))
OBJF95_Main:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95_Main))))
SRCS_Main= $(SRCF90_Main) $(SRCF95_Main)
OBJS_Main= $(OBJF90_Main) $(OBJF95_Main)

SRCF90_LinAlg:=$(wildcard $(SRCDIR_LinAlg)/*.f90)
SRCF95_LinAlg:=$(wildcard $(SRCDIR_LinAlg)/*.F90)
OBJF90_LinAlg:=$(addprefix $(OBJDIR)/, $(patsubst %f90, %o, $(notdir $(SRCF90_LinAlg))))
OBJF95_LinAlg:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95_LinAlg))))
SRCS_LinAlg= $(SRCF90_LinAlg) $(SRCF95_LinAlg)
OBJS_LinAlg= $(OBJF90_LinAlg) $(OBJF95_LinAlg)


SRCS_ALL = $(SRCS) $(SRCS_Main) $(SRCS_LinAlg)
OBJS_ALL = $(OBJS) $(OBJS_Main) $(OBJS_LinAlg)

#$(info $(SRCS_ALL))
#$(info )
#$(info $(OBJS_ALL))

MODOUT=
ifeq ($(strip $(HOST)),other)
  MODOUT=-J$(MODDIR)
endif

ifeq ($(strip $(HOST)),oak)
  MODOUT=-module $(MODDIR)
endif
#--------------------------------------------------
# Rules
#--------------------------------------------------
all: dirs $(TARGET)
$(TARGET): $(OBJS_ALL)
	$(FC) $(FFLAGS) $(DFLAGS) -o $(TARGET).exe $^ $(LFLAGS)
	if test -d $(EXEDIR); then \
		: ; \
	else \
		mkdir -p $(EXEDIR); \
	fi
	mv $(TARGET).exe $(EXEDIR)
	@echo "#####################################################################################"
	@echo "To complete the installation, do 'make install'."
	@echo "Edit '$(PWD)/exe/AtHamil.py' and excecute."
	@echo "#####################################################################################"

$(OBJDIR)/%.o:$(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT)  -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT)  -o $@ -c $<

$(OBJDIR)/%.o:$(SRCDIR_Main)/%.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT)  -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_Main)/%.F90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT)  -o $@ -c $<

$(OBJDIR)/%.o:$(SRCDIR_LinAlg)/%.f90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT)  -o $@ -c $<
$(OBJDIR)/%.o:$(SRCDIR_LinAlg)/%.F90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT)  -o $@ -c $<

dep:
	$(FDEP) $(SRCS_ALL) -b $(OBJDIR)/ > $(DEPDIR)/makefile.d

dirs:
	if test -d $(OBJDIR); then \
		: ; \
	else \
		mkdir $(OBJDIR); \
	fi
	if test -d $(MODDIR); then \
		: ; \
	else \
		mkdir $(MODDIR); \
	fi

install:
	ln -sf $(EXEDIR)/$(TARGET).exe $(INSTLDIR)
	@printf "####################################################################\n"
	@printf " make sure that '$(INSTLDIR)' is included in PATH \n"
	@printf "####################################################################\n"

clean:
	rm -r obj
	rm -r mod
	rm -f $(TARGET).exe
	rm -f $(MODDIR)/*.mod
	rm -f $(OBJS_ALL)
#--------------------------------------------------
-include $(wildcard $(DEPDIR)/*.d)
