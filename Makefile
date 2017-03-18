#COMPILER (PARALLEL)
FC=mpif90
#gfortran
#PRECOMPILATION FLAG (leave blank for serial code)
FPP=

#EXE=GZ_slab
#EXE=GZ_slab_chem_pot
#EXE=GZ_dynamics_slab_ins
#EXE=GZ_equ_temp_china
EXE=GZ_obs_slab
#EXE=GZ_obs_slab_utest
#EXE=GZ_extrema_ratio

DIR=drivers
DIREXE=$(HOME)/.project_bin

.SUFFIXES: .f90

#REVISION SOFTWARE GIT:
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc
#
# OBJS=RK_VIDE.o MATRIX_SPARSE.o AMOEBA.o GZ_VARS_INPUT.o GZ_VARS_GLOBAL.o ELECTRIC_FIELD.o  GZ_AUX_FUNX.o GZ_neqAUX_FUNX.o GZ_LOCAL_FOCK_SPACE.o GZ_VARIATIONAL_BASIS.o GZ_LOCAL_HAMILTONIAN.o GZ_EFFECTIVE_HOPPINGS.o GZ_ENERGY.o GZ_OPTIMIZE.o GZ_DYNAMICS.o GZ_GREENS_FUNCTIONS.o
#
OBJS=amoeba.o RK_VIDE.o   VECTORS.o   BZ_POINTS.o  GUTZWILLER.o GLOBAL.o EQUILIBRIUM.o  slabEQS_OF_MOTION.o obsEQS_OF_MOTION.o  slabDYNAMICS.o obsDYNAMICS.o


#FFLAG +=-fpp -D_$(FPP) ONLY WITH mpif90
LIBDIR=$(HOME)/opt_local
#LIBDIR=/opt/


#GALLIBDIR  = $(LIBDIR)/galahad/objects/mac64.osx.gfo/double
#GALLIBMOD  = $(LIBDIR)/galahad/modules/mac64.osx.gfo/double
GALLIBDIR  = $(LIBDIR)/galahad/objects/pc64.lnx.gfo/double
GALLIBMOD  = $(LIBDIR)/galahad/modules/pc64.lnx.gfo/double


GALLIBS1   = -lgalahad -lgalahad_hsl 
GALLIBS2   = -lgalahad_metis 

MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm


INCARGS =-I$(LIBDIR)/SciFortran/gnu/include -L$(LIBDIR)/SciFortran/gnu/lib 
INCARGS+=-I$(LIBDIR)/DMFTtools/gnu/include -L$(LIBDIR)/DMFTtools/gnu/lib 
INCARGS+=-I$(GALLIBDIR) -L$(GALLIBDIR)
FFLAG += -ffree-line-length-none -cpp $(INCARGS)

#FFLAG+=-O0 -p -g -Wall -fbounds-check -fbacktrace -Wuninitialized



ARGS= -L$(GALLIBDIR) $(GALLIBS1) $(GALLIBS2) -I$(GALLIBMOD) -ldmftt -lscifor  -lfftpack -lminpack  -llapack -lblas -larpack #-lparpack    

all:compile


lib: ed_solver

compile: version $(OBJS)
	@echo " !+------------------------------------------------- "
	@echo " ..................... compile ..................... "
	@echo " !+------------------------------------------------- "
	$(FC) $(FFLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " !+--------------------------------------+! "
	@echo " .................. done .................. "
	@echo " !+--------------------------------------+! "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)

ed_solver:
	@make -C ED_SOLVER/

.f90.o:	
	$(FC) $(FFLAG)  -c $< 

completion:
	sf_lib_completion.sh $(DIR)/$(EXE).f90
	@echo "run: . .bash_completion.d/$(EXE) to add completion for $(EXE) in this shell"

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

all_clean: clean
	@make -C ED_SOLVER/ clean

version:
	@echo $(VER)

