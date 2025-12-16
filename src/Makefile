.SUFFIXES: .f90 .o
#---------------------------------------------------------------------------------
#(1)perpare for 
#---------------------------------------------------------------------------------
FC = mpiifx
CC = g++
LINKER = $(FC)
LOPT_DIR :=LBFGS CGplus

MODULES =Constants.o Smpi_math_module.o Parameters.o Math.o Cellsymmetry_module.o \
	 Struct_module.o Fourier.o Read_module.o Grid_symmetry.o Finite_module.o\
	  Begin_module.o Smearing_module.o  Lapack_module.o Matvec_module.o\
	 Chebyshev_fliter.o lbfgsb.o cgfam.o cgsearch.o blas.o Opt_mdoule.o\
	 IKS_module.o

BASE = Main.o
#BASEF = dcsrch.o dcstep.o

#lib or other include file path
FFTW_DIR=/home/xq/Softwares/fftw-3.3.9/FFTW
LIB_PATH_ROOT=${FFTW_DIR}/lib
LIB_MKL_ROOT=/opt/intel/oneapi/mkl/2024.2/
LIB_PATH=-L ${LIB_PATH_ROOT}
# LIB_MKL_PATH=-L ${LIB_MKL_ROOT}
#----------------------------------------------------------------------------------
# (2) PLATFORM INDEPENDENT VARIABLES.
#----------------------------------------------------------------------------------

GOAL = iks.x

# containing object files for serial version
WORKDIR = ../bin

FFLAGS = -O3 -r8 -i4 -I${LIB_PATH_ROOT}/include -cpp -DMPI -module ${WORKDIR} -qmkl=parallel
#FFLAGS = -O3 -cpp  -module ${WORKDIR}
#parallrl swetich
#debug
#FFLAGS = -g -CA -CU -CB -traceback -check all -cpp -DMPI  -I${LIB_PATH_ROOT}/include -module ${WORKDIR}
LIBS = -I${LIB_PATH_ROOT}/bin -lfftw3
#LIBS = -lfftw3  -larpack_ifort -lxc
#MKLLIBS= -I${LIB_PATH_ROOT}/bin -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -liomp5 -lpthread
#MKLLIBS= -I${LIB_PATH_ROOT}/bin -qmkl

MKL_PATH = $(MKLROOT)/lib/intel64
BLAS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
BLACS = -lmkl_blacs_intelmpi_lp64 
#$(info $(MKL_PATH))
SCALAPACK = $(BLACS) -lmkl_scalapack_lp64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64
MKLLIBS = -L$(MKL_PATH) $(SCALAPACK) $(BLAS) 

#ALL = $(MODULES) $(BASE) $(BASEF) $(RUN) $(IO) $(CMOD)
vpath %.f ${LOPT_DIR}

ALLS = $(MODULES) $(BASE) $(BASEF)

$(WORKDIR)/%.o : %.f
	$(FC) -c $(FFLAGS) $< -o $@

$(WORKDIR)/%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

$(WORKDIR)/%.o : %.c
	$(CC) -c $(FFLAGS) $< -o $@

# Absolute path of the object files (with working dir attached)
OBJS = $(patsubst %.o, $(WORKDIR)/%.o, $(ALLS))

default :
	@ if [ ! -d $(WORKDIR) ]; then mkdir $(WORKDIR); fi
	@ chmod a+rwX $(WORKDIR)
	@ $(MAKE) $(GOAL)

# Rules
$(GOAL) : $(OBJS)
	$(LINKER) $(FFLAGS) -o $(GOAL) $(OBJS) ${LIB_PATH} $(LIBS) $(MKLLIBS)

cleanall : clean

all : default

.PHONY : cleanall all tags

clean :
	@ if [ -d $(WORKDIR) ]; then $(RM) -r $(WORKDIR); fi
	@ if [ -f $(GOAL) ] ; then $(RM) $(GOAL); fi

tags :
	@ ctags *.f90
