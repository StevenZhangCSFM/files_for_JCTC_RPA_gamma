# Choose whether to use MKL or ScaLAPACK or neither
# Warning: don't turn on USE_MKL and USE_SCALAPACK simultaneously
# Set USE_MKL = 1 to compile with MKL BLAS, LAPACK, and ScaLAPACK
# Set USE_MKL = 0 otherwise
USE_MKL       = 1
# Set USE_SCALAPACK = 1 to compile with non-MKL BLAS, LAPACK, and ScaLAPACK
# Set USE_SCALAPACK = 0 to compile with non-MKL BLAS and LAPACK only
USE_SCALAPACK = 0
# Set USE_DP_SUBEIG = 1 to use SPARC rather than ScaLAPACK routines for matrix data distribution
# (USE_DP_SUBEIG = 1 is required if both USE_MKL = 0 and USE_SCALAPACK = 0)
# Set USE_DP_SUBEIG = 0 to use ScaLAPACK rather than SPARC routines
USE_DP_SUBEIG = 0
# Set USE_FFTW = 1 to use FFTW for fast Fourier transform in vdWDF. Don't open it together with USE_MKL
USE_FFTW      = 0
# Set DEBUG_MODE = 1 to run with debug mode and print debug output
DEBUG_MODE    = 1


# Enable SIMD vectorization for complex stencil routines
# CAUTION: for some compilers this results in wrong results! Use for intel/19.0.3 or later versions
ENABLE_SIMD_COMPLEX = 1

# Specify the path MKLROOT if it's not already set to compile with MKL, e.g,
# MKLROOT = /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl

# Specify the path to ScaLAPACK, LAPACK and BLAS if necessary, and 
# add to LDFLAGS. Note that sometimes LDFLAGS already contains the default 
# path to these libraries, or the libraries are located in the default search
# path. In those cases, the following is not needed.
# SCALAPACKROOT = /nv/hp27/qxu78/data/scalapack-2.0.2
# LDFLAGS += -L$(SCALAPACKROOT)
# LAPACKROOT = /usr/local/pacerepov1/lapack/3.6.0
# LDFLAGS += -L$(LAPACKROOT)
# BLASROOT = /usr/lib64
# LDFLAGS += -L$(BLASROOT)

CPPFLAGS = -Iinclude -I../../src/include -I../../src/xc/exx/include -I../../src/xc/vdW/d3/include -I../../src/xc/vdW/vdWDF/include -I../../src/xc/mgga/include -I../../src/linearScaling/include -I../../src/ofdft/ -I../../src/ddbp/ -I../../src/extfpmd/ -I../../src/cyclix/include  -I../../src/cyclix/mlff -I../../src/highT/mlff -I../../src/mlff/
LDLIBS   = -lrt

ifeq ($(USE_MKL), 1)
CPPFLAGS += -m64 -I${MKLROOT}/include -DUSE_MKL
LDFLAGS   = -L${MKLROOT}/lib/intel64
LDLIBS   += -Wl,-rpath=${MKLROOT}/lib/intel64,--no-as-needed -lmkl_scalapack_lp64 -lmkl_cdft_core -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl
else ifeq ($(USE_SCALAPACK), 1)
CPPFLAGS += -DUSE_SCALAPACK
LDLIBS   += -lscalapack
endif

# if MKL is not used, link BLAS and LAPACK directly
ifeq ($(USE_MKL), 0)
# if you are using BLAS instead of OpenBLAS, change -lopenblas to -lblas
# and move it to after -llapack
LDLIBS += -lopenblas -lm
#LDLIBS += -llapacke -llapack -lblas -lm
endif

# if FFTW is used, link fftw
ifeq ($(USE_FFTW), 1)
CPPFLAGS += -DUSE_FFTW
LDLIBS += -lfftw3_mpi -lfftw3 # fftw needs to be loaded before compile
endif

# To use domain parallelization + LAPACK for solving sobspace eigen problem
ifeq ($(USE_DP_SUBEIG), 1)
CPPFLAGS += -DUSE_DP_SUBEIG
endif

# to compile with DEBUG mode
ifeq ($(DEBUG_MODE), 1)
CPPFLAGS += -Wall -g -DDEBUG # -Wno-format-truncation -Wno-stringop-truncation
endif

# to enable SIMD for complex stencil routines
ifeq ($(ENABLE_SIMD_COMPLEX), 1)
CPPFLAGS += -DENABLE_SIMD_COMPLEX
endif

# for old Intel compiler, use -qopenmp instead of -fopenmp. ICC 17 and later also accepts -fopenmp. 
CFLAGS = -std=gnu99 -O3 -fopenmp

SPARC_OBJSC = ../../src/initialization.o ../../src/readfiles.o ../../src/atomdata.o ../../src/parallelization.o ../../src/relax.o ../../src/tools.o ../../src/md.o      \
        ../../src/electrostatics.o ../../src/electronicGroundState.o ../../src/electronDensity.o ../../src/orbitalElecDensInit.o           \
        ../../src/occupation.o ../../src/gradVecRoutines.o ../../src/gradVecRoutinesKpt.o ../../src/nlocVecRoutines.o     \
        ../../src/hamiltonianVecRoutines.o ../../src/lapVecRoutines.o ../../src/lapVecRoutinesKpt.o \
        ../../src/linearSolver.o ../../src/mixing.o ../../src/exchangeCorrelation.o ../../src/eigenSolver.o ../../src/eigenSolverKpt.o ../../src/energy.o      \
        ../../src/forces.o ../../src/stress.o ../../src/pressure.o ../../src/finalization.o ../../src/spinOrbitCoupling.o ../../src/printing.o \
        ../../src/linearAlgebra.o ../../src/kroneckerLaplacian.o\
        ../../src/xc/vdW/d3/d3correction.o ../../src/xc/vdW/d3/d3findR0ab.o ../../src/xc/vdW/d3/d3copyC6.o                       \
        ../../src/xc/vdW/d3/d3initialization.o ../../src/xc/vdW/d3/d3finalization.o ../../src/xc/vdW/d3/d3forceStress.o          \
        ../../src/xc/vdW/vdWDF/vdWDFinitialization.o ../../src/xc/vdW/vdWDF/vdWDFfinalization.o                        \
        ../../src/xc/vdW/vdWDF/vdWDFnonlinearCorre.o ../../src/xc/vdW/vdWDF/vdWDFreadKernel.o                 \
        ../../src/xc/vdW/vdWDF/vdWDFstress.o ../../src/xc/vdW/vdWDF/vdWDFparallelization.o                        \
        ../../src/xc/mgga/mGGAscan.o ../../src/xc/mgga/mGGArscan.o ../../src/xc/mgga/mGGAr2scan.o ../../src/xc/mgga/mGGAhamiltonianTerm.o  \
        ../../src/xc/mgga/mGGAinitialization.o ../../src/xc/mgga/mGGAfinalization.o ../../src/xc/mgga/mGGAtauTransferTauVxc.o  \
        ../../src/xc/mgga/mGGAstress.o \
        ../../src/xc/exx/exactExchange.o ../../src/xc/exx/exactExchangeKpt.o ../../src/xc/exx/exactExchangeInitialization.o      \
        ../../src/xc/exx/exactExchangeFinalization.o ../../src/xc/exx/exactExchangeStress.o                            \
        ../../src/xc/exx/exactExchangePressure.o ../../src/xc/exx/exactExchangeEnergyDensity.o \
        ../../src/linearScaling/sq.o ../../src/linearScaling/sqInitialization.o ../../src/linearScaling/sqFinalization.o ../../src/linearScaling/sqEnergy.o\
        ../../src/linearScaling/sqDensity.o ../../src/linearScaling/sqNlocVecRoutines.o ../../src/linearScaling/sqParallelization.o ../../src/linearScaling/sqProperties.o\
        ../../src/highT/spdft.o ../../src/highT/sqHighTExactExchange.o ../../src/highT/sqHighT.o ../../src/highT/sqHighTDensity.o ../../src/highT/sqHighTInitialization.o\
        ../../src/ddbp/ddbp_tools.o ../../src/ddbp/ddbp_paral.o ../../src/ddbp/ddbp_basis.o ../../src/ddbp/ddbp_hamiltonian.o \
        ../../src/ddbp/ddbp_nloc.o ../../src/ddbp/ddbp_matvec.o ../../src/ddbp/ddbp_eigensolver.o ../../src/ddbp/ddbp.o ../../src/ddbp/ddbp_exactexchange.o\
        ../../src/cyclix/cyclix_tools.o ../../src/cyclix/cyclix_forces.o ../../src/cyclix/cyclix_gradVec.o ../../src/cyclix/cyclix_lapVec.o \
        ../../src/cyclix/cyclix_stress.o \
        ../../src/sq3/sq3.o ../../src/cs/cs.o ../../src/extfpmd/extFPMD.o ../../src/ofdft/ofdft.o \
        ../../src/mlff/covariance_matrix.o ../../src/mlff/regression.o ../../src/mlff/spherical_harmonics.o ../../src/mlff/soap_descriptor.o ../../src/mlff/sparc_mlff_interface.o \
        ../../src/mlff/besselfunctions.o ../../src/mlff/bessel_NR.o ../../src/mlff/nrutils.o ../../src/mlff/sparsification.o ../../src/mlff/tools_mlff.o ../../src/mlff/mlff_read_write.o ../../src/mlff/gmp_deriv_descriptor.o \
        ../../src/mlff/helper_gmp.o ../../src/mlff/solid_harmonics.o ../../src/mlff/descriptor.o ../../src/highT/mlff/internal_energy_model.o ../../src/cyclix/mlff/cyclix_mlff_tools.o \
        ../../src/mlff/hnl_soap.o \

OBJSC = main.o initialization_SPARC.o initialization_RPA.o generateKgrid.o parallelization_RPA.o parallelization_SPARC.o \
	restoreElectronicGroundState.o prepare_PQ_operators.o collectOrbitals.o cheFSI.o cheFSISQ.o eigenSolverGamma_RPA.o eigenSolverKpt_RPA.o \
        nuChi0VecRoutines.o sternheimerEquation.o sternheimerEquationKpt.o linearSolvers.o electrostatics_RPA.o printResult.o \
        finalization_RPA.o test_Hx_nuChi0_eigSolver.o tools_RPA.o traceSolverSQ.o

LIBBASE = ../lib/rpacalc

override CC=mpicc

all: rpacalc

# Note the implicit rule to compile '.c' files into '.o' files is
# %.o : %.c
# 	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

rpacalc: $(OBJSC)
	$(CC) $(CFLAGS) $(LDFLAGS) $(SPARC_OBJSC) -o $(LIBBASE) $^ $(LDLIBS)

.PHONY: clean
clean:
	rm -f  $(OBJSC) $(LIBBASE)
