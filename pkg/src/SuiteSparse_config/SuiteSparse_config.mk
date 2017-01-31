#===============================================================================
# SuiteSparse_config.mk:  common configuration file for the SuiteSparse
#===============================================================================

# This file contains all configuration settings for all packages in SuiteSparse,
# except for CSparse (which is stand-alone) and the packages in MATLAB_Tools.

SUITESPARSE_VERSION = 4.5.4

#===============================================================================
# Options you can change without editing this file:
#===============================================================================

    # To list the options you can modify at the 'make' command line, type
    # 'make config', which also lists their default values.  You can then
    # change them with 'make OPTION=value'.  For example, to use an INSTALL
    # path of /my/path, and to use your own BLAS and LAPACK libraries, do:
    #
    #   make install INSTALL=/my/path BLAS=-lmyblas LAPACK=-lmylapackgoeshere
    #
    # which will install the package into /my/path/lib and /my/path/include,
    # and use -lmyblas -lmylapackgoes here when building the demo program.

#===============================================================================
# Defaults for any system
#===============================================================================

    #---------------------------------------------------------------------------
    # SuiteSparse root directory
    #---------------------------------------------------------------------------

    # Most Makefiles are in SuiteSparse/Pkg/Lib or SuiteSparse/Pkg/Demo, so
    # the top-level of SuiteSparse is in ../.. unless otherwise specified.
    # This is true for all but the SuiteSparse_config package.
    SUITESPARSE ?= $(realpath $(CURDIR)/../..)

    #---------------------------------------------------------------------------
    # installation location
    #---------------------------------------------------------------------------

    # For "make install" and "make uninstall", the default location is
    # SuiteSparse/lib, SuiteSparse/include, and
    # SuiteSparse/share/doc/suitesparse-x.y.z
    # If you do this:
    #    make install INSTALL=/usr/local
    # then the libraries are installed in /usr/local/lib, include files in
    # /usr/local/include, and documentation in
    # /usr/local/share/doc/suitesparse-x.y.z.
    # You can instead specify the install location of each of these 3 components
    # separately, via (for example):
    #    make install INSTALL_LIB=/yada/mylibs INSTALL_INCLUDE=/yoda/myinc  \
    #                 INSTALL_DOC=/solo/mydox
    # which puts the libraries in /yada/mylibs, include files in /yoda/myinc,
    # and documentation in /solo/mydox.
    INSTALL ?= $(SUITESPARSE)
    INSTALL_LIB ?= $(INSTALL)/lib
    INSTALL_INCLUDE ?= $(INSTALL)/include
    INSTALL_DOC ?= $(INSTALL)/share/doc/suitesparse-$(SUITESPARSE_VERSION)

    #---------------------------------------------------------------------------
    # optimization level
    #---------------------------------------------------------------------------

    OPTIMIZATION ?= -O3

    #---------------------------------------------------------------------------
    # CFLAGS for the C/C++ compiler
    #---------------------------------------------------------------------------

    # The CF macro is used by SuiteSparse Makefiles as a combination of
    # CFLAGS, CPPFLAGS, TARGET_ARCH, and system-dependent settings.
    CF ?= $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(OPTIMIZATION) -fexceptions

    #---------------------------------------------------------------------------
    # OpenMP is used in CHOLMOD
    #---------------------------------------------------------------------------

    # with gcc, enable OpenMP directives via -fopenmp
    # This is not supported on Darwin, so this string is cleared, below.
    CFOPENMP ?= -fopenmp

    #---------------------------------------------------------------------------
    # compiler
    #---------------------------------------------------------------------------

    # By default, look for the Intel compilers.  If present, they are used
    # instead of $(CC), $(CXX), and $(F77).  To disable this feature and
    # use the $(CC), $(CXX), and $(F77) compilers, use 'make AUTOCC=no'

    AUTOCC ?= yes

    ifneq ($(AUTOCC),no)
        ifneq ($(shell which icc 2>/dev/null),)
            # use the Intel icc compiler for C codes, and -qopenmp for OpenMP
            CC = icc -D_GNU_SOURCE
            CXX = $(CC)
            CFOPENMP = -qopenmp -I$(MKLROOT)/include
	    LDFLAGS += -openmp
        endif
        ifneq ($(shell which ifort 2>/dev/null),)
            # use the Intel ifort compiler for Fortran codes
            F77 = ifort
        endif
    endif

    #---------------------------------------------------------------------------
    # code formatting (for Tcov on Linux only)
    #---------------------------------------------------------------------------

    PRETTY ?= grep -v "^\#" | indent -bl -nce -bli0 -i4 -sob -l120

    #---------------------------------------------------------------------------
    # required libraries
    #---------------------------------------------------------------------------

    # SuiteSparse requires the BLAS, LAPACK, and -lm (Math) libraries.
    # It places its shared *.so libraries in SuiteSparse/lib.
    # Linux also requires the -lrt library (see below)
    LDLIBS ?= -lm
    LDFLAGS += -L$(INSTALL_LIB)

    # See http://www.openblas.net for a recent and freely available optimzed
    # BLAS.  LAPACK is at http://www.netlib.org/lapack/ .  You can use the
    # standard Fortran LAPACK along with OpenBLAS to obtain very good
    # performance.  This script can also detect if the Intel MKL BLAS is
    # installed.

    LAPACK ?= -llapack

    ifndef BLAS
        ifdef MKLROOT
            # use the Intel MKL for BLAS and LAPACK
            # using static linking:
            # BLAS = -Wl,--start-group \
            #   $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
            #   $(MKLROOT)/lib/intel64/libmkl_core.a \
            #   $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
            #   -Wl,--end-group -lpthread -lm
            # using dynamic linking:
            BLAS = -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm
            LAPACK =
        else
            # use the OpenBLAS at http://www.openblas.net
            BLAS = -lopenblas
        endif
    endif

    # For ACML, use this instead:
    #   make BLAS='-lacml -lgfortran'

    #---------------------------------------------------------------------------
    # shell commands
    #---------------------------------------------------------------------------

    # ranlib, and ar, for generating libraries.  If you don't need ranlib,
    # just change it to RANLAB = echo
    RANLIB ?= ranlib
    ARCHIVE ?= $(AR) $(ARFLAGS)
    CP ?= cp -f
    MV ?= mv -f

    #---------------------------------------------------------------------------
    # Fortran compiler (not required for 'make' or 'make library')
    #---------------------------------------------------------------------------

    # A Fortran compiler is optional.  Only required for the optional Fortran
    # interfaces to AMD and UMFPACK.  Not needed by 'make' or 'make install'
    F77 ?= gfortran
    F77FLAGS ?= $(FFLAGS) $(OPTIMIZATION)

    #---------------------------------------------------------------------------
    # NVIDIA CUDA configuration for CHOLMOD and SPQR
    #---------------------------------------------------------------------------

    # CUDA is detected automatically, and used if found.  To disable CUDA,
    # use CUDA=no

    ifneq ($(CUDA),no)
        CUDA_PATH = $(shell which nvcc 2>/dev/null | sed "s/\/bin\/nvcc//")
    endif

    ifeq ($(wildcard $(CUDA_PATH)),)
        # CUDA is not present
        CUDA_PATH     =
        GPU_BLAS_PATH =
        GPU_CONFIG    =
        CUDART_LIB    =
        CUBLAS_LIB    =
        CUDA_INC_PATH =
        CUDA_INC      =
        NVCC          = echo
        NVCCFLAGS     =
    else
        # with CUDA for CHOLMOD and SPQR
        GPU_BLAS_PATH = $(CUDA_PATH)
        # GPU_CONFIG must include -DGPU_BLAS to compile SuiteSparse for the
        # GPU.  You can add additional GPU-related flags to it as well.
        # with 4 cores (default):
        GPU_CONFIG    = -DGPU_BLAS
        # For example, to compile CHOLMOD for 10 CPU cores when using the GPU:
        # GPU_CONFIG  = -DGPU_BLAS -DCHOLMOD_OMP_NUM_THREADS=10
        CUDART_LIB    = $(CUDA_PATH)/lib64/libcudart.so
        CUBLAS_LIB    = $(CUDA_PATH)/lib64/libcublas.so
        CUDA_INC_PATH = $(CUDA_PATH)/include/
        CUDA_INC      = -I$(CUDA_INC_PATH)
        NVCC          = $(CUDA_PATH)/bin/nvcc
        NVCCFLAGS     = -Xcompiler -fPIC -O3 \
                            -gencode=arch=compute_30,code=sm_30 \
                            -gencode=arch=compute_35,code=sm_35 \
                            -gencode=arch=compute_50,code=sm_50 \
                            -gencode=arch=compute_50,code=compute_50
    endif

    #---------------------------------------------------------------------------
    # UMFPACK configuration:
    #---------------------------------------------------------------------------

    # Configuration for UMFPACK.  See UMFPACK/Source/umf_config.h for details.
    #
    # -DNBLAS       do not use the BLAS.  UMFPACK will be very slow.
    # -D'LONGBLAS=long' or -DLONGBLAS='long long' defines the integers used by
    #               LAPACK and the BLAS (defaults to 'int')
    # -DNSUNPERF    do not use the Sun Perf. Library on Solaris
    # -DNRECIPROCAL do not multiply by the reciprocal
    # -DNO_DIVIDE_BY_ZERO   do not divide by zero
    # -DNCHOLMOD    do not use CHOLMOD as a ordering method.  If -DNCHOLMOD is
    #               included in UMFPACK_CONFIG, then UMFPACK does not rely on
    #               CHOLMOD, CAMD, CCOLAMD, COLAMD, and METIS.

    UMFPACK_CONFIG ?=

    # For example, uncomment this line to compile UMFPACK without CHOLMOD:
    UMFPACK_CONFIG = -DNCHOLMOD
    # or use 'make UMFPACK_CONFIG=-DNCHOLMOD'

    #---------------------------------------------------------------------------
    # CHOLMOD configuration
    #---------------------------------------------------------------------------

    # CHOLMOD Library Modules, which appear in -lcholmod
    # Core       requires: none
    # Check      requires: Core
    # Cholesky   requires: Core, AMD, COLAMD. optional: Partition, Supernodal
    # MatrixOps  requires: Core
    # Modify     requires: Core
    # Partition  requires: Core, CCOLAMD, METIS.  optional: Cholesky
    # Supernodal requires: Core, BLAS, LAPACK
    #
    # CHOLMOD test/demo Modules (these do not appear in -lcholmod):
    # Tcov       requires: Core, Check, Cholesky, MatrixOps, Modify, Supernodal
    #            optional: Partition
    # Valgrind   same as Tcov
    # Demo       requires: Core, Check, Cholesky, MatrixOps, Supernodal
    #            optional: Partition
    #
    # Configuration flags:
    # -DNCHECK      do not include the Check module.
    # -DNCHOLESKY   do not include the Cholesky module.
    # -DNPARTITION  do not include the Partition module.
    #               also do not include METIS.
    # -DNCAMD       do not use CAMD & CCOLAMD in Parition Module.
    # -DNMATRIXOPS  do not include the MatrixOps module.
    # -DNMODIFY     do not include the Modify module.
    # -DNSUPERNODAL do not include the Supernodal module.
    #
    # -DNPRINT      do not print anything.
    # -D'LONGBLAS=long' or -DLONGBLAS='long long' defines the integers used by
    #               LAPACK and the BLAS (defaults to 'int')
    # -DNSUNPERF    for Solaris only.  If defined, do not use the Sun
    #               Performance Library
    # -DGPU_BLAS    enable the use of the CUDA BLAS

    CHOLMOD_CONFIG ?= $(GPU_CONFIG)

    #---------------------------------------------------------------------------
    # SuiteSparseQR configuration:
    #---------------------------------------------------------------------------

    # The SuiteSparseQR library can be compiled with the following options:
    #
    # -DNPARTITION      do not include the CHOLMOD partition module
    # -DNEXPERT         do not include the functions in SuiteSparseQR_expert.cpp
    # -DHAVE_TBB        enable the use of Intel's Threading Building Blocks
    # -DGPU_BLAS        enable the use of the CUDA BLAS

    SPQR_CONFIG ?= $(GPU_CONFIG)

    # to compile with Intel's TBB, use TBB=-ltbb -DSPQR_CONFIG=-DHAVE_TBB
    TBB ?=
    # TBB = -ltbb -DSPQR_CONFIG=-DHAVE_TBB

    # TODO: this *mk file should auto-detect the presence of Intel's TBB,
    # and set the compiler flags accordingly.

#===============================================================================
# finalize the CF compiler flags
#===============================================================================

    CF += $(CFOPENMP)

#===============================================================================
# internal configuration
#===============================================================================

    # The user should not have to change these definitions, and they are
    # not displayed by 'make config'

    #---------------------------------------------------------------------------
    # for removing files not in the distribution
    #---------------------------------------------------------------------------

    # remove object files, but keep compiled libraries via 'make clean'
    CLEAN = *.o *.obj *.ln *.bb *.bbg *.da *.tcov *.gcov gmon.out *.bak *.d \
        *.gcda *.gcno *.aux *.bbl *.blg *.log *.toc *.dvi *.lof *.lot

    # also remove compiled libraries, via 'make distclean'
    PURGE = *.so* *.a *.dll *.dylib *.dSYM

    # location of TCOV test output
    TCOV_TMP ?= /tmp

#===============================================================================
# Building the shared and static libraries
#===============================================================================

# How to build/install shared and static libraries for Mac and Linux/Unix.
# This assumes that LIBRARY and VERSION have already been defined by the
# Makefile that includes this file.

AR_TARGET = $(LIBRARY).a

#===============================================================================
# display configuration
#===============================================================================

ifeq ($(LIBRARY),)
    # placeholders, for 'make config' in the top-level SuiteSparse
    LIBRARY=PackageNameWillGoHere
    VERSION=x.y.z
    SO_VERSION=x
endif

# 'make config' lists the primary installation options
config:
	@echo ' '
	@echo '----------------------------------------------------------------'
	@echo 'SuiteSparse package compilation options:'
	@echo '----------------------------------------------------------------'
	@echo ' '
	@echo 'SuiteSparse Version:     ' '$(SUITESPARSE_VERSION)'
	@echo 'SuiteSparse top folder:  ' '$(SUITESPARSE)'
	@echo 'Package:                  LIBRARY=        ' '$(LIBRARY)'
	@echo 'Version:                  VERSION=        ' '$(VERSION)'
	@echo 'SO version:               SO_VERSION=     ' '$(SO_VERSION)'
	@echo 'System:                   UNAME=          ' '$(UNAME)'
	@echo 'Install directory:        INSTALL=        ' '$(INSTALL)'
	@echo 'Install libraries in:     INSTALL_LIB=    ' '$(INSTALL_LIB)'
	@echo 'Install include files in: INSTALL_INCLUDE=' '$(INSTALL_INCLUDE)'
	@echo 'Install documentation in: INSTALL_DOC=    ' '$(INSTALL_DOC)'
	@echo 'Optimization level:       OPTIMIZATION=   ' '$(OPTIMIZATION)'
	@echo 'BLAS library:             BLAS=           ' '$(BLAS)'
	@echo 'LAPACK library:           LAPACK=         ' '$(LAPACK)'
	@echo 'Intel TBB library:        TBB=            ' '$(TBB)'
	@echo 'Other libraries:          LDLIBS=         ' '$(LDLIBS)'
	@echo 'static library:           AR_TARGET=      ' '$(AR_TARGET)'
	@echo 'shared library (full):    SO_TARGET=      ' '$(SO_TARGET)'
	@echo 'shared library (main):    SO_MAIN=        ' '$(SO_MAIN)'
	@echo 'shared library (short):   SO_PLAIN=       ' '$(SO_PLAIN)'
	@echo 'shared library options:   SO_OPTS=        ' '$(SO_OPTS)'
	@echo 'shared library name tool: SO_INSTALL_NAME=' '$(SO_INSTALL_NAME)'
	@echo 'ranlib, for static libs:  RANLIB=         ' '$(RANLIB)'
	@echo 'static library command:   ARCHIVE=        ' '$(ARCHIVE)'
	@echo 'copy file:                CP=             ' '$(CP)'
	@echo 'move file:                MV=             ' '$(MV)'
	@echo 'remove file:              RM=             ' '$(RM)'
	@echo 'pretty (for Tcov tests):  PRETTY=         ' '$(PRETTY)'
	@echo 'C compiler:               CC=             ' '$(CC)'
	@echo 'C++ compiler:             CXX=            ' '$(CXX)'
	@echo 'CUDA compiler:            NVCC=           ' '$(NVCC)'
	@echo 'CUDA root directory:      CUDA_PATH=      ' '$(CUDA_PATH)'
	@echo 'OpenMP flags:             CFOPENMP=       ' '$(CFOPENMP)'
	@echo 'C/C++ compiler flags:     CF=             ' '$(CF)'
	@echo 'LD flags:                 LDFLAGS=        ' '$(LDFLAGS)'
	@echo 'Fortran compiler:         F77=            ' '$(F77)'
	@echo 'Fortran flags:            F77FLAGS=       ' '$(F77FLAGS)'
	@echo 'Intel MKL root:           MKLROOT=        ' '$(MKLROOT)'
	@echo 'Auto detect Intel icc:    AUTOCC=         ' '$(AUTOCC)'
	@echo 'UMFPACK config:           UMFPACK_CONFIG= ' '$(UMFPACK_CONFIG)'
	@echo 'CHOLMOD config:           CHOLMOD_CONFIG= ' '$(CHOLMOD_CONFIG)'
	@echo 'SuiteSparseQR config:     SPQR_CONFIG=    ' '$(SPQR_CONFIG)'
	@echo 'CUDA library:             CUDART_LIB=     ' '$(CUDART_LIB)'
	@echo 'CUBLAS library:           CUBLAS_LIB=     ' '$(CUBLAS_LIB)'
	@echo 'METIS and CHOLMOD/Partition configuration:'
	@echo 'Your METIS library:       MY_METIS_LIB=   ' '$(MY_METIS_LIB)'
	@echo 'Your metis.h is in:       MY_METIS_INC=   ' '$(MY_METIS_INC)'
	@echo 'METIS is used via the CHOLMOD/Partition module, configured as follows.'
	@echo 'If the next line has -DNPARTITION then METIS will not be used:'
	@echo 'CHOLMOD Partition config: ' '$(CONFIG_PARTITION)'
	@echo 'CHOLMOD Partition libs:   ' '$(LIB_WITH_PARTITION)'
	@echo 'CHOLMOD Partition include:' '$(I_WITH_PARTITION)'
ifeq ($(TCOV),yes)
	@echo 'TCOV=yes, for extensive testing only (gcc, g++, vanilla BLAS)'
endif

