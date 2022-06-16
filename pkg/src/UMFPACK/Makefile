#-------------------------------------------------------------------------------
# UMFPACK Makefile
#-------------------------------------------------------------------------------

SUITESPARSE ?= $(realpath $(CURDIR)/..)
export SUITESPARSE

default: all

include ../SuiteSparse_config/SuiteSparse_config.mk

# compile all C code (except hb, fortran, and fortran64), incl. AMD and demos
# argument -P of command cd is necessary if Rtools42 is used and 
# the package is installed on the CPB file system (the m:-disk).
all:
	( cd -P Lib    ; $(MAKE) )

clean:
	( cd -P Lib    ; $(MAKE) clean )
