#-------------------------------------------------------------------------------
# UMFPACK Makefile
#-------------------------------------------------------------------------------

SUITESPARSE ?= $(realpath $(CURDIR)/..)
export SUITESPARSE

default: all

include ../SuiteSparse_config/SuiteSparse_config.mk

# compile all C code (except hb, fortran, and fortran64), incl. AMD and demos
all:
	( cd -P Lib    ; $(MAKE) )

clean:
	( cd -P Lib    ; $(MAKE) clean )
