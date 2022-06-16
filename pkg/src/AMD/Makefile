#------------------------------------------------------------------------------
# AMD Makefile
#------------------------------------------------------------------------------
SHELL=/usr/bin/bash
SUITESPARSE ?= $(realpath $(CURDIR)/..)
export SUITESPARSE

default: all

include ../SuiteSparse_config/SuiteSparse_config.mk

demos: all

# Compile all C code.  Do not compile the FORTRAN versions.
# argument -P of command cd is necessary if Rtools42 is used and 
# the package is installed on the CPB file system (the m:-disk).
all:
	( cd -P Lib    ; $(MAKE) )

# clean, and then remove compiled programs and library archives
clean:
	( cd -P Lib    ; $(MAKE) clean )
