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
all:
	( cd Lib    ; $(MAKE) )

# clean, and then remove compiled programs and library archives
clean:
	( cd Lib    ; $(MAKE) clean )
