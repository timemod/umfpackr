#-------------------------------------------------------------------------------
# SuiteSparse_config Makefile
#-------------------------------------------------------------------------------
SHELL=/usr/bin/bash
SUITESPARSE ?= $(realpath $(CURDIR)/..)
export SUITESPARSE

# version of SuiteSparse_config is also version of SuiteSparse meta-package
LIBRARY = SuiteSparse_config
VERSION = 4.5.4
SO_VERSION = 4

default: static

include SuiteSparse_config.mk

ccode: all

all: static

OBJ = SuiteSparse_config.o

SuiteSparse_config.o: SuiteSparse_config.c SuiteSparse_config.h
	$(CC) $(CF) -c SuiteSparse_config.c

static: $(AR_TARGET)

$(AR_TARGET): $(OBJ)
	$(ARCHIVE) $(AR_TARGET) SuiteSparse_config.o
	$(RANLIB) $(AR_TARGET)

clean:
	- $(RM) -r $(CLEAN)


