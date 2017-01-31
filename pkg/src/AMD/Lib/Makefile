#-------------------------------------------------------------------------------
# AMD Lib/Makefile
#-------------------------------------------------------------------------------

LIBRARY = libamd

all: static

include ../../SuiteSparse_config/SuiteSparse_config.mk

C = $(CC) $(CF) -I../Include -I../../SuiteSparse_config

#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------

AMD = amd_aat amd_1 amd_2 amd_dump amd_postorder amd_defaults \
        amd_post_tree  \
	amd_order amd_control amd_info amd_valid amd_preprocess

INC = ../Include/amd.h ../Include/amd_internal.h \
      ../../SuiteSparse_config/SuiteSparse_config.h

#-------------------------------------------------------------------------------
# object files for each version
#-------------------------------------------------------------------------------

AMDI = $(addsuffix .o, $(subst amd_,amd_i_,$(AMD)))
OBJ = $(AMDI)

#-------------------------------------------------------------------------------
# compile each int and long routine (with no real/complex version)
#-------------------------------------------------------------------------------

amd_i_%.o: ../Source/amd_%.c $(INC)
	$(C) -DDINT -c $< -o $@

static: $(AR_TARGET)

$(AR_TARGET): $(OBJ)
	$(ARCHIVE) $@ $^
	- $(RANLIB) $@

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

clean:
	- $(RM) -r $(CLEAN)
