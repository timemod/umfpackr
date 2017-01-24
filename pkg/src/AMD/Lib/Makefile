#-------------------------------------------------------------------------------
# AMD Makefile for compiling on Unix systems (for GNU make only)
#-------------------------------------------------------------------------------
CC = gcc
O = o
LIBAMD = libamd.a
MKLIB = ar -rucs
MKLIB_OUT =

C = $(CC) $(CFLAGS) -I../Include -I../../UFconfig

#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------

AMD = amd_aat amd_1 amd_2 amd_dump amd_postorder amd_post_tree amd_defaults \
	amd_order amd_control amd_info amd_valid amd_preprocess

UFCONFIG = ../../UFconfig/UFconfig.h

INC = ../Include/amd.h ../Include/amd_internal.h $(UFCONFIG)

#-------------------------------------------------------------------------------
# object files for each version
#-------------------------------------------------------------------------------

AMDI = $(addsuffix .$O, $(subst amd_,amd_i_,$(AMD)))
AMDL = $(addsuffix .$O, $(subst amd_,amd_l_,$(AMD)))

#-------------------------------------------------------------------------------
# Create the libamd.a library (C versions only)
#-------------------------------------------------------------------------------

$(LIBAMD): amd_global.$O $(AMDI)
	$(MKLIB) $(MKLIB_OUT) $(LIBAMD) $^

#-------------------------------------------------------------------------------
# compile each int and long routine (with no real/complex version)
#-------------------------------------------------------------------------------

amd_global.$O: ../Source/amd_global.c $(INC)
	$(C) -c $< -o $@

amd_i_%.$O: ../Source/amd_%.c $(INC)
	$(C) -c $<  -o $@

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

clean:
	-rm *.o
	-rm *.a

