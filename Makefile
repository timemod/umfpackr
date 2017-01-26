# This is a gnu makefile with several commands to build, document and test
# the package.  The actual building and installation of the package is achieved
# with the standard R commands R CMD BUOLD and R CMD INSTALL.

PKGDIR=pkg
INSTALL_FLAGS=--no-multiarch --with-keep.source
RCHECKARG=--no-multiarch
R_HOME=$(shell R RHOME)
PKG_CXXFLAGS = `"$(R_HOME)/bin/Rscript" -e "Rcpp:::CxxFlags()"` \
	       -I pkg/src/UFconfig -I pkg/src/UMFPACK/Include -I pkg/src/AMD/Include



# Package name, Version and date from DESCIPTION
PKG=$(shell grep 'Package:' $(PKGDIR)/DESCRIPTION  | cut -d " " -f 2)
PKGTAR=$(PKG)_$(shell grep 'Version' $(PKGDIR)/DESCRIPTION  | cut -d " " -f 2).tar.gz
PKGDATE=$(shell grep 'Date' $(PKGDIR)/DESCRIPTION  | cut -d " " -f 2)
TODAY=$(shell date "+%Y-%m-%d")

OSNAME := $(shell uname | tr A-Z a-z)
ifeq ($(findstring windows, $(OSNAME)), windows)
    OSTYPE = windows
else
    # Linux or MAC OSX
    OSTYPE = unix
endif

help:
	@echo
	@echo "The following targets are available:"
	@echo "   help      - displays this help"
	@echo "   test      - run the tests"
	@echo "   covr      - check package coverage (package covr)"
	@echo "   check     - Run R CMD check $(PKGDIR)"
	@echo "   syntax    - check syntax .cpp files"
	@echo "   document  - run roxygen to generate Rd files and make pdf Reference manual"
	@echo "   mkpkg     - builds source package and checks with --as-cran"
	@echo "   bin       - builds binary package in ./tmp"
	@echo "   install   - install package in .libPaths()[1]"
	@echo "   uninstall - uninstall package from .libPaths()[1]"
	@echo "   clean     - cleans up everything"
	@echo "   flags     - display R config flags and some macros"

# make sure that R's variables are used
# if you don't do this you'll get make's initial values
# gives error doing syntax target
#R_CPPFLAGS=$(shell R CMD config --cppflags)
CC=$(shell R CMD config CC)
CXX=$(shell R CMD config CXX)
CPP_FLAGS=$(shell R CMD config --cppflags)

flags:
	@echo "R_HOME=$(R_HOME)"
	@echo "OSTYPE=$(OSTYPE)"
	@echo "SHELL=$(SHELL)"
	@echo "CPP_FLAGS=$(CPP_FLAGS)"
	@echo "PKG_CXXFLAGS=$(PKG_CXXFLAGS)"
	@echo "PKGDIR=$(PKGDIR)"
	@echo "PKG=$(PKG)"
	@echo "PKGTAR=$(PKGTAR)"
	@echo "PKGDATE=$(PKGDATE)"
	@echo "FC=$(FC)"
	@echo "F77=$(F77)"
	@echo "CC=$(CC)"
	@echo "CPP=$(CPP)"
	@echo "CPP_FLAGS=$(CPP_FLAGS)"
	@R --no-save --quiet --slave -e '.libPaths()'

test:
	Rscript test.R

test_covr:
	Rscript test_covr.R

check: cleanx syntax
	@echo " *** Running R CMD check ***"
	R CMD build $(PKGDIR)
	R CMD check $(RCHECKARG) $(PKGTAR)
	@rm -f  $(PKGTAR)
	@echo "Today                           : $(TODAY)"
	@echo "Checked package description date: $(PKGDATE)"
# 	@Rscript -e 'cat("Installed version date          :",packageDescription("nleqslv", fields="Date"))'
	@echo ""

syntax:
ifneq ($(findstring windows, $(OSNAME)), windows)
	$(CXX) "$(CPP_FLAGS)" $(PKG_CXXFLAGS) -c -fsyntax-only -Wall -pedantic $(PKGDIR)/src/*.c*
else
        # On Windows, the compiler complains about misssig Rcpp.h.
        # I have no idea why (the path to the Rcpp headers is added, see PKG_CXXFLAGS).
        # It appears to have something to do with the make system of Rtools.
	@echo syntax checking does not work yet on Windows
endif

cleanx:
ifneq ($(findstring windows, $(OSNAME)), windows)
# Apple Finder rubbish
	@find . -name '.DS_Store' -delete
	@rm -f $(PKGTAR)
	@rm -fr $(PKG).Rcheck
endif

# build date of package must be at least today
# build source package for submission to CRAN
# after building do a check as CRAN does it
mkpkg: cleanx syntax
	R CMD build $(PKG)
	R CMD check --as-cran $(RCHECKARG) $(PKGTAR)
	@cp -nv $(PKGTAR) archive
	@echo "Today                           : $(TODAY)"
	@echo "Checked package description date: $(PKGDATE)"
# 	@Rscript -e 'cat("Installed version date          :",packageDescription("nleqslv", fields="Date"))'
	@echo ""

bin:
	-@rm -rf tmp
	mkdir tmp
	R CMD INSTALL $(INSTALL_FLAGS) -l ./tmp --build $(PKGDIR)

document:
	-@rm -f $(PKGDIR).pdf
	R -e "roxygen2::update_collate('"$(PKGDIR)"'); devtools::document('"$(PKGDIR)"')"
	R CMD Rd2pdf --batch $(PKGDIR) 2>$(PKGDIR).log

install:
	-@rm -rf tmp
	R CMD INSTALL $(INSTALL_FLAGS) $(PKGDIR)

uninstall:
	R CMD REMOVE $(PKG)

clean:
	rm -fr $(PKGDIR).Rcheck
	rm -fr tmp
	rm -f $(PKGTAR)
	rm -f $(PKGDIR).pdf
	rm -f $(PKGDIR).log
