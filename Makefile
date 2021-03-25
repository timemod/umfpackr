# This is a gnu makefile with several commands to build, document and test
# the package.  The actual building and installation of the package is achieved
# with the standard R commands R CMD BUOLD and R CMD INSTALL.

PKGDIR=pkg
INSTALL_FLAGS=--no-multiarch --with-keep.source

# Package name, Version and date from DESCIPTION
PKG=$(shell grep Package: $(PKGDIR)/DESCRIPTION  | cut -d " " -f 2)
PKGTAR=$(PKG)_$(shell grep Version $(PKGDIR)/DESCRIPTION  | cut -d " " -f 2).tar.gz
OSTYPE=$(shell Rscript -e "cat(.Platform[['OS.type']])")
RCPP_CXXFLAGS = $(shell Rscript -e "Rcpp:::CxxFlags()")
PKG_CXXFLAGS = $(RCPP_CXXFLAGS) -I pkg/src/SuiteSparse_config \
	       -I pkg/src/UMFPACK/Include -I pkg/src/AMD/Include

ifeq ($(OSTYPE), windows)
# for unknown reason R CMD check --as-cran does not work on Windows
RCHECKARG=--no-multiarch
else
RCHECKARG=--no-multiarch --as-cran
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
	@echo OSTYPE=$(OSTYPE)
	@echo CPP_FLAGS=$(CPP_FLAGS)
	@echo RCPP_CXXFLAGS=$(RCPP_CXXFLAGS)
	@echo PKG_CXXFLAGS=$(PKG_CXXFLAGS)
	@echo PKGDIR=$(PKGDIR)
	@echo PKG=$(PKG)
	@echo PKGTAR=$(PKGTAR)
	@echo CC=$(CC)
	@echo CPP=$(CPP)
	@echo CPP_FLAGS=$(CPP_FLAGS)
	@echo libPaths:
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

syntax:
	$(CXX) $(CPP_FLAGS) $(PKG_CXXFLAGS) -c -fsyntax-only -Wall -pedantic $(PKGDIR)/src/*.c*

cleanx:
ifneq ($(OSTYPE), windows)
# Apple Finder rubbish
	@find . -name '.DS_Store' -delete
	@rm -f $(PKGTAR)
	@rm -fr $(PKG).Rcheck
endif

# build date of package must be at least today
# build source package for submission to CRAN
# after building do a check as CRAN does it
mkpkg: cleanx syntax install_deps
	R CMD build $(PKGDIR)
	R CMD check --as-cran $(RCHECKARG) $(PKGTAR)
	@cp -nv $(PKGTAR) archive
	./drat.sh --pkg=$(PKGTAR)


bin: install_deps
	-@rm -rf tmp
	mkdir tmp
	R CMD INSTALL $(INSTALL_FLAGS) -l ./tmp --build $(PKGDIR)

document: install_deps
	-@rm -f $(PKGDIR)/vignettes/umfpackr_refman.pdf
	R -e "devtools::document('"$(PKGDIR)"')"
	R CMD Rd2pdf --batch $(PKGDIR) -o $(PKGDIR)/vignettes/umfpackr_refman.pdf 2>&1 refman.log

install: install_deps
	-@rm -rf tmp
	R CMD INSTALL $(INSTALL_FLAGS) $(PKGDIR)

install_deps:
	R --slave -f install_deps.R

uninstall:
	R CMD REMOVE $(PKG)

clean:
	-rm -fr $(PKGDIR).Rcheck
	-rm -fr tmp
	-rm -f $(PKGTAR)
	-rm -f $(PKGDIR).pdf
	-rm -f $(PKGDIR).log
	-rm -f $(PKGDIR)/src/*.o
	-rm -f $(PKGDIR)/src/*.a
	-rm -f $(PKGDIR)/src/*.so
	-$(MAKE) -C $(PKGDIR)/src/SuiteSparse_config clean
	-$(MAKE) -C $(PKGDIR)/src/AMD clean
	-$(MAKE) -C $(PKGDIR)/src/UMFPACK clean

