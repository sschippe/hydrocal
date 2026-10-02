# SPDX-License-Identifier: MIT

OBJECTS = hydrocal.o radrate.o clebsch.o polari.o osci.o lifetime.o RRDRxsec.o RRphotons.o CoolerFractions.o RRDRratecoef.o matrix.o hydromath.o gaussint.o rydserie.o fele.o DRmodelspe.o readxsec.o RRDRconvolve.o fieldion.o fieldion2.o numerov.o kinema.o Stark.o BetheBloch.o autoso1.o Faddeeva_w.o peakfunctions.o stripping.o TOPbase.o burgess.o Cooler.o CoolerMonteCarlo.o FranckCondon.o storagetime.o PIPE-MonteCarlo.o CoolerDeconvolve.o readadas.o readstrz.o

LIBRARIES = -fopenmp -lm 

GUIOBJECTS = hydrocalGUI.o

vpath %.cxx src
vpath %.h src

GUILIBRARIES = -lfltk

CC = g++
#CC = x86_64-w64-mingw32-gcc -static

CFLAGS = -std=c++17 -Wall -Wno-unused-result -Wno-format -Wno-unused-variable -O3 -fopenmp -I. -Isrc
ccgui: CFLAGS += -I/usr/include/cairo/ -I. -Isrc

# Revision counter, continuing the Subversion numbering (see CMakeLists.txt).
# Read from the REVISION file, which the pre-commit hook in .githooks/ bumps
# (activate it with: git config core.hooksPath .githooks).  If the file is
# missing, reconstruct the value from the history instead.
REVISION_BASE = 2122
REVISION = $(shell tr -d ' \t\r\n' < REVISION 2>/dev/null | \
              grep -E '^[0-9]+$$' || \
              awk -v b=$(REVISION_BASE) -v n="$$(git rev-list --first-parent --count HEAD 2>/dev/null)" \
                  'BEGIN { print (n + 0 > 0) ? b + n - 1 : b }')
GITID = $(shell git describe --tags --always --dirty 2>/dev/null || echo unknown)

PREFIX=         ${HOME}
BIN_PREFIX=	${PREFIX}/bin

ifdef BOOST
	CFLAGS += -DuseBOOST
endif
.SUFFIXES: .o .cxx

.cxx.o:
	$(CC) $(CFLAGS) -c $<


########################### help message (first target)
.PHONY: help
help:
	@echo 'Please type ...'
	@echo 'make hydrocal                 to compile the hydrocal source files'
	@echo 'make gui                      to compile the hydrocal GUI (requires libfltk1.3-dev)'
	@echo 'make install                  to install the hydrocal and hcGUI programs in the ~/bin folder'
	@echo 'make test                     run test with inputs from /doc/examples'
	@echo 'make doc-html                 to generate the html documentation using doxygen'
	@echo 'make doc-pdf                  to generate the pdf documentation using doxygen'
	@echo 'make clean                    to remove all intermediate files'
	@echo 'make distclean                to remove all files that were created during installation'
	@echo 'make dist                     to pack all files of the installation into ../hydrocal.tar.gz'
	@echo '---'
	@echo 'NOTE: CMake build is also available:'
	@echo '  mkdir build && cd build && cmake .. && make'
	@echo '  Or with optional dependencies:'
	@echo '  cmake .. -DUSE_BOOST=ON' 

########################### create header file defining the current git revision
.PHONY: buildinfo
buildinfo:
	@echo "#ifndef HYDROCAL_BUILDINFO_H" > buildinfo.h
	@echo "#define HYDROCAL_BUILDINFO_H" >> buildinfo.h
	@echo "#define HYDROCAL_REVISION \"$(REVISION)\"" >> buildinfo.h
	@echo "#define HYDROCAL_GITID \"$(GITID)\"" >> buildinfo.h
	@echo "#endif" >> buildinfo.h 

########################### select compilation the hydrocal program using boost library
.PHONY: hydrocal-boost
hydrocal-boost:
	make BOOST=yes hydrocal

########################### compile the hydrocal program
.PHONY: hydrocal
hydrocal: buildinfo $(OBJECTS)
	$(CC) -o hydrocal.exe $(OBJECTS) $(LIBRARIES) 

########################### run tests
.PHONY: test
test: hydrocal
	@$(RM) -rf test
	@mkdir -p test
	@echo "Running tests ..."
	@hcins=""; for f in doc/examples/*.hcin; do hcins="$$hcins $$f"; done; \
	passed=0; failed=0; \
	for hcin in $$hcins; do \
	  name=$$(basename $$hcin .hcin); \
	  res="test/.result_$$name"; \
	  ./run_test.sh $$hcin $$res; \
	  read line < $$res; \
	  printf "  %-12s ... %s\n" "$$name" "$$line"; \
	  case "$$line" in ok*) passed=$$((passed+1));; *) failed=$$((failed+1));; esac; \
	  rm -f $$res; \
	done; \
	echo ""; \
	echo "$$passed passed, $$failed failed"

########################### select compilation the hydrocal GUI
.PHONY: gui
gui:
	make ccgui

########################### compile the hydrocal GUI
.PHONY: ccgui
ccgui: buildinfo $(GUIOBJECTS)
	touch src/hydrocalGUI.cxx
	$(CC) -o hydrocalGUI.exe $(GUIOBJECTS) $(GUILIBRARIES) 

########################### delete all temporary files
.PHONY: clean
clean:
	rm -f *.o *.log *.aux *.out *.core core *# *~ src/*.o

########################### delete all files that were created during installation
.PHONY: distclean
distclean: clean
	rm -f *.tar
	rm -f *.exe
	rm -rf doc/html
	rm -rf doc/latex

########################### install the executables
.PHONY: install
install: 
	ln -fs $(PWD)/hydrocal.exe $(BIN_PREFIX)/hydrocal
	ln -fs $(PWD)/hydrocalGUI.exe $(BIN_PREFIX)/hydrocalgui

########################### install the html documentation
.PHONY: doc-html
doc-html:
	@echo "git revision: $(REVISION)" > doc/buildinfo.txt
	doxygen doc/hydrocal_html.dox

########################### install the latex-pdf documentation
.PHONY: doc-pdf
doc-pdf:
	@echo "git revision: $(REVISION)" > doc/buildinfo.txt
	doxygen doc/hydrocal_pdf.dox
	$(MAKE) -C doc/latex/ # invoke makefile created by doxygen
	cp doc/latex/refman.pdf doc/hydrocal_manual.pdf

########################## pack all files of the installation into the compressed tar-file ../hydrocal.tar.gz
.PHONY: dist
dist:	distclean
	tar --transform='s,.,hydrocal,' -cvf hydrocal.tar --exclude-vcs $(@D)
	gzip hydrocal.tar
	mv hydrocal.tar.gz ../

