OBJECTS = hydrocal.o radrate.o clebsch.o polari.o osci.o lifetime.o sigmarr.o RRphotons.o fraction.o ratecoef.o hydromath.o gaussint.o rydserie.o fele.o drmodel.o readxsec.o convolute.o fieldion.o fieldion2.o numerov.o kinema.o Stark.o BetheBloch.o autoso1.o Faddeeva_w.o peakfunctions.o stripping.o TOPbase.o burgess.o CoolerMonteCarlo.o 

EIGEN3OBJECTS = PIPE-MonteCarlo.o

LIBRARIES = -lm 

CC = g++

CFLAGS = -std=c++0x  -Wall -Wno-unused-result -Wno-format -Wno-unused-variable -O3

PREFIX=         $(HOME)/bin

ifdef EIGEN
	CFLAGS += -DEIGEN
	OBJECTS += $(EIGEN3OBJECTS)
endif

.SUFFIXES: .o .cxx

.cxx.o:
	$(CC) $(CFLAGS) -c $*.cxx


########################### help message (first target)
.PHONY: help
help:
	@echo 'Please type ...'
	@echo 'make hydrocal                 to compile the hydrocal source files'
	@echo 'make hydrocal-eigen           to compile hydrocal with including the modules that require libeigen3-dev'
	@echo 'make install                  to install the hydrocal program in the ~/bin folder'
	@echo 'make install-html             to generate the html documentation using doxygen'
	@echo 'make install-pdf              to generate the pdf documentation using doxygen'
	@echo 'make clean                    to remove all intermediate files'
	@echo 'make distclean                to remove all files that were created during installation'
	@echo 'make dist                     to pack all files of the installation into ../hydrocal.tar.gz' 

########################### create header file defing the curretn SVN revision
.PHONY: svnrevision
svnrevision:
	echo "#define SVNrevision \""$(shell svnversion -n .)\" > svnrevision.h 

########################### select compilation of the hydrocal program
.PHONY: hydrocal
hydrocal: cchydrocal

########################### select compilation the hydrocal program with the modules relying on eigen3 
.PHONY: hydrocal-eigen
hydrocal-eigen:
	make EIGEN=yes cchydrocal

########################### compile the hydrocal program
.PHONY: cchydrocal
cchydrocal: svnrevision $(OBJECTS)
	touch hydrocal.cxx
	$(CC) -o hydrocal.exe $(OBJECTS) $(LIBRARIES) 


########################### delete all temporary files
.PHONY: clean
clean:
	rm -f *.o *.log *.aux *.out *.core core *# *~

########################### delete all files that were created during installation
.PHONY: distclean
distclean: clean
	rm -f *.tar
	rm -f *.exe
	rm -rf doc/html
	rm -rf doc/latex

########################### install the executable
.PHONY: install
install: hydrocal
	ln -fs $(PWD)/hydrocal.exe $(PREFIX)/hydrocal

########################### install the html documentation
.PHONY: install-html
install-html:
	echo "SVN revision: "$(shell svnversion -n .) > doc/svnversion.h
	doxygen doc/hydrocal_html.dox

########################### install the latex-pdf documentation
.PHONY: install-pdf
install-pdf:
	echo "SVN revision: "$(shell svnversion -n .) > doc/svnversion.h
	doxygen doc/hydrocal_pdf.dox
	$(MAKE) -C doc/latex/ # invoke makefile created by doxygen
	cp doc/latex/refman.pdf doc/hydrocal_manual.pdf

########################## pack all files of the installation into the compressed tar-file ../hydrocal.tar.gz
.PHONY: dist
dist:	distclean
	tar --transform='s,.,hydrocal,' -cvf hydrocal.tar --exclude-vcs $(@D)
	gzip hydrocal.tar
	mv hydrocal.tar.gz ../
