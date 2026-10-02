#CC = f95
CC = gfortran -w -static-libgfortran # gnu fortran compiler
#CC = pgf95                           # portland fortran compiler
#CC = ifort -fp-model precise          # intel fortran compiler
#CC = g77 -w
PCC = mpif90 -lfmpich -static -O3 -parallel # wrapper compiler producing parallel code
AUTOS=asdeck.f                        # autostructure source code 
PAUTOS=pasdeck.f                      # autostructure source code (parallel version)
AUTOSPP=mdrcs13.f                     # general postprocessor source code (outdated)
AUTOSDR=adasdr.f                      # DR postprocessor source code
AUTOSRR=adasrr.f                      # RR postprocessor source code

note:
	@echo ' '
	@echo '  Define the variables AUTOS, AUTOSPP etc. in the makefile as appropriate.'
	@echo '  All files should be in the directory $(HOME)/autostructure'
	@echo ' '
	@echo '  Usage:_'
	@echo ' '
	@echo '  The following make commands are required for a minimal installation'
	@echo '  make AUTOS        ... to compile $(AUTOS) (the autostructure code)'
	@echo '  make AUTOSDR      ... to compile $(AUTOSDR)   (the DR postprocessor)'
	@echo '  make install-all  ... to make soft links in your ~/bin directory'
	@echo '  make doc          ... to extract detailed information on namelists from sources'
	@echo ' '
	@echo '  The following make commands are optional:'
	@echo '  make AUTOSRR      ... to compile $(AUTOSRR) (the RR postprocessor)'
	@echo '  make AUTOSPP      ... to compile $(AUTOSPP) (an outdated general postprocessor)'
	@echo '  make PAUTOS       ... to compile $(PAUTOS)(a parallelized version of AUTOS)'
	@echo '  make HF           ... to compile Hartree-Fock programs for the generation of external wavefunctions'
	@echo ' '
	@echo '  make clean        ... to delete some intermediate files'
	@echo '  make distclean    ... to delete all files made'
	@echo ' '
	@echo '  for updates check at http://amdpp.phys.strath.ac.uk/autos'
	@echo ' '

.SUFFIXES: .o .f

.f.o:
	$(CC) -c $*.f

AUTOS: $(AUTOS)	
	$(CC) -o AUTOS.exe $(AUTOS)

PAUTOS: $(PAUTOS)	
	$(PCC) -o PAUTOS.exe $(PAUTOS)

AUTOSPP: $(AUTOSPP)
	$(CC) -o AUTOSPP.exe $(AUTOSPP)

AUTOSDR: $(AUTOSDR)
	$(CC) -o AUTOSDR.exe $(AUTOSDR)

AUTOSRR: $(AUTOSRR)
	$(CC) -o AUTOSRR.exe $(AUTOSRR)

HF:
	$(CC) -o AUTOSHF.exe  FroeseFischerHF.f
	$(CC) -o AUTOSHFR.exe HFreformat.f

install-all: 
	@echo ' '
	@echo ' The following commands are defined:'
	@echo ' '
	@ln -fs $(PWD)/autos.sh $(HOME)/bin/autos 
	@echo ' autos:    AUTOSTRUCTURE'
	@ln -fs $(PWD)/autosdr.sh $(HOME)/bin/autosdr
	@echo ' autosdr:  DR postprocessor'
	@ln -fs $(PWD)/autosrr.sh $(HOME)/bin/autosrr
	@echo ' autosrr:  RR postprocessor'
	@ln -fs $(PWD)/autospp.sh $(HOME)/bin/autospp
	@echo ' autospp:  general postprocessor (not any longer maintained)'
	@ln -fs $(PWD)/autoshf.sh $(HOME)/bin/autoshf
	@echo ' autoshf:  Hartree-Fock wavefunctions (Froese-Fischer)'
	@echo '  '
	@echo ' The following commands require the "qsub" batch job system'
	@echo ' '
	@ln -fs $(PWD)/qautos.sh $(HOME)/bin/qautos
	@echo ' qautos:   qsub AUTOSTRUCTURE job'
	@ln -fs $(PWD)/qautos.sh $(HOME)/bin/qautosdr
	@echo ' qautosdr: qsub AUTOSTRUCTURE DR postprocessor job'
	@ln -fs $(PWD)/qpautos.sh $(HOME)/bin/pautos
	@echo ' pautos:   qsub AUTOSTRUCTURE multi-processor job'
	@ln -fs $(PWD)/qautospp.sh $(HOME)/bin/qautospp
	@echo ' qautospp: qsub postprocessor job'
	@echo ' '

.PHONY: doc
doc:
	@awk -f doc/docpp.awk $(AUTOSPP) > doc/AUTOSPPnamelists
	@awk -f doc/docpp.awk $(AUTOSDR) > doc/AUTOSDRnamelists 
	@awk -f doc/docpp.awk $(AUTOSRR) > doc/AUTOSRRnamelists 
	@echo ' '
	@echo ' For documentation see WRITEUP, AUTOSDRnamelists,  AUTOSRRnamelists, and AUTOSPPnamelists,  '
	@echo ' '

.PHONY: clean
clean:
	@rm -f *.o *~ *.mod

.PHONY: disclean
distclean:
	@rm -f *.o *~ *.mod
	@rm -f *.ex
	@rm -f *.ps*
	@echo ' '
	@echo ' all compiled files deleted
	@echo ' '
