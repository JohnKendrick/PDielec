-include ../../Makefile.mk
test:	
	@echo "Calculating absorption for L-Aspartic acid from CASTEP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for L-Aspartic acid from CASTEP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv

pdgui:	
	@echo "PDGui Calculating absorption for L-Aspartic acid from CASTEP...."
	@../../../pdgui castep phonon.castep results.xlsx -exit -script ./script.py
	@../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:	
	@echo "PDGui regenerating reference files for L-Aspartic acid from CASTEP...."
	../../../pdgui castep phonon.castep results.xlsx -exit -script ./script.py
	@mv results.xlsx results.ref.xlsx

