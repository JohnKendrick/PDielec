-include ../../Makefile.mk
pdgui:	
	@echo "PDGui Calculating absorption for MgO from CASTEP...."
	@python castep phonon.castep results.xlsx -script script.py > pdgui.out
	@../../../checkexcel results.ref.xlsx results.xlsx

test:	
	@echo "Calculating absorption for MgO from CASTEP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for MgO from CASTEP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv

