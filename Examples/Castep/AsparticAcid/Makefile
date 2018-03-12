-include ../../Makefile.mk
pdgui:	
	@echo "Calculating absorption for L-Aspartic acid from CASTEP...."
	@bash pdgui.sh > pdgui.out
	@../../../checkexcel results.ref.xlsx results.xlsx

test:	
	@echo "Calculating absorption for L-Aspartic acid from CASTEP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for L-Aspartic acid from CASTEP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv
