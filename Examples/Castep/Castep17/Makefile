-include ../../Makefile.mk
pdgui:	
	@echo "Calculating CASTEP17 test...."
	@bash pdgui.sh > pdgui.out
	@../../../checkexcel results.ref.xlsx results.xlsx

test:	
	@echo "Calculating CASTEP17 test...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating CASTEP17 test...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv
