-include ../../Makefile.mk
pdgui:	
	@echo "PDGui Calculation for absorption of BaTiO3 from Abinit...."
	@../../../pdgui abinit BaTiO3.out results.xlsx -script script.py
	@../../../checkexcel results.ref.xlsx results.xlsx

test:	
	@echo "Calculating absorption for BaTiO3 from Abinit...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for BaTiO3 from Abinit...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv
