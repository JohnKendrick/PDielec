-include ../../Makefile.mk
pdgui:	
	@echo "PDGui Calculation for absorption of BaTiO3 from Abinit...."
	@../../../pdgui abinit BaTiO3.out results.xlsx -script script.py -exit
	@../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:	
	@echo "PDGui regeneration for absorption of BaTiO3 from Abinit...."
	@../../../pdgui abinit BaTiO3.out results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:	
	@echo "PDGui viewing for absorption of BaTiO3 from Abinit...."
	@../../../pdgui abinit BaTiO3.out -script script.py

test:	
	@echo "Calculating absorption for BaTiO3 from Abinit...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for BaTiO3 from Abinit...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv

