-include ../../Makefile.mk


pdgui:	
	@echo "PDGui Calculating size effects for BaTiO3 from Abinit...."
	@../../../pdgui abinit BaTiO3.out results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for size effects BaTiO3 from Abinit...."
	@../../../pdgui abinit BaTiO3.out results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Regenerating reference files for size effects BaTiO3 from Abinit...."
	@../../../pdgui abinit BaTiO3.out -script script.py

test:	
	@echo "Calculating size effects for BaTiO3 from Abinit...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for size effects BaTiO3 from Abinit...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv
