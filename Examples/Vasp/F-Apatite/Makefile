-include ../../Makefile.mk
pdgui:	
	@echo "PDGui Calculating absorption for F-Apatite from VASP...."
	@../../../pdgui vasp OUTCAR results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for F-Apataite from CASTEP...."
	@../../../pdgui vasp OUTCAR results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing reference files for F-Apataite from CASTEP...."
	@../../../pdgui vasp OUTCAR -script script.py

test:	
	@echo "Calculating absorption for F-Apatite from VASP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for F-Apataite from CASTEP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv
