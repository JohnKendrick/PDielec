-include ../../Makefile.mk
pdgui:	
	@echo "PDGui Calculating absorption for MgO from CASTEP...."
	@../../../pdgui castep phonon.castep results.xlsx -script script.py -exit
	@../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:	
	@echo "PDGui Calculating absorption for MgO from CASTEP...."
	@../../../pdgui castep phonon.castep results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:	
	@echo "PDGui Viewing absorption for MgO from CASTEP...."
	@../../../pdgui castep phonon.castep -script script.py

test:	
	@echo "Calculating absorption for MgO from CASTEP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for MgO from CASTEP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv

