-include ../../Makefile.mk
pdgui:	
	@echo "PDGui Calculating size effects for MgO from CASTEP...."
	@../../../pdgui castep phonon.castep results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for size effects in MgO from CASTEP...."
	@../../../pdgui castep phonon.castep results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing reference files for size effects in MgO from CASTEP...."
	@../../../pdgui castep phonon.castep -script script.py

ttest:	
	@echo "Calculating size effects for MgO from CASTEP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for size effects in MgO from CASTEP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv

