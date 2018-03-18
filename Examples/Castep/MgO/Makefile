-include ../../Makefile.mk
test:	
	@echo "Calculating absorption for MgO from CASTEP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for MgO from CASTEP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv

pdgui:	
	@echo "PDGui Calculating absorption for MgO from CASTEP...."
	@../../../pdgui castep phonon.castep results.xlsx -script script.py -exit 2> /dev/null 1> pdgui.out
	@../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:	
	@echo "PDGui Calculating absorption for MgO from CASTEP...."
	@../../../pdgui castep phonon.castep results.xlsx -script script.py -exit 2> /dev/null 1> pdgui.out
	@mv results.xlsx results.ref.xlsx

