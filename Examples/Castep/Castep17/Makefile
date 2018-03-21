-include ../../Makefile.mk
pdgui:	
	@echo "PDGui calculating CASTEP17 test...."
	@../../../pdgui castep phonon.castep results.xlsx -exit -script ./script.py
	@../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:	
	@echo "PDGui Regenerating CASTEP17 test...."
	@../../../pdgui castep phonon.castep results.xlsx -exit -script ./script.py
	@mv results.xlsx results.ref.xlsx

pdgui-view:	
	@echo "PDGui Viewing CASTEP17 test...."
	@../../../pdgui castep phonon.castep -script ./script.py

test:	
	@echo "Calculating CASTEP17 test...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating CASTEP17 test...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv

