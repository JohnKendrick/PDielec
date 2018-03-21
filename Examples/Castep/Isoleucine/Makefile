-include ../../Makefile.mk
pdgui:	
	@echo "PDGui on Isoleucine from CASTEP...."
	@../../../pdgui.py castep phonon.castep results.xlsx -script script.py -exit
	@../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:	
	@echo "PDGui regenerating Isoleucine from CASTEP...."
	@../../../pdgui.py castep phonon.castep results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:	
	@echo "PDGui Viewing Isoleucine from CASTEP...."
	@../../../pdgui.py castep phonon.castep -script script.py

test:	
	@echo "Analysing molecular modes for Isoleucine from CASTEP...."
	@bash command.sh ${PYTHON_EXE} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating molecular modes for Isoleucine from CASTEP...."
	@bash command.sh ${PYTHON_EXE} > command.ref.out
	@mv command.csv command.ref.csv

