-include ../../Makefile.mk
test:	
	@echo "Analysing molecular modes for Isoleucine from CASTEP...."
	@bash command.sh ${PYTHON_EXE} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating molecular modes for Isoleucine from CASTEP...."
	@bash command.sh ${PYTHON_EXE} > command.ref.out
	@mv command.csv command.ref.csv

pdgui:	
	@echo "PDGui on Isoleucine from CASTEP...."
	@../../../pdgui.py castep phonon.castep results.xlsx -script script.py -exit 2> /dev/null 1> pdgui.out
	@../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:	
	@echo "PDGui regenerating Isoleucine from CASTEP...."
	@../../../pdgui.py castep phonon.castep results.xlsx -script script.py -exit 2> /dev/null 1> pdgui.out
	@mv results.xlsx results.ref.xlsx
