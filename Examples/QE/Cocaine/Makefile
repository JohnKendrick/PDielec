-include ../../Makefile.mk
pdgui:	
	@echo "PDGui Calculating absorption for cocaine salt from Quantum Espresso...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../pdgui qe phonon.dynG results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for cocaine salt from Quantum Espresso...."
	@../../../pdgui qe phonon.dynG results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Regenerating reference files for cocaine salt from Quantum Espresso...."
	@../../../pdgui qe phonon.dynG -script script.py 

test:	
	@echo "Calculating absorption for cocaine salt from Quantum Espresso...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for cocaine salt from Quantum Espresso...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv
