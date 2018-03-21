-include ../../Makefile.mk
pdgui:	
	@echo "PDGui Calculating absorption for Na2(SO4)2 from Quantum Espresso...."
	@../../../pdgui qe Na2SO42.dynG results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for Na2(SO4)2 from Quantum Espresso...."
	@../../../pdgui qe Na2SO42.dynG results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Regenerating reference files for Na2(SO4)2 from Quantum Espresso...."
	@../../../pdgui qe Na2SO42.dynG -script script.py

test:	
	@echo "Calculating absorption for Na2(SO4)2 from Quantum Espresso...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for Na2(SO4)2 from Quantum Espresso...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv
