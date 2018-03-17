-include ../../Makefile.mk
pdgui:	
	@echo "PDGui Calculation of absorption in Na2(SO4)2 from Abinit...."
	@../../../pdgui abinit Na2SO42.out results.xlsx -script script.py -exit
	@../../../checkexcel results.ref.xlsx results.xlsx

test:	
	@echo "Calculating absorption for Na2(SO4)2 from Abinit...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for Na2(SO4)2 from Abinit...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv
