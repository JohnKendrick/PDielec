-include ../../Makefile.mk
pdgui:	
	@echo "PDGui Calculating absorption for Na2(SO4)2 from Crystal...."
	@../../../pdgui crystal Na2SO42.out results.xlsx -script script.py -exit
	@../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for Na2(SO4)2 from Crystal...."
	@../../../pdgui crystal Na2SO42.out results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

test:	
	@echo "Calculating absorption for Na2(SO4)2 from Crystal...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for Na2(SO4)2 from Crystal...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv
