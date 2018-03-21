-include ../../../Makefile.mk
pdgui:	
	@echo "PDGui Calculating absorption for ZnO (CPHF) from Crystal...."
	@../../../../pdgui crystal ZnO_CPHF.out results.xlsx -script script.py -exit
	@../../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for ZnO (CPHF) from Crystal ...."
	@../../../../pdgui crystal ZnO_CPHF.out results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing reference files for ZnO (CPHF) from Crystal ...."
	@../../../../pdgui crystal ZnO_CPHF.out -script script.py 

test:	
	@echo "Calculating absorption for ZnO (CPHF) from Crystal...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../../checkcsv command.ref.out command.out
	@../../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for ZnO (CPHF) from Crystal ...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv
