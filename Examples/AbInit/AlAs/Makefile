-include ../../Makefile.mk
pdgui:
	@echo "PDGui calculation for absorption of AlAs from Abinit...."
	@../../../pdgui abinit AlAs.out results.xlsx -script script.py -exit
	@../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for AlAs from Abinit...."
	@../../../pdgui abinit AlAs.out results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "Viewing pdgui reference files for AlAs from Abinit...."
	@../../../pdgui abinit AlAs.out -script script.py 

test:	
	@echo "Calculating absorption for AlAs from Abinit...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for AlAs from Abinit...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv

