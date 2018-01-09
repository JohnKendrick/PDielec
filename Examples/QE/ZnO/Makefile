-include ../../Makefile.mk
test:	
	@echo "Calculating absorption for ZnO from Quantum Espresso...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for ZnO from Quantum Espresso...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv
