-include ../../Makefile.mk
test:	
	@echo "Calculating size effects for MgO from CASTEP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for size effects in MgO from CASTEP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv

