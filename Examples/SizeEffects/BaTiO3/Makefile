-include ../../Makefile.mk
test:	
	@echo "Calculating size effects for BaTiO3 from Abinit...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for size effects BaTiO3 from Abinit...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv
