-include ../../Makefile.mk
test:	
	@echo "Testing preader for QE...."
	@bash command.sh ${PYTHON_EXE} ${PREADER_OPTIONS} > command.csv
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files from QE...."
	@bash command.sh ${PYTHON_EXE} ${PREADER_OPTIONS} > command.ref.csv
