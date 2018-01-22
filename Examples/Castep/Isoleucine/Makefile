-include ../../Makefile.mk
test:	
	@echo "Analysing molecular modes for Isoleucine from CASTEP...."
	@bash command.sh ${PYTHON_EXE} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating molecular modes for Isoleucine from CASTEP...."
	@bash command.sh ${PYTHON_EXE} > command.ref.out
	@mv command.csv command.ref.csv
