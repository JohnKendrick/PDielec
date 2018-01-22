-include ../../Makefile.mk
test:	
	@echo "Analysing phonon modes for BaTiO3 from Abinit...."
	@bash command.sh ${PYTHON_EXE} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating phonon mode files for BaTiO3 from Abinit...."
	@bash command.sh ${PYTHON_EXE} > command.ref.out
	@mv command.csv command.ref.csv
