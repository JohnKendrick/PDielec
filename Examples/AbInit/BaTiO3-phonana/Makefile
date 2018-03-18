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

pdgui:	
	@echo "PDgui analysing phonon modes for BaTiO3 from Abinit...."
	@../../../pdgui abinit BaTiO3.out results.xlsx -exit
	@../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:
	@echo "Regenerating phonon mode files for BaTiO3 from Abinit...."
	@../../../pdgui abinit BaTiO3.out results.xlsx -exit
	@mv results.xlsx results.ref.xlsx
