-include ../../Makefile.mk
pdgui:	
	@echo "PDGui Calculating absorption for Leucine from Crystal...."
	@../../../pdgui crystal LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out results.xlsx -script script.py -exit 
	@../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for Leucine from Crystal...."
	@../../../pdgui crystal LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out results.xlsx -script script.py -exit 
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing reference files for Leucine from Crystal...."
	@../../../pdgui crystal LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out -script script.py

test:	
	@echo "Calculating absorption for Leucine from Crystal...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for Leucine from Crystal...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv

