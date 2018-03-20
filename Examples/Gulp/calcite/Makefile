-include ../../Makefile.mk
pdgui:	
	@echo "PDGui Calculating absorption for Calcite from GULP...."
	@../../../pdgui gulp calcite.gout results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for Calcite from GULP...."
	@../../../pdgui gulp calcite.gout results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx
test:	
	@echo "Calculating absorption for Calcite from GULP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for Calcite from GULP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv
