-include ../../Makefile.mk
pdgui:	
	@echo "PDGui Calculating absorption for Na2(SO4)2 from VASP...."
	@../../../pdgui vasp OUTCAR results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for Na2(SO4)2 from VASP...."
	@../../../pdgui vasp OUTCAR results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing reference files for Na2(SO4)2 from VASP...."
	@../../../pdgui vasp OUTCAR -script script.py

test:	
	@echo "Calculating absorption for Na2(SO4)2 from VASP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.out
	@../../../checkcsv command.ref.out command.out
	@../../../checkcsv command.ref.csv command.csv

regenerate:
	@echo "Regenerating reference files for Na2(SO4)2 from VASP...."
	@bash command.sh ${PYTHON_EXE} ${PDIELEC_OPTIONS} > command.ref.out
	@mv command.csv command.ref.csv
