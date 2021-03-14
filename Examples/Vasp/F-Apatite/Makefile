pdgui:	
	@echo "PDGui Calculating absorption for F-Apatite from VASP...."
	@../../../pdgui -nosplash vasp OUTCAR results.xlsx -script script.py -exit
	@../../../vibanalysis --autosel vasp OUTCAR
	@../../../checkexcel results.xlsx results.ref.xlsx
	@diff OUTCAR.nma OUTCAR.nma.ref

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for F-Apataite from CASTEP...."
	@../../../pdgui -nosplash vasp OUTCAR results.xlsx -script script.py -exit
	@../../../vibanalysis --autosel vasp OUTCAR
	@mv results.xlsx results.ref.xlsx
	@mv OUTCAR.nma OUTCAR.nma.ref

pdgui-view:
	@echo "PDGui Viewing reference files for F-Apataite from CASTEP...."
	@../../../pdgui vasp OUTCAR -script script.py

