pdgui:	
	@echo "PDGui Calculating absorption for ZnO from VASP...."
	@../../../pdgui -nosplash vasp OUTCAR results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for ZnO from VASP...."
	@../../../pdgui -nosplash vasp OUTCAR results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing reference files for ZnO from VASP...."
	@../../../pdgui vasp OUTCAR -script script.py

