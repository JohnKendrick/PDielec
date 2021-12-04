pdgui:	
	@echo "PDGui Calculating size effects for ZnO from VASP...."
	@../../../pdgui -nosplash vasp OUTCAR results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

benchmark:	
	@echo "PDGui benchmark size effects for ZnO from VASP...."
	@../../../pdgui -nosplash vasp OUTCAR -script script.py -exit

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for size effects in ZnO from VASP...."
	@../../../pdgui -nosplash vasp OUTCAR results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Regenerating reference files for size effects in ZnO from VASP...."
	@../../../pdgui vasp OUTCAR -script script.py 

