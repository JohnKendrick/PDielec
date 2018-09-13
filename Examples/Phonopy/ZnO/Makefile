pdgui:	
	@echo "PDGui Calculating absorption for ZnO from PHONOPY...."
	@../../../pdgui -nosplash phonopy OUTCAR results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for ZnO from PHONOPY...."
	@../../../pdgui -nosplash phonopy OUTCAR results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Regenerating reference files for ZnO from PHONOPY...."
	@../../../pdgui phonopy OUTCAR -script script.py 

