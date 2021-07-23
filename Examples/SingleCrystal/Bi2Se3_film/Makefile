pdgui:	
	@echo "PDGui Calculating Single Crystal for Bi2Se3 VASP calculation...."
	@../../../pdgui -nosplash vasp OUTCAR results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating Single Crystal for Bi2Se3 VASP calculation...."
	@../../../pdgui -nosplash vasp OUTCAR results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing Single Crystal for Bi2Se3 VASP calculation...."
	@../../../pdgui vasp OUTCAR -script script.py

