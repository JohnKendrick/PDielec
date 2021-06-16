pdgui:	
	@echo "PDGui Calculating Single Crystal for forsterite experiment file...."
	@../../../pdgui -nosplash experiment forsterite.exp results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating Single Crystal for forsterite experiment file...."
	@../../../pdgui -nosplash experiment forsterite.exp results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing Single Crystal for forsterite experiment file...."
	@../../../pdgui vasp OUTCAR -script script.py
	@../../../pdgui experiment forsterite.exp -script script.py

