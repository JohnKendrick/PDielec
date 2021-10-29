pdgui:	
	@echo "PDGui Calculating Single Crystal for forsterite experiment file...."
	@../../../pdgui -nosplash experiment forsterite.exp results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

benchmark:	
	@echo "PDGui benchmark Single Crystal for forsterite experiment file...."
	@../../../pdgui -nosplash experiment forsterite.exp -script script.py -exit

pdgui-regenerate:
	@echo "PDGui Regenerating Single Crystal for forsterite experiment file...."
	@../../../pdgui -nosplash experiment forsterite.exp results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing Single Crystal for forsterite experiment file...."
	@../../../pdgui experiment forsterite.exp -script script.py

