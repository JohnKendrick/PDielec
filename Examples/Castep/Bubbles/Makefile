-include ../../Makefile.mk
pdgui:	
	@echo "PDGui Calculating bubbles for MgO from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -script script.py -exit
	@../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:	
	@echo "PDGui Regenerating absorption for MgO from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:	
	@echo "PDGui Viewing absorption for MgO from CASTEP...."
	@../../../pdgui castep phonon.castep -script script.py

