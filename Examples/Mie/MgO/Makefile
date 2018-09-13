pdgui:	
	@echo "PDGui Calculating Mie absorption for MgO from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for Mie absorption in MgO from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:	
	@echo "PDGui Calculating Mie absorption for MgO from CASTEP...."
	@../../../pdgui castep phonon.castep -script script.py 

