-include ../../Makefile.mk
pdgui:	
	@echo "PDGui Calculating absorption for MgO from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -script script.py -exit
	@../../../checkexcel results.ref.xlsx results.xlsx

benchmark:	
	@echo "PDGui benchmark absorption for MgO from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep -script script.py -exit

pdgui-regenerate:	
	@echo "PDGui Regenerating absorption for MgO from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:	
	@echo "PDGui Viewing absorption for MgO from CASTEP...."
	@../../../pdgui castep phonon.castep -script script.py

