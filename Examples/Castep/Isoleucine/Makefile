pdgui:	
	@echo "PDGui on Isoleucine from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -script script.py -exit
	@../../../checkexcel results.ref.xlsx results.xlsx

benchmark:	
	@echo "PDGui benchmark Isoleucine from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep -script script.py -exit

pdgui-regenerate:	
	@echo "PDGui regenerating Isoleucine from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:	
	@echo "PDGui Viewing Isoleucine from CASTEP...."
	@../../../pdgui castep phonon.castep -script script.py

