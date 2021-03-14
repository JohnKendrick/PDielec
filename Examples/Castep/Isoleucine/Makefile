pdgui:	
	@echo "PDGui on Isoleucine from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -script script.py -exit
	@../../../vibanalysis --autosel castep phonon.castep 
	@../../../checkexcel results.ref.xlsx results.xlsx
	@diff phonon.nma phonon.nma.ref

pdgui-regenerate:	
	@echo "PDGui regenerating Isoleucine from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -script script.py -exit
	@../../../vibanalysis --autosel castep phonon.castep 
	@mv results.xlsx results.ref.xlsx
	@mv phonon.nma phonon.nma.ref

pdgui-view:	
	@echo "PDGui Viewing Isoleucine from CASTEP...."
	@../../../pdgui castep phonon.castep -script script.py

