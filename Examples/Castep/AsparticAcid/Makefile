pdgui:	
	@echo "PDGui Calculating absorption for L-Aspartic acid from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -exit -script ./script.py
	@../../../vibanalysis --autosel phonon.castep
	@../../../checkexcel results.ref.xlsx results.xlsx
	@diff phonon.nma phonon.nma.ref

pdgui-regenerate:	
	@echo "PDGui regenerating reference files for L-Aspartic acid from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -exit -script ./script.py
	@../../../vibanalysis --autosel phonon.castep
	@mv results.xlsx results.ref.xlsx
	@mv phonon.nma phonon.nma.ref

pdgui-view:	
	@echo "PDGui Viewing reference files for L-Aspartic acid from CASTEP...."
	../../../pdgui castep phonon.castep -script ./script.py

