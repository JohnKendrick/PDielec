pdgui:	
	@echo "PDGui Calculating absorption for L-Aspartic acid from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -exit -script ./script.py
	@../../../checkexcel results.ref.xlsx results.xlsx

benchmark:	
	@echo "PDGui benchmark absorption for L-Aspartic acid from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep -exit -script ./script.py

pdgui-regenerate:	
	@echo "PDGui regenerating reference files for L-Aspartic acid from CASTEP...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -exit -script ./script.py
	@mv results.xlsx results.ref.xlsx

pdgui-view:	
	@echo "PDGui Viewing reference files for L-Aspartic acid from CASTEP...."
	../../../pdgui castep phonon.castep -script ./script.py

