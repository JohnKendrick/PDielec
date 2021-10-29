pdgui:	
	@echo "PDGui calculating CASTEP17 test...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -exit -script ./script.py
	@../../../checkexcel results.ref.xlsx results.xlsx

benchmark:	
	@echo "PDGui benchmark CASTEP17 test...."
	@../../../pdgui -nosplash castep phonon.castep -exit -script ./script.py

pdgui-regenerate:	
	@echo "PDGui Regenerating CASTEP17 test...."
	@../../../pdgui -nosplash castep phonon.castep results.xlsx -exit -script ./script.py
	@mv results.xlsx results.ref.xlsx

pdgui-view:	
	@echo "PDGui Viewing CASTEP17 test...."
	@../../../pdgui castep phonon.castep -script ./script.py

