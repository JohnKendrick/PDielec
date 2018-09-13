pdgui:	
	@echo "PDGui Calculating absorption for cocaine salt from Quantum Espresso...."
	@../../../pdgui -nosplash qe phonon.dynG results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for cocaine salt from Quantum Espresso...."
	@../../../pdgui -nosplash qe phonon.dynG results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Regenerating reference files for cocaine salt from Quantum Espresso...."
	@../../../pdgui qe phonon.dynG -script script.py 

