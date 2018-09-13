pdgui:	
	@echo "PDGui Calculating absorption for ZnO from Quantum Espresso...."
	@../../../pdgui -nosplash qe zno.ph.dynG results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for ZnO from Quantum Espresso...."
	@../../../pdgui -nosplash qe zno.ph.dynG results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing reference files for ZnO from Quantum Espresso...."
	@../../../pdgui qe zno.ph.dynG -script script.py 

