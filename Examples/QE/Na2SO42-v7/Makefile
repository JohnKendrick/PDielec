pdgui:	
	@echo "PDGui Calculating absorption for Na2(SO4)2 from Quantum Espresso...."
	@../../../pdgui -nosplash qe Na2SO42.dynG results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

benchmark:	
	@echo "PDGui benchmark absorption for Na2(SO4)2 from Quantum Espresso...."
	@../../../pdgui -nosplash qe Na2SO42.dynG -script script.py -exit

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for Na2(SO4)2 from Quantum Espresso...."
	@../../../pdgui -nosplash qe Na2SO42.dynG results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Regenerating reference files for Na2(SO4)2 from Quantum Espresso...."
	@../../../pdgui qe Na2SO42.dynG -script script.py

