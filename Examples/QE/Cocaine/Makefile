pdgui:	
	@echo "PDGui Calculating absorption for cocaine salt from Quantum Espresso...."
	@../../../pdgui -nosplash qe phonon.dynG results.xlsx -script script.py -exit
	@../../../vibanalysis --autosel --vmp --novmard qe phonon.dynG 
	@../../../checkexcel results.xlsx results.ref.xlsx
	@diff phonon.nma phonon.nma.ref

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for cocaine salt from Quantum Espresso...."
	@../../../pdgui -nosplash qe phonon.dynG results.xlsx -script script.py -exit
	@../../../vibanalysis --autosel --vmp --novmard qe phonon.dynG 
	@mv results.xlsx results.ref.xlsx
	@mv phonon.nma phonon.nma.ref

pdgui-view:
	@echo "PDGui Regenerating reference files for cocaine salt from Quantum Espresso...."
	@../../../pdgui qe phonon.dynG -script script.py 

