vibanalysis:	
	@echo "Vibanalysis calculating absorption for cocaine salt from Quantum Espresso...."
	@../../../vibanalysis --autosel --vmp --novmard qe phonon.dynG 
	@diff phonon.nma phonon.nma.ref

vibanalysis-regenerate:
	@echo "Vibanalysis Regenerating reference files for cocaine salt from Quantum Espresso...."
	@../../../vibanalysis --autosel --vmp --novmard qe phonon.dynG 
	@mv phonon.nma phonon.nma.ref
