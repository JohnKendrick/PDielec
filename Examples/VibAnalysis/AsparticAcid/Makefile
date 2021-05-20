vibanalysis:	
	@echo "Vibanalysis calculation for L-Aspartic acid from CASTEP...."
	@../../../vibanalysis --autosel phonon.castep
	@diff phonon.nma phonon.nma.ref

vibanalysis-regenerate:	
	@echo "Vibanalysis regenerating reference files for L-Aspartic acid from CASTEP...."
	@../../../vibanalysis --autosel phonon.castep
	@mv phonon.nma phonon.nma.ref

