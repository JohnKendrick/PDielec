vibanalysis:	
	@echo "Vibanalysis on Isoleucine from CASTEP...."
	@../../../vibanalysis --autosel castep phonon.castep 
	@diff phonon.nma phonon.nma.ref

vibanalysis-regenerate:	
	@echo "Vibanalysis regenerating Isoleucine from CASTEP...."
	@../../../vibanalysis --autosel castep phonon.castep 
	@mv phonon.nma phonon.nma.ref
