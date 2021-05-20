vibanalysis:	
	@echo "Vibanalysis Calculating absorption for F-Apatite from VASP...."
	@../../../vibanalysis --autosel vasp OUTCAR
	@diff OUTCAR.nma OUTCAR.nma.ref

vibanalysis-regenerate:
	@echo "Vibanalysis Regenerating reference files for F-Apataite from CASTEP...."
	@../../../vibanalysis --autosel vasp OUTCAR
	@mv OUTCAR.nma OUTCAR.nma.ref
