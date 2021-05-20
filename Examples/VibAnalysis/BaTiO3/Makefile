vibanalysis:	
	@echo "Vibanalysis calculation for absorption of BaTiO3 from Abinit...."
	@../../../vibanalysis BaTiO3.out
	@diff BaTiO3.nma.ref BaTiO3.nma

vibanalysis-regenerate:	
	@echo "vibanalysis regeneration for absorption of BaTiO3 from Abinit...."
	@../../../vibanalysis BaTiO3.out
	@mv BaTiO3.nma BaTiO3.nma.ref
