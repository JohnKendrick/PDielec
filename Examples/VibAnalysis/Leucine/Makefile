vibanalysis:	
	@echo "Vibanalysis Calculating absorption for Leucine from Crystal...."
	@../../../vibanalysis --autosel crystal LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out 
	@diff LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.nma LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.nma.ref

vibanalysis-regenerate:
	@echo "Vibanalysis Regenerating reference files for Leucine from Crystal...."
	@../../../vibanalysis crystal --autosel LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out 
	@mv LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.nma LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.nma.ref
