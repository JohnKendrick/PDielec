pdgui:	
	@echo "PDGui Calculating absorption for Leucine from Crystal...."
	@../../../pdgui -nosplash crystal LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out results.xlsx -script script.py -exit 
	@../../../vibanalysis --autosel crystal LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out 
	@../../../checkexcel results.ref.xlsx results.xlsx
	@diff LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.nma LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.nma.ref

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for Leucine from Crystal...."
	@../../../pdgui -nosplash crystal LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out results.xlsx -script script.py -exit 
	@../../../vibanalysis crystal --autosel LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out 
	@mv results.xlsx results.ref.xlsx
	@mv LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.nma LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.nma.ref

pdgui-view:
	@echo "PDGui Viewing reference files for Leucine from Crystal...."
	@../../../pdgui crystal LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out -script script.py
