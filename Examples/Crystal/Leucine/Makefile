pdgui:	
	@echo "PDGui Calculating absorption for Leucine from Crystal...."
	@../../../pdgui -nosplash crystal LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out results.xlsx -script script.py -exit 
	@../../../checkexcel results.ref.xlsx results.xlsx

benchmark:	
	@echo "PDGui benchmark absorption for Leucine from Crystal...."
	@../../../pdgui -nosplash crystal LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out -script script.py -exit 

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for Leucine from Crystal...."
	@../../../pdgui -nosplash crystal LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out results.xlsx -script script.py -exit 
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing reference files for Leucine from Crystal...."
	@../../../pdgui crystal LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out -script script.py
