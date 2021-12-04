pdgui:	
	@echo "PDGui Calculating absorption for ZnO (CPHF) from Crystal...."
	@../../../../pdgui -nosplash crystal ZnO_CPHF.out results.xlsx -script script.py -exit
	@../../../../checkexcel results.ref.xlsx results.xlsx

benchmark:	
	@echo "PDGui benchmark absorption for ZnO (CPHF) from Crystal...."
	@../../../../pdgui -nosplash crystal ZnO_CPHF.out -script script.py -exit

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for ZnO (CPHF) from Crystal ...."
	@../../../../pdgui -nosplash crystal ZnO_CPHF.out results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing reference files for ZnO (CPHF) from Crystal ...."
	@../../../../pdgui crystal ZnO_CPHF.out -script script.py 

