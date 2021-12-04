pdgui:	
	@echo "PDGui calculation absorption for ZnO (default settings) from Crystal...."
	@../../../../pdgui -nosplash crystal ZnO_default.out results.xlsx -script script.py -exit
	@../../../../checkexcel results.ref.xlsx results.xlsx

benchmark:	
	@echo "PDGui benchmark absorption for ZnO (default settings) from Crystal...."
	@../../../../pdgui -nosplash crystal ZnO_default.out -script script.py -exit

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for ZnO (default settings) from Crystal ...."
	@../../../../pdgui -nosplash crystal ZnO_default.out results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing reference files for ZnO (default settings) from Crystal ...."
	@../../../../pdgui crystal ZnO_default.out -script script.py 

