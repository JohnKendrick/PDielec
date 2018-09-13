pdgui:	
	@echo "PDGui Calculating absorption for ZnO (NOECKART settings) from Crystal...."
	@../../../../pdgui -nosplash crystal ZnO_noeckart.out results.xlsx -script script.py -exit
	@../../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for ZnO (NOECKART settings) from Crystal ...."
	@../../../../pdgui -nosplash crystal ZnO_noeckart.out results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing reference files for ZnO (NOECKART settings) from Crystal ...."
	@../../../../pdgui crystal ZnO_noeckart.out -script script.py 
