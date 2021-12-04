pdgui:	
	@echo "PDGui Calculating absorption for Na2(SO4)2 from Crystal...."
	@../../../pdgui -nosplash crystal Na2SO42.out results.xlsx -script script.py -exit
	@../../../checkexcel results.ref.xlsx results.xlsx

benchmark:	
	@echo "PDGui benchmark absorption for Na2(SO4)2 from Crystal...."
	@../../../pdgui -nosplash crystal Na2SO42.out -script script.py -exit

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for Na2(SO4)2 from Crystal...."
	@../../../pdgui -nosplash crystal Na2SO42.out results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Viewing reference files for Na2(SO4)2 from Crystal...."
	@../../../pdgui crystal Na2SO42.out -script script.py

