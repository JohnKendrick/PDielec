pdgui:	
	@echo "PDGui Calculation of absorption in Na2(SO4)2 from Abinit...."
	@../../../pdgui -nosplash abinit Na2SO42.out results.xlsx -script script.py -exit
	@../../../checkexcel results.ref.xlsx results.xlsx

benchmark:	
	@echo "PDGui benchmark of absorption in Na2(SO4)2 from Abinit...."
	@../../../pdgui -nosplash abinit Na2SO42.out -script script.py -exit

pdgui-regenerate:	
	@echo "PDGui Regenerating of absorption in Na2(SO4)2 from Abinit...."
	@../../../pdgui -nosplash abinit Na2SO42.out results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:	
	@echo "PDGui Viewing of absorption in Na2(SO4)2 from Abinit...."
	@../../../pdgui abinit Na2SO42.out -script script.py

