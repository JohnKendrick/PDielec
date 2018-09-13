pdgui:
	@echo "PDGui calculation for absorption of AlAs from Abinit...."
	@../../../pdgui -nosplash abinit AlAs.out results.xlsx -script script.py -exit
	@../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for AlAs from Abinit...."
	@../../../pdgui -nosplash abinit AlAs.out results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "Viewing pdgui reference files for AlAs from Abinit...."
	@../../../pdgui abinit AlAs.out -script script.py 

