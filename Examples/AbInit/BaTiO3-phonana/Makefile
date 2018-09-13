pdgui:
	@echo "PDGui calculation for absorption of AlAs from Abinit...."
	@../../../pdgui -nosplash abinit BaTiO3.out results.xlsx -script script.py -exit
	@../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for BaTiO3 from Abinit...."
	@../../../pdgui -nosplash abinit BaTiO3.out results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "Viewing pdgui reference files for BaTiO3 from Abinit...."
	@../../../pdgui abinit BaTiO3.out -script script.py 

