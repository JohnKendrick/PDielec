pdgui:	
	@echo "PDGui Calculation for absorption of BaTiO3 from Abinit...."
	@../../../pdgui -nosplash abinit BaTiO3.out results.xlsx -script script.py -exit
	@../../../checkexcel results.ref.xlsx results.xlsx

pdgui-regenerate:	
	@echo "PDGui regeneration for absorption of BaTiO3 from Abinit...."
	@../../../pdgui -nosplash abinit BaTiO3.out results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:	
	@echo "PDGui viewing for absorption of BaTiO3 from Abinit...."
	@../../../pdgui abinit BaTiO3.out -script script.py

