pdgui:	
	@echo "PDGui Calculation for absorption of BaTiO3 from Abinit...."
	@../../../pdgui -nosplash abinit BaTiO3.out results.xlsx -script script.py -exit
	@../../../vibanalysis BaTiO3.out
	@../../../checkexcel results.ref.xlsx results.xlsx
	@diff BaTiO3.nma.ref BaTiO3.nma

pdgui-regenerate:	
	@echo "PDGui regeneration for absorption of BaTiO3 from Abinit...."
	@../../../pdgui -nosplash abinit BaTiO3.out results.xlsx -script script.py -exit
	@../../../vibanalysis BaTiO3.out
	@mv results.xlsx results.ref.xlsx
	@mv BaTiO3.nma BaTiO3.nma.ref

pdgui-view:	
	@echo "PDGui viewing for absorption of BaTiO3 from Abinit...."
	@../../../pdgui abinit BaTiO3.out -script script.py

