pdgui:	
	@echo "PDGui Calculating size effects for BaTiO3 from Abinit...."
	@../../../pdgui -nosplash abinit BaTiO3.out results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

benchmark:	
	@echo "PDGui benchmark size effects for BaTiO3 from Abinit...."
	@../../../pdgui -nosplash abinit BaTiO3.out -script script.py -exit

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for size effects BaTiO3 from Abinit...."
	@../../../pdgui -nosplash abinit BaTiO3.out results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "PDGui Regenerating reference files for size effects BaTiO3 from Abinit...."
	@../../../pdgui abinit BaTiO3.out -script script.py

