pdgui:	
	@echo "Calculating absorption for Na2(SO4)2 from PHONOPY...."
	@../../../pdgui -nosplash phonopy OUTCAR.born results.xlsx -script script.py -exit
	@../../../checkexcel results.xlsx results.ref.xlsx

pdgui-regenerate:
	@echo "PDGui Regenerating reference files for Na2(SO4)2 from PHONOPY...."
	@../../../pdgui -nosplash phonopy OUTCAR.born results.xlsx -script script.py -exit
	@mv results.xlsx results.ref.xlsx

pdgui-view:
	@echo "Regenerating reference files for Na2(SO4)2 from PHONOPY...."
	@../../../pdgui phonopy OUTCAR.born -script script.py

